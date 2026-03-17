# How Clann Calculates RF Distances Using Hashing

## Overview

The Robinson-Foulds (RF) distance between two trees counts the number of
bipartitions (splits) that appear in one tree but not the other. Clann
represents and compares bipartitions using 64-bit integer hashes rather than
explicit taxon-set strings. This makes comparison fast (integer equality
instead of set operations) and allows the final sorted-array intersection to
run in linear time.

---

## Step 1 — Assign a random weight to each taxon (done once at load time)

When source trees are loaded, every taxon index `i` gets a unique 64-bit
random-looking value using the **splitmix64** finaliser:

```
for each taxon i:
    x = i + 1
    x += 0x9e3779b97f4a7c15          # golden-ratio increment
    x = (x XOR (x >> 30)) * 0xbf58476d1ce4e5b9
    x = (x XOR (x >> 27)) * 0x94d049bb133111eb
    taxon_hash_vals[i] = x XOR (x >> 31)
```

**Why:** splitmix64 spreads bits uniformly across the full 64-bit range, so
that XOR-ing subsets of these values gives distinct results with overwhelming
probability. Two different subsets of taxa will almost certainly produce
different XOR sums.

---

## Step 2 — Represent a bipartition as a single 64-bit number

A bipartition divides the taxon set into two complementary sides, say **A**
and **B**. The hash of side A is:

```
hash_A = XOR of taxon_hash_vals[t] for every taxon t in side A
```

Because `A union B = all taxa in this tree`, we have:

```
hash_B = total_hash XOR hash_A
         where total_hash = XOR of taxon_hash_vals[t] for every taxon in the tree
```

The same bipartition can be written as A|B or B|A — it is unrooted, so the
two sides are interchangeable. To get a single **canonical** value that is the
same regardless of which side you are looking at:

```
canonical_hash = min(hash_A, hash_B)
```

Taking the smaller of the two values makes the hash orientation-independent.

**Trivial bipartitions** (a single taxon vs everything else) are excluded.
They are always shared between any two trees on the same taxa, so they carry
no RF information. In practice `canonical_hash == 0` is used as the sentinel
for trivial bipartitions and those entries are skipped.

---

## Step 3 — Extract all bipartitions from a Newick string with a stack

`collect_biparts_newick` reads the integer-indexed Newick string with a single
left-to-right pass and uses a small stack to accumulate XOR sums:

```
stack = [0]          # one slot per bracket depth; depth 0 = root
depth = 0

for each character in the Newick string:

    '('  →  push a new 0 onto the stack; depth++

    taxon_index  →  stack[depth] XOR= taxon_hash_vals[taxon_index]
                    # accumulate this leaf's weight into the current subtree

    ':'  →  skip everything until the next ',' or ')' or ';'
            # ignore branch lengths

    ')'  →  child_sum = stack[depth]; depth--
            complement = total_hash XOR child_sum
            canonical  = min(child_sum, complement)
            if canonical != 0:
                add canonical to output list   # record this bipartition
            stack[depth] XOR= child_sum
            # propagate the child subtree's XOR sum up to the parent level

output list is sorted ascending
return (output list, count)
```

**Why a stack works:** when you close a `)`  the top of the stack holds the
XOR of all taxon weights inside that subtree. That XOR sum *is* the bipartition
hash for one side. XOR-ing it into the parent level propagates the contribution
upward so that the parent's eventual XOR sum covers all taxa in the parent's
subtree. At the end of every `)`  you have captured exactly one non-trivial
internal edge as a canonical bipartition hash.

---

## Step 4 — Precompute source-tree bipartitions once before the search

Before the heuristic search begins, `rf_precompute_fund_biparts` calls
`collect_biparts_newick` on every source tree and stores the results in
`fund_bipart_sets[]`:

```
for each source tree i:
    if tree is excluded: skip

    total_hash = XOR of taxon_hash_vals[t] for every taxon present in tree i

    copy fundamentals[i] to a temporary string
    unroot the temporary string (remove bifurcating root if present)

    fund_bipart_sets[i].hashes = collect_biparts_newick(tmp, total_hash)
    fund_bipart_sets[i].count  = number of bipartitions found
```

This is done once and reused for every candidate supertree evaluated during
the search. For bootstrap, it is redone after each resample because the source
trees change.

---

## Step 5 — Extract bipartitions from the pruned supertree

For each source tree `i`, `compare_trees_rf` prunes the supertree down to
only the taxa present in source tree `i`, then extracts its bipartitions:

```
for each source tree i:

    prune_tree:   tag only the branches covering taxa present in tree i
    shrink_tree:  collapse internal nodes with fewer than 2 active children
                  (this handles partial taxon overlap cleanly)

    print_pruned_tree → pruned_nwk   (integer-indexed Newick string)
    while unroottree(pruned_nwk):    (normalise away any bifurcating root)

    total_hash = XOR of taxon_hash_vals[t] for taxa present in tree i
    super_bp[] = collect_biparts_newick(pruned_nwk, total_hash)
    super_cnt  = number of bipartitions in the pruned supertree
```

Pruning ensures the comparison is fair: you only compare the supertree on the
exact same taxon set as the source tree.

---

## Step 6 — Count shared bipartitions with a merge-walk

Both `super_bp[]` and `fund_bipart_sets[i].hashes[]` are sorted arrays of
`uint64_t`. Counting shared elements is a standard merge-walk (like merge sort
but only counting matches):

```
i = 0, j = 0, shared = 0
while i < super_cnt AND j < gene_cnt:
    if super_bp[i] == gene_bp[j]:   shared++; i++; j++
    elif super_bp[i] < gene_bp[j]:  i++
    else:                           j++
```

This runs in O(super_cnt + gene_cnt) — linear in the number of bipartitions.

---

## Step 7 — Compute the normalised RF distance

```
raw_RF    = super_cnt + gene_cnt - 2 * shared
            # bipartitions only in supertree  +  bipartitions only in gene tree
            # (symmetric difference)

max_RF    = 2 * (ntaxa_i - 3)
            # maximum possible RF for an unrooted binary tree on ntaxa_i leaves

normalised_RF = raw_RF / max_RF      (skipped if ntaxa_i < 4)

score_for_tree_i = normalised_RF * tree_weight[i]
```

`max_RF` is `2*(n-3)` because an unrooted binary tree on *n* leaves has
exactly `n-3` internal edges (bipartitions), so the two trees together have at
most `2*(n-3)` distinct bipartitions. Dividing by this makes the score 0
(identical) to 1 (maximally different), and comparable across source trees of
different sizes.

The **total score** returned by `compare_trees_rf` is the sum of
`normalised_RF * weight` over all source trees. The search minimises this sum.

---

## How ML scoring differs

`compare_trees_ml` is identical up to the normalisation step. Instead of
normalising, it uses the raw RF distance *d* in an exponential likelihood
model (Steel & Rodrigo 2008):

```
-ln L_i = beta * d_i               (mlscale=paper or lnl)
-ln L_i = beta * d_i * log10(e)    (mlscale=lust, Akanni et al. 2014)
```

The total score is the sum of `-ln L_i` across source trees, which is
minimised (equivalent to maximising the likelihood).

---

## Worked example — 4 taxa, 1 supertree vs 2 gene trees

### Setup

Four taxa with integer indices:

```
A = 0,  B = 1,  C = 2,  D = 3
```

For clarity the taxon hashes are kept small (real splitmix64 values are
64-bit numbers, but the arithmetic is identical):

```
h(A) = 5   (binary 0101)
h(B) = 3   (binary 0011)
h(C) = 6   (binary 0110)
h(D) = 9   (binary 1001)
```

All four taxa are present in every tree, so:

```
total_hash = h(A) XOR h(B) XOR h(C) XOR h(D)
           = 5 XOR 3 XOR 6 XOR 9
           = 6 XOR 6 XOR 9          (5 XOR 3 = 6)
           = 0 XOR 9                (6 XOR 6 = 0)
           = 9
```

Three trees to compare:

```
Supertree  S  :  ((A,B),(C,D))   — the AB|CD split
Gene tree G1  :  ((A,B),(C,D))   — same as supertree
Gene tree G2  :  ((A,C),(B,D))   — the AC|BD split
```

---

### Step A — Precompute gene-tree bipartitions

An unrooted binary tree on 4 leaves has exactly **1** non-trivial
bipartition (n − 3 = 1).  Both gene trees are first unrooted:

`((A,B),(C,D))` has a bifurcating root (two children) →
unroottree collapses it to `(A,B,(C,D))`.

`((A,C),(B,D))` similarly becomes `(A,C,(B,D))`.

#### Gene tree G1 — stack parse of `(0,1,(2,3));`

```
stack = [0],  depth = 0

'('   → depth=1,  stack = [0, 0]
'0'   → stack[1] ^= h(A)=5      stack = [0, 5]
','   → skip
'1'   → stack[1] ^= h(B)=3      stack = [0, 5^3=6]
','   → skip
'('   → depth=2,  stack = [0, 6, 0]
'2'   → stack[2] ^= h(C)=6      stack = [0, 6, 6]
','   → skip
'3'   → stack[2] ^= h(D)=9      stack = [0, 6, 6^9=15]

')'   → child_sum = 15,  depth=1
        complement = total_hash XOR child_sum = 9 XOR 15 = 6
        canonical  = min(15, 6) = 6   ← record this bipartition
        stack[1] ^= 15  →  stack = [0, 6^15=9]

')'   → child_sum = 9,  depth=0
        complement = 9 XOR 9 = 0
        canonical  = min(9, 0) = 0    ← SKIP  (trivial, sentinel value)
        stack[0] ^= 9  →  stack = [9]

';'   → done
```

`fund_bipart_sets[G1] = [6]`   (one bipartition, hash = 6)

**What does hash 6 mean?**  6 = h(A) XOR h(B) = 5 XOR 3.  It is the XOR
sum of the *AB side* of the split.  The other side C,D gives
h(C) XOR h(D) = 6 XOR 9 = 15.  canonical = min(6, 15) = 6.  So hash 6
uniquely represents the bipartition **AB | CD**.

#### Gene tree G2 — stack parse of `(0,2,(1,3));`

```
stack = [0],  depth = 0

'('   → depth=1,  stack = [0, 0]
'0'   → stack[1] ^= h(A)=5      stack = [0, 5]
','   → skip
'2'   → stack[1] ^= h(C)=6      stack = [0, 5^6=3]
','   → skip
'('   → depth=2,  stack = [0, 3, 0]
'1'   → stack[2] ^= h(B)=3      stack = [0, 3, 3]
','   → skip
'3'   → stack[2] ^= h(D)=9      stack = [0, 3, 3^9=10]

')'   → child_sum = 10,  depth=1
        complement = 9 XOR 10 = 3
        canonical  = min(10, 3) = 3   ← record this bipartition
        stack[1] ^= 10  →  stack = [0, 3^10=9]

')'   → child_sum = 9,  depth=0
        complement = 9 XOR 9 = 0
        canonical  = min(9, 0) = 0    ← SKIP
        stack[0] ^= 9  →  stack = [9]

';'   → done
```

`fund_bipart_sets[G2] = [3]`   (one bipartition, hash = 3)

**What does hash 3 mean?**  3 = h(A) XOR h(C) = 5 XOR 6.  This is the
XOR sum of the *AC side*.  The BD side gives 3 XOR 9 = 10.
canonical = min(3, 10) = 3.  So hash 3 represents **AC | BD**.

---

### Step B — Extract supertree bipartitions (pruned to each gene tree's taxa)

Both gene trees contain all four taxa, so no pruning is needed.
The supertree `((A,B),(C,D))` unroots to `(A,B,(C,D))` — the same string as
G1.  The stack parse is identical:

```
super_bp = [6]     (the AB|CD bipartition)
super_cnt = 1
```

---

### Step C — Count shared bipartitions

Both arrays are sorted.  The merge-walk compares them element by element.

#### Supertree vs G1

```
super_bp  = [6]
gene_bp   = [6]

i=0, j=0:  super_bp[0]=6  ==  gene_bp[0]=6  →  shared++
           shared = 1,  i=1,  j=1   (both exhausted)
```

#### Supertree vs G2

```
super_bp  = [6]
gene_bp   = [3]

i=0, j=0:  super_bp[0]=6  >  gene_bp[0]=3  →  j++
           j=1  (gene_bp exhausted, loop ends)
           shared = 0
```

---

### Step D — Compute RF distances and final score

For a 4-taxon unrooted binary tree: `max_RF = 2*(4−3) = 2`.

```
vs G1:
    raw_RF        = super_cnt + gene_cnt − 2*shared
                  = 1 + 1 − 2*1 = 0
    normalised_RF = 0 / 2 = 0.0     ← identical topology, perfect score

vs G2:
    raw_RF        = 1 + 1 − 2*0 = 2
    normalised_RF = 2 / 2 = 1.0     ← maximally different
                                       (the only possible split disagrees)

Total score (uniform weights) = 0.0 + 1.0 = 1.0
```

---

### Step E — Intuition check

There are only three possible unrooted topologies for 4 taxa:

```
AB|CD   (hash 6)   ← supertree and G1
AC|BD   (hash 3)   ← G2
AD|BC   (hash ?)   ← not used here
```

The supertree shares its one bipartition with G1 (RF = 0) and shares
nothing with G2 (RF = 1).  The total score 1.0 reflects one disagreement
out of two comparisons.  A supertree with topology AC|BD would score
1.0 + 0.0 = 1.0 — the same total, because the two gene trees conflict
and no single topology can satisfy both.

---

## Summary diagram

```
taxa loaded
    │
    └─ splitmix64 → taxon_hash_vals[0..n-1]   (one 64-bit value per taxon)

rf_precompute_fund_biparts   (once before search, or once per bootstrap replicate)
    │
    for each source tree i:
        total_hash = XOR of hash_vals for taxa in tree i
        stack-parse Newick → sorted list of canonical bipart hashes
        store in fund_bipart_sets[i]

compare_trees_rf   (called for every candidate supertree)
    │
    for each source tree i:
        prune + shrink supertree to taxa of tree i
        stack-parse pruned Newick → sorted super_bp[]
        merge-walk super_bp[] ∩ fund_bipart_sets[i].hashes → shared count
        RF_i = (super_cnt + gene_cnt - 2*shared) / max_RF
        score += RF_i * weight_i
    │
    return total score   (search minimises this)
```
