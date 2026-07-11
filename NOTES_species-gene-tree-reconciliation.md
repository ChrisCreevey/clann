# Species/Gene-tree Reconciliation and Rooting in Clann

The `recon` criterion scores a candidate supertree by **reconciling** every
source gene tree against it — mapping each gene tree onto the species tree and
counting the **duplications** and **losses** needed to explain it. The score is
`duplication_weight × duplications + loss_weight × losses`, summed over all gene
trees; the heuristic search (`hs criterion=recon`) minimises it, and the
`reconstruct` command reports and illustrates it.

Reconciliation only makes sense on **rooted** trees, but the input gene trees —
and the candidate supertrees the search generates — are **unrooted**. So every
reconciliation has to pick a rooting, and that choice turned out to dominate both
the accuracy and the speed of the whole criterion. This note explains the rooting
problem, the `mindup` solution that is now the default, the polytomy-resolution
choice underneath it, how `hs` and `reconstruct` are kept consistent, and the
approaches we tried and rejected along the way.

---

## 1. The rooting problem

Scoring one candidate supertree means, for every gene tree, choosing:

- a **gene-tree rooting** — where to root the unrooted gene tree (≈ 2N−3
  possible branches), and
- a **species-tree rooting** — where to root the candidate supertree (also
  ≈ 2N−3 branches).

The reconciliation cost is the minimum over these choices. Trying **every**
combination (`numspeciesrootings=all numgenerootings=all`) is exact but costs
O(N²) reconciliations *per candidate*, and a search visits thousands of
candidates — so exhaustive rooting is far too slow for a normal search.

The old default avoided that cost by **randomly sampling** just 2 species × 2
gene rootings (`2×2`). That was a bad trade:

- **It badly overestimates the score.** On the tutorial dataset a tree whose
  true (exhaustive) cost is **17** scored **29–47** under `2×2`, run to run.
- **The error is noise, not a constant offset.** Because it's random, two
  scorings of the *same* tree disagree by more than the real difference between
  a good and a bad tree. During a hill-climb that **misranks candidates**, so
  the search follows the wrong gradient and settles on genuinely worse trees
  (true cost 18–23 instead of 17).

So the search was being steered by rooting noise, not by the biology.

---

## 2. The `mindup` solution (the default)

The key insight: **you don't need to score every rooting to find the best one.**
Duplications are cheap to count and they *constrain* where the good rootings are.

For each gene tree we:

1. **Count duplications at every rooting — cheaply.** A linear-time
   reroot algorithm gives the duplication count of all rootings in one pass,
   without running the expensive loss reconstruction (`add_losses`).
2. **Keep only the minimum-duplication rootings.** These are the candidates for
   the true optimum.
3. **Score the full dup+loss cost only on those**, and take the best.

The same idea is applied to the species-tree dimension: a cheap duplications-only
pass finds the minimum-duplication species rootings, and the full cost is
computed only at those.

**Why this works.** The rooting that minimises duplications+losses almost always
also minimises duplications (measured: 88–93% of the time it is *exactly* a
minimum-duplication rooting, and 99.99% within one duplication). So restricting
to those rootings loses almost nothing, while skipping the expensive loss
computation on the ~60% of rootings that can't be optimal.

**What it delivers (measured on `examples/tutorial_multicopy.ph`):**

| Rooting mode | True quality of tree the search finds | Time (one search) |
|---|---|---|
| `2×2` random (old default) | 18–23 — misses the optimum | < 1 s |
| **`mindup` (new default)** | **17 — the true optimum, every run** | **~5–7 s** |
| `all` (exhaustive) | 17 — the true optimum | > 110 s (didn't finish) |

So `mindup` gives **exhaustive-quality search outcomes deterministically, more
than 15× faster than exhaustive**, and it replaces the noise of random sampling
with a principled, repeatable choice.

**One honest caveat.** The *gene*-rooting restriction is very accurate (rank
correlation 0.974 vs exhaustive). The *species*-rooting restriction is a weaker
approximation (0.806) — it can misjudge some trees badly. But those trees are
always *already-bad* candidates (which the search rejects anyway); every
near-optimal candidate is scored exactly. So it doesn't mislead the search where
it actually matters. If you want the strongest guarantee, use
`numspeciesrootings=all numgenerootings=mindup` (exhaustive species, fast gene).

---

## 3. Quick reference — the rooting settings

Both settings apply to `hs criterion=recon`. `reconstruct` always roots
exhaustively (see §5) and ignores them.

| Value | Meaning | When to use |
|---|---|---|
| `mindup` *(default)* | Restrict scoring to the minimum-duplication rootings. Near-exhaustive accuracy, deterministic, fast. | The default; the right choice for almost all analyses. |
| `all` | Try every rooting (exhaustive). Exact but O(N²) — slow. | Small/medium datasets, or to produce a guaranteed-exact reference. |
| `<integer>` | Sample that many **random** rootings. Fast but noisy. | Legacy behaviour; not recommended (this is what caused the mis-ranking). |

```
hs criterion=recon                                   # mindup / mindup (default)
hs criterion=recon numspeciesrootings=all numgenerootings=mindup   # strongest gene guarantee
hs criterion=recon numspeciesrootings=all numgenerootings=all      # fully exhaustive (slow, exact)
```

Related fix worth knowing: `numgenerootings=all` used to be silently broken (a
copy-paste bug meant it did *one random* gene rooting instead of all). It now
genuinely tries every gene rooting.

---

## 4. Polytomy resolution — the soft, minimum-duplication choice

Gene trees often contain **polytomies** (multifurcations — nodes with more than
two children, usually from collapsed low-support branches). Reconciliation needs
a binary tree, so a polytomy must be resolved into some binary shape first, and
*different* resolutions give *different* duplication/loss counts.

Clann used to resolve polytomies by a species-tree **path-distance** heuristic,
which could force duplications that a smarter resolution would avoid. We replaced
that with a **soft, minimum-duplication resolution**: a polytomy is treated as
uncertainty about the true binary shape, and it's resolved into the shape with
the **fewest duplications**.

- The rule is a greedy one — repeatedly join the pair of subtrees whose
  common ancestor is a *speciation* (not a duplication), deepest first.
- We proved this attains the true minimum: it was validated against exhaustive
  enumeration of all resolutions plus an independent counting method over
  **600,000 random test cases** with zero mismatches.
- On binary gene trees it does nothing, so those results are unchanged.

This makes the reconciliation score *more accurate* (fewer spurious
duplications) and, because it needs the species tree's parent/depth tables
anyway, it also let us drop an old per-tree O(species²) distance-matrix step.

---

## 5. Keeping `hs` and `reconstruct` consistent

A common workflow is: find the best supertree with `hs`, then run
`reconstruct ... showrecon` to *illustrate* the duplications and losses on that
tree (ASCII drawing, per-branch tables, NHX output). For that to be meaningful,
`reconstruct`'s numbers must match the score `hs` optimised.

They didn't, originally — the same tree could score 17 in `hs` but 30+ in
`reconstruct`. Both use the identical duplication/loss counter, but `reconstruct`
prepared the trees differently:

1. **It wrapped both the gene tree and the species tree in an extra "root"
   node** before reconciling. That artificial node maps to the same species node
   as its single child, so it was miscounted as a spurious duplication (plus
   losses) on every gene tree.
2. **Its rooting-selection pass didn't resolve polytomies** (while `hs` and
   `reconstruct`'s own illustration pass did).

Both are now fixed: the wrappers are removed and polytomies are resolved in both
passes, so `reconstruct`'s total equals the `hs` score for the best tree, and the
per-tree numbers sum to it. The illustration outputs (ASCII, per-branch `.recon`
table, NHX) all still work.

**Note on rooting:** `reconstruct` still roots **exhaustively** (every species ×
gene rooting), whereas `hs` uses `mindup`. For the *best* tree they agree exactly
(mindup is exact at the optimum). If you deliberately reconstruct a clearly
*suboptimal* tree, its exhaustive total can come in slightly below the `mindup`
score for that tree — expected, and only for non-optimal trees.

---

## 6. Things we tried and rejected

Recording these so they aren't re-attempted:

- **Two-phase search (explore on duplications only, then finish on dup+loss).**
  The idea was to explore tree space cheaply using duplications alone, then
  polish with the full objective. It was **~4× slower with no accuracy gain**:
  the duplication-only optimum sits in a *different basin* from the dup+loss
  optimum (it scored 20 under dup+loss — a trap a local search can't escape), so
  the finish phase needed a full multi-restart search anyway, making the cheap
  first phase pure overhead.

- **Making exhaustive rooting the default.** Accurate but too slow to be a
  blanket default (O(N²) per candidate). The obvious constant-factor speedups
  only bought ~3% — the real cost is the loss reconstruction, not the setup —
  so exhaustive-by-default would make searches unusable. `mindup` is the answer
  instead.

- **A closed-form linear-time *loss* count.** Duplications have a clean
  linear-time all-rootings algorithm; losses do not, in Clann's hand-rolled loss
  model (the correction for duplications is not local). Speciation edges follow
  a simple depth formula, but multi-copy overlaps don't — so there's no exact
  linear loss DP to lean on. This is *why* `mindup` filters by duplications (which
  are cheap and linear) rather than trying to compute losses for all rootings.

---

## 6b. Two loss models: `legacy` (default) and `standard`

Clann's original ("legacy") loss count *reconstructs* the missing species
subtrees (`add_losses`/`construct_tree` build the implied lost lineages with
`make_taxon`, then `count_losses` counts the maximal fully-lost subtrees). It is
self-consistent (identical gene/species tree → 0) but idiosyncratic in two ways
that make its numbers **not comparable to other reconciliation tools**:

- it counts overlapping multi-copy losses as a **union** (shared region once)
  rather than the standard per-lineage **sum**; and
- it charges losses at the unrooted forest top per-sibling, unlike the textbook
  rooted model.

`lossmodel=standard` (opt-in on both `hs` and `reconstruct`; `legacy` stays the
default) instead scores the **textbook LCA-mapping duplication-loss model** used
by NOTUNG/DupTree/ete3: for each internal gene node mapped to species node `T`
with children mapped to depth `d_i`,
`losses += Σ_i(d_i − depth(T)) − (#children at speciations)`, duplication iff a
child maps to `T`; the gene tree's root (the implicit parent of the top-level
forest siblings) is charged too. It runs in a single **O(N) pass with no
allocation**, so it is also **faster** than the legacy reconstruction (~1.6–1.8×
on the tutorials). It was **cross-validated exactly against ete3** (0 mismatches
on dups *and* losses over 1700+ reconciliations, single- and multi-copy, at every
rooting the search visits).

The one implementation subtlety: Clann's rooted species tree is stored as a
2-component *forest* with the true root implicit, and the LCA sentinel reuses the
`xnum` tag of a real top component — so a root-spanning gene node and a node
mapping to that real component share a tag. `node_comp()` resolves this by
deriving each gene node's forest component from its leaves (spanning → the
implicit root at depth −1), independent of the collided tag.

**Illustration (`reconstruct`).** The arithmetic scorer does no tree
reconstruction, so `reconstruct`'s drawing/NHX/`.recon` outputs can't be read off
it. `annotate_standard()` builds a standard-model reconciled tree just for the
illustration: it wraps both trees in explicit roots (so the gene-tree-root
duplication is counted) and inserts losses **per edge, per copy** — `std_mirror`
walks the species path `T→M(child)` inserting each off-path species subtree as a
lost lineage, and duplication copies each descend the whole `T` subtree
independently. This is exactly the fix for the legacy under-count: legacy
`add_losses` reconstructs a node's subtree once with all children present (so a
duplication's overlapping losses **merge**, a union), whereas per-edge insertion
never merges (a sum). Validated: `count_losses()` on the annotated tree and the
duplication marks equal the arithmetic standard score, 0 mismatches over 80,000
reconciliations. So under `lossmodel=standard` the ASCII/NHX/`.recon` events match
the reported score, exactly as legacy's do (both count a lost *clade* as one loss
even though it draws as several `*LOST` leaves).

**Caveats.** Switching **changes the scores** (e.g. `tutorial_multicopy.ph`: best
tree = 30 under `standard` vs 17 under `legacy`) and can change which supertree
the search selects, so the two models are not comparable — keep `legacy` to
reproduce earlier results. And `mindup` is a *looser* proxy under `standard`
(gene/species min-duplication rootings track the min-dup+loss rooting less
tightly than in the legacy model — on the tutorial `hs` finds 32 while exhaustive
`reconstruct` finds 30), so pair `standard` with `numspeciesrootings=all
numgenerootings=all` when you need the exact standard optimum. This section is
also *why* §6's "no clean linear loss form" caveat is specifically about an
**all-rootings** DP for Clann's *legacy* model — the standard *per-rooting* loss
is a clean O(N) formula, which is exactly what `lossmodel=standard` uses.

---

## 7. Practical recommendations

- **Just use the defaults.** `hs criterion=recon` now gives accurate,
  deterministic, fast reconciliation search out of the box.
- **For maximum gene-rooting confidence** on an important result, run
  `numspeciesrootings=all numgenerootings=mindup`.
- **For a guaranteed-exact reference** on a small/medium dataset, run
  `numspeciesrootings=all numgenerootings=all` — and expect it to be slow.
- **To illustrate a result**, `reconstruct speciestree=memory showrecon=yes`
  after the search now reports numbers consistent with the `hs` score.

---

## 8. Open / related items

- The **species-rooting `mindup` approximation** is the weakest link (§2); if it
  ever proves too coarse on real data, the fallback is exhaustive species rooting
  (`numspeciesrootings=all`), and a stronger species-rooting method would be the
  natural next improvement.
- Run-to-run variation in `reconstruct` on multicopy data comes from
  **paralog-representative selection**, not rooting (see
  `NOTES_tobeaddressed.md`, issue #5).
- `hs start=<file>` now does proper multi-restart (honours `nreps`), so you can
  seed a recon search from your own starting tree and it will reach the optimum.
