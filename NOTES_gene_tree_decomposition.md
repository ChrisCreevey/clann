# Gene-tree Decomposition into Ortholog Subtrees

`decomposegenetrees` (and its load-time counterpart, `exe ... autodecompose=yes`)
splits multicopy gene family trees into maximal single-copy-ish **ortholog
subtree fragments**, so that supertree criteria which assume single-copy input
(quartet-based criteria like `qfit`, and any future coalescent criterion) can
use signal from gene families that currently contribute nothing at all. This
note explains why that's needed, the three-stage algorithm, the design
decisions behind its defaults, and how the two commands fit into Clann's
existing multicopy-handling machinery.

---

## 1. Why decompose, instead of just filtering or pruning?

Clann already has two ways of dealing with multicopy gene family trees:

- **`criterion=recon`** uses *every* tree, scoring duplication/loss cost
  directly — no filtering needed, but this is a different objective from the
  quartet-based criteria.
- **The single-copy filter** (active automatically for `hs`/`nj`/`alltrees`
  when `criterion != recon`) excludes a multicopy tree from the search
  entirely. This is all-or-nothing: one in-paralog anywhere in a 50-taxon gene
  family tree discards the whole tree's signal.
- **`autoprunemono`/`prunemonophylies`** do better, but only for one specific
  case: *monophyletic* same-species clades (in-paralogs sitting together as
  sisters). They collapse those to one representative, promoting the tree to
  single-copy if that's all that was wrong with it. They cannot do anything
  for a tree where the duplicate copies of a species are **not** sisters —
  e.g. gene-tree error, or a genuine ancient duplication, has scattered them
  to different parts of the tree.

`decomposegenetrees`/`autodecompose` close that last gap. Where
`autoprunemono` can only collapse, this also **cuts**: it maps the gene tree
onto a guide species tree, finds well-supported duplication nodes, and splits
the family tree apart at genuine out-paralog boundaries — so a family that
`autoprunemono` cannot touch at all can still yield one or more usable,
single-copy ortholog fragments.

### Terminology

- **Orthologs** — genes related by a speciation event.
- **Paralogs** — genes related by a duplication event.
  - **In-paralogs** (relative to duplication node D) — copies that arose
    *after* D, redundant within one descendant lineage.
  - **Out-paralogs** (relative to D) — the two (or more) clades *on either
    side* of D. These must never be treated as orthologous to each other.
- Decomposition means: cut at out-paralog splits, collapse in-paralog
  redundancy, and do neither where the evidence is too weak or the resulting
  pieces would be uninformative.

---

## 2. The three-stage algorithm

### Stage 1 — Monophyly collapse (no guide tree needed)

Every source tree first goes through the same collapse logic as
`prunemonophylies`: species-specific (monophyletic in-paralog) clades are
identified topologically and pruned to one representative, using
`identify_species_specific_clades()` (`prune.c`) — this needs no reference
species tree at all, since "these leaves are all the same species and form a
clade" is a fact about the gene tree alone.

```
for each source tree:
    collapse every monophyletic same-species clade to one representative
    if the tree is now fully single-copy:
        -> one fragment, weight 1.0 (done, skip Stage 2)
    else:
        -> still has non-monophyletic duplicate species, hand off to Stage 2
```

`prune_monophylies()` (the `prunemonophylies` command) and
`decomposegenetrees`'s Stage 1 share one function,
`collapse_monophyly_in_tree()` (`prune.c`), so this is exactly the same
collapse a user already gets from `prunemonophylies`/`autoprunemono` — Stage 1
is not new logic, it's reuse of it.

### Stage 2 — LCA-mapping decomposition (needs a guide tree)

For a tree still multicopy after Stage 1, Clann's existing gene-tree/
species-tree reconciliation machinery (`label_gene_tree()`/
`reconstruct_map()`, used unmodified by `criterion=recon` and `reconstruct`)
is reused to find duplication nodes:

```
best_rooting = None
for each possible rooting of the gene tree:
    map the gene tree onto the guide species tree (LCA mapping)
    tag duplication nodes: a node is a duplication iff one of its
        children maps to the same species-tree node as the node itself
    count duplications
    keep this rooting if it has the fewest duplications seen so far
```

Duplication calls are rooting-dependent, so the rooting that minimises the
duplication count is used — the same "try every rooting, keep the best" idea
`get_recon_score()` already uses for `criterion=recon`, just scored on
duplication count alone rather than the full duplication+loss cost (losses
aren't needed for a cut/keep decision).

The tagged tree is then walked top-down. At each duplication node `D`:

```
walk(D):
    if D is not a duplication, or the call at D isn't well-supported:
        treat D as an ordinary node, recurse through it unchanged

    else if every leaf under D maps to a single species (pure in-paralog):
        collapse D's whole clade to one representative

    else:
        # D genuinely separates out-paralogs -- check each side independently
        for each child side of D:
            evaluate against minfragtaxa/minfragspecies

        if both/all sides pass the floor:
            CUT: each child becomes an independent fragment,
                 recurse into each one separately (repeat this same walk)

        else if both/all sides fail the floor:
            collapse the whole clade under D to one representative

        else (mixed):
            keep recursing into the side(s) that passed
            prune the failing side(s) to a representative and merge
                it into the passing side -- UNLESS that representative's
                species already survives on the passing side, in which
                case drop it instead of merging in a duplicate
```

A subtree whose root and every retained descendant is a speciation node (no
further duplication tags) is a finished ortholog fragment — recursion stops
there.

**The "mixed" branch's collision check** (drop rather than merge a
representative whose species is already kept on the passing side) exists
because a merge is not always safe just because it clears the informativeness
floor: merging a picked representative into an already-resolved sibling can
silently reintroduce the exact multicopy-ness the whole procedure exists to
remove, if that representative's species happens to already be present
elsewhere in the passing side (from a nested duplication decision made
earlier in the recursion). Since the entire point of `decomposegenetrees` is
to hand back single-copy fragments, this case is resolved conservatively —
drop the representative rather than ever emit a fragment with a duplicate
species in it.

As a final safety net, every fragment is swept once more immediately before
it's written out: if more than one surviving leaf for the same species still
made it through (possible when the colliding leaf sits entirely outside any
node that got flagged as a duplication — the walk above only makes decisions
*at* duplication nodes), only the first one encountered is kept.

#### The support gate

`node_well_supported(D, dupsupport)` decides whether a duplication call is
trustworthy enough to act on, using whatever confidence data is attached to
the branch above `D` in the gene tree's Newick string:

| Data present on the branch above D | Gate |
|---|---|
| A numeric label before `:` (bootstrap-style support) | trust the call iff `label ≥ dupsupport` (labels > 1.0 are treated as percentages and divided by 100) |
| Only a branch length (`:0.0123`, no label) | trust the call iff `branchlength > 0` — a literal zero/near-zero branch is treated as noise, `dupsupport` isn't used for this test |
| Nothing at all (no label, no branch length) | **trust the call** |

The "nothing at all ⇒ trust it" default deliberately mirrors an existing
Clann convention: `autoweight=bootstrap` already documents "trees with no
bootstrap labels are assigned weight 1.0" — absence of confidence data is
treated as full confidence, not as ambiguity, everywhere else in this
codebase. The alternative (treating "no data" as "not well-supported") would
make the cut logic essentially unusable on any tree that doesn't carry
bootstrap values or branch lengths, which describes a lot of real gene-family
input, including this feature's own worked example below.

`dupsupport` defaults to `0.5`.

#### The informativeness floor

`minfragtaxa` (default `4`) and `minfragspecies` (default `4`) gate whether a
candidate fragment is worth keeping on its own, or should be pruned down to a
representative and merged back into its sibling instead. `minfragtaxa=4`
reuses the same floor `qfit` already applies internally (`scoring.c`, quartet
scoring is undefined below 4 taxa) — exposing it as an option here rather
than adding a second hardcoded literal for the same constraint.

`minfragspecies` defaults to `4` as well, **not** the smaller value you might
expect for a secondary check — the reasoning is worth spelling out, because
it's not obvious from the option's one-line description alone. A fragment can
satisfy `minfragtaxa=4` with, say, 4 leaves from only 2 species (two
in-paralog pairs that Stage 1's topology-only monophyly test happened to
miss). Since Clann's delimiter mode already collapses every gene copy of a
species to the same taxon identity for any species-tree-level comparison
(`qfit` and friends), such a fragment is exactly as phylogenetically
uninformative — a degenerate quartet — as the raw `<4`-taxa case `minfragtaxa`
was already there to prevent. Setting `minfragspecies=4` closes that gap: the
floor is really "4 distinct species", and `minfragtaxa` alone doesn't
guarantee that.

### Stage 3 — Weighting

A family tree that survives Stage 1/2 as `k` fragments has each fragment
weighted `1/k`, so a family that duplicated heavily and produced many small
fragments doesn't out-vote a family that never duplicated and contributes one
fragment at weight 1.0. This is deliberately the simplest possible scheme —
a size- or informativeness-proportional weighting is a plausible future
refinement, not implemented here.

---

## 3. Guide-tree resolution

Stage 2 needs a species tree to map gene trees onto. Resolution order
(`resolve_guide_tree()`, `reconcile.c`), matching `reconstruct`'s existing
`speciestree=` option for consistency:

```
1. speciestree=<file> supplied
       -> read it, hard-error if its taxon set doesn't exactly match
          the loaded data (both directions: nothing missing, nothing extra)
2. speciestree=memory, or nothing supplied but a supertree is already
   in memory (trees_in_memory > 0)
       -> use retained_supers[0]
3. Nothing available
       -> generate one automatically with nj()
```

Step 3 is the one deliberate difference from `reconstruct`'s own behaviour,
which hard-errors if nothing is available. `decomposegenetrees` should never
be fatal for lack of a species tree, since the whole point is that it can run
standalone with zero setup — `nj()` (an existing, already-used function) is a
perfectly serviceable fallback.

---

## 4. Two commands, two destructiveness levels

Mirroring the existing `prunemonophylies`/`autoprunemono` split exactly:

| | `decomposegenetrees` | `exe ... autodecompose=yes` |
|---|---|---|
| When it runs | On demand, standalone command | At load time, as part of `exe` |
| Destructive? | No — only writes `<filename>` + `<filename>_info.txt` | Yes — commits fragments into `fundamentals[]` immediately |
| Adopting the result | `exe <filename>` afterwards | Already live |

### Why the destructive path reloads through `exe`, not a hand-rolled update

Decomposition changes `Total_fund_trees` (one family tree becomes `k`
fragments), which invalidates every piece of per-tree derived state computed
at load time — `presence_of_taxa[][]`, the RF bipartition cache, the `dfit`
distance matrices, `sourcetree_scores[]`, `tree_weights[]` sizing, taxon
hashes, and so on. Re-deriving all of that by hand in a second code path,
parallel to what `execute_command()` already does for a freshly loaded file,
would duplicate logic Clann already trusts to do this correctly. So
`autodecompose_apply()` instead: writes the fragments to
`autodecomposed_fragments.txt` using the same `<filename>`/`<filename>_info.txt`
machinery `decomposegenetrees` uses, then calls `execute_command()` on that
file — the exact same registration path a fresh `exe <file>` takes. This is
the only place in the implementation that reuses Clann's own file-loading
parser rather than working with the in-memory tree structures directly.

### Getting a per-fragment weight into that file

Clann's plain-Newick ("Phylip format") reader already understands a bracketed
float placed **before** a tree's opening `(`, parsed straight into
`tree_weights[]`:

```
[0.333333](A,(B,C));[fragment_name]
```

— so Stage 3's `1/k` weights ride through the existing loader for free, no
parser changes needed. There is one quirk to work around: the code that runs
*before* the very first tree in a file unconditionally treats any leading
`[...]` as a discardable comment, so **the first tree in a file can never
pick up a weight via this annotation** — it silently keeps its default of
1.0. `build_decompose_output_text()` works around this by always writing a
naturally weight-1.0 fragment first when one exists (a Stage-1 passthrough,
or any family that didn't need splitting) — true for realistic datasets,
since most families in a run are single-copy already. For `autodecompose`
specifically this quirk doesn't matter in practice anyway: after the reload,
the correct per-fragment weight is re-applied explicitly in memory from
Stage 3's own record, matched by fragment name rather than relying on the
file annotation surviving intact.

### Reconciliation still needs the real, undecomposed trees

`reconstruct` (and `criterion=recon`) exist specifically to report
duplication/loss counts — reconciling against already-decomposed,
already-single-copy fragments would silently report near-zero duplications
everywhere, which defeats the purpose of running them at all. `autodecompose`
therefore keeps a full snapshot of the pristine (pre-decomposition) gene-tree
set — as resolved full-name Newick text, taken immediately before Stage 1
runs — independent of `fundamentals[]`'s own lifetime, so it survives the
`execute_command()` reload that follows. `reconstruct` checks a
`decompose_active` flag (mirroring `autoprunemono_active`'s existing
precedent exactly) and, when set, temporarily reloads that pristine snapshot
in place of the live decomposed pool for the duration of the reconciliation,
then restores the decomposed pool afterward — printing an equivalent message
to the existing `"(Using original unpruned trees for reconstruct)"` one.

No separate "undo" command is needed: decomposition never modifies the file
the user originally loaded from, so `exe <original file>` — a plain reload —
always fully restores the pre-decomposition state from scratch, in any
session.

---

## 5. Worked example — `examples/tutorial_multicopy.ph`, tree 6

```
((((Human.alpha,Chimp),(Human.beta,Gorilla)),Orangutan),Macaque);
```

`Human.alpha` and `Human.beta` are *not* sisters — one pairs with Chimp, the
other with Gorilla — so Stage 1's monophyly test doesn't fire at all; this
tree reaches Stage 2 unchanged. Running

```
decomposegenetrees minfragtaxa=2 minfragspecies=1
```

(a looser floor than the default, chosen only because this fixture's
families are too small — 6-9 leaves — to clear the real `minfragtaxa=4
minfragspecies=4` defaults) produces:

```
Tree # 6 [  ]: cut at duplication node -> 2 new fragment(s)
```

and three surviving fragments in `decomposedtrees.txt`:

```
[0.333333]((Human.beta,Gorilla));[tree6_frag1]
[0.333333]((Orangutan,Macaque));[tree6_frag2]
[0.333333]((Human.alpha,Chimp));[tree6_frag3]
```

Each is single-copy, each is weighted `1/3` so the family's total
contribution still sums to 1, and together they salvage real phylogenetic
signal (the Human/Chimp and Human/Gorilla speciation relationships) from a
tree that the existing single-copy filter would otherwise have discarded
completely.

Tree 7 of the same fixture —

```
(((Human.a,(Chimp,Gorilla)),(Human.b,(Orangutan,Macaque))),(Mouse,Rat));
```

— demonstrates the mixed cut/merge branch: it cuts once, then on one side
prunes an uninformative 1-taxon sibling clade (`Human.a`, below even the
loosened `minfragtaxa=2` floor here) and merges its representative into the
passing side, while the other spun-off side is dropped entirely for still
being too small:

```
Tree # 7 [  ]: cut at duplication node -> 2 new fragment(s)
Tree # 7 [  ]: pruned uninformative sibling clade under duplication node
    (< minfragtaxa=2/minfragspecies=1) -> kept Human.b, merged into sibling
Tree # 7 [  ]: dropped a 1-leaf fragment (below minfragtaxa=2)
```

---

## 6. Summary diagram

```
exe <file>                                    exe <file> autodecompose=yes
    │                                                  │
    ▼                                                  ▼
fundamentals[] loaded normally               fundamentals[] loaded normally
                                                         │
decomposegenetrees [options]                  Stage 0: snapshot pristine
    │                                              trees (pre_decompose_*)
    ├─ resolve_guide_tree()                          │
    │                                          resolve_guide_tree()
    ├─ Stage 1: collapse_monophyly_in_tree()          │
    │      per source tree                    Stage 1 + Stage 2 + Stage 3
    │                                              (same functions as the
    ├─ Stage 2: decompose_gene_tree_stage2()          standalone command)
    │      per still-multicopy tree                   │
    │                                          write autodecomposed_fragments.txt
    ├─ Stage 3: 1/k weighting                          │
    │                                          execute_command(fragments file)
    ├─ write <filename> + <filename>_info.txt     — the SAME registration path
    │      (fundamentals[] untouched)               a fresh `exe` uses —
    │                                              presence_of_taxa[][], RF
    └─ print "exe <filename> to adopt"             bipartition cache, dfit
                                                    matrices, tree_weights[]
                                                    sizing all rebuilt for free
                                                         │
                                              decompose_active = 1
                                              fragments live in fundamentals[]
                                                         │
                                    ┌────────────────────┴───────────────────┐
                                    ▼                                        ▼
                          hs / nj / alltrees                          reconstruct
                    use the decomposed fragment pool          swaps to pre_decompose_*
                                                                for real dup/loss counts,
                                                                then swaps back
```

---

## 7. Explicitly out of scope

- The deep-coalescence/MDC criterion itself — this feature only produces
  clean ortholog-subtree input for existing criteria (`qfit` today, `mdc` if
  it's built later); no scoring code is touched.
- Size- or informativeness-proportional fragment weighting — ship equal-split
  (`1/k`) first.
- Iterative refinement (decompose → build supertree → re-decompose against
  the improved guide tree → repeat) — worth revisiting once the single-pass
  version has been used on real data.

See [USER_MANUAL.md](USER_MANUAL.md) (`Autodecompose`/`decomposegenetrees`
sections) for the full command reference, and
[TUTORIAL.md](TUTORIAL.md) (Part 5) for a hands-on walkthrough.
