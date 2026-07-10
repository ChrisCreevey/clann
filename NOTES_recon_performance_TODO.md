# DONE: Speed up the recon-criterion search (linear-DP dup scans)

**Status:** implemented (2026-07-10). Both O(gene²)-with-copies duplication-scan
loops in `get_recon_score` were replaced with a single `lr_root_counts()`
linear-DP call each (§2 below). Verified identical to the old scan (0 mismatches
over 153,218 gene-mindup and 619,665 species-mindup rooting comparisons across
`tutorial_multicopy`, `tutorial_single`, `tutorial_trees`, `bsweight_demo`),
correctness targets met (`hs` reaches 17 and `reconstruct` reports 17 on
`tutorial_multicopy`), and measured wall-clock speedups (5 reps, nthreads=1):

| Dataset | Before | After | Speedup |
|---|---|---|---|
| `tutorial_multicopy.ph` (8 trees) | 11.8 s | 3.0 s | ~4.0× |
| `tutorial_single.ph` (28 trees) | 40.3 s | 7.3 s | ~5.5× |

The larger dataset benefits more, as expected: the species pre-pass amplifier
drops from O(species × genes × gene²) to O(species × genes × gene).

**Note on reproduction:** the correct command sequence is `set criterion=recon`
then `hs` — `hs criterion=recon` does *not* parse the criterion (heuristic_search
reads the global `criterion`, which only `set` updates), so it silently runs as
DFIT and the recon path (and its single-copy-filter skip) never executes.

The remainder of this brief is retained for context.

---

## 1. The finding (already assessed — don't redo this part)

A user reported the machine running hot during `hs criterion=recon`. The obvious
hypothesis was excessive `malloc`/`free` (which had caused heat before). **It is
NOT `malloc` this time.** Evidence:

- One candidate scoring allocates **31,775 tree nodes** (`make_taxon`) plus 528
  small arrays — huge *counts*, on a 9-taxon dataset.
- But a sampling profile of the reconciliation shows `malloc`/`free`
  (`szone_malloc_should_clear`, `free_tiny`) is only **~4%** of time. The hot
  functions are `duplicate_tree`, `dismantle_tree`, `add_losses`, and
  `make_taxon`'s field-init — i.e. **per-rooting tree copying and
  reconciliation**, not allocation.
- A node-recycling pool (reuse freed nodes instead of malloc/free) was
  implemented and measured: **~0% on the search, ~3.5% on a 30-taxon
  reconstruct.** It was reverted (marginal, and immediate reuse risks surfacing
  latent use-after-free right before a merge).

**Conclusion:** the cost is the reconciliation *compute*, specifically the
minimum-duplication rooting scans that reroot the gene tree to **every branch**
and **rebuild the whole tree from a copy each time** (O(gene²) per gene tree),
repeated for every species rooting in the species pre-pass.

## 2. The proposed optimisation

Replace the reroot-to-every-branch duplication scans with the **linear-time
all-rootings duplication DP** that already exists in this file for `decompose`:

- `lr_root_counts(struct taxon *gene_top, int n, int xnum, int *out)`
  — `reconcile.c:1074`. It fills `out[j]` = duplication count when the gene tree
  is rooted at branch `j`, for **all** `n` rootings, in **O(gene) total**, from a
  single fixed tree (two passes computing `mdown`/`mup` LCA maps) — **no
  rerooting, no per-rooting `duplicate_tree`/`dismantle_tree`**. Helper functions:
  `lr_dupsum` (1009), `lr_dp_down` (1061). It's already used by decompose at
  `reconcile.c:1176` and was validated there (0 mismatches vs `tag_duplications_only`
  over 80 trees, all rootings).

The recon path currently computes the same per-rooting dup counts the slow way in
two places inside `get_recon_score` (`reconcile.c:2440`):

1. **Gene-rooting Pass 1** (`mindup_mode`): a `for` loop that `reroot_tree`s to
   each branch, calls `recon_dups_only()` (`reconcile.c:2434` = `label_gene_tree_rec`
   + `tag_duplications_only`), then rebuilds the gene tree from `copy` and
   `number_tree`s it — per rooting. Search for `mindup_mode` and the `dupc[]`
   array. This is the O(gene²)-with-copies loop to replace with one
   `lr_root_counts()` call.
2. **Species pre-pass** (`spec_mindup_mode`, the block at `reconcile.c:2501`):
   for **every** species rooting it re-runs that same per-gene gene-rooting dup
   scan to build `spec_dupc[m]`. This is the biggest amplifier — replace its
   inner gene-rooting scan with `lr_root_counts()` too.

`add_losses` (the full dup+loss scoring, Pass 2) still needs the actual rerooted
tree, but it runs only on the ~40% minimum-duplication rootings, so it stays as
is — the win is removing the *dup-scan* rerooting/copying.

## 3. What to watch out for / validate

- **Do the counts match?** `lr_root_counts` must produce the same per-rooting
  duplication counts as the current `recon_dups_only()` path *in the recon
  context*. The recon gene trees have been through `do_resolve_tricotomies`
  (soft min-dup polytomy resolution) and are unrooted (top is a sibling forest —
  which `lr_root_counts` models with its "virtual centre"). Add a temporary
  verification: for each gene tree, compute `dupc[]` both ways and assert equal,
  over the tutorial and a larger set, before deleting the old loop (this mirrors
  how the decompose linear DP was validated).
- **`g_lr_*` globals are not `threadprivate`.** They're fine for recon (which is
  single-threaded — `nthreads` is disabled for recon) and for decompose
  (single-threaded). If recon is ever parallelised, make them threadprivate.
- **Species tree state:** `lr_root_counts` needs `build_species_partag()` done
  and `xnum` = the species root sentinel. In `get_recon_score` these are already
  set per species rooting (`sp_xnum`, the `build_species_partag` hoist).
- **Correctness target:** on `examples/tutorial_multicopy.ph`, `hs criterion=recon`
  must still reach the true optimum **17**, and `reconstruct` on the best tree
  must still report **17** (they were aligned — see
  `NOTES_species-gene-tree-reconciliation.md`).

## 4. How to reproduce the measurements

```
# hot workload (stochastic — average several runs):
hs criterion=recon   on examples/tutorial_multicopy.ph, nreps=5 nthreads=1

# deterministic heavy workload:
reconstruct speciestree=<fixed tree> on a ~30-taxon synthetic set
```
Profile with macOS `sample <pid> 3` while a search runs to see the hot functions.
Instrument `make_taxon` (there's a `malloc_check`/`count_now` counter already) to
count node allocations per `get_recon_score` call.

## 5. Context files

- `NOTES_species-gene-tree-reconciliation.md` — the rooting/mindup design and why
  it filters by duplications (which is exactly what makes the linear DP
  applicable here).
- `NOTES_gene_tree_decomposition.md` — where the linear DP is already in use.
- The reconciliation cost is `dup_weight·dups + loss_weight·losses`; losses have
  **no** clean linear all-rootings form in Clann's hand-rolled model (that dead
  end is documented in the reconciliation NOTES), which is *why* only the
  duplication scan is being optimised.
