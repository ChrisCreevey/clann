# Clann — Issues To Be Addressed

Running log of problems discovered incidentally while working on other tasks.
Each entry is something worth a dedicated investigation in a future session —
they are **not** blocking the current work, but should not be forgotten.

Format: newest issues can go at the top or bottom; keep the discovery date and
enough context (files, symptoms, how it was found) to pick it up cold.

---

## 1. `start=<file>` heuristic search does not actually search

**Discovered:** 2026-07-03 (while building fixed-start test cases for ILS)

**Symptom:** Running the heuristic search with a user-supplied starting tree
file, e.g.

```
exe synthA.ph
set criterion=ml
hs nreps=1 nthreads=1 swap=spr start=start_0.ph
```

produces **no `rep 1/1 ... converged` progress line** and returns a score that
is essentially the starting tree's own score (observed `lnL = -154.0` on a
16-taxon concordant dataset whose true optimum is `0`). By contrast,
`start=random` and `start=nj` on the same data run the full search and reach
`0` (or good local optima). So the `start=<file>` path (`start == 1` in
`heuristic_search`, `treecompare2.c`) appears to read the tree but then skip or
short-circuit the actual `do_search()` hill-climb.

**Where to look:** `heuristic_search()` in `treecompare2.c` — the `start == 1`
branch that opens `userfile` / reads `userfilename`, and how it feeds the
user's tree into `do_search()` compared with the `start == 0` (random) and
`start == 2` (NJ) branches. Suspect the user-file tree is scored/printed but
not passed through the same `do_search()` loop, or `nreps`/loop wiring differs.

**Impact:** Users cannot reliably start a heuristic search from their own tree
(e.g. to refine a specific hypothesis). Also blocks using fixed starting trees
as a clean, deterministic test harness for search-algorithm changes.

---

## 2. Random seed does not make runs deterministic

**Discovered:** 2026-07-03 (while trying to build a paired ILS on/off test)

**Symptom:** With an explicit `seed=<n>` set, repeated runs of the *same*
command give *different* results. Example: `hs ... start=random sample=1
ils=15 seed=3` over six runs produced `0, 0, 0, -8, -4, 0`. The variation is in
the **random starting-tree generation**, not the search: with a deterministic
start (`start=nj seed=3`) the search is fully reproducible (`0,0,0,0`).

**Likely cause:** `random_star_decom()` (and possibly other spots) use the
non-thread-safe global `rand()`, which is not seeded from the user `seed=`
value. `seed=` is wired into `thread_seed` (used by `rand_r()`) in some paths
(`treecompare2.c` ~lines 1046, 3222, 4508) but `rand()` itself is never
`srand()`-ed from it. So anything using `rand()` is effectively seeded
nondeterministically (or from a fixed default, depending on libc) and is not
controlled by `seed=`.

**Where to look:** all `rand()` call sites (grep `rand()` in `treecompare2.c` —
e.g. `random_star_decom`, and lines ~1709, 1928, 3522, 6678, 6800, 6875), and
where/whether `srand()` is ever called. Decide on one RNG: either route
everything through `rand_r(&thread_seed)` seeded from `seed=`, or `srand(seed)`
at startup. Thread-safety matters under OpenMP — prefer `rand_r`/per-thread
seeds derived deterministically from `seed=` and the thread id.

**Impact:** `seed=` is advertised as the way to get reproducible runs, but it
does not fully work. This blocks reproducibility for users and makes
paired/deterministic testing of search changes impossible (had to fall back to
statistical comparison over many runs).

---

## 3. Search-time scoring is slow on large trees (incremental cache is dead)

**Discovered:** 2026-07-03 (building a large test dataset for the ILS evaluation)

**Symptom:** A single heuristic-search replicate on a 40-taxon dataset (30
source trees) takes ~30 s; a 26-taxon dataset ~2–6 s. The cost scales steeply
with taxon count, making large real datasets slow.

**Cause / background:** `evaluate_candidate`/`probe_candidate` were changed
(this session) to score every criterion with the full-recompute path
(`compare_trees_*(FALSE)`), because the cached path (`TRUE`) returned scores
inconsistent with the authoritative recompute. That fixed correctness but
removed all incremental-scoring benefit: every candidate re-scores **all**
source trees from scratch, even though a single SPR/TBR move only changes the
induced topology of the few source trees whose taxa span the moved subtree.

The cache was *meant* to exploit this: `compare_trees_ml`/`_rf`/etc. skip
re-scoring source tree `i` when the SPR move didn't touch its taxa, gated by
`presenceof_SPRtaxa[]`. But the string-native SPR/TBR driver (`spr_new3`/
`tbr_new2`) never sets `presenceof_SPRtaxa[]` (it rebuilds `tree_top` from
Newick per candidate), so the cache either mis-fired (wrong scores — the bug we
fixed) or, if made to fire, would be unsound.

**Proper fix (future):** make incremental scoring correct *and* fast — after
producing a candidate, mark which taxa moved (set `presenceof_SPRtaxa[]` from
the prune/graft in `spr_new3`/`tbr_new2`) so `compare_trees_*` can soundly reuse
cached per-source-tree scores for untouched trees, then re-enable the `TRUE`
path. Alternatively, keep a per-source-tree cache keyed on the induced
topology hash. Either restores the O(moved) instead of O(all) scaling.

**Where to look:** `compare_trees_ml`/`compare_trees_rf`/`compare_trees_sfit`/
`compare_trees_qfit` in `scoring.c` (the `spr && here1||here2` reuse branch and
`presenceof_SPRtaxa[]`); `spr_new3`/`tbr_new2` and `evaluate_candidate`/
`probe_candidate` in `treecompare2.c` (where the move is applied and scored).

**Impact:** Correctness is fine now, but large-tree searches are much slower
than they could be. This also limited how hard a dataset we could use when
evaluating ILS (very hard datasets became too slow to run many replicates).

---

## 4. Poor parallel scaling of the multi-replicate heuristic search

**Discovered:** 2026-07-03 (while comparing ILS vs multi-rep under parallelism)

**Symptom:** `hs nreps=16` on a 16-core machine gave `real 53.5s` vs
`user 227s` — an effective speed-up of only ~4.2x, not ~16x. Most of the
machine sits idle even though 16 replicates are dispatched to 16 threads.

**Likely cause:** load imbalance across replicates. Replicates that converge
quickly finish early and their thread goes idle, while a few long-running
replicates (bad starting trees, many swaps to converge) dominate the wall time.
With a static `#pragma omp for` schedule over `nreps` this leaves cores unused.

**Where to look:** the OpenMP parallel region over the replicate loop in
`heuristic_search()` (`treecompare2.c`) — the `#pragma omp` schedule clause and
how replicates are assigned to threads. Options: `schedule(dynamic,1)` so idle
threads pick up new replicates; or decouple replicate count from thread count
and keep a work queue of "starts to optimise" that any idle thread drains; or,
if `nreps < nthreads`, split idle threads onto independent searches. Also
consider that ILS replicates have very different lengths, worsening imbalance.

**Impact:** This is arguably a *bigger practical lever than any search-strategy
change*: recovering the lost ~12x of a 16-core machine would dwarf the marginal
differences between ILS and extra replicates. Worth fixing before investing
further in which search variant is nominally "best per core".

---
