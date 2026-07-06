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

## 5. Decompose paralog-representative selection is non-deterministic

**Discovered:** 2026-07-06 (while validating the linear-time reconciliation)

**Symptom:** Running `autodecompose=yes` on a multicopy dataset twice with the
**same binary and same input** produces different fragment files. The
difference is *which* paralog copy is kept as a species' representative in a
multicopy tree — e.g. `tree1_frag1` comes out as
`...5693.XP_817042.1...` on one run and `...5693.XP_821347.1...` (or
`5691.EAN80177`) on the next. It even shifts a few per-tree fragment *counts*
(observed 538 vs 536 total fragment tags between two runs of the committed
binary on the same 80-tree input). The decomposition *structure* (which
`[treeN_fragM]` tags exist) is stable modulo those few shifts; the leaf labels
within a fragment are what vary.

**How found:** Comparing linear-DP vs pre-linear decompose output looked like a
regression until two runs of the *same* pre-linear binary were diffed and
showed the identical variation — so it is pre-existing, not from the linear DP.

**Where to look:** the multicopy paralog-pruning / representative-choice logic
in the decompose path — `decompose_walk` and the "keep one copy per species"
selection in `reconcile.c` (the duplication-node collapse that picks a `rep`),
and `prune.c`. Suspect iteration/pointer-order dependence or an unseeded
`rand()` (see issue #2 — the global-`rand()` seeding problem could be the same
root cause). Determine whether representative choice should be made
deterministic (e.g. lexicographically smallest name, or longest sequence) so
runs are reproducible.

**Impact:** Decomposition (and any supertree built from it) is not reproducible
run-to-run. Blocks byte-for-byte regression testing of the decompose path and
surprises users who expect identical output from identical input.

---

## 6. Post-decompose SIGILL/segfault drawing very large trees

**Discovered:** 2026-07-06 (running the 121k-tree / 477-taxon set end-to-end)

**Symptom:** On the big gene-tree file (`/Users/chriscreevey/Downloads/2759_trees.ph`,
~296 MB, 121,501 trees) with `autodecompose=yes autoprunemono=yes; set
criterion=ml; hs`, the run reaches and passes the decomposition
("Single-copy filter: using 813 single-copy trees...") and then dies with
`illegal hardware instruction` (exit 132) or a segfault (exit 139), sometimes
exiting 1 with empty output — while **drawing a very large ASCII tree** in the
post-decompose output/supertree-search phase. Non-deterministic across runs.

**How found:** Chasing whether the linear DP caused it — the committed
pre-linear binary reaches the *same* post-decompose drawing phase and churns
there just as long (>7 min) before dying, with none of the linear-DP code
involved. So it is pre-existing and unrelated to the decompose change; it is
simply reached now because decomposition itself got fast.

**Ruled out:** out-of-range gene-tree tags (a bounds guard in the linear DP
never fired); stack overflow (`ulimit -s 65500` did not help). ASan is
unavailable for gcc on macOS, which hampered pinpointing.

**Where to look:** the ASCII tree-drawing / tree-reload path and its fixed-size
`char` buffers (`TREE_LENGTH`-class), the same family as the overflow fixes made
earlier this session (the `returntree`/`print_*_tree` buffers). A ~2 MB Newick
tree very likely overruns a `TREE_LENGTH` (or similar) buffer during
drawing/printing, tripping the fortified `__chk_fail_overflow` (SIGILL) or
corrupting the heap (SIGSEGV). Size the relevant buffers from the actual tree
string length (as done for the autoprune buffers) rather than a fixed constant.

**Impact:** The large-scale ML `hs` workflow the user wants cannot complete on
the full dataset, even though the memory and decompose-speed problems around it
are now solved.

---
