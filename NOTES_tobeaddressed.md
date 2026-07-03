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
