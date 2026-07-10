# Heuristic Search Options in Clann

A guide to the options that control the `hs` (heuristic supertree) search, what
each one does in plain terms, what we learned about whether it helps, and which
optimality criteria each applies to.

All options are set on the `hs` command line, e.g.

```
exe gene_trees.ph
set criterion=ml
hs nreps=10 vns=yes ils=10
```

The banner printed at the top of each `hs` run echoes the active settings, so you
can always confirm what you are running.

---

## 1. The problem these options address

A heuristic tree search climbs from a starting tree by making small
rearrangements (moves) that improve the score, until no move helps — a **local
optimum**. The trouble is that the score landscape is **rugged**: there are many
local optima, and the score is a coarse (Robinson-Foulds–based) measure, so a
tree only a few branches away from the true optimum can score badly. A single
climb therefore often stops short of the **global optimum** (the best possible
supertree).

Because the true tree is unknown in a real analysis, **the only practical
evidence that you have found the global optimum is that independent searches keep
converging to the same tree.** So the options below are ultimately about one
thing: making each search reliable enough that replicates agree.

---

## 2. Quick reference

| Option | Values | Default | One-line purpose |
|---|---|---|---|
| `vns` | `yes` / `no` | **yes** | Variable Neighborhood Descent: escalate NNI→SPR→TBR |
| `ils` | integer (0 = off) | `10` | Iterated Local Search: perturb + re-optimise to escape optima |
| `ilsstrength` | integer ≥ 1 | `3` | Size of each ILS perturbation (grows adaptively when stuck) |
| `ilsguided` | `yes` / `no` | `no` | Aim ILS kicks at conflicted taxa (opt-in; no measured benefit) |
| `plateau` | `yes` / `no` | `yes` | Allow equal-score sideways moves to cross flat regions |
| `maxplateau` | integer (0 = off) | auto `2×N` | Cap on consecutive sideways moves |
| `maxexhaustive` | number (0 = never) | `1 000 000` | Auto-switch tiny trees to a guaranteed exhaustive search |
| `smoothsearch` | `yes` / `no` | `no` | ML only: smoother search objective (opt-in; no measured benefit) |
| `maxskips` | integer (0 = off) | auto `N³` | Stop a replicate after this many repeat-topology visits |
| `strategy` | `first` / `best` | `first` | Take the first improving move, or survey and take the best |

`N` = number of taxa in the search.

---

## 3. What each option does — and what we learned

### `vns` — Variable Neighborhood Descent  *(the main win; default on)*
Instead of using one kind of move, VNS uses a hierarchy: it searches with cheap
**NNI** moves to convergence; if stuck, it escalates to **SPR**; if still stuck,
to **TBR** (the largest neighborhood); and whenever a bigger move finds an
improvement it drops back to NNI. The result is a tree that is optimal under *all
three* move sets.

**Why it helps:** many local optima that trap SPR can be escaped by a single TBR
move — TBR can cross rearrangements the smaller moves cannot. This is the one
change we tested that gave a large, statistically significant improvement.

**What we measured** (single search reaching the known optimum on concordant
test data): 20 taxa 49% → 83%; 30 taxa 36% → 89%; and the worst-case result
improved dramatically (e.g. −126 → −12). Its cost is extra runtime (~2–3× a
single-operator climb), which is why it is the default only because reliability
matters more than speed — it is never *worse* in quality (its neighborhood
includes NNI/SPR/TBR), only slower. Turn it off with `vns no` for speed.

### `ils` — Iterated Local Search  *(default on, 10 kicks)*
When a climb converges, ILS **perturbs** the best tree with a few random SPR
moves and re-optimises, keeping the better of the two. It repeats until `ils`
consecutive perturbations yield no improvement. Set `ils 0` to disable.

**Why it helps:** the perturbation jumps to a nearby but different region, so the
next climb can find a different — possibly deeper — optimum. It reliably removes
the very worst outcomes (it shrinks the "deep trap" tail).

**What we learned:** ILS clearly beats a plain single climb, but on *easy* data,
at a matched compute budget, simply running more independent replicates reaches
the optimum slightly more often (cheap restarts saturate faster). Its real value
is on hard problems and in combination with VNS, where together they give the
best reliability.

### `ilsstrength` — perturbation size  *(default 3, adaptive)*
The number of random SPR moves in each ILS kick. It is **adaptive**: it grows by
one each time a kick fails to improve (up to 4× the base) and resets on any
improvement — small nudges when progress is steady, bigger jumps to break out of
a persistent trap.

### `ilsguided` — conflict-directed kicks  *(opt-in; default off)*
Instead of perturbing random branches, aim the kick at the taxa that sit in
source trees the current supertree fails to display (the taxa "in conflict").

**What we learned:** despite being well-motivated, this did **not** beat random
kicks (65% vs 62%, not significant). On a deceptive landscape the locally
conflicted taxa do not reliably point toward the global optimum, and the
re-optimisation after the kick does the real work regardless of where it lands.
Kept as opt-in; the underlying per-taxon conflict signal is also useful on its
own as a **rogue-taxon diagnostic** (which taxa are hardest to place).

### `plateau` / `maxplateau` — cross flat regions  *(default on, budget 2×N)*
The RF-based scores are full of large **plateaus** (many different trees with the
exact same score). A strict "only move if strictly better" climber stops at the
edge of a plateau even when crossing it would open a path to a better tree.
`plateau yes` lets the search take equal-score sideways moves; `maxplateau` caps
how many in a row so it cannot wander a plateau forever (a per-search visited-set
prevents it revisiting the same trees).

**What we learned:** clearly helps the **dfit** criterion; roughly neutral for
ML. Low-risk, so left on by default.

### `maxexhaustive` — guaranteed answer for tiny trees  *(default 1,000,000)*
If the total number of possible supertrees is at or below this limit, `hs`
automatically runs an **exhaustive** search (scoring every tree) instead of the
heuristic — faster on tiny problems and *guaranteed* to find the global optimum.
Set `maxexhaustive 0` to always use the heuristic. Applies to the tree space
*after* Clann collapses paralogues, so e.g. a 9-effective-taxon problem
(≈135,135 trees) auto-routes here.

### `smoothsearch` — smoother ML objective  *(ML only; opt-in; default off)*
The RF score gives a branch all-or-nothing credit (present or absent), which is
why the landscape is flat and rugged. `smoothsearch yes` scores candidates during
the ML search with a **transfer-distance** term (partial credit for "almost
right" branches) as a tie-breaker on RF plateaus, then reports the true ML score.

**What we learned:** did **not** significantly improve reaching the optimum
(52% vs 47%, not significant) and costs ~2× the scoring time. An earlier variant
that *replaced* RF with transfer distance was actively worse (it rewarded
near-misses). Kept as opt-in for other datasets/criteria.

### `maxskips` — convergence stopping  *(auto N³)*
A replicate stops once it has re-visited already-seen topologies `maxskips` times
in a row (a sign it has converged). Auto-scaled to N³ so bigger trees get more
budget. `maxskips 0` disables it.

### `strategy` — first vs best improvement  *(default first)*
`first` commits to the first improving move found; `best` surveys the whole
neighborhood and takes the single best move. `best` is more thorough per step but
slower; `first` was as reliable in our tests.

---

## 4. A correctness fix underneath all of this

Independently of the options above, we fixed two scoring bugs that were the real
cause of the original symptom ("even on data with no conflict, the ML search
sometimes can't reach the zero-conflict optimum"):

1. **Search-time scoring consistency.** The SPR/TBR search was scoring
   sfit/qfit/rf/ml candidates with a stale cached shortcut that disagreed with
   the authoritative scoring used by `usertrees`/`alltrees`. It now recomputes
   fully during search (as `dfit` always did), so every code path agrees.
2. **Trivial-bipartition over-count.** A source tree stored with redundant outer
   parentheses (a spurious degree-2 root) was crediting an extra, meaningless
   split, making an otherwise-displayed tree look like a conflict. Bipartition
   counting now skips these, so a truly compatible dataset scores 0 everywhere.

These apply to all the bipartition/distance criteria (dfit, sfit, qfit, rf, ml).

---

## 5. Which criteria can use each option

Criteria: **dfit** (Most Similar Supertree, `criterion=dfit`), **sfit** (maximum
split fit), **qfit** (maximum quartet fit), **rf** (Robinson-Foulds), **ml**
(maximum likelihood supertree), **recon** (duplication/loss reconciliation).
MRP (`criterion=mrp`) and AVCON (average consensus) do **not** use Clann's
SPR/TBR search engine at all (MRP hands off to PAUP*; AVCON builds a
neighbour-joining tree from average-consensus distances), so none of these search
options apply to them.

| Option | dfit | sfit | qfit | rf | ml | recon | MRP / AVCON |
|---|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| `vns` (VND local search) | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✗ |
| `ils` + `ilsstrength` (adaptive) | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✗ |
| `plateau` / `maxplateau` | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✗ |
| `maxskips` / `strategy` | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✗ |
| `maxexhaustive` (exhaustive routing) | ✅ | ✅ | ✅ | ✅ | ✅ | ✗ | ✗ |
| `ilsguided` (conflict-directed kicks) | ✅ | ✅ | ✅ | ✅ | ✅ | ✗ | ✗ |
| `smoothsearch` (transfer surrogate) | ✗ | ✗ | ✗ | ✗ | ✅ | ✗ | ✗ |
| scoring-consistency + bipartition fix | ✅ | ✅ | ✅ | ✅ | ✅ | — | — |

Notes on the exceptions:
- **recon** is excluded from `maxexhaustive` (the exhaustive scorer doesn't
  handle reconciliation) and from `ilsguided` (its per-source residual doesn't
  map to per-taxon conflict the same way), but it *does* use VNS, ILS, plateau.
- **smoothsearch** is ML-only today; the transfer-distance idea would extend
  naturally to **rf**, less naturally to sfit/qfit — not yet done, because it has
  not shown a clear win on any criterion.

---

## 6. What we learned overall

- **The big win is VNS** (escalating neighborhood size). It is the only
  search-strategy change that produced a large, significant, repeatable
  improvement, and it directly improves the thing that matters when the true tree
  is unknown: how often independent searches converge to the same tree.
- **Random restarts are hard to beat per unit compute.** Directed kicks
  (`ilsguided`) and a smoother objective (`smoothsearch`) were well-motivated but
  did not help — on a deceptive landscape, local information doesn't reliably
  point to the global optimum. Neighborhood *size* (VNS) was the lever that
  worked, not kick *direction* or objective *smoothing*.
- **Correctness first.** Part of the original problem was a scoring
  inconsistency, not the search. With that fixed, a compatible dataset now scores
  0 consistently across every code path.

### Practical recommendation
For a careful analysis where you want confidence in the result, use the defaults
(`vns yes`, `ils 10`) and **run several replicates** (`nreps`), then check how
many independent replicates ended at the same best score — that agreement is your
evidence of having found the global optimum. Use `vns no` only when you need
speed and are willing to trade reliability.
