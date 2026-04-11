# Clann User Manual

**Version 5.0.0**
Copyright (C) 2003–2026 Chris Creevey <chris.creevey@gmail.com>
http://www.creeveylab.org

Distributed under the GNU General Public License v2 or later.
See COPYING.txt for full licence terms.

---

## Table of Contents

1. [Overview](#1-overview)
2. [Installation and Startup](#2-installation-and-startup)
3. [Input File Formats](#3-input-file-formats)
4. [Key Concepts](#4-key-concepts)
   - 4.1 [Optimality Criteria](#41-optimality-criteria)
   - 4.2 [Delimiter Mode and Species Names](#42-delimiter-mode-and-species-names)
   - 4.3 [Single-Copy and Multicopy Gene Families](#43-single-copy-and-multicopy-gene-families)
5. [Command Reference](#5-command-reference)
   - [execute (exe)](#execute-exe)
   - [set](#set)
   - [hs](#hs--hsearch)
   - [nj](#nj)
   - [alltrees](#alltrees)
   - [bootstrap](#bootstrap--boot)
   - [usertrees](#usertrees)
   - [consensus](#consensus)
   - [reconstruct](#reconstruct)
   - [savetrees](#savetrees)
   - [showtrees](#showtrees)
   - [excludetrees](#excludetrees)
   - [includetrees](#includetrees)
   - [deletetaxa](#deletetaxa)
   - [restoretaxa](#restoretaxa)
   - [randomisetrees](#randomisetrees)
   - [rfdists](#rfdists)
   - [generatetrees](#generatetrees)
   - [yaptp](#yaptp)
   - [Autoprunemono](#autoprunemono)
   - [prunemonophylies](#prunemonophylies)
   - [sprdists](#sprdists)
   - [log](#log)
   - [tips](#tips)
   - 5b. [Tree-space landscape analysis](#5b-tree-space-landscape-analysis)
6. [Worked Examples](#6-worked-examples)
7. [Output Files Reference](#7-output-files-reference)
8. [Tips and Troubleshooting](#8-tips-and-troubleshooting)

---

## 1. Overview

Clann constructs supertrees from collections of gene trees and provides tools to explore the underlying phylogenomic signal. It implements several supertree optimality criteria, bootstrap support estimation, consensus methods, and gene-tree reconciliation. It can also calculate Robinson-Foulds and SPR distances between trees, perform randomisation tests, and manage large collections of gene trees.

---

## 2. Installation and Startup

### Direct command-line usage (recommended for pipelines)

Clann can be called directly from the shell with a command, an input file, and options:

```
clann <command> <treefile> [key=value ...]
```

**Examples:**
```bash
clann hs trees.ph
clann hs trees.ph criterion=ml nreps=5 nthreads=4
clann hs trees.ph --criterion=ml --nreps=5       # GNU-style flags also work
clann alltrees trees.ph criterion=rf
clann usertrees source.ph candidates.ph criterion=ml tests=yes nboot=1000
clann consensus trees.ph
clann hs --help                                   # per-command help
clann --help                                      # general help
```

**Available commands:**

| Command | Description |
|---------|-------------|
| `hs` / `hsearch` | Heuristic supertree search |
| `alltrees` | Exhaustive search (small datasets) |
| `usertrees` | Score / test user-supplied topologies |
| `consensus` | Consensus tree from source trees |
| `nj` | Neighbour-joining supertree |

**Global options** (applied before the command):

| Option | Description | Default |
|--------|-------------|---------|
| `criterion=<c>` | Scoring criterion: `dfit`, `ml`, `rf`, `sfit`, `qfit`, `avcon` | `dfit` |
| `nthreads=<n>` | Threads for parallel search | `1` |
| `mlbeta=<f>` | β parameter for ML criterion | `1.0` |
| `mlscale=<s>` | ML score display: `lnl`, `paper`, `lust` | `lnl` |
| `seed=<n>` | Random seed | (random) |

Command-specific options (e.g. `nreps`, `swap`, `tests`, `nboot`) are passed directly to the command — the same options available in interactive mode.

For `usertrees`, the first file is the source trees and the second file is the candidate topologies.

---

### Legacy command-line syntax

The original flag-based syntax is still fully supported:

```
clann [-lnh] [-c commands_file] [tree_file]
```

| Flag | Description |
|------|-------------|
| `-l` | Turn on logging; screen output is written to `clann.log` |
| `-n` | Non-interactive (batch) mode — requires a commands file (`-c`) or a Nexus clann-block in the tree file |
| `-c <file>` | Execute commands from `<file>` (one command per line; `#` starts a comment) |
| `-h` | Print usage and exit |
| `[tree_file]` | Optional Nexus or Phylip tree file to load on startup |

### Interactive mode

Without any arguments, Clann presents a `clann>` prompt. Type commands interactively. Append `?` to any command to see its options:

```
clann> hs ?
clann> reconstruct ?
```

### Batch / script mode

For complex workflows involving multiple commands, use a commands file:

```
clann -n -c commands.txt trees.ph
```

Or embed commands in a Nexus clann-block:

```
#NEXUS
begin clann;
  exe trees.ph;
  set criterion=dfit;
  hs nreps=10;
end;
```

---

## 3. Input File Formats

Clann auto-detects the file format when loading trees.

### Phylip / Newick format

Plain-text file, one tree per line, each ending with `;`:

```
((Human,Mouse),(Apple,Lemon),(Horse,Donkey));
(Human,(Mouse,(Apple,(Lemon,Horse))));
```

Optional features on each line:
- **Branch lengths**: `(Human:0.12,Mouse:0.08)` — preserved and used by average-consensus methods
- **Tree weight** (appended in brackets): `(Human,(Mouse,Apple));[2.5]`
- **Tree name** (appended in brackets): `(Human,(Mouse,Apple));[RAxML_tree1]`

Tree weights default to `1.0` if not supplied. Where weights are present they scale each gene tree's contribution to the total score for all criteria. For `criterion=ml` this has a direct probabilistic interpretation: the weight acts as a **per-tree β**, so a gene tree with weight 2.0 has twice the likelihood sharpness of a tree with weight 1.0 — it penalises RF disagreement more strongly. This makes tree weights a natural way to encode confidence in individual gene trees (e.g. derived from bootstrap support or sequence length).

### Nexus format

Detected by `#NEXUS` or `#` as the first non-whitespace character. Trees are read from a standard `trees` block. Tree weights stored in PAUP\* `[&W value]` annotations (including fractional form `[&W 1/2]`) are read automatically.

### Taxon name conventions

By default, Clann operates in **delimiter mode** and treats `'.'` as a separator between the species name and a gene-copy suffix. For example, `Human.1` and `Human.2` are both recognised as species `Human`. See [Section 4.2](#42-delimiter-mode-and-species-names) for details.

---

## 4. Key Concepts

### 4.1 Optimality Criteria

Set with `set criterion=<value>` before running any search command.

| Keyword | Full name | Description |
|---------|-----------|-------------|
| `dfit` | Distance Fit (DFIT) | Minimises the sum of squared differences between supertree and source-tree path-length distances. **Default.** |
| `sfit` | Splits Fit (SFIT) | Maximises the number of source-tree splits present in the supertree. |
| `qfit` | Quartet Fit (QFIT) | Maximises the number of source-tree quartet topologies consistent with the supertree. |
| `mrp` | Matrix Representation Parsimony (MRP) | Parsimony analysis of a matrix encoding of source-tree splits. Requires PAUP\*. |
| `avcon` | Average Consensus (AVCON) | Average-consensus distance matrix approach. Requires PAUP\*. |
| `recon` | Reconstruction / DL (RECON) | Minimises the total weighted duplication and loss cost of reconciling all source trees against the supertree. Uses **all** source trees including multicopy families. |
| `rf` | Robinson-Foulds (RF) | Minimises the normalised Robinson-Foulds distance summed across all source trees. RF distance counts the number of bipartitions present in one tree but not the other, normalised to [0, 1] per gene tree. |
| `ml` | Maximum Likelihood (ML) | Maximum-likelihood supertree criterion based on the exponential model of Steel & Rodrigo (2008). The probability of observing each gene tree given the supertree is modelled as P(G_i \| T) ∝ e^(−β·d_i), where d_i is the RF distance. The score reported is the total log-likelihood lnL = −β·Σd_i (negative, higher is better). Controlled by `mlbeta` and `mlscale`. If per-tree weights are supplied in the input file, each weight acts as a per-tree β multiplier (see [Section 3](#3-input-file-formats)), providing a principled way to encode confidence in individual gene trees. |

> **Note:** `mrp` and `avcon` require an external installation of [PAUP\*](https://paup.phylosolutions.com/).

> **Note:** `rf` and `ml` are newer criteria that are still under active testing. Use `dfit` for production analyses.

---

### 4.2 Delimiter Mode and Species Names

Clann's **delimiter mode** is active by default. It extracts the species name from gene-copy names by truncating at the first occurrence of the delimiter character (default `'.'`):

| Gene-copy name | Extracted species name |
|----------------|----------------------|
| `Human.1` | `Human` |
| `Mouse.2a` | `Mouse` |
| `E_coli` | `E_coli` (no dot → name unchanged) |
| `A.3.1` | `A` (truncates at first dot) |

This means a file containing both `(A,(B,(C,D)));` and `(A.1,B.1,(C.1,D.1));` loads as 4 unique species: A, B, C, D.

**To disable delimiter mode** and use full raw taxon names:
```
exe mydata.ph maxnamelen=full
```

**To change the delimiter character:**
```
exe mydata.ph delimiter_char=_
```

Delimiter mode is safe for datasets without dots in names — taxa without the delimiter character use their full name unchanged.

---

### 4.3 Single-Copy and Multicopy Gene Families

When delimiter mode is active, Clann automatically classifies each source tree:

- **Single-copy**: every species appears exactly once — suitable for all supertree methods
- **Multicopy**: at least one species appears more than once — indicates gene duplication

This affects how supertree search commands behave:

| Command | Behaviour |
|---------|-----------|
| `hs` (criterion ≠ recon) | Scores only single-copy trees; multicopy trees reserved for `reconstruct` |
| `hs criterion=recon` | Scores **all** trees; multicopy trees contribute via DL reconciliation |
| `nj` | Always uses only single-copy trees for distance matrix construction |
| `alltrees` | Uses only single-copy trees for scoring |
| `reconstruct` | Always uses **all** trees (single-copy and multicopy) |

A message is printed at the start of any filtered search:
```
Single-copy filter: using 12 single-copy trees for supertree search (4 multicopy trees reserved for 'reconstruct').
```

All trees are fully restored after the search — other commands such as `showtrees`, `savetrees`, and `reconstruct` always see the complete set.

---

## 5. Command Reference

---

### execute (exe)

Load source trees from a file. Must be run before any analysis command.

```
exe <filename> [options]
```

| Option | Values | Default | Description |
|--------|--------|---------|-------------|
| `maxnamelen` | `<integer>`, `delimited`, `full` | `delimited` | Maximum characters per taxon name. `delimited` extracts names before the delimiter character. `full` uses complete names as-is. An integer caps name length at N characters. |
| `delimiter_char` | `<character>` | `.` | Character used to split gene-copy names into species names. Only active when delimiter mode is on. |
| `summary` | `short`, `long` | `long` | Controls how much detail is printed about loaded trees. |
| `autoprunemono` | `yes`, `no` | `no` | At load time, prune monophyletic same-species clades from multicopy trees. Trees that become single-copy after pruning are promoted into the supertree search pool. Trees that remain multicopy after pruning (genuine deep paralogs) stay in the multicopy pool for `reconstruct`. Original unpruned trees are preserved for `reconstruct`. See [Autoprunemono](#autoprunemono) below. |

**Examples:**
```
exe trees.ph
exe trees.ph maxnamelen=full
exe trees.ph delimiter_char=_ summary=short
exe mydata.ph autoprunemono=yes
```

---

### set

Set global parameters that affect subsequent analyses.

```
set <parameter>=<value>
```

| Parameter | Values | Description |
|-----------|--------|-------------|
| `criterion` | `dfit`, `sfit`, `qfit`, `mrp`, `avcon`, `recon`, `rf`, `ml` | Optimality criterion for supertree reconstruction. See [Section 4.1](#41-optimality-criteria). |
| `seed` | `<integer>` | Random number seed for reproducibility. Default is based on system time and process ID. |
| `mlbeta` | `<float > 0>` | Global slope parameter β for the ML exponential model (Steel & Rodrigo 2008). Controls how steeply the likelihood decays with increasing RF distance — larger β penalises disagreement more strongly and implies higher confidence in the gene trees. Default: `1.0`. If per-tree weights are supplied in the input file, the effective per-tree slope is β × w_i, so weights are a natural way to encode differential confidence across gene trees. Only relevant when `criterion=ml`. |
| `mlscale` | `paper`, `lust`, `lnl` | Scoring convention for the ML criterion. `lnl` (default) reports lnL = −β·Σd_i, matching the sign convention of standard ML tools (negative, higher is better). `paper` reports β·Σd_i directly as in Steel & Rodrigo (2008) (positive, lower is better). `lust` applies an additional log₁₀(e) factor to match the original L.U.st Python tool of Akanni *et al.* (2014) exactly. |

**Examples:**
```
set criterion=dfit
set criterion=rf
set criterion=ml
set mlbeta=2.0
set mlscale=lnl
set seed=12345
```

---

### hs / hsearch

Heuristic search for the best-scoring supertree using SPR (subtree pruning and regrafting) branch swapping. The available options depend on the current criterion. When compiled with OpenMP, multiple independent replicates run in parallel.

```
hs [options]
```

#### Options for dfit, sfit, qfit, rf, ml (criteria 0, 2, 3, 6, 7)

| Option | Values | Default | Description |
|--------|--------|---------|-------------|
| `sample` | `<integer>` | 10,000 | Number of random starting trees evaluated before selecting the best `nreps` as starting points for heuristic search. Larger values find better starting trees at the cost of time. |
| `nreps` | `<integer>` | 10 | Number of independent heuristic search replicates. |
| `swap` | `nni`, `spr`, `tbr` | `spr` | Tree rearrangement algorithm: NNI (nearest-neighbour interchange), SPR (subtree pruning and regrafting), or TBR (tree bisection and reconnection). SPR and TBR generally explore more of tree space than NNI. |
| `nsteps` | `<integer>` | 5 | Number of improvement steps per replicate before stopping. |
| `start` | `nj`, `random`, `<filename>` | `nj` | Starting tree. `nj` uses a neighbour-joining tree; `random` uses the best tree from random sampling; a filename loads a user-provided starting tree. |
| `maxswaps` | `<integer>` | 1,000,000 | Maximum number of branch swaps per replicate. |
| `savetrees` | `<filename>` | `Heuristic_result.txt` | Output file for the best supertree(s) found. |
| `nthreads` | `<integer>` | all CPUs | Number of OpenMP threads. Each thread runs an independent search replicate in parallel. Not available for `criterion=recon`. |
| `maxskips` | `<integer>` | auto (2N²) | Stop a replicate after this many consecutive already-visited SPR moves. Set to `0` to disable. Default auto-scales to 2×(number of taxa)². Prevents long runs that have converged. |
| `progress` | `<integer>` | 5 | *(OpenMP only)* How often (in seconds) to print a best-so-far status line when running multiple threads. Set to `0` to print every time a new global best is found. Has no effect when running single-threaded. |
| `droprep` | `<float>` | 0 (disabled) | *(OpenMP only)* Abandon a replicate early if its current best score is more than this fraction above the current global best across all threads (e.g. `0.1` = 10%). The per-rep completion line will show `droprep` as the stop reason. The freed thread immediately starts the next queued replicate. Set to `0` to disable. Has no effect when running single-threaded. |
| `visitedtrees` | `<filename>` | *(disabled)* | Record every unique topology visited during the search to a tab-separated file (columns: `newick`, `score`, `visit_count`). Accumulated across all replicates and threads. See [Tree-space landscape analysis](#tree-space-landscape-analysis) for post-processing. |
| `autoprunemono` | — | — | Set via the `exe` command (`exe myfile.ph autoprunemono=yes`), not `hs` directly. Prunes monophyletic same-species clades from multicopy trees at load time so more source trees contribute to the search. See [Autoprunemono](#autoprunemono). |

**Weight options (criterion-specific):**

| Criterion | Option | Values | Description |
|-----------|--------|--------|-------------|
| `dfit` | `weight` | `equal`, `comparisons` | Whether source trees are weighted equally or by the number of taxon-pair comparisons they contribute. |
| `sfit` | `weight` | `equal`, `splits` | Equal weighting or weighting by number of splits. |
| `qfit` | `weight` | `equal`, `taxa`, `quartets` | Equal, taxa-weighted, or quartet-count-weighted. |

**Histogram options (dfit only):**

| Option | Values | Default | Description |
|--------|--------|---------|-------------|
| `drawhistogram` | `yes`, `no` | `no` | Generate a histogram of scores across the random sample. |
| `nbins` | `<integer>` | 20 | Number of bins for the histogram. |
| `histogramfile` | `<filename>` | `Heuristic_histogram.txt` | Output file for histogram data. |

#### Options for ml criterion (criterion 7) only

| Option | Values | Default | Description |
|--------|--------|---------|-------------|
| `mlbeta` | `<float > 0>` | 1.0 | Slope parameter β for the exponential likelihood model P(G_i\|T) ∝ e^(−β·d_i). Larger values make the likelihood sharper and penalise RF distance more strongly. Can also be set via `set mlbeta=<value>`. |
| `mlscale` | `paper`, `lust`, `lnl` | `lnl` | Scoring convention. `lnl` reports lnL = −β·Σd_i (negative, standard ML convention). `paper` reports the Steel & Rodrigo (2008) positive score directly. `lust` uses Akanni *et al.* (2014) log₁₀ scaling. See `set mlscale` for details. |
| `mleta` | `<float ≥ 0>` | `0.0` | **[Experimental]** Tree-size scaling exponent η. Each source tree's RF distance is divided by k_i^η before scoring, where k_i is the number of internal splits in that tree. η = 0 recovers the original Steel & Rodrigo (2008) model; η = 1 fully normalises by split count; η > 1 actively down-weights large trees beyond normalisation. Can be estimated from data using `mlscores eta=auto`. See note below. |

> **mleta — experimental tree-size scaling exponent.** The Steel model implicitly treats each split disagreement as an independent event with equal cost, meaning source trees with more splits contribute proportionally more to WD and therefore dominate β̂ and the search objective. The `mleta` parameter addresses this by scaling each tree's RF distance by k_i^η before scoring, derived from the probability model P(d_i | S, β, η) = (β/k_i^η) · exp(−(β/k_i^η) · d_i). This gives log L(β, η) = n·log(β) − η·Σ log(k_i) − β·Σ w_i·d_i/k_i^η. The term −η·Σ log(k_i) is a natural penalty from the model's normalisation constant that prevents η from growing without bound, unlike ad hoc weighting schemes. η = 0 is the original Steel (2008) model; η = 1 gives equal per-tree contribution regardless of size; η > 1 actively down-weights large trees (appropriate if large trees are less reliable due to ILS or estimation error). The optimal η can be estimated jointly with β using `mlscores eta=auto`, which performs a 1-D grid search since β profiles analytically at each fixed η. If the optimal supertree topology changes substantially between η = 0 and the estimated η*, tree-size heterogeneity is a confounding factor in your dataset.
>
> **Note on naming:** Clann's η is a *tree-size scaling exponent* and should not be confused with the normalising constant α used in Steel & Rodrigo (2008) or the L.U.St paper (Akanni *et al.* 2014). Those papers use α to denote the partition function Z_T — a completely different quantity. Clann uses η to avoid this naming clash.

#### Options for recon criterion (criterion 5) only

| Option | Values | Default | Description |
|--------|--------|---------|-------------|
| `duplications` | `<float>` | 1.0 | Cost weight applied to each duplication event. |
| `losses` | `<float>` | 1.0 | Cost weight applied to each gene loss event. |
| `numspeciesrootings` | `<integer>`, `all` | 2 | Number of species-tree rootings to trial when computing the reconciliation score. `all` tries every possible rooting (slow but exact). |
| `numgenerootings` | `<integer>`, `all` | 2 | Number of gene-tree rootings to trial per species-tree rooting. |

#### Options for mrp criterion (criterion 1)

| Option | Values | Default | Description |
|--------|--------|---------|-------------|
| `analysis` | `parsimony`, `nj` | `parsimony` | Whether to use parsimony or NJ for the MRP matrix. |
| `weighted` | `yes`, `no` | `no` | Weight characters in the MRP matrix. |
| `swap` | `nni`, `spr`, `tbr` | `tbr` | Branch-swap algorithm. |
| `addseq` | `simple`, `closest`, `asis`, `random`, `furthest` | `random` | Taxon-addition sequence for parsimony. |
| `nreps` | `<integer>` | 10 | Number of parsimony replicates. |
| `savetrees` | `<filename>` | `MRP.tree` | Output file. |

**Examples:**
```
hs
hs nreps=20 sample=50000 swap=spr
hs nreps=5 savetrees=my_supertree.ph
set criterion=rf
hs nreps=10
set criterion=ml
set mlbeta=1.0
hs nreps=10 nthreads=8
set criterion=recon
hs nreps=10 duplications=1.0 losses=0.5 numspeciesrootings=5
hs nreps=20 nthreads=8 visitedtrees=landscape.tsv
```

---

### nj

Build a single supertree by the neighbour-joining method. Pairwise taxon distances are derived from the source trees using the average-consensus approach, with missing pairwise distances estimated from the `missing` option.

Multicopy gene trees are always excluded from the distance matrix calculation when delimiter mode is active.

```
nj [options]
```

| Option | Values | Default | Description |
|--------|--------|---------|-------------|
| `missing` | `4point`, `ultrametric` | `4point` | Method for estimating pairwise distances for taxa pairs not co-occurring in any source tree. `4point` uses the four-point condition; `ultrametric` assumes an ultrametric tree. |
| `savetrees` | `<filename>` | `NJ-tree.ph` | Output file for the NJ supertree. |

**Examples:**
```
nj
nj missing=ultrametric savetrees=my_nj.ph
```

---

### alltrees

Exhaustively evaluate every possible supertree topology within a specified range. Practical only for small numbers of taxa (≤ ~8). Multicopy gene trees are excluded from scoring when delimiter mode is active. Supported for criteria: `dfit`, `sfit`, `qfit`, `rf`, and `ml`.

```
alltrees [options]
```

| Option | Values | Default | Description |
|--------|--------|---------|-------------|
| `range` | `<start> - <end>` | all | Restrict evaluation to supertree topologies in this index range. |
| `savetrees` | `<filename>` | `top_alltrees.txt` | Output file for the best-scoring tree(s). |
| `create` | `yes`, `no` | `no` | Write all evaluated trees to `alltrees.ph`. |
| `nbest` | `<integer>` | — | Keep only the top N scoring trees. |
| `keep` | `<float>` | — | Retain trees scoring at least this value. |

**Weight options** are the same as for `hs` (criterion-dependent, see above).

> `alltrees` with `criterion=rf` or `criterion=ml` is the recommended way to verify that `hs` has found the global optimum on small datasets, since it scores every possible topology exactly.

**Examples:**
```
alltrees
alltrees savetrees=best_topology.ph
set criterion=rf
alltrees
set criterion=ml
alltrees
```

---

### bootstrap / boot

Perform a bootstrap supertree analysis. Source trees are resampled with replacement for each replicate, a heuristic search is run on the bootstrap replicate, and the bootstrap support values are mapped onto a consensus tree.

```
bootstrap [options]
```

| Option | Values | Default | Description |
|--------|--------|---------|-------------|
| `nreps` | `<integer>` | 100 | Number of bootstrap replicates. |
| `hsreps` | `<integer>` | 10 | Number of heuristic search replicates per bootstrap replicate. |
| `sample` | `<integer>` | 10,000 | Random starting trees per heuristic replicate. |
| `swap` | `nni`, `spr`, `tbr`, `all` | `spr` | Branch-swap algorithm. |
| `start` | `random`, `<filename>` | `random` | Starting tree for each heuristic search. |
| `nsteps` | `<integer>` | 5 | Steps per replicate. |
| `maxswaps` | `<integer>` | 1,000,000 | Maximum swaps per replicate. |
| `nthreads` | `<integer>` | all CPUs | Number of OpenMP threads. Each thread runs an independent bootstrap replicate in parallel. |
| `progress` | `<integer>` | 5 | *(OpenMP only)* How often (in seconds) to print a best-so-far status line when running multiple threads. Set to `0` to print every time a new global best is found. |
| `droprep` | `<float>` | 0 (disabled) | *(OpenMP only)* Abandon a bootstrap heuristic search replicate early if its score is more than this fraction above the current global best. See `hs droprep` for details. |
| `treefile` | `<filename>` | `bootstrap.txt` | Output file for all bootstrap trees. |
| `consensus` | `strict`, `majrule`, `minor`, `<float 0–1>` | `majrule` | Consensus method: strict (1.0), majority-rule (0.5), or a custom threshold (e.g. `0.75` for 75%). |
| `consensusfile` | `<filename>` | `consensus.ph` | Output file for the consensus tree with bootstrap support. |
| `autoprunemono` | — | — | Set via the `exe` command (`exe myfile.ph autoprunemono=yes`), not `boot` directly. Prunes monophyletic same-species clades from multicopy trees at load time. See [Autoprunemono](#autoprunemono). |

**Weight and criterion-specific options** are the same as for `hs`.

**Examples:**
```
bootstrap
bootstrap nreps=500 hsreps=5 consensus=majrule
bootstrap nreps=100 treefile=bs_trees.ph consensusfile=bs_consensus.ph
```

---

### usertrees

Score one or more user-provided supertree topologies against the source trees in memory using the current criterion. Useful for testing specific hypotheses, scoring reference topologies, or statistical comparison of candidate trees under `criterion=ml`.

```
usertrees <filename> [options]
```

| Option | Values | Default | Description |
|--------|--------|---------|-------------|
| `outfile` | `<filename>` | `Usertrees_result.txt` | Output file for scores. |
| `printsourcescores` | `yes`, `no` | `no` | Also print the score of each individual source tree against the best user tree. |
| `tests` | `yes`, `no` | `no` | Run ML topology tests after scoring. Only valid with `criterion=ml`. |
| `nboot` | `<integer>` | `1000` | Bootstrap replicates for the SH test. |
| `testsfile` | `<filename>` | `mltest_results.txt` | File for per-gene-tree δ breakdown table. |
| `normcorrect` | *(flag)* | off | Apply the Bryant & Steel (2008) normalising constant correction to all reported lnL values. Only active with `criterion=ml` and `mlscale=lnl`. See below. |

**Weight options** are criterion-dependent (same as `hs`).

#### ML topology tests (`tests=yes`)

When `criterion=ml` is active and `tests=yes` is specified, CLANN runs three standard topology tests comparing each candidate tree to the best-scoring tree (T1):

- **Winning Sites (Steel & Rodrigo 2008):** Counts gene trees that favour T1 (n+) versus the candidate (n−). Reports an exact two-sided binomial p-value.
- **Kishino–Hasegawa (KH) test (Kishino & Hasegawa 1989):** Parametric test. Computes the per-gene-tree log-likelihood difference δ_i, estimates its standard error under the null, and reports a one-sided normal p-value.
- **Shimodaira–Hasegawa (SH) test (Shimodaira & Hasegawa 1999):** Bootstrap version of the KH test that corrects for the selection of the best tree. Resamples gene trees with replacement (`nboot` replicates) after centring the differences.

The AU (approximately unbiased) test is not yet implemented.

Each gene tree contributes one independent observation (analogous to a site in sequence-based tests). Only gene trees that are informative for at least one topology (score > 0) are counted. A warning is issued when fewer than 4 informative gene trees are available.

**Typical workflow:**
```
exe examples/tutorial_trees.ph
set criterion=ml
alltrees savetrees=all.ph create=yes   # or produce candidate topologies any other way
usertrees all.ph tests=yes
usertrees all.ph tests=yes nboot=5000 testsfile=detailed_results.txt
```

**Output example:**
```
ML topology tests  (T1 = tree 1,  lnL = -43.0000,  beta = 1.00)
=============================================================================
  Tree      lnL         WinSites p    KH z / p         SH p
-----------------------------------------------------------------------------
  T2        lnL=-44.000  p=0.2380      z=1.42 p=0.078 ns  p=0.1430 ns
  T3        lnL=-47.000  p=0.0312      z=2.30 p=0.011 *   p=0.0280 *
-----------------------------------------------------------------------------
  * p<alpha=0.05; ** p<0.01; *** p<0.001; ns=not significant
```

A tab-separated table with full per-gene-tree details is written to `testsfile`.

**References:**
- Steel, M. & Rodrigo, A. (2008) *Maximum likelihood supertrees.* Syst. Biol. 57:243–250.
- Kishino, H. & Hasegawa, M. (1989) *Evaluation of the maximum likelihood estimate of the evolutionary tree topologies from DNA sequence data.* J. Mol. Evol. 29:170–179.
- Shimodaira, H. & Hasegawa, M. (1999) *Multiple comparisons of log-likelihoods with applications to phylogenetic inference.* Mol. Biol. Evol. 16:1114–1116.

#### Normalising constant correction (`normcorrect`)

The Steel & Rodrigo (2008) model requires a normalising constant Z_T so that probabilities sum to 1 over all possible source trees:

    Z_T = Σ_{T'} exp(−β · d(T', T|X_i))

The original Steel & Rodrigo (2008) paper omitted this constant when deriving the ML condition (noted by Bryant & Steel 2008), arguing it was approximately equal across candidate supertrees. Bryant & Steel (2008) showed that Z_T does vary with tree shape — specifically with the number of *cherries* (pairs of adjacent leaves) in T|X_i — and derived an efficient O(n⁵) algorithm to compute it exactly.

The corrected per-source-tree lnL contribution is:

    lnL_i = w_i · [log(β) + log(w_i) − β · d_i − log(Z_{T|X_i})]

Without `normcorrect`, `log(Z_{T|X_i})` is always zero (ignored). With `normcorrect`, Clann subtracts this term using a truncated large-β expansion:

    log Z_T ≈ log(1 + b₂ε + b₄(c_T)ε²)

where ε = e^{−2β}, b₂ = 2(n−3) is the number of trees at RF distance 2 (same for all binary trees), and b₄(c_T) = 4·C(n−3,2) + 6·(n−6+c_T) depends on the cherry count c_T of T|X_i. This approximation is accurate for β > ~1.5.

**When does it matter?** Bryant & Steel (2008) showed that ignoring Z_T changes the *ranking* of candidate supertrees only when |log Z_{T1} − log Z_{T2}| ≥ 2β. The topology-dependent part of log Z is driven entirely by the cherry count and enters only at second order (ε²), so:

| Dataset | Effect on rankings | Effect on absolute lnL |
|---|---|---|
| Source trees n ≤ 20, any β | Negligible | Minor |
| Source trees n = 33, β ≈ 1.9 | Negligible (Δ ≈ 0.05 ≪ 2β = 3.8) | ~0.83 nats per tree |
| Source trees n ≥ 50, β ≈ 1.5 | Potentially significant | Large |

The correction is most valuable when comparing absolute lnL values between datasets or between different methods. It does **not** affect KH/SH test statistics (differences in lnL cancel the topology-independent part of log Z), but makes the absolute lnL more interpretable and model-correct.

```
usertrees candidates.ph normcorrect
usertrees candidates.ph tests=yes normcorrect
```

**Reference:**
- Bryant, D. & Steel, M. (2008) *Computing the distribution of a tree metric.* arXiv:0810.0868.

**Examples:**
```
usertrees mytopology.ph
usertrees candidates.ph outfile=scores.txt printsourcescores=yes
usertrees candidates.ph tests=yes nboot=1000
usertrees candidates.ph tests=yes normcorrect          # correct absolute lnL
```

---

### consensus

Calculate a consensus tree from the source trees currently in memory (or from supertrees/bootstrap trees if specified).

```
consensus [options]
```

| Option | Values | Default | Description |
|--------|--------|---------|-------------|
| `method` | `strict`, `majrule`, `minor`, `<float 0–1>` | `minor` | Consensus type. `strict` retains only splits present in all trees; `majrule` retains splits in >50% of trees; `minor` retains any split; a value between 0 and 1 sets a custom threshold. |
| `filename` | `<filename>` | `consensus.ph` | Output file for the consensus tree. |
| `data` | `source`, `supertrees`, `bootstraps` | `source` | Which set of trees to build consensus from. |
| `guidetree` | `<filename>` | — | Optional file containing a guide tree to constrain the consensus. |

**Examples:**
```
consensus
consensus method=majrule filename=majority_rule.ph
consensus method=0.75
consensus data=bootstraps filename=boot_consensus.ph
```

---

### reconstruct

Carry out gene-tree reconciliation by mapping each source gene tree onto a species tree and counting the minimum number of duplications and losses required. Can be used to identify paralogs, orthologs, and to study gene-family evolution.

`reconstruct` always uses **all** source trees in memory, including multicopy gene families.

```
reconstruct [options]
```

| Option | Values | Default | Description |
|--------|--------|---------|-------------|
| `speciestree` | `memory`, `first`, `<filename>` | `memory` | Source of the species tree. `memory` uses the supertree most recently produced by `hs`, `nj`, or `alltrees`. `first` uses the first source tree (legacy behaviour). A filename reads the species tree from a file. |
| `dups` / `duplications` | `<float>` | 1.0 | Cost applied to each duplication event. |
| `losses` | `<float>` | 1.0 | Cost applied to each gene loss event. |
| `basescore` | `<float>` | 1.0 | Score used when resolving polytomies in gene trees. |
| `showrecon` | `yes`, `no` | `no` | Print detailed reconciliation information for each node. |
| `printfiles` | `yes`, `no` | `yes` | Generate detailed output files (see [Section 7](#7-output-files-reference)). |
| nhxfile | `<filename>` | none | Prints nhx-formatted file of resulting reconstructions for all source trees. This allows viewing of the reconstructions (inclusding indications of duplications, losses, etc) with compatible tools. |

**Typical workflow:**

```
exe mydata.ph
hs nreps=10
reconstruct speciestree=memory nhxfile=my_reconstructions.nhx
```

**Using an external species tree:**

```
reconstruct speciestree=my_species_tree.ph dups=2.0 losses=1.0
```

This command can also be used to view the gene-tree reconstructions from a heuristic search of treespace using the *recon* criterion:



```
exe mydata.ph
set criterion=recon
hs nreps=1
reconstruct speciestree=memory nhxfile=my_reconstructions.nhx
```



---

### savetrees

Save a subset of source trees currently in memory to a file in Phylip format.

```
savetrees [filter options]
```

| Option | Values | Default | Description |
|--------|--------|---------|-------------|
| `filename` | `<filename>` | `savedtrees.txt` | Output file. |
| `range` | `<start> - <end>` | all | Save only trees with indices in this range (1-based). |
| `size` | `equalto <N>`, `lessthan <N>`, `greaterthan <N>` | — | Filter by number of taxa. |
| `namecontains` | `<string>` | — | Save only trees whose name contains this string. |
| `containstaxa` | `<string>` | — | Save only trees that contain a taxon with this name. |
| `score` | `<min> - <max>` | — | Save only trees whose score falls in this range (requires a previous scoring step). |

**Examples:**
```
savetrees filename=subset.ph range=1-50
savetrees size=greaterthan 10 filename=large_trees.ph
savetrees containstaxa=Human
```

---

### showtrees

Display source trees in ASCII format on screen and/or save them to a file.

```
showtrees [filter options]
```

| Option | Values | Default | Description |
|--------|--------|---------|-------------|
| `display` | `yes`, `no` | `yes` | Whether to print trees to the screen. |
| `savetrees` | `yes`, `no` | `no` | Whether to save the displayed trees to a file. |
| `filename` | `<filename>` | `showtrees.txt` | Output file (when `savetrees=yes`). |

Filter options are the same as `savetrees` (`range`, `size`, `namecontains`, `containstaxa`, `score`).

**Examples:**
```
showtrees
showtrees range=1-5
showtrees namecontains=RAxML savetrees=yes filename=raxml_trees.txt
```

---

### excludetrees

Flag source trees so they are excluded from all subsequent analyses (`hs`, `boot`, `nj`, `rfdists`, etc.). The trees remain in memory and can be restored with `includetrees`. Both `excludetrees` and `includetrees` accept the same filter criteria.

```
excludetrees [filter options]
```

| Option | Values | Description |
|--------|--------|-------------|
| `singlecopy` | (flag, no value) | Exclude all single-copy trees. |
| `multicopy` | (flag, no value) | Exclude all multicopy trees. |
| `range` | `<start> <end>` | Exclude trees in the given index range (1-based). |
| `size` | `equalto \| lessthan \| greaterthan <n>` | Exclude trees by taxon count. |
| `namecontains` | `<string>` | Exclude trees whose name contains the given string. |
| `containstaxa` | `<taxon>` | Exclude trees containing the specified taxon. |
| `score` | `<min> <max>` | Exclude trees whose score falls in the given range. |

**Examples:**
```
excludetrees range 1 10
excludetrees size lessthan 5
excludetrees multicopy
```

---

### includetrees

Restore previously excluded source trees so they are included in subsequent analyses. Accepts the same filter criteria as `excludetrees`; with no options, restores all excluded trees.

```
includetrees [filter options]
```

| Option | Values | Description |
|--------|--------|-------------|
| `range` | `<start> <end>` | Restore trees in the given index range (1-based). |
| `size` | `equalto \| lessthan \| greaterthan <n>` | Restore trees by taxon count. |
| `namecontains` | `<string>` | Restore trees whose name contains the given string. |
| `containstaxa` | `<taxon>` | Restore trees containing the specified taxon. |
| `score` | `<min> <max>` | Restore trees whose score falls in the given range. |

**Examples:**
```
includetrees
includetrees range 1 10
includetrees namecontains RAxML
```

---

### deletetaxa

Remove specified taxa from all source trees in memory, pruning branches while preserving the remaining topology and branch lengths. Trees that fall below the minimum taxon count after pruning are removed entirely.

A full snapshot of the original trees is saved automatically before the operation. Use `restoretaxa` immediately afterwards to undo it. Only one snapshot is held at a time — a second `deletetaxa` overwrites the previous snapshot.

```
deletetaxa <taxon1> [<taxon2> ...] [options]
```

| Option | Values | Default | Description |
|--------|--------|---------|-------------|
| `mintaxa` | `<integer>` | 4 | Trees with fewer than this many taxa after pruning are deleted. |

**Examples:**
```
deletetaxa Human
deletetaxa Human Mouse mintaxa=5
```

---

### restoretaxa

Restore the full set of source trees that existed before the most recent `deletetaxa` operation. The snapshot includes all trees, their names, and their weights. Calling `restoretaxa` clears the snapshot — a subsequent `restoretaxa` without a preceding `deletetaxa` will report an error.

```
restoretaxa
```

No options.

**Example:**
```
deletetaxa Outgroup
hs
restoretaxa
```

---

### randomisetrees

Randomise all source trees in memory while preserving the taxon composition of each tree. The resulting trees have the same set of taxa per tree but random topologies. Useful as null-model datasets.

```
randomisetrees
```

No options.

---

### rfdists

Calculate all pairwise Robinson-Foulds (RF) topological distances between the source trees in memory and write them to a file.

```
rfdists [options]
```

| Option | Values | Default | Description |
|--------|--------|---------|-------------|
| `filename` | `<filename>` | `robinson_foulds.txt` | Output file. |
| `output` | `matrix`, `vector` | `matrix` | Format of the output: full distance matrix or a flat list. |
| `missing` | `none`, `4point`, `ultrametric` | `none` | How to handle taxon pairs absent from some trees. |

**Examples:**
```
rfdists
rfdists filename=RF_distances.txt output=matrix
```

---

### generatetrees

Generate random supertrees, score them against the source trees, and produce a histogram of the score distribution. Used to assess the significance of supertree scores relative to a null distribution.

```
generatetrees [options]
```

| Option | Values | Default | Description |
|--------|--------|---------|-------------|
| `method` | `equiprobable`, `markovian` | `equiprobable` | Randomisation method. `equiprobable` samples topologies with equal probability; `markovian` uses a Markov chain process. |
| `ntrees` | `<integer>`, `all` | 100 | Number of random trees to generate and score. |
| `nbins` | `<integer>` | 20 | Number of bins in the output histogram. |
| `outfile` | `<filename>` | `histogram.txt` | Output file for histogram data. |
| `sourcedata` | `real`, `randomised`, `ideal` | `real` | Data to score against: the real source trees, randomised source trees, or an ideal dataset derived from a specified supertree. |
| `supertree` | `memory`, `<filename>` | `memory` | Supertree to use when `sourcedata=ideal`. |
| `savescores` | `yes`, `no` | `no` | Save individual scores to a separate file. |
| `savesourcetrees` | `yes`, `no` | `no` | Save the randomised source trees. |

**Examples:**
```
generatetrees ntrees=1000 outfile=null_distribution.txt
generatetrees sourcedata=ideal supertree=memory ntrees=500
```

---

### yaptp

**Yet Another Permutation Tail Probability** test. Assesses whether the source trees contain more phylogenetic signal than expected by chance. Repeatedly randomises the source trees, performs a heuristic search on each randomisation, and compares the resulting score distribution to the real-data score.

```
yaptp [options]
```

| Option | Values | Default | Description |
|--------|--------|---------|-------------|
| `method` | `equiprobable`, `markovian` | `equiprobable` | Randomisation method. |
| `nreps` | `<integer>` | 100 | Number of randomisation replicates. |
| `hsreps` | `<integer>` | 10 | Heuristic search replicates per randomisation. |
| `sample` | `<integer>` | 10,000 | Random starting trees per heuristic search. |
| `search` | `nni`, `spr`, `all` | `spr` | Branch-swap method. |
| `nsteps` | `<integer>` | 5 | Steps per heuristic replicate. |
| `maxswaps` | `<integer>` | 1,000,000 | Maximum swaps per heuristic search. |
| `treefile` | `<filename>` | `yaptp.ph` | Output file. |

**Weight and criterion-specific options** (recon, dfit, sfit, qfit) are the same as for `hs`.

**Examples:**
```
yaptp nreps=100 hsreps=5
yaptp method=markovian nreps=200 treefile=yaptp_result.ph
```

---

### Autoprunemono

The `autoprunemono=yes` option on the `exe` command provides an automated, in-memory version of `prunemonophylies`. It prunes monophyletic same-species clades at load time so that more source trees can contribute to the supertree search.

**Behaviour:**
- For each multicopy tree, monophyletic same-species clades are pruned to a single representative (chosen at random).
- If the tree becomes single-copy after pruning, it is **promoted** to the supertree search pool (`hs`, `nj`, `alltrees`).
- If the tree is still multicopy after pruning (genuine deep paralogs with non-monophyletic copies), it remains in the multicopy pool and is used only by `reconstruct`.
- The **original unpruned trees are preserved** internally for `reconstruct`, which always sees the full data.

**Example output:**
```
Autoprunemono: pruned monophyletic same-species clades in 47 multicopy trees.
  Promoted to single-copy pool:   38 trees
  Still multicopy after pruning:  9 trees (retained for reconstruct)
  Original (unpruned) trees stored for reconstruct.
```

**Workflow:**
```
clann> exe mydata.ph autoprunemono=yes
clann> hs nreps=10              # uses single-copy trees + promoted trees
clann> reconstruct speciestree memory   # uses all original trees
```

> **Note:** The original unpruned trees are only preserved within the current session. If you save and reload the pruned trees (e.g. using `prunemonophylies`), `reconstruct` will see only the pruned versions.

---

### prunemonophylies

For each source tree, identify clades that consist entirely of gene copies from the same species (monophyletic inparalogs) and prune them down to a single representative. Useful for pre-processing multicopy gene trees before supertree analysis with methods that require single-copy input.

> **Tip:** For an automated in-session version of this workflow, use `exe mydata.ph autoprunemono=yes` — it prunes at load time, promotes qualifying trees to the supertree pool, and preserves originals for `reconstruct`. See [Autoprunemono](#autoprunemono).

```
prunemonophylies [options]
```

| Option | Values | Default | Description |
|--------|--------|---------|-------------|
| `filename` | `<filename>` | `prunedtrees.txt` | Output file for the pruned trees. |
| `selection` | `random`, `length` | `random` | How to choose the retained representative from a monophyletic group. `random` picks one at random. `length` picks the longest sequence, using the value encoded in the gene-copy name after the species and an additional delimiter (e.g. `Species.length.XXXXX`). |

**Examples:**
```
prunemonophylies
prunemonophylies selection=length filename=pruned.ph
```

---

### sprdists

Estimate the SPR (Subtree Pruning and Regrafting) distance between: the real source trees and a supertree, an ideal dataset derived from the supertree, and a randomised dataset. Used to assess how well the supertree explains the source data.

```
sprdists [options]
```

| Option | Values | Default | Description |
|--------|--------|---------|-------------|
| `supertree` | `create`, `memory`, `<filename>` | `create` | Supertree to use. `create` generates one; `memory` uses the in-memory supertree; a filename reads from file. |
| `dorandomisation` | `yes`, `no` | `yes` | Whether to also compute SPR distances to a randomised dataset. |
| `outfile` | `<filename>` | `SPRdistances.txt` | Output file. |

**Examples:**
```
sprdists
sprdists supertree=memory outfile=spr_results.txt
sprdists dorandomisation=no
```

---

### log

Control logging of all screen output to a file.

```
log status=<on|off> [file=<filename>]
```

| Option | Values | Default | Description |
|--------|--------|---------|-------------|
| `status` | `on`, `off` | `off` | Enable or disable logging. |
| `file` | `<filename>` | `clann.log` | Logfile name. Can only be set when turning logging on. |

**Examples:**
```
log status=on
log status=on file=my_analysis.log
log status=off
```

---

### tips

Display usage hints.

```
tips [number=<1-10>]
```

Without options, shows a random tip. Specify `number=N` to display a specific tip.

---

## 5b. Tree-space landscape analysis

The `visitedtrees=<filename>` option on `hs` writes a tab-separated file
recording every unique topology encountered during the search (across all
replicates and threads), with three columns:

| Column | Description |
|--------|-------------|
| `newick` | Named-taxon unrooted Newick string |
| `score` | Criterion score on first visit (lower = better for dfit/RF; more negative = better for ML) |
| `visit_count` | Total times this topology was proposed across all replicates and threads |

The intended use is post-hoc exploration of the score landscape in tree space:
compute pairwise SPR distances between the visited topologies, embed with MDS,
and visualise to identify local and global optima and understand convergence.

For the full workflow — SPR distance computation, Python visualisation scripts
(MDS scatter plot, score histogram, neighbour graph, 3D surface), identifying
local optima, and troubleshooting — see
[NOTES_treespace_landscape.md](NOTES_treespace_landscape.md).

---

## 6. Worked Examples

### Example 1: Basic supertree from single-copy genes

```
exe my_gene_trees.ph
hs nreps=10 sample=10000 swap=spr
```

### Example 2: Full workflow with mixed single-copy and multicopy dataset

```
exe my_dataset.ph
# Delimiter mode is on by default — species names extracted from gene-copy names
# hs will automatically use only single-copy trees for scoring
hs nreps=10 sample=10000
# Reconstruct uses ALL trees, including multicopy families
reconstruct speciestree=memory dups=1.0 losses=1.0
```

### Example 3: Supertree with DL reconciliation criterion (uses all trees)

```
exe my_dataset.ph
set criterion=recon
hs nreps=10 sample=5000 numspeciesrootings=5 numgenetries=5
reconstruct speciestree=memory
```

### Example 4: Bootstrap analysis

```
exe trees.ph
set criterion=dfit
hs nreps=10
bootstrap nreps=100 hsreps=5 consensus=majrule consensusfile=bs_support.ph
```

### Example 5: Testing significance of phylogenetic signal

```
exe trees.ph
hs nreps=10
yaptp nreps=100 hsreps=5
```

### Example 6: Batch mode with a commands file

`commands.txt`:
```
set criterion=dfit
hs nreps=20 sample=50000 swap=spr savetrees=result.ph
bootstrap nreps=200 hsreps=10 consensusfile=bootstrap_consensus.ph
reconstruct speciestree=memory
```

Run:
```
clann -n -c commands.txt trees.ph
```

### Example 7: Tree-space landscape recording

```
exe my_gene_trees.ph
set criterion=dfit
hs nreps=20 nthreads=8 visitedtrees=landscape.tsv
```

Post-process with Python (see [Tree-space landscape analysis](#5b-tree-space-landscape-analysis)):
```bash
tail -n +2 landscape.tsv | cut -f1 > visited.ph
rspr visited.ph > spr_matrix.txt
python3 landscape_plot.py
```

### Example 8: Working with datasets that have dots in species names

```
# Disable delimiter mode if your species names contain dots (e.g. E.coli, B.subtilis)
exe my_prokaryote_trees.ph maxnamelen=full
hs nreps=10
```

---

## 7. Output Files Reference

| File | Produced by | Contents |
|------|-------------|----------|
| `Heuristic_result.txt` | `hs` | Best supertree(s) in Phylip format with scores |
| `bootstrap.txt` | `bootstrap` | All bootstrap supertrees |
| `consensus.ph` | `bootstrap`, `consensus` | Consensus tree with support values |
| `NJ-tree.ph` | `nj` | Neighbour-joining supertree |
| `top_alltrees.txt` | `alltrees` | Best-scoring supertree(s) from exhaustive search |
| `Usertrees_result.txt` | `usertrees` | Scores for each user-provided topology |
| `mltest_results.txt` | `usertrees tests=yes` | Per-gene-tree δ table for ML topology tests (criterion=ml only) |
| `<filename>` | `hs visitedtrees=<filename>` | Tab-separated landscape file: all unique topologies visited during the search with their scores and visit counts (see [Tree-space landscape analysis](#5b-tree-space-landscape-analysis)) |
| `robinson_foulds.txt` | `rfdists` | Pairwise RF distances between source trees |
| `SPRdistances.txt` | `sprdists` | SPR distance analysis results |
| `prunedtrees.txt` | `prunemonophylies` | Pruned source trees |
| `clann.log` | `log` | Complete session log |
| `supertree.ps` | `hs`, `alltrees` | PostScript rendering of the supertree |
| `<input>.recon` | `reconstruct` | Reconciliation mapping for each gene tree |
| `<input>.recon.dist` | `reconstruct` | Distribution of duplications and losses per gene tree |
| `<input>.onetoone` | `reconstruct` | Putative one-to-one orthologs |
| `<input>.strict.onetoone` | `reconstruct` | Strict one-to-one orthologs |
| `<input>.gene.births` | `reconstruct` | Gene-tree birth-point assignments |
| `<input>.species.decendents` | `reconstruct` | Species-tree descendant information |

---

## 8. Tips and Troubleshooting

**Q: My `hs` runs but finds a very poor supertree.**
Increase `nreps` and `sample`. Also try `swap=tbr` for more thorough tree-space exploration. Check that you have loaded enough informative source trees.

**Q: Taxon names appear truncated or merged unexpectedly.**
Check whether your taxon names contain dots. If they do and you do not want delimiter splitting, add `maxnamelen=full` to your `exe` command.

**Q: I have multicopy gene families but only single-copy trees are used in `hs`.**
This is expected behaviour when delimiter mode is active and criterion is not `recon`. To use all trees (including multicopy), set `criterion=recon`:
```
set criterion=recon
hs nreps=10
```

**Q: Clann says "number of single copy trees: 0".**
All trees are multicopy. Options:
- Use `exe mydata.ph autoprunemono=yes` to automatically prune monophyletic same-species clades at load time — trees that become single-copy after pruning are promoted into the supertree search pool.
- Set `criterion=recon` to use all multicopy trees directly via DL reconciliation.
- Run `prunemonophylies`, save the output, reload the pruned file, then run `hs`.

**Q: `bootstrap` or `hs` with `mrp` or `avcon` fails with an error about PAUP\*.**
These criteria require PAUP\* to be installed and accessible on your PATH. Install PAUP\* and ensure the `paup` executable is findable, or switch to a different criterion (e.g., `dfit`).

**Q: How do I reproduce a previous analysis exactly?**
Use `set seed=<N>` with a fixed integer before running any search command. Record the seed value in your log file (`log status=on`).

**Q: `reconstruct` gives different scores each run.**
By default, only a small number of gene-tree and species-tree rootings are sampled. Increase `numspeciesrootings` and `numgenerootings` (or set to `all`) for deterministic results:
```
reconstruct speciestree=memory numspeciesrootings=all numgenerootings=all
```

**Q: How do I run Clann non-interactively on a cluster or in a pipeline?**
The simplest approach is the direct CLI:
```bash
clann hs trees.ph criterion=ml nthreads=8 > clann_output.txt 2>&1
```
For complex multi-step workflows, use a commands file:
```bash
clann -n -c commands.txt trees.ph > clann_output.txt 2>&1
```
Both approaches exit automatically when done and return a shell exit code.

**Q: What does "Investigating phylogenetic information through supertree analyses" actually mean?**
See the primary Clann publications (Creevey & McInerney 2005; Creevey et al. 2004) cited in the README.

---

*For further information, consult the [Clann publication](https://academic.oup.com/bioinformatics/article/21/3/390/238167) or contact chris.creevey@gmail.com.*
