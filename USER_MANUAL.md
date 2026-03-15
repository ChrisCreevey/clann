# Clann User Manual

**Version 4.3.0**
Copyright Chris Creevey 2003–2026
http://www.creeveylab.org

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
   - [deletetrees](#deletetrees)
   - [deletetaxa](#deletetaxa)
   - [randomisetrees](#randomisetrees)
   - [rfdists](#rfdists)
   - [generatetrees](#generatetrees)
   - [yaptp](#yaptp)
   - [prunemonophylies](#prunemonophylies)
   - [sprdists](#sprdists)
   - [log](#log)
   - [tips](#tips)
6. [Worked Examples](#6-worked-examples)
7. [Output Files Reference](#7-output-files-reference)
8. [Tips and Troubleshooting](#8-tips-and-troubleshooting)

---

## 1. Overview

Clann constructs supertrees from collections of gene trees and provides tools to explore the underlying phylogenomic signal. It implements several supertree optimality criteria, bootstrap support estimation, consensus methods, and gene-tree reconciliation. It can also calculate Robinson-Foulds and SPR distances between trees, perform randomisation tests, and manage large collections of gene trees.

---

## 2. Installation and Startup

### Command-line syntax

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

Without `-n`, Clann presents a `clann>` prompt. Type commands interactively. Append `?` to any command to see its options:

```
clann> hs ?
clann> reconstruct ?
```

### Batch mode

Create a plain-text commands file (one command per line, or separated by `;`) and run:

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

### Nexus format

Detected by `#NEXUS` or `#` as the first non-whitespace character. Trees are read from a standard `trees` block.

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
| `ml` | Maximum Likelihood (ML) | Maximum-likelihood supertree criterion based on the L.U.st exponential model (Akanni *et al.* 2013). The probability of observing each gene tree given the supertree is modelled as P(G_i \| T) ∝ e^(−β·d_i), where d_i is the RF distance. The score reported is the total log-likelihood lnL = −β·Σd_i (negative, higher is better). Controlled by `mlbeta` and `mlscale`. |

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

**Examples:**
```
exe trees.ph
exe trees.ph maxnamelen=full
exe trees.ph delimiter_char=_ summary=short
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
| `mlbeta` | `<float > 0>` | Slope parameter β for the ML exponential model. Controls how steeply the likelihood decays with increasing RF distance. Default: `1.0`. Only relevant when `criterion=ml`. |
| `mlscale` | `paper`, `lust`, `lnl` | Scoring convention for the ML criterion. `lnl` (default) reports lnL = −β·Σd_i, matching the sign convention of standard ML tools (negative, higher is better). `paper` reports β·Σd_i directly (positive, lower is better). `lust` applies an additional log₁₀(e) factor to match the original L.U.st Python tool output exactly. |

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
| `maxskips` | `<integer>` | 1,000 | Stop a replicate after this many consecutive already-visited SPR moves. Set to `0` to disable. Prevents long runs that have converged. |

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
| `mlscale` | `paper`, `lust`, `lnl` | `lnl` | Scoring convention. `lnl` reports lnL = −β·Σd_i (negative, standard ML convention). `paper` and `lust` report positive minimisation scores. See `set mlscale` for details. Can also be set via `set mlscale=<value>`. |

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
| `treefile` | `<filename>` | `bootstrap.txt` | Output file for all bootstrap trees. |
| `consensus` | `strict`, `majrule`, `minor`, `<float 0–1>` | `majrule` | Consensus method: strict (1.0), majority-rule (0.5), or a custom threshold (e.g. `0.75` for 75%). |
| `consensusfile` | `<filename>` | `consensus.ph` | Output file for the consensus tree with bootstrap support. |

**Weight and criterion-specific options** are the same as for `hs`.

**Examples:**
```
bootstrap
bootstrap nreps=500 hsreps=5 consensus=majrule
bootstrap nreps=100 treefile=bs_trees.ph consensusfile=bs_consensus.ph
```

---

### usertrees

Score one or more user-provided supertree topologies against the source trees in memory using the current criterion. Useful for testing specific hypotheses.

```
usertrees <filename> [options]
```

| Option | Values | Default | Description |
|--------|--------|---------|-------------|
| `outfile` | `<filename>` | `Usertrees_result.txt` | Output file for scores. |
| `printsourcescores` | `yes`, `no` | `no` | Also print the score of each individual source tree against the best user tree. |

**Weight options** are criterion-dependent (same as `hs`).

**Examples:**
```
usertrees mytopology.ph
usertrees candidates.ph outfile=scores.txt printsourcescores=yes
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

**Typical workflow:**
```
exe mydata.ph
hs nreps=10
reconstruct speciestree=memory
```

**Using an external species tree:**
```
reconstruct speciestree=my_species_tree.ph dups=2.0 losses=1.0
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

### deletetrees

Permanently remove source trees from memory based on filter criteria. This operation cannot be undone within a session.

```
deletetrees [filter options]
```

| Option | Values | Description |
|--------|--------|-------------|
| `singlecopy` | (flag, no value) | Delete all single-copy trees. |
| `multicopy` | (flag, no value) | Delete all multicopy trees. |

Filter options are the same as `savetrees` (`range`, `size`, `namecontains`, `containstaxa`, `score`).

**Examples:**
```
deletetrees range=1-10
deletetrees size=lessthan 5
deletetrees multicopy
```

---

### deletetaxa

Remove specified taxa from all source trees in memory, pruning branches while preserving the remaining topology and branch lengths. Trees that fall below the minimum taxon count after pruning are removed entirely.

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

### prunemonophylies

For each source tree, identify clades that consist entirely of gene copies from the same species (monophyletic inparalogs) and prune them down to a single representative. Useful for pre-processing multicopy gene trees before supertree analysis with methods that require single-copy input.

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

### Example 7: Working with datasets that have dots in species names

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
All trees are multicopy. Either set `criterion=recon` to use them via DL reconciliation, or use `prunemonophylies` to reduce each multicopy tree to a single copy per species before analysis.

**Q: `bootstrap` or `hs` with `mrp` or `avcon` fails with an error about PAUP\*.**
These criteria require PAUP\* to be installed and accessible on your PATH. Install PAUP\* and ensure the `paup` executable is findable, or switch to a different criterion (e.g., `dfit`).

**Q: How do I reproduce a previous analysis exactly?**
Use `set seed=<N>` with a fixed integer before running any search command. Record the seed value in your log file (`log status=on`).

**Q: `reconstruct` gives different scores each run.**
By default, only a small number of gene-tree and species-tree rootings are sampled. Increase `numspeciesrootings` and `numgenerootings` (or set to `all`) for deterministic results:
```
reconstruct speciestree=memory numspeciesrootings=all numgenerootings=all
```

**Q: How do I run Clann non-interactively on a cluster?**
Use batch mode with a commands file:
```
clann -n -c commands.txt trees.ph > clann_output.txt 2>&1
```

**Q: What does "Investigating phylogenetic information through supertree analyses" actually mean?**
See the primary Clann publications (Creevey & McInerney 2005; Creevey et al. 2004) cited in the README.

---

*For further information, consult the [Clann publication](https://academic.oup.com/bioinformatics/article/21/3/390/238167) or contact chris.creevey@gmail.com.*
