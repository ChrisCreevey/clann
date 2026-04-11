
	*********************************************************************
	*                                                                   *
	*                             Clann                                 *
	* Investigating phylogenetic information through supertree analyses *
	*                                                                   *
	*                  lab: http://www.creeveylab.org                   *
	*                  email: chris.creevey@gmail.com                   *
	*                                                                   *
	*           Copyright (C) Chris Creevey 2003-2026                   *
	*     SPDX-License-Identifier: GPL-2.0-or-later                     *
	*  Free software under GNU GPL v2 or later — see COPYING.txt        *
	*                                                                   *
	*          HINT: Type "help" to see all available commands          *
	*********************************************************************

# Aim

To construct supertrees and explore the underlying phylogenomic information from gene trees. Clann implements several supertree optimality criteria — including DFIT (most-similar supertree), Robinson-Foulds (RF), and a maximum likelihood criterion based on Steel & Rodrigo (2008) — with bootstrapping support for all criteria. It can calculate consensus trees and provides tools to manage sets of gene trees: pruning taxa, filtering by a variety of criteria, automatically pruning lineage-specific duplicates, and more. A gene-tree reconciliation approach is also implemented to utilise information from multicopy gene families for species tree estimation. Clann can calculate RF and SPR distances between sets of trees, create randomised versions of trees, and perform a PTP test for informativeness (the YAPTP test). Clann has been continually developed since 2003 and suggestions for new features and tools are welcome.


# Referencing Clann

Clann has been published in Bioinformatics in 2005 under the following title:

Creevey C.J. and McInerney J.O. 2005 Clann: investigating phylogenetic information through supertree analyses. Bioinformatics 21(3): 390-2. [Link](https://academic.oup.com/bioinformatics/article/21/3/390/238167)

The Bootstrapping and YAPTP methods and the DFIT (most similar supertree algorithm) have all been described in the paper:

Creevey C.J., Fitzpatrick, D.A., Philip, G.A., Kinsella, R.J., O'Connell M.J., Travers, S.A, Wilkinson M. and McInerney J.O. 2004 Does a tree-like phylogeny only exist at the tips in the prokaryotes? Proceedings of the Royal Society London, B series: Biological Sciences 271(1557): 2551-8. [Link](http://rspb.royalsocietypublishing.org/content/271/1557/2551)

Either or both of these publications should be cited if you use Clann in published work.

For a list of papers that have cited Clann in their work see the [Clann google scholar page](https://scholar.google.co.uk/citations?view_op=view_citation&hl=en&citation_for_view=7JkjEd4AAAAJ:UeHWp8X0CEIC)

# What's new in v5

- **ML supertree criterion** — maximum likelihood criterion based on Steel & Rodrigo (2008), scoring candidate supertrees via the sum of exponential-decay likelihoods over RF distances to source trees. Supports Bryant & Steel (2008) normalising constant correction (`normcorrect`), including an exact subset-enumeration method (`normcorrect=exact`) for small trees.
- **RF criterion** — Robinson-Foulds bipartition distance as an optimality criterion for `hs`, `bootstrap`, `alltrees`, and `usertrees`.
- **ML topology tests** — Winning Sites, Kishino-Hasegawa (KH), and Shimodaira-Hasegawa (SH) tests for comparing candidate supertrees in `usertrees`.
- **Per-source-tree score matrix** — `usertrees scorematrix=yes` outputs a tab-delimited table of per-source-tree scores against each candidate topology, including tree size and weight columns, for identifying decisive source trees.
- **OpenMP parallelism** — `hs` and `bootstrap` run in parallel across CPU cores by default (`nthreads=` to control). Parallel replicates with periodic best-score reporting.
- **Improved heuristic search** — best-improvement SPR strategy, convergence detection via skip-streak, bipartition-hash visited-set to avoid redundant scoring, and `start=memory` to resume from a previous result.
- **Source tree weighting** — `autoweight=bootstrap` weights source trees by bipartition support; `autoweight=clan` weights by compatibility with user-defined clans; per-tree weights reported in `showtrees`.
- **Interactive CLI improvements** — tab completion for commands, options and values; coloured prompt; Ctrl+R reverse history search; persistent command history (`~/.clannhistory`).
- **Code refactoring** — split from a monolithic source into modular files (`scoring.c`, `tree_io.c`, `topology.c`, `tree_ops.c`, `consensus.c`, `viz.c`, `utils.c`).

# Documentation

A complete reference for every command and its options is provided in **[USER_MANUAL.md](USER_MANUAL.md)**. It covers:

- All command-line flags and startup options
- Input file formats (Phylip/Newick and Nexus)
- All optimality criteria and when to use each
- Delimiter mode and single-copy/multicopy gene family handling
- Full option tables for every command (`hs`, `nj`, `alltrees`, `bootstrap`, `reconstruct`, and more)
- Worked examples for common workflows
- Output file reference
- Troubleshooting

A tutorial walking through the new features in clann is provided in **[TUTORIAL.md](TUTORIAL.md)**.

# Installation

## Package managers (recommended)

**Conda (Bioconda):**

```bash
conda install -c bioconda -c conda-forge clann
```

Or with `mamba` for faster solves:

```bash
mamba install -c bioconda -c conda-forge clann
```

**Homebrew (macOS / Linux):**

```bash
brew tap ChrisCreevey/clann
brew install clann
```

Or, once accepted into [homebrew-bio](https://github.com/brewsci/homebrew-bio):

```bash
brew tap brewsci/bio
brew install clann
```

See **[PACKAGING.md](PACKAGING.md)** for detailed submission and update instructions for both channels.

---

## Building from source

### Prerequisites

| Dependency | Purpose | Install |
|---|---|---|
| GCC ≥ 4.9 | Compiler | `apt install build-essential` / `brew install gcc` |
| GNU readline | Interactive prompt (tab completion, history) | `apt install libreadline-dev` / `brew install readline` |
| GCC with OpenMP | Parallel heuristic search | included in Linux GCC; `brew install gcc` on macOS |

On **macOS**, Apple's Clang does not ship an OpenMP runtime. `configure` automatically detects Apple Clang and switches to a Homebrew GCC (`gcc-15` … `gcc-11`) if one is installed. If none is found, Clann builds single-threaded and prints:

```
To enable multi-threading: brew install gcc  (then re-run ./configure)
```

### Build steps

```bash
./configure
make
sudo make install    # optional — installs clann to /usr/local/bin
```

To keep the binary in the source directory instead, skip `make install` and run `./clann` directly.

See the **[INSTALL](INSTALL)** file for full prerequisites, macOS notes, and troubleshooting (permission errors, automake version mismatches, etc.).

---

# Usage

## Direct command-line (recommended)

Clann can be called directly from the shell — no interactive session required:

```bash
clann hs trees.ph
clann hs trees.ph criterion=ml nreps=5 nthreads=4
clann hs trees.ph --criterion=ml --nreps=5          # GNU-style flags also work
clann alltrees trees.ph criterion=rf
clann usertrees source.ph candidates.ph criterion=ml tests=yes nboot=1000
clann consensus trees.ph
clann --help                 # general help
clann hs --help              # per-command help
```

**Available commands:** `hs` (heuristic search), `alltrees` (exhaustive), `usertrees` (score/test candidate topologies), `consensus`, `nj` (neighbour-joining).

**Global options** (prefix with `--` or use native `key=value`):

| Option | Values | Default |
|--------|--------|---------|
| `criterion` | `dfit`, `ml`, `rf`, `sfit`, `qfit`, `avcon` | `dfit` |
| `nthreads` | integer | all CPUs |
| `mlbeta` | float | `1.0` |
| `mleta` | float | `0.0` |
| `mlscale` | `lnl`, `paper`, `lust` | `lnl` |
| `seed` | integer | (random) |

Command-specific options (`nreps`, `swap`, `savetrees`, `tests`, `nboot`, etc.) are identical to their interactive equivalents — see `clann <command> --help` or the [USER_MANUAL.md](USER_MANUAL.md).

This interface is designed for use in shell scripts, Snakemake rules, Nextflow processes, and cluster job submission. Clann exits automatically and returns a standard shell exit code.

---

## Interactive mode

Typing `clann` with no arguments starts the interactive REPL:

```
clann> exe trees.ph
clann> set criterion=ml
clann> hs nreps=10
clann> quit
```

The interactive prompt supports **tab completion** for commands and options, a **coloured prompt**, and **Ctrl+R** reverse history search. Command history is saved across sessions in `~/.clannhistory`.

Append `?` to any command to see its options: `hs ?`, `usertrees ?`, etc.

---

## Batch / script mode (legacy)

For complex multi-step workflows, use a commands file:

```bash
clann -n -c commands.txt trees.ph
```

Or embed commands in a Nexus clann-block inside the tree file:

```
#NEXUS
begin clann;
  exe trees.ph;
  set criterion=ml;
  hs nreps=10;
end;
```

Full flag reference: `clann -h`


## Working with Single-Copy and Multicopy Gene Families

Clann handles datasets that contain a mixture of single-copy and multicopy gene families automatically. This is particularly relevant for phylogenomics datasets where some gene families are present as a single copy per species (single-copy orthologs) and others have undergone duplications (multicopy gene families).

### Species name extraction (delimiter mode)

By default, Clann extracts species names from gene-copy names using `'.'` as a delimiter. For example, gene-copy names like `Human.1`, `Human.2` are both recognised as belonging to species `Human`. This means you can load a file containing both single-copy and multicopy gene trees without any additional setup:

```
clann> exe mydata.ph
Delimiter mode active: species names extracted before '.' (use maxnamelen=full to disable)
    Number of unique taxa: 4
```

If your taxon names contain dots but you do **not** want this splitting behaviour (e.g. names like `E.coli` should remain as-is), disable delimiter mode at load time:

```
clann> exe mydata.ph maxnamelen=full
Full taxon names in use (delimiter mode disabled).
    Number of unique taxa: 17
```

You can also change the delimiter character from `'.'` to another character using `delimiter_char`:

```
clann> exe mydata.ph delimiter_char=_
```

### Supertree searches with mixed datasets (hs, nj, alltrees)

When delimiter mode is active and your dataset contains multicopy gene trees, the supertree search commands (`hs`, `nj`, `alltrees`) automatically restrict scoring to **single-copy trees only**. Multicopy trees are reserved for use with `reconstruct`. An informational message is printed at the start of any search:

```
clann> hs
Single-copy filter: using 12 single-copy trees for supertree search (4 multicopy trees reserved for 'reconstruct').
```

After the search completes, all trees (including multicopy) are fully restored in memory for use by any other command. This filtering is completely transparent — `reconstruct`, `showtrees`, `savetrees` and all other commands continue to see the full set of source trees.

If you disable delimiter mode (`maxnamelen=full`), no filtering occurs and all trees contribute to supertree scoring.

### Automatic pruning of monophyletic clades (autoprunemono)

Many multicopy trees are "almost" single-copy: they contain lineage-specific duplications that form monophyletic same-species clades. The `autoprunemono=yes` option prunes these clades at load time so that qualifying trees can contribute to the supertree search:

```
clann> exe mydata.ph autoprunemono=yes

Autoprunemono: pruned monophyletic same-species clades in 47 multicopy trees.
  Promoted to single-copy pool:   38 trees
  Still multicopy after pruning:  9 trees (retained for reconstruct)
  Original (unpruned) trees stored for reconstruct.
```

- **Promoted trees** (single-copy after pruning): join `hs`, `nj`, `alltrees`
- **Still-multicopy trees** (genuine deep paralogs): remain in the multicopy pool for `reconstruct`
- `reconstruct` always uses the **original unpruned trees** — no information is discarded

### Gene-tree reconciliation with multicopy families (reconstruct)

The `reconstruct` command performs duplication/loss reconciliation of source trees against a species tree. Unlike the supertree search methods, `reconstruct` uses **all** source trees — both single-copy and multicopy — making it the appropriate tool for extracting information from multicopy gene families.

```
clann> exe mydata.ph
clann> hs
clann> reconstruct speciestree memory
```

In the example above, `hs` builds a supertree from single-copy trees only, then `reconstruct` reconciles all source trees (including multicopy families) against that supertree.

See `reconstruct ?` for full options including supplying an external species tree file.


## Available Commands:


*The following commands are always available:*

	execute (exe)	- Read in a file of source trees
	help		- Display this message
	quit		- Quit Clann
	set		- Set global parameters such as optimality criterion for carrying reconstructing a supertree
	!		- Run a shell session, while preserving the current Clann session (type 'exit' to return)
	tips		- Show tips and hints for better use of Clann
	log		- Control logging of screen output to a log file

*The following commands are only available when there are source trees in memory:*

**Supertree reconstruction:**

	hs		- Heuristic search for the best supertree (SPR moves, parallelised with OpenMP)
			  Supports best-improvement strategy, convergence detection, and visited-topology deduplication
			  (multicopy trees auto-excluded when delimiter mode is active)
	bootstrap	- Bootstrap supertree analysis; parallelised across replicates (all criteria supported)
	nj		- Construct a neighbour-joining supertree
			  (multicopy trees auto-excluded when delimiter mode is active)
	alltrees	- Exhaustively search all possible supertrees
			  (multicopy trees auto-excluded when delimiter mode is active)
	usertrees	- Score candidate supertrees against source trees; supports ML topology tests
			  (Winning Sites, KH, SH), per-source-tree score matrix output (scorematrix=yes),
			  and Bryant & Steel normalising constant correction (normcorrect / normcorrect=exact)
	consensus	- Calculate a consensus tree of all trees containing all taxa

**Source tree selection and modification:**

	showtrees	- Visualise selected source trees in ASCII format (also can save selected trees to file)
	deletetrees	- Specify source trees to delete from memory (based on a variety of criteria)
			  Includes excludetrees= / includetrees= filtering by tree index range
	deletetaxa	- Specify taxa to delete from all source trees in memory (i.e. prune from the trees while preserving branch lengths)
	randomisetrees	- Randomises the source trees in memory, while preserving taxa composition in each tree

**Miscellaneous calculations:**

	rfdists		- Calculate Robinson-Foulds distances between all source trees
	generatetrees	- Generate random supertrees & assess against source trees in memory
	yaptp		- "Yet another permutation-tail-probability" test - performs a randomisation test

**Experimental Options:**

	reconstruct	- Carry out a gene-tree reconciliation (source trees against a species tree)
			  Uses all trees including multicopy families
	prunemonophylies - Prunes clades which consist of multiple sequences from the same species, to a single representative
	sprdists	- Carry out estimation of SPR distances of real data versus ideal and randomised versions of the data



In interactive mode, type a command followed by `?` to see its options: `hs ?`, `usertrees ?`, etc.
From the shell, use `clann <command> --help` for the same information.

Full descriptions of all commands and options are available in **[USER_MANUAL.md](USER_MANUAL.md)**.
