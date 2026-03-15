
	*********************************************************************
	*                                                                   *
	*                             Clann                                 *
	* Investigating phylogenetic information through supertree analyses *
	*                                                                   *
	*                  lab: http://www.creeveylab.org                   *
	*                  email: chris.creevey@gmail.com                   *
	*                                                                   *
	*                 Copyright Chris Creevey 2003-2026                 *
	*                                                                   *
	*          HINT: Type "help" to see all available commands          *
	*********************************************************************

# Aim

To construct supertrees and explore the underlying phylogenomic information from gene trees. Clann implements several well known supertree methods (including the ability to apply bootstrapping to them) can calculate consensus trees and provides methods to manage sets of gene trees, pruning selected taxa, filtering based on a number of criteria, automatically pruning species-specific duplicates and more. Also implemented is a gene-tree reconciliation approach to try to utilise genetic information from more than just single-copy gene-trees, this can be used as a criterion for estimating species trees. Feed back on these approaches are welcome! Finally, Clann can also calculate Robinson-Foulds (RF) and SPR distances between sets of trees, create randomised versions of sets of trees and implementes a PTP test for informativeness (the YAPTP test) of a set of trees provided. Clann has been continually developed since 2003 and suggestions for new features and tools are welcome.


# Referencing Clann

Clann has been published in Bioinformatics in 2005 under the following title:

Creevey C.J. and McInerney J.O. 2005 Clann: investigating phylogenetic information through supertree analyses. Bioinformatics 21(3): 390-2. [Link](https://academic.oup.com/bioinformatics/article/21/3/390/238167)

The Bootstrapping and YAPTP methods and the DFIT (most similar supertree algorithm) have all been described in the paper:

Creevey C.J., Fitzpatrick, D.A., Philip, G.A., Kinsella, R.J., O'Connell M.J., Travers, S.A, Wilkinson M. and McInerney J.O. 2004 Does a tree-like phylogeny only exist at the tips in the prokaryotes? Proceedings of the Royal Society London, B series: Biological Sciences 271(1557): 2551-8. [Link](http://rspb.royalsocietypublishing.org/content/271/1557/2551)

Either or both of these publications should be cited if you use Clann in published work.

For a list of papers that have cited Clann in their work see the [Clann google scholar page](https://scholar.google.co.uk/citations?view_op=view_citation&hl=en&citation_for_view=7JkjEd4AAAAJ:UeHWp8X0CEIC)

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

# Usage

Usage:
```"clann -lnh [-c commands file] [tree file]"```

	Where [tree file] is an optional Nexus or Phylip formatted file of phylogenetic trees
	-l turn on logging of screen output to file "clann.log"
	-n turns off interactive mode - requires commands to be provided in a nexus 'clann block' or with '-c'
	-c <file name> specifies a file with commands to be executed (each separated by ';')
	-h prints this message


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

	hs		- Carry out a heuristic search for the best supertree using the criterion selected
			  (multicopy trees auto-excluded when delimiter mode is active)
	bootstrap	- Carry out a bootstrap supertree analysis using the criterion selected
	nj		- Construct a neighbour-joining supertree
			  (multicopy trees auto-excluded when delimiter mode is active)
	alltrees	- Exhaustively search all possible supertrees
			  (multicopy trees auto-excluded when delimiter mode is active)
	usertrees	- Assess user-defined supertrees (from separate file), to find the best scoring
	consensus	- Calculate a consensus tree of all trees containing all taxa

**Source tree selection and modification:**

	showtrees	- Visualise selected source trees in ASCII format (also can save selected trees to file)
	deletetrees	- Specify source trees to delete from memory (based on a variety of criteria)
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



Type a command followed by '?' in interactive mode to get information on the options available i.e.: `exe ?`

Full descriptions of all commands and options are available in **[USER_MANUAL.md](USER_MANUAL.md)**.
