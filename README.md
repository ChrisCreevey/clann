
	*********************************************************************
	*                                                                   *
	*                             Clann                                 *
	* Investigating phylogenetic information through supertree analyses *
	*                                                                   *
	*                  lab: http://www.creeveylab.org                   *
	*                  email: chris.creevey@gmail.com                   *
	*                                                                   *
	*                 Copyright Chris Creevey 2003-2018                 *
	*                                                                   *
	*          HINT: Type "help" to see all available commands          *
	*********************************************************************

# Aim

To construct supertrees and explore the underlying phylogenomic information from gene trees. Clann implements several well known supertree methods (including the ability to apply bootstrapping to them) can calculate consensus trees and provides methods to manage sets of gene trees, pruning selected taxa, filtering based on a number of criteria, automatically pruning species-specific duplicates and more. Also implemented is a gene-tree reconciliation approach to try to utilise genetic information from more than just single-copy gene-trees, this can be used as a criterion for estimating species trees. Feed back on these approaches are welcome! Finally, Clann can also calculate Robinson-Doulds (RF) and SPR distances between sets of trees, create randomised versions of sets of trees and implementes a PTP test for informativeness (the YAPTP test) of a set of trees provided. Clann has been continually developed since 2003 and suggestions for new features and tools are welcome.


# Referencing Clann

Clann has been published in Bioinformatics in 2005 under the following title:

Creevey C.J. and McInerney J.O. 2005 Clann: investigating phylogenetic information through supertree analyses. Bioinformatics 21(3): 390-2. [Link](https://academic.oup.com/bioinformatics/article/21/3/390/238167)

The Bootstrapping and YAPTP methods and the DFIT (most similar supertree algorithm) have all been described in the paper:

Creevey C.J., Fitzpatrick, D.A., Philip, G.A., Kinsella, R.J., O’Connell M.J., Travers, S.A, Wilkinson M. and McInerney J.O. 2004 Does a tree-like phylogeny only exist at the tips in the prokaryotes? Proceedings of the Royal Society London, B series: Biological Sciences 271(1557): 2551-8. [Link](http://rspb.royalsocietypublishing.org/content/271/1557/2551)

Either or both of these publications should be cited if you use Clann in published work. 

For a list of papers that have cited Clann in their work see the [Clann google scholar page](https://scholar.google.co.uk/citations?view_op=view_citation&hl=en&citation_for_view=7JkjEd4AAAAJ:UeHWp8X0CEIC)

# Usage

Usage: 
```"clann -lnh [-c commands file] [tree file]"```

	Where [tree file] is an optional Nexus or Phylip formatted file of phylogenetic trees
	-l turn on logging of screen output to file "clann.log"
	-n turns off interactive mode - requires commands to be provided in a nexus 'clann block' or with '-c'
	-c <file name> specifies a file with commands to be executed (each seperated by ';')
	-h prints this message


## Available Commands:


*The following commands are always available:*

	execute		- Read in a file of source trees
	help		- Display this message
	quit		- Quit Clann
	set		- Set global parameters such as optimality criterion for carrying reconstructing a supertree
	!		- Run a shell session, while preserving the current Clann session (type 'exit' to return)
	tips		- Show tips and hints for better use of Clann
	log		- Control logging of screen output to a log file

*The following commands are only available when there are source trees in memory:*

**Supertree reconstruction:**

	hs		- Carry out a heuristic search for the best supertree usign the criterion selected
	bootstrap	- Carry out a bootstrap supertree analysis using the criterion selected
	nj		- Construct a neighbour-joining supertree
	alltrees	- Exhaustively search all possible supertrees
	usertrees	- Assess user-defined supertrees (from seperate file), to find the best scoring
	consensus	- Calculate a consensus tree of all trees containing all taxa

**Source tree selection and modification:**

	showtrees	- Visualise selected source trees in ASCII format (also can save selected trees to file)
	deletetrees	- Specify source trees to delete from memory (based on a variety of criteria)
	deletetaxa	- Specify taxa to delete from all source trees in memory (i.e. prune from the trees while preserving branch lengths)
	randomisetrees	- Randomises the source trees in memory, while preserving taxa composition in each tree

**Miscellaneous calculations:**

	rfdists		- Calculate Robinson-Foulds distances between all source trees
	generatetrees	- Generate random supertrees & assess  against source trees in memory
	yaptp		- "Yet another permutation-tail-probability" test - performs a randomisation test

**Experimental Options:**

	reconstruct	- Carry out a gene-tree reconciliation (source trees against a species tree)
	prunemonophylies - Prunes clades which consist of multiple sequences from the same species, to a single representative
	sprdists	- Carry out estimation of SPR distances of real data versus ideal and randomised versions of the data



Type a command followed by '?' in interactive mode to get information on the options available i.e.: "exe ?"
Full descriptions of the commands are available in the manual


