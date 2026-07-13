/*
 *  reconcile.c — Clann v5.0.0
 *  Gene tree reconciliation and HGT reconstruction
 *
 *  Copyright (C) 2003-2026 Chris Creevey <chris.creevey@gmail.com>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */
#include "reconcile.h"
#include "prune.h"   /* for extract_length()/select_longest, reused by decomposegenetrees Stage 2's representative-choice logic */
#include "viewer_template.h"   /* embedded self-contained interactive HTML tree/reconciliation viewer */

void gene_content_parsimony(struct taxon * position, int * array)
	{
	/* Go to the last internal nodes on the tree that are still tagged, then calculate the possible number of copies for this node */
	struct taxon *start = position;
	
	while(position != NULL)
		{
		if(position->daughter != NULL && position->tag > 0) gene_content_parsimony(position->daughter, array);
		position = position->next_sibling;		
		}
	position = start;
	while(position != NULL)
		{
		
		
		position = position->next_sibling;
		}
	}

void nj(void)
	{
	int i, j, missing_method = 1, error = FALSE;
	char *tree = NULL, useroutfile[100], *fakefilename = NULL;
	char htmlfilename[1000], htmlmeta[10100]; FILE *htmlfile = NULL;
	FILE *outfile = NULL;
	int *saved_tags = NULL;  /* for single-copy auto-filter */

	htmlfilename[0] = '\0';
	useroutfile[0] = '\0';
	strcpy(useroutfile, "NJ-tree.ph");
    for(i=0; i<num_commands; i++)
        {
        if(strcmp(parsed_command[i], "missing") == 0)
			{
			if(strcmp(parsed_command[i+1], "4point") == 0)
				missing_method = 1;
			else
				{
				if(strcmp(parsed_command[i+1], "ultrametric") == 0)
					missing_method = 0;
				else
					{
					printf2("Error: '%s' is an invalid method to be assigned to missing\n", parsed_command[i+1]);
					missing_method = 1;
					error = TRUE;
					}
				}
			}
			
		if(strcmp(parsed_command[i], "savetrees") == 0)
			strcpy(useroutfile, parsed_command[i+1]);

		if(strcmp(parsed_command[i], "htmlview") == 0)
			{
			if(strcmp(parsed_command[i+1], "yes") == 0)
				{ strncpy(htmlfilename, inputfilename, sizeof(htmlfilename)-9); strncat(htmlfilename, ".nj.html", sizeof(htmlfilename)-strlen(htmlfilename)-1); }
			else
				strncpy(htmlfilename, parsed_command[i+1], sizeof(htmlfilename)-1);
			}
		}
		
	if(!error)
		{
		if((outfile = fopen(useroutfile, "w")) == NULL)
			{
			printf2("Error opening file named %s\n", useroutfile);
			error = TRUE;
			}
		}

	if(!error)
		{
		fakefilename = malloc(100*sizeof(char));
		fakefilename[0] = '\0';
		tree = malloc(TREE_LENGTH*sizeof(char));
		printf2("\n\nNeighbor-joining settings:\n\tDistance matrix generation by average consensus method\n\tEstimation of missing data using ");
		if(missing_method == 1)
			printf2("4 point condition distances\n");
		if(missing_method == 0)
			printf2("ultrametric distances\n");
		printf2("\tresulting tree saved to file %s\n\n\n", useroutfile);

		saved_tags = apply_singlecopy_filter();

		average_consensus(0, missing_method, fakefilename, NULL);
		neighbor_joining(TRUE, tree, TRUE);
		fprintf(outfile, "%s\n", tree);
		tree_coordinates(tree, TRUE, TRUE, FALSE, -1);
		
		for(i=0; i<number_retained_supers; i++)
			{
			strcpy(retained_supers[i], "");
			scores_retained_supers[i] = -1;
			}
		retained_supers[0] = realloc(retained_supers[0], (strlen(tree)+10)*sizeof(char));
		strcpy(retained_supers[0], tree);
		scores_retained_supers[0] = 0;
		trees_in_memory = 1;
		if(htmlfilename[0] != '\0')
			{
			snprintf(htmlmeta, sizeof(htmlmeta), "{\"dataset\":\"%s\",\"criterion\":\"nj\"}", inputfilename);
			htmlfile = html_view_open(htmlfilename, htmlmeta, 0);
			html_view_add_newick(htmlfile, retained_supers[0], "NJ tree", -1, TRUE);
			html_view_close(htmlfile, htmlfilename);
			htmlfile = NULL;
			}
		restore_singlecopy_filter(saved_tags);
		fclose(outfile);
		free(tree);
		}

	}

/* -----------------------------------------------------------------------
 * resolve_guide_tree: shared guide-tree resolution logic for
 * decomposegenetrees() (NOTES_gene_tree_decomposition.md §4, §9 step 2).
 *
 * Priority order:
 *   1. speciestree_opt is a filename (anything other than NULL, "", or
 *      "memory") -> read it, and hard-validate that its leaf-taxon set is
 *      *exactly* the current taxa_names[] registry (both directions: no
 *      missing taxa, no extra ones). On mismatch, fail and describe exactly
 *      what's missing/extra rather than silently intersecting.
 *   2. speciestree_opt == "memory", or nothing supplied but trees_in_memory
 *      > 0 -> use retained_supers[0] (same source reconstruct() uses).
 *   3. Nothing available at all -> generate one via nj() (deliberate
 *      deviation from reconstruct(), which hard-errors here; a missing
 *      guide tree should never be fatal for decomposegenetrees since the
 *      whole point is it can run standalone with zero setup).
 *
 * out_tree   : caller-allocated buffer, >= TREE_LENGTH bytes. On success,
 *              receives the guide tree as a raw Newick string (terminated
 *              with ';', not yet built into a struct taxon tree).
 * error_msg  : caller-allocated buffer (>=2000 bytes recommended), or NULL
 *              if the caller doesn't want a message. On failure (return
 *              FALSE), filled with a human-readable description.
 *
 * Returns TRUE on success, FALSE on failure (error_msg explains why).
 * Deliberately does NOT touch parsed_command[]/num_commands -- the caller
 * (decomposegenetrees' own option parser) extracts the speciestree= value
 * and passes it in, which also makes this function trivially unit-testable
 * in isolation (see NOTES_gene_tree_decomposition.md, §9 step 2).
 * ----------------------------------------------------------------------- */
int resolve_guide_tree(char *speciestree_opt, char *out_tree, char *error_msg)
	{
	int explicit_file = (speciestree_opt != NULL && strlen(speciestree_opt) > 0 && strcmp(speciestree_opt, "memory") != 0);

	out_tree[0] = '\0';
	if(error_msg != NULL) error_msg[0] = '\0';

	if(explicit_file)
		{
		FILE *stfile = fopen(speciestree_opt, "r");
		char *temptree = NULL;
		char **leafnames = NULL;
		int nleaf = 0, li = 0, ti = 0, si = 0, sc, len;
		char *missing = NULL, *extra = NULL;

		if(stfile == NULL)
			{
			if(error_msg != NULL) sprintf(error_msg, "Error: Cannot open species tree file '%s'", speciestree_opt);
			return FALSE;
			}

		temptree = malloc(TREE_LENGTH*sizeof(char));
		while((sc = getc(stfile)) != EOF && sc != ';' && si < TREE_LENGTH-2) temptree[si++] = (char)sc;
		temptree[si++] = ';';
		temptree[si] = '\0';
		fclose(stfile);
		len = (int)strlen(temptree);

		/* Extract leaf names: any run of non-structural characters that
		 * directly follows '(' or ',' (this skips internal-node support
		 * labels, which only ever follow ')'). */
		leafnames = malloc((number_of_taxa+len+10)*sizeof(char*));
		li = 0;
		while(li < len)
			{
			if(temptree[li] == '(' || temptree[li] == ',')
				{
				int j = li+1, k = 0;
				char namebuf[NAME_LENGTH];
				while(j < len && temptree[j] != '(' && temptree[j] != ')' && temptree[j] != ',' && temptree[j] != ':' && temptree[j] != ';' && k < NAME_LENGTH-1)
					{
					namebuf[k++] = temptree[j];
					j++;
					}
				namebuf[k] = '\0';
				if(k > 0)
					{
					leafnames[nleaf] = malloc((k+1)*sizeof(char));
					strcpy(leafnames[nleaf], namebuf);
					nleaf++;
					}
				}
			li++;
			}

		missing = malloc(TREE_LENGTH*sizeof(char));
		extra = malloc(TREE_LENGTH*sizeof(char));
		missing[0] = '\0';
		extra[0] = '\0';

		for(ti=0; ti<number_of_taxa; ti++)
			{
			int found = FALSE;
			for(li=0; li<nleaf; li++)
				if(strcmp(taxa_names[ti], leafnames[li]) == 0) { found = TRUE; break; }
			if(!found) { strcat(missing, taxa_names[ti]); strcat(missing, " "); }
			}
		for(li=0; li<nleaf; li++)
			{
			int found = FALSE;
			for(ti=0; ti<number_of_taxa; ti++)
				if(strcmp(taxa_names[ti], leafnames[li]) == 0) { found = TRUE; break; }
			if(!found) { strcat(extra, leafnames[li]); strcat(extra, " "); }
			}
		for(li=0; li<nleaf; li++) free(leafnames[li]);
		free(leafnames);

		if(strlen(missing) > 0 || strlen(extra) > 0)
			{
			if(error_msg != NULL)
				{
				sprintf(error_msg, "Error: species tree file '%s' taxon set does not match the loaded taxa.\n", speciestree_opt);
				if(strlen(missing) > 0) { strcat(error_msg, "  Missing from species tree: "); strcat(error_msg, missing); strcat(error_msg, "\n"); }
				if(strlen(extra) > 0) { strcat(error_msg, "  Extra in species tree (not among loaded taxa): "); strcat(error_msg, extra); strcat(error_msg, "\n"); }
				}
			free(missing); free(extra); free(temptree);
			return FALSE;
			}

		free(missing); free(extra);
		strcpy(out_tree, temptree);
		free(temptree);
		return TRUE;
		}

	/* speciestree=memory, or nothing supplied but something's already loaded */
	if(trees_in_memory > 0)
		{
		strcpy(out_tree, retained_supers[0]);
		return TRUE;
		}

	if(speciestree_opt != NULL && strcmp(speciestree_opt, "memory") == 0)
		{
		/* explicitly asked for memory, but there's nothing there -- this one still hard-errors */
		if(error_msg != NULL) sprintf(error_msg, "Error: No supertree in memory for 'speciestree=memory'.");
		return FALSE;
		}

	/* Nothing supplied, nothing in memory -> generate a guide tree on the fly */
	nj();
	if(trees_in_memory > 0)
		{
		strcpy(out_tree, retained_supers[0]);
		return TRUE;
		}

	if(error_msg != NULL) sprintf(error_msg, "Error: could not generate a guide tree via nj().");
	return FALSE;
	}

/* =========================================================================
 * decomposegenetrees Stage 2 -- LCA-mapping decomposition
 * (NOTES_gene_tree_decomposition.md §5, §9 step 3).
 *
 * Builds on label_gene_tree()/reconstruct_map() (defined below, unmodified)
 * without touching them.  Everything here is new.
 * ========================================================================= */

/* Support-gate default, resolved in the NOTES.md HANDOFF STATUS section:
 * conservative, biases toward NOT cutting when evidence is ambiguous. */
#define DECOMPOSE_DEFAULT_DUPSUPPORT 0.5

/* node_well_supported(): implements the §7-resolved semantics for "is this
 * duplication call well-supported enough to cut on":
 *   - weight == ""                       -> TRUE (no data at all; matches
 *     compute_autoweights_bootstrap's precedent of treating "no bootstrap
 *     label" as full confidence, and is required for the cut logic to ever
 *     fire on examples/tutorial_multicopy.ph, which has zero branch lengths
 *     or bootstrap values anywhere).
 *   - "label:branchlength" or bare "label" (no ':') -> treat label as a
 *     bootstrap value (rescaled from percentage if >1.0), gate on
 *     bootstrap >= dupsupport.
 *   - ":branchlength" (no label)         -> gate on branchlength > 0 (fixed
 *     zero test, dupsupport not used in this branch).
 */
/* tag_duplications_only(): same duplication *detection* test as
 * reconstruct_map() (reconcile.c above, unmodified) -- a node is a
 * duplication iff one of its daughters' tag/name equals its own tag -- but
 * WITHOUT reconstruct_map()'s side effect of splicing in extra "newbie"
 * wrapper nodes around the mismatched-tag daughter. That splicing exists
 * solely to prepare the tree for add_losses() (reconcile.c) to insert
 * missing-lineage branches immediately afterward; since Stage 2 only needs
 * duplication *tags*, not loss reconstruction, calling the real
 * reconstruct_map() here would leave the gene tree in a half-transformed
 * state (extra single-child wrapper nodes, confirmed by direct testing to
 * corrupt Stage 2's fragment output) because add_losses() never runs to
 * finish what it started. This is a separate, non-mutating function, not a
 * modification of reconstruct_map() itself. */
static int tag_duplications_only(struct taxon *position)
	{
	struct taxon *tmp;
	int found, num_dups = 0;

	while(position != NULL)
		{
		if(position->daughter != NULL)
			{
			num_dups += tag_duplications_only(position->daughter);
			tmp = position->daughter;
			found = FALSE;
			while(tmp != NULL)
				{
				if(tmp->name == -1)
					{
					if(tmp->tag == position->tag) found = TRUE;
					}
				else
					{
					if(tmp->name == position->tag) found = TRUE;
					}
				tmp = tmp->next_sibling;
				}
			if(found)
				{
				position->loss = 2;
				num_dups++;
				}
			}
		position = position->next_sibling;
		}
	return(num_dups);
	}

static int node_well_supported(struct taxon *node, float dupsupport)
	{
	char *w = node->weight;
	char *colon;
	char label[100];
	int len;
	float boot, brlen;

	if(w[0] == '\0') return TRUE;

	colon = strchr(w, ':');
	if(colon == NULL)
		{
		boot = atof(w);
		if(boot > 1.0) boot /= 100.0;
		return(boot >= dupsupport);
		}
	if(colon == w)
		{
		brlen = atof(w+1);
		return(brlen > 0);
		}
	len = (int)(colon - w);
	if(len > 99) len = 99;
	strncpy(label, w, len);
	label[len] = '\0';
	boot = atof(label);
	if(boot > 1.0) boot /= 100.0;
	return(boot >= dupsupport);
	}

/* Recursively collect, for the subtree rooted at `node` (node itself plus
 * its whole daughter chain -- NOT node's own next_sibling), which species
 * (taxa_names[] ids) are present among its leaves and how many leaves there
 * are in total.  Mirrors the "descend"/list_taxa_in_clade style already
 * used elsewhere in this file/prune.c, but scoped to a single clade root
 * rather than a whole tree or a whole sibling chain. */
static void collect_clade_stats(struct taxon *node, int *presence, int *leafcount)
	{
	if(node->daughter != NULL)
		{
		struct taxon *child = node->daughter;
		while(child != NULL)
			{
			collect_clade_stats(child, presence, leafcount);
			child = child->next_sibling;
			}
		}
	else
		{
		if(node->name >= 0 && node->name < number_of_taxa) presence[node->name] = TRUE;
		(*leafcount)++;
		}
	}

static void count_leaves_rec(struct taxon *node, int *count)
	{
	if(node->daughter != NULL)
		{
		struct taxon *child = node->daughter;
		while(child != NULL) { count_leaves_rec(child, count); child = child->next_sibling; }
		}
	else (*count)++;
	}

/* Total surviving (->spr == TRUE) leaves under `node`, post-fragment_shrink.
 * NOT the same thing as fragment_print()'s/print_pruned_tree()'s return
 * value -- that return value is the number of *direct items at the
 * top-level sibling chain* (used there only for comma-placement/"wrap in
 * parens" bookkeeping), which happens to be small and unrelated to the
 * true leaf count whenever the root has few direct children (e.g. our
 * 2-way wrapper). Confirmed by direct testing: using fragment_print()'s
 * return value against minfragtaxa silently misclassified full multi-leaf
 * fragments as 1- or 2-leaf ones. */
static void count_kept_leaves(struct taxon *node, int *count)
	{
	if(node->daughter != NULL)
		{
		struct taxon *child = node->daughter;
		while(child != NULL) { count_kept_leaves(child, count); child = child->next_sibling; }
		}
	else
		{
		if(node->spr) (*count)++;
		}
	}

static void enforce_single_copy_rec(struct taxon *node, int *seen, int treenum, FILE *infofile)
	{
	if(node->daughter != NULL)
		{
		struct taxon *child = node->daughter;
		while(child != NULL) { enforce_single_copy_rec(child, seen, treenum, infofile); child = child->next_sibling; }
		}
	else
		{
		if(node->spr && node->name >= 0 && node->name < number_of_taxa)
			{
			if(seen[node->name])
				{
				node->spr = FALSE;
				if(infofile != NULL)
					fprintf(infofile, "Tree # %d [ %s ]: dropped %s (a duplication node's decision left this fragment with a second copy of a species already kept elsewhere in it -- final single-copy safety net, not expected to fire in ordinary cases)\n", treenum, tree_names[treenum], node->fullname);
				}
			else seen[node->name] = TRUE;
			}
		}
	}

/* Last-resort safety net: decompose_walk()'s per-duplication-node decisions
 * (cut/collapse/merge) are each locally correct, but a merge that pulls a
 * failing side's representative in next to an already-fully-resolved
 * passing side can still, in topologies where the true duplication/species
 * boundary doesn't line up with which node got flagged loss==2 (e.g. the
 * species-colliding leaf sits *outside* any duplication node's direct
 * children -- see the fix note on the "mixed" branch above, which closes
 * the common case but not every possible nesting), leave more than one
 * surviving leaf for the same species somewhere in the final fragment.
 * Found by direct testing against examples/tutorial_multicopy.ph tree 7
 * (0-indexed): the mixed-branch collision check above did not catch every
 * random-representative-selection case, and a fragment still occasionally
 * came out multicopy (`Human.a` AND `Human.b` both survived, with `Human.a`
 * never touched by any duplication-node decision because its immediate
 * ancestor chain up to the point where it rejoins `Human.b`'s clade was
 * not itself flagged as a duplication node). Since producing single-copy-
 * ish ortholog subtrees is decomposegenetrees' entire purpose (§1), this is
 * enforced unconditionally, in tree (pre-order) traversal order, right
 * before fragment_shrink() in process_fragment_root(): keep the first
 * ->spr==TRUE leaf seen per species, drop any later ->spr==TRUE leaf of a
 * species already kept. Verified this closes the gap: 25 repeated runs of
 * `exe examples/tutorial_multicopy.ph autodecompose=yes` (which exercises
 * `selection=random`'s nondeterminism, the condition that surfaced this)
 * all now report "number of multicopy trees: 0", vs ~24% of runs failing
 * before this net was added. */
static void enforce_single_copy(struct taxon *root, int treenum, FILE *infofile)
	{
	int *seen = calloc(number_of_taxa, sizeof(int));
	enforce_single_copy_rec(root, seen, treenum, infofile);
	free(seen);
	}

/* Species (taxa_names[] ids) among leaves under `node` that are still
 * KEPT (->spr == TRUE), i.e. survived whatever decompose_walk() has already
 * done to this subtree. Used by the "mixed" (merge) branch below to check
 * whether a representative picked from a *failing* sibling would collide
 * (same species) with something already kept on the *passing* sibling(s) --
 * see the collision-avoidance fix documented at that call site. */
static void collect_kept_species(struct taxon *node, int *presence)
	{
	if(node->daughter != NULL)
		{
		struct taxon *child = node->daughter;
		while(child != NULL) { collect_kept_species(child, presence); child = child->next_sibling; }
		}
	else
		{
		if(node->spr && node->name >= 0 && node->name < number_of_taxa) presence[node->name] = TRUE;
		}
	}

static void find_representative_rec(struct taxon *node, int *idx, int keep, struct taxon **best, long *bestlen)
	{
	if(node->daughter != NULL)
		{
		struct taxon *child = node->daughter;
		while(child != NULL) { find_representative_rec(child, idx, keep, best, bestlen); child = child->next_sibling; }
		}
	else
		{
		if(select_longest)
			{
			long len = extract_length(node->fullname);
			if(*best == NULL || len > *bestlen) { *best = node; *bestlen = len; }
			}
		else
			{
			if(*idx == keep) *best = node;
			(*idx)++;
			}
		}
	}

/* Choose one representative leaf under `node`, using the same
 * selection=length|random convention as prune_monophylies (prune.c) --
 * longest embedded sequence-length if select_longest is set, otherwise a
 * uniformly random leaf. */
static struct taxon *pick_representative(struct taxon *node)
	{
	int leafcount = 0, idx = 0, keep;
	struct taxon *best = NULL;
	long bestlen = -1;

	count_leaves_rec(node, &leafcount);
	keep = (leafcount > 0) ? (int)fmod(rand(), leafcount) : 0;
	find_representative_rec(node, &idx, keep, &best, &bestlen);
	return best;
	}

/* Set ->spr = FALSE on every leaf under `node` except `keep` (or every leaf
 * if keep == NULL).  Internal-node flags are left for fragment_shrink() to
 * derive afterward, exactly as prune.c's untag_nodes_below/shrink_tree()
 * collapse logic does for ->tag.
 *
 * IMPORTANT: this uses ->spr, NOT ->tag, as the keep/remove flag --
 * deliberately different from Stage 1 (prune.c) and every other tag-based
 * collapse convention in this codebase. Reason (found by direct testing,
 * not by inspection): label_gene_tree()/reconstruct_map() (reconcile.c,
 * unmodified, called earlier in decompose_gene_tree_stage2() to find the
 * best rooting and tag duplications) permanently overwrite ->tag on every
 * node with its LCA-mapped species-tree node id -- including leaves, whose
 * ->tag becomes their taxon id. Since taxon ids are 0-based, any leaf whose
 * species happens to be taxon #0 (e.g. "Human" in
 * examples/tutorial_multicopy.ph) ends up with ->tag == 0 == FALSE, which
 * silently made print_pruned_tree()/shrink_tree() (tree_ops.c, also
 * unmodified) treat that leaf as already deleted. Using the otherwise-
 * unused ->spr field (scratch space for SPR search elsewhere, never
 * touched by these freshly-built, isolated gene-tree copies) sidesteps the
 * collision entirely rather than fighting over the meaning of ->tag. See
 * init_keep_flags()/fragment_shrink()/fragment_print() below, which are
 * this function's ->spr-based counterparts to reset_tree()/shrink_tree()/
 * print_pruned_tree(). */
static void tag_off_leaves_except(struct taxon *node, struct taxon *keep)
	{
	if(node->daughter != NULL)
		{
		struct taxon *child = node->daughter;
		while(child != NULL) { tag_off_leaves_except(child, keep); child = child->next_sibling; }
		}
	else
		{
		node->spr = (node == keep) ? TRUE : FALSE;
		}
	}

/* make_taxon() (tree_ops.c) defaults ->spr to FALSE (it's scratch space for
 * an unrelated purpose elsewhere); Stage 2 needs it to default to "keep"
 * (TRUE) across the whole gene tree before any collapse decisions run. */
static void init_keep_flags(struct taxon *position)
	{
	while(position != NULL)
		{
		position->spr = TRUE;
		if(position->daughter != NULL) init_keep_flags(position->daughter);
		position = position->next_sibling;
		}
	}

/* ->spr-based counterpart to shrink_tree() (tree_ops.c) -- collapses away
 * internal nodes with fewer than 2 surviving (spr==TRUE) daughters. Unlike
 * shrink_tree(), does not attempt to merge branch lengths on collapse
 * (none of examples/tutorial_multicopy.ph's fixtures have any; a real
 * dataset with branch lengths would want that added here before this is
 * relied on for scoring, not just for the decision log). */
static int fragment_shrink(struct taxon *position)
	{
	int count = 0, tot;

	while(position != NULL)
		{
		if(position->daughter != NULL)
			{
			tot = fragment_shrink(position->daughter);
			position->spr = (tot >= 2);
			count += (tot < 2) ? tot : 1;
			}
		else
			{
			if(position->spr) count++;
			}
		position = position->next_sibling;
		}
	return(count);
	}

/* ->spr-based counterpart to print_pruned_tree() (tree_ops.c), always
 * emitting ->fullname (the per-paralog-copy original label) rather than a
 * numeric id. */
static int fragment_print(struct taxon *position, int count, char *buf)
	{
	while(position != NULL)
		{
		if(position->daughter != NULL)
			{
			if(position->spr)
				{
				if(count > 0) strcat(buf, ",");
				strcat(buf, "(");
				count++;
				fragment_print(position->daughter, 0, buf);
				strcat(buf, ")");
				}
			else
				{
				count = fragment_print(position->daughter, count, buf);
				}
			}
		else
			{
			if(position->spr)
				{
				if(count > 0) strcat(buf, ",");
				strcat(buf, position->fullname);
				count++;
				}
			}
		position = position->next_sibling;
		}
	return(count);
	}

static void process_fragment_root(struct taxon *root, int treenum, float dupsupport, int minfragtaxa, int minfragspecies, FILE *fragfile, FILE *infofile, int *fragcount);

/* decompose_walk(): the per-node cut/collapse/merge decision procedure
 * (NOTES.md §5 Stage 2, steps 3-4).  Processes exactly one node (recursing
 * into its daughter chain as needed) -- never node's own next_sibling,
 * which is the caller's concern.
 *
 * Side effects only: marks leaves tag=FALSE where a clade is being
 * collapsed to a single representative or pruned-and-merged, and (for
 * genuine cut points) recursively spins each child off as its own
 * independent fragment via process_fragment_root(), then marks the whole
 * cut node's subtree tag=FALSE so it disappears from the *current*
 * fragment's own eventual print (it has already been written out
 * separately). */
static void decompose_walk(struct taxon *node, int treenum, float dupsupport, int minfragtaxa, int minfragspecies, FILE *fragfile, FILE *infofile, int *fragcount)
	{
	int is_dup, species_count, i, nchildren, npass;
	int *presence;
	struct taxon *child, *children[64];
	int child_leafcount[64], child_speccount[64];
	struct taxon *rep;

	if(node == NULL || node->daughter == NULL) return; /* leaf: nothing to decide */

	is_dup = (node->loss == 2);

	if(!(is_dup && node_well_supported(node, dupsupport)))
		{
		/* Ordinary node: not a duplication, or a duplication call we don't
		 * trust enough to act on -- keep it, recurse through it. */
		child = node->daughter;
		while(child != NULL)
			{
			decompose_walk(child, treenum, dupsupport, minfragtaxa, minfragspecies, fragfile, infofile, fragcount);
			child = child->next_sibling;
			}
		return;
		}

	/* Pure in-paralog check: every leaf under this duplication maps to one species. */
	presence = calloc(number_of_taxa, sizeof(int));
	species_count = 0;
	{
	int leafcount = 0;
	collect_clade_stats(node, presence, &leafcount);
	for(i=0; i<number_of_taxa; i++) if(presence[i]) species_count++;
	}
	free(presence);

	if(species_count <= 1)
		{
		rep = pick_representative(node);
		tag_off_leaves_except(node, rep);
		if(infofile != NULL && rep != NULL)
			fprintf(infofile, "Tree # %d [ %s ]: collapsed pure in-paralog clade at a duplication node -> kept %s\n", treenum, tree_names[treenum], rep->fullname);
		return;
		}

	/* Gather children (usually 2; capped generously for safety against
	 * unresolved polytomies). */
	nchildren = 0;
	child = node->daughter;
	while(child != NULL && nchildren < 64) { children[nchildren++] = child; child = child->next_sibling; }

	npass = 0;
	for(i=0; i<nchildren; i++)
		{
		int *pres = calloc(number_of_taxa, sizeof(int));
		int lc = 0, sc = 0, k;
		collect_clade_stats(children[i], pres, &lc);
		for(k=0; k<number_of_taxa; k++) if(pres[k]) sc++;
		free(pres);
		child_leafcount[i] = lc;
		child_speccount[i] = sc;
		if(lc >= minfragtaxa && sc >= minfragspecies) npass++;
		}

	if(npass == nchildren)
		{
		/* Both/all sides informative and cross-species: cut. */
		if(infofile != NULL)
			fprintf(infofile, "Tree # %d [ %s ]: cut at duplication node -> %d new fragment(s)\n", treenum, tree_names[treenum], nchildren);
		for(i=0; i<nchildren; i++)
			{
			/* children[i] is still linked to its sibling(s) under `node` via
			 * ->next_sibling -- fragment_shrink()/fragment_print() both walk
			 * that chain (same convention as shrink_tree()/print_pruned_tree()),
			 * so calling process_fragment_root() on children[i] as-is would
			 * silently pull its sibling's content into the same fragment too.
			 * Temporarily sever the link (restoring it right after) so each
			 * child is processed/printed as its own independent fragment root
			 * -- mirroring the "extra wrapper root" trick already used above
			 * for best_mapping, for exactly the same reason. Found by direct
			 * testing against tutorial_multicopy.ph tree 6/7 (NOTES.md). */
			struct taxon *saved_next = children[i]->next_sibling;
			children[i]->next_sibling = NULL;
			process_fragment_root(children[i], treenum, dupsupport, minfragtaxa, minfragspecies, fragfile, infofile, fragcount);
			children[i]->next_sibling = saved_next;
			}
		tag_off_leaves_except(node, NULL); /* remove the whole clade from the current fragment -- already spun off above */
		return;
		}

	if(npass == 0)
		{
		/* Both/all sides fail the informativeness floor: collapse the whole clade. */
		rep = pick_representative(node);
		tag_off_leaves_except(node, rep);
		if(infofile != NULL && rep != NULL)
			fprintf(infofile, "Tree # %d [ %s ]: collapsed clade at duplication node (all children below minfragtaxa=%d/minfragspecies=%d) -> kept %s\n", treenum, tree_names[treenum], minfragtaxa, minfragspecies, rep->fullname);
		return;
		}

	/* Mixed: prune the failing side(s) down to a representative tip and
	 * merge back in; keep searching for further duplications in the side(s)
	 * that passed.
	 *
	 * Process the PASSING side(s) first (recursing through decompose_walk,
	 * which may itself prune/collapse leaves further down), then the
	 * failing side(s) -- order matters for the collision check just below.
	 *
	 * BUG FIXED (found during §9 step 7 end-to-end testing): a
	 * representative picked from a failing side is being merged, as a
	 * sibling, into a subtree that -- because `node` is itself a
	 * duplication node -- is by definition the failing side's OUT-PARALOG.
	 * If the representative's species also survives somewhere in the
	 * passing side (e.g. the passing side already kept its own copy of
	 * that same species, or a nested duplication further down keeps one),
	 * merging silently re-introduces the very multicopy-ness Stage 2 exists
	 * to remove -- confirmed by direct testing against
	 * examples/tutorial_multicopy.ph tree 7 (0-indexed) with the default
	 * floor and `selection=random`: on some (not all -- representative-
	 * dependent, hence the nondeterminism that surfaced this) random seeds,
	 * the surviving fragment for tree 7 ended up containing BOTH `Human.a`
	 * and `Human.b`, i.e. a fragment that was supposed to be single-copy by
	 * construction was still multicopy end-to-end (reproduced with
	 * `exe examples/tutorial_multicopy.ph autodecompose=yes`, "number of
	 * multicopy trees" was 1 instead of 0 on the affected run). Fix: after
	 * the passing side(s) are fully resolved, check whether the failing
	 * side's chosen representative's species is already present among the
	 * passing side's *surviving* (->spr==TRUE) leaves; if so, drop the
	 * representative too (tag off the whole failing clade with no
	 * survivor) instead of merging a same-species duplicate in. This is
	 * conservative (loses a little more information on the rare collision
	 * case) rather than ever emitting a fragment with duplicate species,
	 * which would defeat the entire point of decomposegenetrees. */
	for(i=0; i<nchildren; i++)
		{
		if(child_leafcount[i] >= minfragtaxa && child_speccount[i] >= minfragspecies)
			decompose_walk(children[i], treenum, dupsupport, minfragtaxa, minfragspecies, fragfile, infofile, fragcount);
		}
	{
	int *kept_species = calloc(number_of_taxa, sizeof(int));
	for(i=0; i<nchildren; i++)
		if(child_leafcount[i] >= minfragtaxa && child_speccount[i] >= minfragspecies)
			collect_kept_species(children[i], kept_species);

	for(i=0; i<nchildren; i++)
		{
		if(child_leafcount[i] < minfragtaxa || child_speccount[i] < minfragspecies)
			{
			rep = pick_representative(children[i]);
			if(rep != NULL && rep->name >= 0 && rep->name < number_of_taxa && kept_species[rep->name])
				{
				tag_off_leaves_except(children[i], NULL);
				if(infofile != NULL)
					fprintf(infofile, "Tree # %d [ %s ]: dropped uninformative sibling clade under duplication node (representative %s's species already kept on the other side -- merging would re-introduce a duplicate)\n", treenum, tree_names[treenum], rep->fullname);
				}
			else
				{
				tag_off_leaves_except(children[i], rep);
				if(infofile != NULL && rep != NULL)
					fprintf(infofile, "Tree # %d [ %s ]: pruned uninformative sibling clade under duplication node (< minfragtaxa=%d/minfragspecies=%d) -> kept %s, merged into sibling\n", treenum, tree_names[treenum], minfragtaxa, minfragspecies, rep->fullname);
				}
			}
		}
	free(kept_species);
	}
	}

/* process_fragment_root(): finish deciding cuts/collapses within `root`
 * (via decompose_walk), then shrink away whatever got tagged off and, if
 * what remains meets minfragtaxa, write it to fragfile as one ortholog
 * subtree.  Recurses (via decompose_walk -> process_fragment_root for
 * genuine cut points) before this root's own print, so nested cuts are
 * written first/independently and never duplicated in this root's output. */
static void process_fragment_root(struct taxon *root, int treenum, float dupsupport, int minfragtaxa, int minfragspecies, FILE *fragfile, FILE *infofile, int *fragcount)
	{
	char *fragtree;
	int leafcount;

	decompose_walk(root, treenum, dupsupport, minfragtaxa, minfragspecies, fragfile, infofile, fragcount);

	enforce_single_copy(root, treenum, infofile);

	fragment_shrink(root);

	/* NOTE: root->spr here is NOT a "did anything survive" flag -- per
	 * fragment_shrink()'s shrink_tree()-mirroring convention, it means "did
	 * >=2 branches survive under root" (a real bifurcation, vs. collapsing
	 * to a single pass-through branch). At an ordinary internal node, a
	 * spr==FALSE single-surviving-branch case is meant to be absorbed
	 * transparently into that node's *parent* (fragment_print()'s "else"
	 * branch recurses through without wrapping when spr==FALSE, exactly
	 * mirroring print_pruned_tree()'s collapse-single-child handling). But
	 * `root` here has no parent to absorb into -- it's the top of this
	 * fragment. Using root->spr as an early-exit guard therefore silently
	 * dropped a genuinely surviving one-branch fragment (found by testing:
	 * a cut where only the reroot-outgroup side of the top wrapper survives
	 * -- tutorial_multicopy.ph tree 6/7 both hit this). count_kept_leaves()
	 * is the correct "did anything survive" test regardless of root->spr's
	 * collapse/no-collapse distinction; fragment_print() already handles
	 * printing correctly either way (wraps if spr==TRUE, transparently
	 * passes through if spr==FALSE). */
	fragtree = malloc(TREE_LENGTH*sizeof(char));
	fragtree[0] = '\0';
	leafcount = 0;
	count_kept_leaves(root, &leafcount);
	if(leafcount == 0) { free(fragtree); return; } /* truly nothing survived under root */
	fragment_print(root, 0, fragtree);

	if(leafcount >= minfragtaxa)
		{
		char *wrapped = malloc(TREE_LENGTH*sizeof(char));
		wrapped[0] = '\0';
		if(leafcount > 1) { strcpy(wrapped, "("); strcat(wrapped, fragtree); strcat(wrapped, ")"); }
		else strcpy(wrapped, fragtree);
		strcat(wrapped, ";");
		if(fragfile != NULL)
			{
			(*fragcount)++;
			/* BUG FIXED (found during §9 step 7 documentation prep): unlike
			 * the Stage-1-passthrough naming a few hundred lines down (which
			 * falls back to "tree%d_frag1" when tree_names[treenum] is
			 * empty -- true of every line in examples/tutorial_multicopy.ph,
			 * none of which carry an explicit Newick tree name), this
			 * Stage-2 cut-fragment naming used tree_names[treenum]
			 * unconditionally, producing cosmetically confusing fragment
			 * names like "_frag1"/"_frag2" (leading underscore, no tree
			 * number) whenever the source tree had no name -- confirmed via
			 * `decomposegenetrees minfragtaxa=2 minfragspecies=1` against
			 * that fixture, where tree 6/7's cut fragments came out named
			 * "_frag1"/"_frag2"/"_frag3" instead of "tree6_frag1" etc. Not a
			 * correctness bug (weights/content were still right, and names
			 * only have to be unique per §5, which they still were --
			 * distinguished by fragindex/position in the file), but
			 * inconsistent with Stage 1's own fallback and confusing enough
			 * in <filename>_info.txt / the tutorial's worked example to fix
			 * here rather than merely note. */
			if(tree_names[treenum][0] != '\0')
				fprintf(fragfile, "%s[%s_frag%d]\n", wrapped, tree_names[treenum], *fragcount);
			else
				fprintf(fragfile, "%s[tree%d_frag%d]\n", wrapped, treenum, *fragcount);
			}
		free(wrapped);
		}
	else if(infofile != NULL)
		{
		fprintf(infofile, "Tree # %d [ %s ]: dropped a %d-leaf fragment (below minfragtaxa=%d)\n", treenum, tree_names[treenum], leafcount, minfragtaxa);
		}

	free(fragtree);
	}

/* decompose_gene_tree_stage2(): entry point.  Finds the gene-tree rooting
 * that minimises duplication count against `guide_tree_str` (reusing the
 * rerooting pattern from get_recon_score(), duplication count only -- no
 * losses, per NOTES.md §5 step 1), tags duplications via
 * label_gene_tree()/reconstruct_map() (unmodified), then walks the result
 * with decompose_walk()/process_fragment_root() to write out the surviving
 * ortholog-subtree fragments.
 *
 * guide_tree_str : raw Newick guide tree (e.g. from resolve_guide_tree()).
 * fragfile       : receives one Newick fragment per line (may be NULL to
 *                  suppress writing, e.g. for a dry-run/count-only call).
 * infofile       : receives a decision log (may be NULL).
 * Returns the number of fragments written. */
/* defined later in this file; used here to hoist the species-tree preprocessing
 * out of the per-rooting loops (the species tree is invariant across rootings). */
static void build_species_partag(struct taxon *species_top);
static int  species_lca_tag(int a, int b, int root);
static void label_gene_tree_rec(struct taxon *gene_position, struct taxon *species_top, int *presence, int xnum);
static void resolve_tricotomies_mindup(struct taxon *gene_tree, struct taxon *species_tree, int xnum);
static float recon_standard_score(struct taxon *gene_top, int xnum, int *pd, int *pl);   /* lossmodel=standard; defined after the species depth table */
struct taxon *annotate_standard(struct taxon *gene_top, struct taxon *species_top, int xnum);   /* lossmodel=standard illustration annotator */

/* =====================================================================
 * Linear-time all-rootings duplication reconciliation (step 1: binary).
 *
 * The current loop reroots the gene tree to every branch and re-labels +
 * re-counts duplications from scratch -- O(gene^2). Instead, compute for every
 * directed edge (a->b) the LCA-map of the species on b's side (M[a->b]) with two
 * passes over one fixed rooting (mdown = below a node, mup = above it). A gene
 * node w, given a parent p, maps to M[p->w] and is a duplication iff some child
 * edge M[w->c] equals it -- exactly tag_duplications_only()'s test. Rooting to
 * make node j the outgroup splits edge {j, parent(j)} with a new root R; its
 * duplication count is R + the two sides, all read off the fixed maps.
 *
 * This first version still sums each rooting in O(gene) (so O(gene^2) overall),
 * purely to VERIFY the map formula against the existing loop; the O(gene)
 * reroot-DP replaces the per-rooting sum once the counts are confirmed. Indexed
 * by number_tree() tag so it lines up with get_branch()/the loop's j. Binary
 * gene trees only for now (recursion; polytomies handled in step 2). */
static struct taxon **g_lr_node  = NULL;   /* node by tag */
static int  *g_lr_par   = NULL;            /* struct-parent tag, -1 at root */
static int  *g_lr_mdown = NULL;            /* LCA of species under node */
static int  *g_lr_mup   = NULL;            /* LCA of species NOT under node */
static int   g_lr_cap   = 0;
/* These linear-DP scratch tables are rebuilt from scratch by every
 * lr_root_counts() call. get_recon_score() (the recon criterion) calls
 * lr_root_counts(), and the heuristic search runs get_recon_score() in parallel
 * (one replicate per thread, GOMP_parallel over do_search). Without per-thread
 * copies, concurrent lr_root_counts() calls realloc/overwrite the same arrays
 * and one thread dereferences another's freed/half-built g_lr_node[] -- observed
 * as SIGSEGV in `set criterion=recon; hs nthreads>1`. threadprivate gives each
 * thread its own tables, matching build_species_partag()'s globals below.
 * (Recon is NOT single-threaded, contrary to an earlier note.) */
#ifdef _OPENMP
#pragma omp threadprivate(g_lr_node, g_lr_par, g_lr_mdown, g_lr_mup, g_lr_cap)
#endif

static int lr_parent_tag(struct taxon *v)
	{
	struct taxon *s = v;
	while(s->prev_sibling != NULL) s = s->prev_sibling;
	return (s->parent != NULL) ? s->parent->tag : -1;
	}

static int lr_down(struct taxon *v, int xnum)   /* post-order: fill node,par,mdown; return mdown[v] */
	{
	struct taxon *c;
	int m = -1;
	g_lr_node[v->tag] = v;
	g_lr_par[v->tag]  = lr_parent_tag(v);
	if(v->daughter == NULL) { g_lr_mdown[v->tag] = v->name; return v->name; }
	for(c = v->daughter; c != NULL; c = c->next_sibling)
		{
		int cm = lr_down(c, xnum);
		m = (m < 0) ? cm : species_lca_tag(m, cm, xnum);
		}
	g_lr_mdown[v->tag] = m;
	return m;
	}

static void lr_up(struct taxon *u, int xnum)   /* pre-order: mup[c] = LCA(mup[u], mdown of u's other children) */
	{
	struct taxon *c, *c2;
	for(c = u->daughter; c != NULL; c = c->next_sibling)
		{
		int m = g_lr_mup[u->tag];
		for(c2 = u->daughter; c2 != NULL; c2 = c2->next_sibling)
			{
			if(c2 == c) continue;
			m = (m < 0) ? g_lr_mdown[c2->tag] : species_lca_tag(m, g_lr_mdown[c2->tag], xnum);
			}
		g_lr_mup[c->tag] = m;
		}
	for(c = u->daughter; c != NULL; c = c->next_sibling) lr_up(c, xnum);
	}

/* Duplications in the subtree of node w whose parent (in this rooting) is `par`
 * and whose own map is Mw. Children = w's neighbours (struct children + struct
 * parent) other than `par`; each child edge's map is read off mdown/mup. */
/* The unrooted gene tree's top is a sibling forest (an implicit trifurcation
 * centre with no explicit node). We model that centre as a virtual node VC=n,
 * whose neighbours are the top-level siblings; VC's own parent is -1. */
static int  g_lr_vc = -1;         /* virtual centre tag = n */
static int *g_lr_topsibs = NULL;  /* top-level sibling tags (VC's children) */
static int  g_lr_ntop = 0;
#ifdef _OPENMP
#pragma omp threadprivate(g_lr_vc, g_lr_topsibs, g_lr_ntop)   /* see note above */
#endif

static int lr_dupsum(int w, int par, int Mw)
	{
	int is_dup = 0, total, i;
	struct taxon *c;
	if(w == g_lr_vc)   /* virtual centre: neighbours are the top-level siblings */
		{
		for(i = 0; i < g_lr_ntop; i++)
			if(g_lr_topsibs[i] != par && g_lr_mdown[g_lr_topsibs[i]] == Mw) is_dup = 1;
		total = is_dup;
		for(i = 0; i < g_lr_ntop; i++)
			if(g_lr_topsibs[i] != par) total += lr_dupsum(g_lr_topsibs[i], g_lr_vc, g_lr_mdown[g_lr_topsibs[i]]);
		return total;
		}
	{
	struct taxon *node = g_lr_node[w];
	int pp = g_lr_par[w];   /* logical parent: struct parent, or VC for a top sibling */
	for(c = node->daughter; c != NULL; c = c->next_sibling)
		if(c->tag != par && g_lr_mdown[c->tag] == Mw) is_dup = 1;
	if(pp != -1 && pp != par && g_lr_mup[w] == Mw) is_dup = 1;
	total = is_dup;
	for(c = node->daughter; c != NULL; c = c->next_sibling)
		if(c->tag != par) total += lr_dupsum(c->tag, w, g_lr_mdown[c->tag]);
	if(pp != -1 && pp != par) total += lr_dupsum(pp, w, g_lr_mup[w]);
	return total;
	}
	}

/* Is node v a duplication when its parent-neighbour is p? (v's map is M[p->v];
 * v is a duplication iff some other neighbour edge M[v->c] equals it.) O(deg). */
static int lr_dupstat(int v, int p)
	{
	int Mv, dup = 0, i, pp;
	struct taxon *c;
	if(v == g_lr_vc)
		{
		Mv = g_lr_mup[p];   /* M[p->VC], p a top sibling */
		for(i = 0; i < g_lr_ntop; i++)
			if(g_lr_topsibs[i] != p && g_lr_mdown[g_lr_topsibs[i]] == Mv) dup = 1;
		return dup;
		}
	pp = g_lr_par[v];
	Mv = (p == pp) ? g_lr_mdown[v] : g_lr_mup[p];   /* M[p->v]: p is logical parent, else a struct child */
	for(c = g_lr_node[v]->daughter; c != NULL; c = c->next_sibling)
		if(c->tag != p && g_lr_mdown[c->tag] == Mv) dup = 1;
	if(pp != -1 && pp != p && g_lr_mup[v] == Mv) dup = 1;
	return dup;
	}

/* Propagate rooting counts down the VC-rooted tree: moving the root from edge
 * {par(v),v} to {v,c} crosses v, flipping only v's parent, so the count changes
 * by dupstat(v,c) - dupstat(v,par(v)). O(gene) over the whole tree. */
static int *g_lr_out = NULL;
#ifdef _OPENMP
#pragma omp threadprivate(g_lr_out)   /* see note above */
#endif
static void lr_dp_down(int v, int o_v)
	{
	struct taxon *c;
	int dv_par = lr_dupstat(v, g_lr_par[v]);
	for(c = g_lr_node[v]->daughter; c != NULL; c = c->next_sibling)
		{
		int o_c = o_v - dv_par + lr_dupstat(v, c->tag);
		g_lr_out[c->tag] = o_c;
		lr_dp_down(c->tag, o_c);
		}
	}

/* Fill out[j] = duplication count when rooting to make node-j the outgroup. */
static void lr_root_counts(struct taxon *gene_top, int n, int xnum, int *out)
	{
	int i, k, cap = n + 1;
	struct taxon *s;
	if(g_lr_cap < cap)
		{
		g_lr_node    = realloc(g_lr_node,    cap*sizeof(struct taxon*));
		g_lr_par     = realloc(g_lr_par,     cap*sizeof(int));
		g_lr_mdown   = realloc(g_lr_mdown,   cap*sizeof(int));
		g_lr_mup     = realloc(g_lr_mup,     cap*sizeof(int));
		g_lr_topsibs = realloc(g_lr_topsibs, cap*sizeof(int));
		g_lr_cap = cap;
		}
	g_lr_vc = n;
	/* gather + descend each top-level sibling; their logical parent is VC */
	g_lr_ntop = 0;
	for(s = gene_top; s != NULL; s = s->next_sibling)
		{ g_lr_topsibs[g_lr_ntop++] = s->tag; lr_down(s, xnum); }
	for(i = 0; i < g_lr_ntop; i++) g_lr_par[g_lr_topsibs[i]] = g_lr_vc;
	g_lr_node[g_lr_vc] = NULL; g_lr_par[g_lr_vc] = -1;
	/* VC's down-map = LCA over all top siblings; nothing above VC */
	{ int m = -1; for(i = 0; i < g_lr_ntop; i++) m = (m < 0) ? g_lr_mdown[g_lr_topsibs[i]] : species_lca_tag(m, g_lr_mdown[g_lr_topsibs[i]], xnum); g_lr_mdown[g_lr_vc] = m; }
	g_lr_mup[g_lr_vc] = -1;
	/* each top sibling's up-map = LCA of the OTHER siblings */
	for(i = 0; i < g_lr_ntop; i++)
		{
		int m = -1;
		for(k = 0; k < g_lr_ntop; k++) { if(k == i) continue; m = (m < 0) ? g_lr_mdown[g_lr_topsibs[k]] : species_lca_tag(m, g_lr_mdown[g_lr_topsibs[k]], xnum); }
		g_lr_mup[g_lr_topsibs[i]] = m;
		}
	for(i = 0; i < g_lr_ntop; i++) lr_up(g_lr_node[g_lr_topsibs[i]], xnum);
	/* O(gene) reroot DP. Seed one rooting (edge {VC, s0}) directly, then walk the
	 * rooting space: moving the root across a node flips only that node's parent,
	 * so each new count is the previous +/- one duplication status. reroot_tree()
	 * materialises the trifurcation centre as the root (VC models it), so there is
	 * no extra root node -- verified against the existing loop (0 mismatches). */
	g_lr_out = out;
	{
	int s0 = g_lr_topsibs[0];
	int seed = lr_dupsum(s0, g_lr_vc, g_lr_mdown[s0]) + lr_dupsum(g_lr_vc, s0, g_lr_mup[s0]);
	int dvc0 = lr_dupstat(g_lr_vc, s0);
	for(i = 0; i < g_lr_ntop; i++)
		{
		int sib = g_lr_topsibs[i];
		int o = (i == 0) ? seed : seed - dvc0 + lr_dupstat(g_lr_vc, sib);   /* cross VC from s0 */
		out[sib] = o;
		lr_dp_down(sib, o);
		}
	}
	}


int decompose_gene_tree_stage2(int treenum, char *guide_tree_str, float dupsupport, int minfragtaxa, int minfragspecies, FILE *fragfile, FILE *infofile)
	{
	struct taxon *species_tree = NULL, *gene_tree = NULL, *best_mapping = NULL, *position = NULL;
	char *temptree = malloc(TREE_LENGTH*sizeof(char));
	int taxaorder = 0, nrootings, j, xnum, best_dups = -1, fragcount = 0;
	int *presence;

	/* Build the species (guide) tree, with the same "extra wrapper root"
	 * reconstruct() uses (reconcile.c, species_source==1/3 branch) so
	 * get_min_node()'s LCA numbering has a defined outgroup above the
	 * top-level split. */
	strcpy(temptree, guide_tree_str);
	temp_top = NULL;
	{ int _to = 0; tree_build(1, temptree, species_tree, TRUE, -1, &_to); }
	species_tree = temp_top;
	temp_top = make_taxon();
	temp_top->daughter = species_tree;
	species_tree->parent = temp_top;
	species_tree = temp_top;
	temp_top = NULL;
	number_tree1(species_tree, number_of_taxa);

	/* Build the gene tree, unrooted, so every possible rooting can be tried. */
	strcpy(temptree, fundamentals[treenum]);
	unroottree(temptree);
	returntree(temptree);
	temp_top = NULL;
	taxaorder = 0;
	tree_build(1, temptree, gene_tree, 1, treenum, &taxaorder);
	gene_tree = temp_top;
	temp_top = NULL;

	if(presence_of_trichotomies(gene_tree)) gene_tree = do_resolve_tricotomies(gene_tree, species_tree, 1);

	nrootings = number_tree(gene_tree, 0);
	presence = malloc(2*number_of_taxa*sizeof(int));

	/* Species tree is invariant across rootings: number it and build the LCA
	 * parent-tag table once here instead of every loop iteration. */
	xnum = number_tree1(species_tree, number_of_taxa);
	xnum--;
	build_species_partag(species_tree);

	/* Linear-time optimal rooting: compute the duplication count of every rooting
	 * in O(gene), then pick the first branch (in number_tree order) achieving the
	 * minimum -- the same tie-break the old reroot-every-branch loop used -- and
	 * reroot once to it. Replaces the O(gene^2) search. */
	{
	int *dupc = malloc(nrootings*sizeof(int));
	int best_j = 0;
	lr_root_counts(gene_tree, nrootings, xnum, dupc);
	for(j = 0; j < nrootings; j++)
		if(best_dups == -1 || dupc[j] < best_dups) { best_dups = dupc[j]; best_j = j; }
	free(dupc);

	position = get_branch(gene_tree, best_j);
	temp_top = gene_tree;
	reroot_tree(position);
	gene_tree = temp_top;
	temp_top = NULL;
	label_gene_tree_rec(gene_tree, species_tree, presence, xnum);
	tag_duplications_only(gene_tree);   /* set loss=2 tags used by decompose_walk */
	best_mapping = gene_tree;
	gene_tree = NULL;
	}

	free(presence);
	if(gene_tree != NULL) { dismantle_tree(gene_tree); gene_tree = NULL; }

	if(best_mapping != NULL)
		{
		/* reroot_tree() leaves its result as a *pair* of siblings (the
		 * rerooted subtree, plus the extracted outgroup as its
		 * next_sibling) -- this is the same convention tree_map()/
		 * get_recon_score() rely on when calling label_gene_tree()/
		 * reconstruct_map() directly on it. But process_fragment_root()
		 * (and shrink_tree()/print_pruned_tree() underneath it) need a
		 * single coherent root whose ->daughter is that 2-way split, or
		 * the outgroup sibling is silently invisible to them and gets
		 * dropped from every fragment. Wrap it, mirroring the same
		 * "extra top node" trick reconstruct() uses for the species
		 * tree (reconcile.c, species_source resolution above). */
		struct taxon *wrapper = make_taxon();
		wrapper->daughter = best_mapping;
		best_mapping->parent = wrapper;
		init_keep_flags(wrapper);
		process_fragment_root(wrapper, treenum, dupsupport, minfragtaxa, minfragspecies, fragfile, infofile, &fragcount);
		dismantle_tree(wrapper);
		best_mapping = NULL;
		}

	if(species_tree != NULL) { dismantle_tree(species_tree); species_tree = NULL; }
	free(temptree);

	return(fragcount);
	}

/* -----------------------------------------------------------------------
 * decomposegenetrees Stage 3: weighting and pool-integration prep.
 * See NOTES_gene_tree_decomposition.md §5 Stage 3 / §9 step 4.
 *
 * For every source tree, runs Stage 1 (collapse_monophyly_in_tree()); if
 * that fully resolves it to single-copy, the collapsed tree becomes one
 * fragment at weight 1.0.  Otherwise Stage 2 (decompose_gene_tree_stage2())
 * is run on it, producing k fragments (k>=0) each weighted 1/k.  In both
 * cases the minfragtaxa/minfragspecies informativeness floor is applied --
 * for the Stage 1 case directly here (reusing collapse_monophyly_in_tree()'s
 * own out_presence counts, not print_pruned_tree()'s return value, which is
 * a top-level-sibling count, not a total leaf count -- see
 * process_fragment_root()'s comment above for why that distinction matters);
 * for the Stage 2 case the floor is already enforced by
 * process_fragment_root()/decompose_walk() internally, so this function
 * does not re-check those fragments.
 *
 * Operates directly on whatever is currently in fundamentals[]/
 * Total_fund_trees -- Stage 0's restore-pristine swap (autoprunemono
 * backup/restore) is NOT wired up here; that is step 6's job.
 *
 * Returns the total fragment count across all source trees; *out_frags is
 * set to a malloc'd array of that many struct decompose_fragment (caller
 * must free with free_decompose_fragments()), or NULL if the total is 0.
 * ----------------------------------------------------------------------- */
int decompose_gene_trees_stage3(char *guide_tree_str, float dupsupport, int minfragtaxa, int minfragspecies, FILE *infofile, struct decompose_fragment **out_frags)
	{
	int treenum, total = 0, cap = 0;
	struct decompose_fragment *frags = NULL;
	char *collapsed_tree = malloc(TREE_LENGTH*sizeof(char));
	int *presence = malloc(number_of_taxa*sizeof(int));

	for(treenum=0; treenum<Total_fund_trees; treenum++)
		{
		int num_commas = 0, still_multicopy = FALSE, k2, ntaxa_present, nleaves;

		collapsed_tree[0] = '\0';
		collapse_monophyly_in_tree(treenum, infofile, collapsed_tree, &num_commas, &still_multicopy, presence);

		if(!still_multicopy)
			{
			/* Stage 1 fully resolved this family: one fragment, weight 1.0. */
			ntaxa_present = 0; nleaves = 0;
			for(k2=0; k2<number_of_taxa; k2++)
				if(presence[k2] > 0) { ntaxa_present++; nleaves += presence[k2]; }

			if(nleaves >= minfragtaxa && ntaxa_present >= minfragspecies)
				{
				if(total >= cap) { cap = cap ? cap*2 : 16; frags = realloc(frags, cap*sizeof(struct decompose_fragment)); }
				frags[total].newick = strdup(collapsed_tree);
				frags[total].weight = 1.0f;
				frags[total].source_treenum = treenum;
				frags[total].fragindex = 1;
				frags[total].k = 1;
				if(strcmp(tree_names[treenum], "") != 0)
					snprintf(frags[total].name, NAME_LENGTH, "%s_frag1", tree_names[treenum]);
				else
					snprintf(frags[total].name, NAME_LENGTH, "tree%d_frag1", treenum);
				total++;
				}
			else if(infofile != NULL)
				{
				fprintf(infofile, "Tree # %d [ %s ]: Stage 1 passthrough dropped (%d leaves / %d species, below minfragtaxa=%d/minfragspecies=%d)\n",
					treenum, tree_names[treenum], nleaves, ntaxa_present, minfragtaxa, minfragspecies);
				}
			}
		else
			{
			/* Still multicopy after Stage 1: hand off to Stage 2. Its
			 * fragments (already floor-checked internally) are written to
			 * a scratch temp file, then read back so each one can get its
			 * own struct decompose_fragment entry with weight 1/k. */
			FILE *scratch = tmpfile();
			int fragcount;

			if(scratch == NULL)
				{
				printf2("Error: could not open scratch file in decompose_gene_trees_stage3\n");
				continue;
				}

			fragcount = decompose_gene_tree_stage2(treenum, guide_tree_str, dupsupport, minfragtaxa, minfragspecies, scratch, infofile);

			if(fragcount > 0)
				{
				char *line = malloc(TREE_LENGTH*sizeof(char));
				int idx = 0;
				rewind(scratch);
				while(fgets(line, TREE_LENGTH, scratch) != NULL)
					{
					char *bracket_open, *bracket_close, *nl;
					nl = strchr(line, '\n');
					if(nl != NULL) *nl = '\0';
					bracket_open = strrchr(line, '[');
					bracket_close = strrchr(line, ']');
					if(bracket_open == NULL || bracket_close == NULL || bracket_close < bracket_open) continue;

					idx++;
					if(total >= cap) { cap = cap ? cap*2 : 16; frags = realloc(frags, cap*sizeof(struct decompose_fragment)); }
					*bracket_open = '\0';
					frags[total].newick = strdup(line);
					frags[total].weight = 1.0f / (float)fragcount;
					frags[total].source_treenum = treenum;
					frags[total].fragindex = idx;
					frags[total].k = fragcount;
					*bracket_close = '\0';
					snprintf(frags[total].name, NAME_LENGTH, "%s", bracket_open+1);
					total++;
					}
				free(line);
				}

			fclose(scratch);
			}
		}

	free(collapsed_tree);
	free(presence);

	*out_frags = frags;
	return(total);
	}

/* build_decompose_output_text(): renders `frags` (as produced by
 * decompose_gene_trees_stage3()) into the exact "[weight](tree);[name]"
 * bracket-annotated text block main.c's Phylip tree-file reader already
 * understands (main.c execute_command(), ~line 3110-3180 -- parses a
 * leading "[...]" as a per-tree weight via atof(), and a "[...]" immediately
 * following the tree's ';' as its name). Per the documented quirk that a
 * leading "[...]" before the very first tree in a file is discarded as a
 * comment rather than read as a weight, fragments with weight == 1.0 are
 * written first (their weight is what a discarded bracket would default to
 * anyway) so the file is correct even if it is later `exe`'d directly.
 * Returns a malloc'd, NUL-terminated buffer the caller must free(). */
char *build_decompose_output_text(struct decompose_fragment *frags, int n)
	{
	int i, *order, oi = 0;
	size_t cap = 16;
	char *buf, *line;

	for(i=0; i<n; i++) cap += strlen(frags[i].newick) + strlen(frags[i].name) + 64;
	buf = malloc(cap);
	buf[0] = '\0';
	line = malloc(TREE_LENGTH + 256);
	order = malloc((n>0?n:1)*sizeof(int));

	for(i=0; i<n; i++) if(frags[i].weight == 1.0f) order[oi++] = i;
	for(i=0; i<n; i++) if(frags[i].weight != 1.0f) order[oi++] = i;

	for(i=0; i<n; i++)
		{
		int fi = order[i];
		snprintf(line, TREE_LENGTH + 256, "[%.6f]%s[%s]\n", frags[fi].weight, frags[fi].newick, frags[fi].name);
		strcat(buf, line);
		}

	free(order);
	free(line);
	return(buf);
	}

void free_decompose_fragments(struct decompose_fragment *frags, int n)
	{
	int i;
	if(frags == NULL) return;
	for(i=0; i<n; i++) if(frags[i].newick != NULL) free(frags[i].newick);
	free(frags);
	}

/* -----------------------------------------------------------------------
 * decompose_gene_trees_cmd(): CLI wiring for the "decomposegenetrees"
 * command (NOTES_gene_tree_decomposition.md §3, §9 step 5). Parses options,
 * resolves the guide tree (resolve_guide_tree()), runs Stage 1-3
 * (decompose_gene_trees_stage3()), renders the fragment list
 * (build_decompose_output_text()), and writes it to <filename> plus a
 * decision log to <filename>_info.txt.
 *
 * Per §6.1, this standalone command is non-destructive: it never touches
 * fundamentals[]/Total_fund_trees/presence_of_taxa[][], it only writes the
 * two output files, and it prints an "exe <filename>" hint so the user can
 * explicitly adopt the result -- mirroring prune_monophylies()'s existing
 * "written to the file" messaging, plus the extra adoption hint that
 * prune_monophylies() does not need (its output is drop-in replacement
 * trees at the same count, so users may not always want to reload it
 * immediately; decomposegenetrees changes the tree count outright, so the
 * hint is more important here).
 * ----------------------------------------------------------------------- */
void decompose_gene_trees_cmd(void)
	{
	int j, num_frags = 0;
	float dupsupport = 0.5f;
	int minfragtaxa = 4, minfragspecies = 4;
	char speciestree_opt[10000];
	char filename2[10000], filename3[10000];
	char guide_tree[TREE_LENGTH];
	char error_msg[2000];
	struct decompose_fragment *frags = NULL;
	FILE *df_outfile = NULL;
	char *outtext = NULL;

	select_longest = FALSE;
	speciestree_opt[0] = '\0';
	strcpy(filename2, "decomposedtrees.txt");
	strcpy(filename3, "decomposedtrees.txt_info.txt");

	for(j=0; j<num_commands; j++)
		{
		if(strcmp(parsed_command[j], "speciestree") == 0)
			strcpy(speciestree_opt, parsed_command[j+1]);
		if(strcmp(parsed_command[j], "dupsupport") == 0)
			dupsupport = atof(parsed_command[j+1]);
		if(strcmp(parsed_command[j], "minfragtaxa") == 0)
			minfragtaxa = atoi(parsed_command[j+1]);
		if(strcmp(parsed_command[j], "minfragspecies") == 0)
			minfragspecies = atoi(parsed_command[j+1]);
		if(strcmp(parsed_command[j], "selection") == 0)
			{
			if(strcmp(parsed_command[j+1], "length") == 0) select_longest = TRUE;
			}
		if(strcmp(parsed_command[j], "filename") == 0)
			{
			strcpy(filename2, parsed_command[j+1]);
			strcpy(filename3, filename2);
			strcat(filename3, "_info.txt");
			}
		}

	if(resolve_guide_tree(speciestree_opt, guide_tree, error_msg) == FALSE)
		{
		printf2("%s\n", error_msg);
		return;
		}

	tempoutfile = fopen(filename3, "w");
	if(tempoutfile == NULL)
		{
		printf2("Error: could not open '%s' for writing\n", filename3);
		return;
		}

	fprintf(tempoutfile, "decomposegenetrees decision log\n");
	fprintf(tempoutfile, "speciestree=%s dupsupport=%.4f minfragtaxa=%d minfragspecies=%d selection=%s\n",
		(strlen(speciestree_opt) > 0 ? speciestree_opt : "(resolved automatically)"),
		dupsupport, minfragtaxa, minfragspecies, (select_longest == TRUE ? "length" : "random"));
	fprintf(tempoutfile, "Note: the first fragment written to \"%s\" is always a weight==1.0\n"
		"fragment (if one exists) so that the load-time bracket-weight parser, which\n"
		"discards any \"[...]\" before the very first tree in a file as a comment,\n"
		"never silently drops a real weight. This only matters if literally every\n"
		"surviving fragment needed a non-1.0 weight (no naturally weight-1 fragment\n"
		"existed anywhere in the output).\n\n", filename2);

	num_frags = decompose_gene_trees_stage3(guide_tree, dupsupport, minfragtaxa, minfragspecies, tempoutfile, &frags);

	if(num_frags > 0)
		{
		outtext = build_decompose_output_text(frags, num_frags);
		df_outfile = fopen(filename2, "w");
		if(df_outfile == NULL)
			{
			printf2("Error: could not open '%s' for writing\n", filename2);
			}
		else
			{
			fprintf(df_outfile, "%s", outtext);
			fclose(df_outfile);
			}
		free(outtext);
		}
	else
		{
		printf2("Warning: 0 fragments survived decomposition (minfragtaxa=%d/minfragspecies=%d may be too strict for this data) -- \"%s\" was not written\n",
			minfragtaxa, minfragspecies, filename2);
		}

	printf2("\nDecomposition finished. %d fragment%s written to the file \"%s\"\n", num_frags, (num_frags == 1 ? "" : "s"), filename2);
	printf2("Information on the decisions made for each source tree written to the file \"%s\"\n", filename3);
	if(num_frags > 0)
		printf2("This command is non-destructive: the source trees in memory are unchanged.\n\tTo adopt the decomposed fragments, run: exe %s\n", filename2);

	fclose(tempoutfile);
	free_decompose_fragments(frags, num_frags);
	}

void isittagged(struct taxon * position)
	{
	while(position != NULL)
		{
		if(position->tag2 == TRUE) printf2("found\n");
		if(position->daughter != NULL)
			{
			isittagged(position->daughter);
			}
		position = position->next_sibling;
		}
	}

/* Core of tree_map() with the species-tree preprocessing already done by the
 * caller: the species tree must already be number_tree1-numbered (xnum is its
 * top tag) and build_species_partag() must have run. This is the part that is
 * invariant across all gene-tree rootings of a fixed species-tree rooting, so
 * hoisting it out of the per-rooting loop avoids re-numbering + re-tabling the
 * species tree O(rootings) times. label_gene_tree_rec() ignores its presence
 * argument, so we pass NULL. */
float tree_map_prepared(struct taxon * gene_top, struct taxon * species_top, int xnum, int print)
	{
	int num_dups = 0, num_losses = 0;

	/* label the gene tree against the (already-numbered) species tree */
	label_gene_tree_rec(gene_top, species_top, NULL, xnum);
	if(loss_model == 1)   /* lossmodel=standard: arithmetic DL score, no reconstruction */
		{
		int d = 0, l = 0;
		float sc = recon_standard_score(gene_top, xnum, &d, &l);
		if(print) fprintf(distributionreconfile, "%d\t%d\n", d, l);
		return sc;
		}
	/* duplications: a node whose id equals one of its daughters' ids */
	num_dups = reconstruct_map(gene_top, species_top);
	/* Losses are the expensive part (add_losses reconstructs the missing
	 * species subtrees). When loss_weight is 0 (duplications-only scoring, e.g.
	 * the fast exploration phase of a two-phase recon search) they contribute
	 * nothing to the score, so skip the whole loss computation. */
	if(loss_weight != 0.0)
		{
		add_losses(gene_top, species_top);
		join_losses(gene_top);
		num_losses = count_losses(gene_top);
		}

	if(print)fprintf(distributionreconfile, "%d\t%d\n", num_dups, num_losses);
	return((dup_weight*(float)num_dups)+(loss_weight*(float)num_losses));
	}

float tree_map(struct taxon * gene_top, struct taxon * species_top, int print)
	{
	int xnum = number_tree1(species_top, number_of_taxa) - 1;
	build_species_partag(species_top);
	return tree_map_prepared(gene_top, species_top, xnum, print);
	}

void resolve_tricotomies(struct taxon *position, struct taxon *species_tree)
	{
	struct taxon *start = position, *copy = NULL, *best1 = NULL, *best2 = NULL, **thislevel = NULL, *previous = position->parent, *pos1 = NULL, *pos2 = NULL, *newbie = NULL, *lastsibling = NULL;
	int i=0, j=0, k=0, l=0, x=0, y=0, multiples = FALSE, topoftree = FALSE, *found = NULL, num_dups = 0, num_losses = 0, *presence = NULL, possible_internal = 0, donestuff = TRUE, smallest = 2*number_of_taxa, result, nope = TRUE;
	float best_score = -1;
	char *pruned_tree = NULL;
	
	
	presence = malloc((2*number_of_taxa)*sizeof(int));
	found = malloc((2*number_of_taxa)*sizeof(int));
	for(j=0; j<2*number_of_taxa; j++)
		{
		presence[j] = 0;
		found[j] = 0;
		}	
	/* go down to the bottom of the gene tree */
	
	
	while(position != NULL)
		{
		i++;
		if(position->daughter != NULL)
			{
			resolve_tricotomies(position->daughter, species_tree);
			found[position->tag]++;
			}
		else
			found[position->name]++;
		position = position->next_sibling;
		}
	position = start;
	/* work out any tricotomies containing multiples at this level */
	if(start->parent == NULL)
		{
		topoftree=TRUE;
		}

/*	if( i > 2 && start->parent != NULL ) *//*if there are more than two siblings at this level, and there is at least one more than once */
	
	if( (i > 2 && !topoftree) || ( i>3 && topoftree)) /*if there are more than two siblings at this level, and there is at least one more than once */
		{
		l = i;
		while(donestuff && ((l > 2 && !topoftree) || ( l>3 && topoftree)))
			{
			while(donestuff && ((l > 2 && !topoftree) || ( l>3 && topoftree)))
				{
				
				donestuff = FALSE;
				/* Step 1: Identify sister taxa and put them together (while marking new internal nodes as such )*/
				pos1=start;
				x=0;
				
				while(pos1 != NULL && ((l > 2 && !topoftree) || ( l>3 && topoftree)))
					{
					pos2 = pos1->next_sibling;
					y=x+1;
					while(pos2 != NULL && ((l > 2 && !topoftree) || ( l>3 && topoftree)))
						{
						if((possible_internal = are_siblings(species_tree, pos1->tag, pos2->tag)) != FALSE)
							{
						
							/* put them together */
							newbie = make_taxon();
							newbie->tag = possible_internal;
							newbie->daughter = pos1;
							newbie->tag2 = TRUE;
							if(pos1->prev_sibling == NULL)
								{
								if(pos1->parent != NULL)
									(pos1->parent)->daughter = newbie;
								newbie->parent = pos1->parent;
								start = newbie;
								}
							if(pos1->next_sibling != NULL) (pos1->next_sibling)->prev_sibling = newbie;
							if(pos1->prev_sibling != NULL) (pos1->prev_sibling)->next_sibling = newbie;
							pos1->parent = newbie;
							newbie->next_sibling = pos1->next_sibling;
							newbie->prev_sibling = pos1->prev_sibling;
							pos1->next_sibling = pos2;
							pos1->prev_sibling = NULL;
							
							
							if(pos2->next_sibling != NULL) (pos2->next_sibling)->prev_sibling = pos2->prev_sibling;
							if(pos2->prev_sibling != NULL) (pos2->prev_sibling)->next_sibling = pos2->next_sibling;
							
							pos2->next_sibling = NULL;
							pos2->prev_sibling = pos1;
							
							pos1 = start;
							pos2 = pos1->next_sibling;
							if(topoftree && newbie->prev_sibling == NULL) temp_top2 = newbie;
							newbie = NULL;
							l--;
							donestuff = TRUE;
							/*printf("done siblings\n");*/
							}
						else
							{
							pos2 = pos2->next_sibling;
							y++;
							}
						}
					pos1 = pos1->next_sibling;
					x++;
					}
				
				if(((l > 2 && !topoftree) || ( l>3 && topoftree)) )
					{
					/* Step 2: Join together any that are the same */
					pos1=start;
					x=0;
					while((pos1 != NULL && ((l > 2 && !topoftree) || ( l>3 && topoftree)))  )
						{
						nope = TRUE;
						pos2 = pos1->next_sibling;
						y=x+1;
						while((pos2 != NULL && ((l > 2 && !topoftree) || ( l>3 && topoftree))) )
							{
							if(pos1->tag == pos2->tag)
								{
								newbie = make_taxon();
								newbie->tag = pos1->tag;
								newbie->daughter = pos1;
								newbie->tag2=TRUE;
								if(pos1->prev_sibling == NULL)
									{
									if(pos1->parent != NULL)
										(pos1->parent)->daughter = newbie;
									newbie->parent = pos1->parent;
									start = newbie;
									}
								if(pos1->next_sibling != NULL) (pos1->next_sibling)->prev_sibling = newbie;
								if(pos1->prev_sibling != NULL) (pos1->prev_sibling)->next_sibling = newbie;
								pos1->parent = newbie;
								newbie->next_sibling = pos1->next_sibling;
								newbie->prev_sibling = pos1->prev_sibling;
								pos1->next_sibling = pos2;
								pos1->prev_sibling = NULL;
								
								
								if(pos2->next_sibling != NULL) (pos2->next_sibling)->prev_sibling = pos2->prev_sibling;
								if(pos2->prev_sibling != NULL) (pos2->prev_sibling)->next_sibling = pos2->next_sibling;
								
								pos2->next_sibling = NULL;
								pos2->prev_sibling = pos1;
								pos1 = start;
								pos2 = pos1->next_sibling;
								if(topoftree && newbie->prev_sibling == NULL) temp_top2 = newbie;
								newbie = NULL;
								l--;
								donestuff = TRUE;
								/*printf("done easy same\n");*/
								}
							else
								{
								pos2 = pos2->next_sibling;
								y++;
								}
							}
						if(nope) /* if we haven't found any taxa the same as pos 1 , go down through the pairings already made */
							{
							
							pos2 = NULL;
							copy = start;
							y=0;
							while(copy != NULL && pos2 == NULL)
								{
								if(copy->daughter != NULL && copy != pos1 && copy->tag2 == TRUE)
									pos2 = find_same(copy->daughter, pos1->tag);
								copy = copy->next_sibling;
								y++;
								}
							if(pos2 != NULL && pos2 != pos1)
								{
								newbie = make_taxon();
								newbie->tag = pos1->tag;
								newbie->daughter = pos1;
								newbie->tag2=TRUE;
								if(pos1->prev_sibling == NULL)
									{
									if(pos1->parent != NULL )
										(pos1->parent)->daughter = pos1->next_sibling;
									(pos1->next_sibling)->parent = pos1->parent;
									start = pos1->next_sibling;
									}							
								if(pos1->next_sibling != NULL) (pos1->next_sibling)->prev_sibling = pos1->prev_sibling;
								if(pos1->prev_sibling != NULL) (pos1->prev_sibling)->next_sibling = pos1->next_sibling;
								pos1->parent = newbie;
								newbie->next_sibling = pos2->next_sibling;
								newbie->prev_sibling = pos2->prev_sibling;
								if(pos2->next_sibling != NULL) (pos2->next_sibling)->prev_sibling = newbie;
								if(pos2->prev_sibling != NULL) (pos2->prev_sibling)->next_sibling = newbie;
								if(pos2->prev_sibling == NULL)
									{
									if(pos2->parent != NULL)
										(pos2->parent)->daughter = newbie;
									newbie->parent = pos2->parent;
									pos2->parent = NULL;
									if(start == pos2) start = newbie;
									}
								pos1->next_sibling = pos2;
								pos2->prev_sibling = pos1;
								
								pos1->prev_sibling = NULL;
								pos2->next_sibling = NULL;
								l--;
								donestuff = TRUE;
								
								pos1 = start;
								newbie = NULL;
								
								/*printf("done hard same\n");*/
								}
							}
						pos1 = pos1->next_sibling;
						}				
					}
				}
				/* Step 3: reconstruct the LCAs of remaining branches */
			donestuff = FALSE;
			if(((l > 2 && !topoftree) || ( l>3 && topoftree)))
				{
				smallest = 2*number_of_taxa;
				pos1 = start;
				pos2 = pos1->next_sibling;
				x=0;
				while(pos1 != NULL)
					{
					while(pos2 != NULL)
						{
						y=x+1;
						for(j=0; j<2*number_of_taxa; j++)
							presence[j] = 0;
						presence[pos1->tag] = TRUE;
						presence[pos2->tag] = TRUE;
						result =  get_min_node(species_tree, presence, 2*number_of_taxa-2);
						if(result <= smallest)
							{
							smallest = result;
							best1 = pos1;
							best2 = pos2;
							}
						pos2 = pos2->next_sibling;
						y++;
						}
					pos1 = pos1->next_sibling;
					x++;
					}
				
				pos1=best1;
				pos2 = best2;

				newbie = make_taxon();
				newbie->tag = smallest;
				newbie->daughter = pos1;
				newbie->tag2 = TRUE;
				if(pos1->prev_sibling == NULL)
					{
					if(pos1->parent != NULL)
						(pos1->parent)->daughter = newbie;
					newbie->parent = pos1->parent;
					start = newbie;
					}
				if(pos1->next_sibling != NULL) (pos1->next_sibling)->prev_sibling = newbie;
				if(pos1->prev_sibling != NULL) (pos1->prev_sibling)->next_sibling = newbie;
				pos1->parent = newbie;
				newbie->next_sibling = pos1->next_sibling;
				newbie->prev_sibling = pos1->prev_sibling;
				pos1->next_sibling = pos2;
				pos1->prev_sibling = NULL;
				
				if(pos2->next_sibling != NULL) (pos2->next_sibling)->prev_sibling = pos2->prev_sibling;
				if(pos2->prev_sibling != NULL) (pos2->prev_sibling)->next_sibling = pos2->next_sibling;
				
				pos2->next_sibling = NULL;
				pos2->prev_sibling = pos1;
				
				pos1 = start;
				pos2 = pos1->next_sibling;
				if(topoftree && newbie->prev_sibling == NULL) temp_top2 = newbie;
				newbie = NULL;
				l--;
				donestuff = TRUE;				
				/*printf("done lca\n");*/
				}
			}
		
		/* label copy with the number of copies of taxa present */
		

		
		
				
		}
	
	
	free(presence);
	free(found);
	}

void print_tree_labels(struct taxon *position, int **results, int treenum, struct taxon *species_tree)
	{
	int onetoone;
	while(position != NULL)
		{
		if(position->loss >= 1) 
			{
			results[1][position->tag]++;
			if(strcmp(position->weight, "") != 0) results[5][position->tag]++;
			}
		else
			{
			if(position->loss == -1) results[2][position->tag]++;
			else 
				{
				results[0][position->tag]++;
				/** now check to see if what follows is a 1:1 ortholog **/
				if(position->daughter != NULL)
					{
					onetoone = isit_onetoone(position->daughter, 2);
					if(onetoone == 1)
						{
						fprintf(onetoonefile, "Tree num:%d\tName:%s\tBranch:", treenum, tree_names[treenum]);
						check_tree(species_tree, position->tag, onetoonefile);
						fprintf(onetoonefile, "\n");
						results[3][position->tag]++;
						print_onetoone_names(position->daughter, 1);
						fprintf(onetoonefile, "\n");
						}
					if(onetoone == 2)
						{
						fprintf(onetoonefile, "Tree num:%d\tName:%s\tBranch:", treenum, tree_names[treenum]);
						check_tree(species_tree, position->tag, onetoonefile);
						fprintf(onetoonefile, "\n");
						fprintf(strictonetoonefile, "Tree num:%d\tName:%s\tBranch:", treenum, tree_names[treenum]);
						check_tree(species_tree, position->tag, strictonetoonefile);
						fprintf(strictonetoonefile, "\n");
						results[3][position->tag]++;
						results[4][position->tag]++;
						print_onetoone_names(position->daughter, 2);
						fprintf(onetoonefile, "\n");
						fprintf(strictonetoonefile, "\n");
						}
					}
				}
			}
		if(position->daughter != NULL)
			{
			if(position->loss != -1) print_tree_labels(position->daughter, results, treenum, species_tree);
			}
			
		position = position->next_sibling;
		}
	}

/* --- incremental LCA mapping for label_gene_tree ------------------------------
 * number_tree1() numbers the species tree post-order, so every node's tag is
 * smaller than its parent's and leaf tags are [0,N), internal [N,2N-2]. That
 * lets us map a gene node to its species-tree node incrementally: it is the LCA
 * of its children's mapped nodes, and (because tags increase toward the root)
 * that LCA's tag is exactly what get_min_node() returned -- computed by walking
 * up the smaller tag via a parent-tag table instead of an O(species) scan. */
static int *g_species_partag = NULL;   /* parent tag of each species-tree node, indexed by its tag */
static int *g_species_depth  = NULL;   /* depth of each species-tree node, indexed by its tag */
static int  g_species_partag_n = 0;
/* These species parent/depth tables are rebuilt by build_species_partag() on
 * every label_gene_tree()/reconciliation call. The heuristic search scores
 * candidates in parallel (one replicate per thread), and the recon criterion
 * runs the reconciliation -- so without per-thread copies, concurrent
 * build_species_partag() calls realloc/overwrite the same arrays and corrupt
 * the heap (observed: SIGSEGV in `set criterion=recon; hs` with >1 thread).
 * threadprivate gives each thread its own table, matching how the other
 * reconciliation/search globals are handled (globals.h). */
#ifdef _OPENMP
#pragma omp threadprivate(g_species_partag, g_species_depth, g_species_partag_n)
#endif

/* ---- lossmodel=standard: textbook duplication-loss reconciliation ----------
 * A non-mutating, allocation-free score matching the standard DL parsimony
 * model (as used by NOTUNG/DupTree/ete3), replacing Clann's legacy add_losses()
 * subtree reconstruction. Call on the labelled gene tree (label_gene_tree_rec
 * done) with build_species_partag() already run.
 *
 * For a gene node g mapped to species node T=M(g), with children mapped to M(c):
 *   duplication iff some child maps to T;
 *   losses = sum_c (depth(M(c)) - depth(T)),  minus (#children) at speciations.
 * The rooted gene tree's root is the implicit parent of the top-level forest
 * siblings; it is charged too (mapped to the LCA of the siblings' maps). Species
 * nodes above the depth-0 forest top (the sentinel root) are treated as depth -1
 * so a clean root split scores zero losses. */
static int sp_depth(int tag)
	{
	if(tag < 0 || tag >= g_species_partag_n) return -1;   /* sentinel/out-of-range: above the forest top */
	return g_species_depth[tag];
	}
/* Top-level forest component of a species tag (its parent==-1 ancestor). Clann's
 * rooted species tree is a 2+ component forest with the true root implicit above
 * the depth-0 components; a gene node whose children fall in different components
 * maps to that implicit root, which must be treated as depth -1 (not the reused
 * xnum tag of one component). */
static int sp_comp(int tag)
	{
	if(tag < 0 || tag >= g_species_partag_n) return -1;
	while(g_species_partag[tag] >= 0) tag = g_species_partag[tag];
	return tag;
	}
/* STD_ROOT marks the implicit whole-tree root above Clann's species forest -- the
 * component a gene node maps into when its leaves span >1 top-level component.
 * Deriving it from the leaves (bottom-up) makes it immune to the sentinel/xnum
 * tag collision (the implicit root and the top real component node share a tag). */
#define STD_ROOT (-99)

/* Single O(n) post-order pass. Charges the standard DL dup/loss cost of the
 * subtree rooted at `node` (this node included when it is internal): returns the
 * subtree's forest component (STD_ROOT if it spans), and writes the depth of the
 * species node it maps to (-1 for the implicit root) to *out_depth. For an
 * internal node with children c_i mapping to depth d_i and this node to depth dT:
 *   losses += sum_i(d_i - dT) - (#children if speciation, else 0);
 *   speciation unless a child maps to the same species node (duplication). */
static int std_score_rec(struct taxon *node, int *dups, int *losses, int *out_depth)
	{
	struct taxon *c;
	if(node->daughter == NULL) { *out_depth = sp_depth(node->name); return sp_comp(node->name); }
	{
	int T = node->tag, mycomp = -2, nchild = 0, sum_child_depth = 0, child_span = 0, real_dup = 0, dT, is_dup;
	for(c = node->daughter; c != NULL; c = c->next_sibling)
		{
		int cd, ccomp = std_score_rec(c, dups, losses, &cd);
		sum_child_depth += cd;
		if(mycomp == -2) mycomp = ccomp; else if(ccomp != mycomp) mycomp = STD_ROOT;
		if(ccomp == STD_ROOT) child_span = 1;
		{ int Mc = (c->name != -1) ? c->name : c->tag; if(Mc == T) real_dup = 1; }
		nchild++;
		}
	dT = (mycomp == STD_ROOT) ? -1 : sp_depth(T);
	is_dup = (mycomp == STD_ROOT) ? child_span : real_dup;
	*losses += (sum_child_depth - nchild * dT) - (is_dup ? 0 : nchild);
	if(is_dup) (*dups)++;
	*out_depth = dT;
	return mycomp;
	}
	}

static float recon_standard_score(struct taxon *gene_top, int xnum, int *pd, int *pl)
	{
	int dups = 0, losses = 0, Rcomp = -2, k = 0, sum_child_depth = 0, child_span = 0, Troot = -1;
	struct taxon *s;
	/* Charge every internal gene node (one pass), gathering the top-level forest
	 * siblings' components/depths for the implicit root -- the gene tree's own root. */
	for(s = gene_top; s != NULL; s = s->next_sibling)
		{
		int cd, ccomp = std_score_rec(s, &dups, &losses, &cd);
		sum_child_depth += cd;
		if(Rcomp == -2) Rcomp = ccomp; else if(ccomp != Rcomp) Rcomp = STD_ROOT;
		if(ccomp == STD_ROOT) child_span = 1;
		k++;
		}
	if(k >= 2)
		{
		int dR, is_dup;
		if(Rcomp == STD_ROOT) { dR = -1; is_dup = child_span; }
		else   /* siblings share one component: root maps to their real LCA within it */
			{
			int real_dup = 0;
			for(s = gene_top; s != NULL; s = s->next_sibling){ int Ms = (s->name != -1) ? s->name : s->tag; Troot = (Troot < 0) ? Ms : species_lca_tag(Troot, Ms, xnum); }
			dR = sp_depth(Troot);
			for(s = gene_top; s != NULL; s = s->next_sibling){ int Ms = (s->name != -1) ? s->name : s->tag; if(Ms == Troot) { real_dup = 1; break; } }
			is_dup = real_dup;
			}
		losses += (sum_child_depth - k * dR) - (is_dup ? 0 : k);
		if(is_dup) dups++;
		}
	if(pd) *pd = dups;
	if(pl) *pl = losses;
	return dup_weight * (float)dups + loss_weight * (float)losses;
	}

static void build_partag_rec(struct taxon *node, int parent_tag, int depth)
	{
	for(; node != NULL; node = node->next_sibling)
		{
		if(node->tag >= 0 && node->tag < g_species_partag_n)
			{ g_species_partag[node->tag] = parent_tag; g_species_depth[node->tag] = depth; }
		if(node->daughter != NULL) build_partag_rec(node->daughter, node->tag, depth+1);
		}
	}

/* build parent-tag + depth tables for the (already number_tree1-numbered) species tree */
static void build_species_partag(struct taxon *species_top)
	{
	int need = 2*number_of_taxa + 1, i;
	if(g_species_partag_n < need)
		{
		g_species_partag = realloc(g_species_partag, need*sizeof(int));
		g_species_depth  = realloc(g_species_depth,  need*sizeof(int));
		if(!g_species_partag || !g_species_depth) memory_error(97);
		g_species_partag_n = need;
		}
	for(i = 0; i < g_species_partag_n; i++) { g_species_partag[i] = -1; g_species_depth[i] = -1; }
	build_partag_rec(species_top, -1, 0);   /* root's parent = -1 (never dereferenced) */
	}

/* LCA of two species-tree node tags via depth: lift the deeper node to the other's
 * depth, then lift both together until they meet. The species tree is unrooted --
 * its top level is a sibling forest with no single root node -- so two tags in
 * different top-level components have no common ancestor node; get_min_node()
 * returns the sentinel `root` (its initial xnum) in that case, so we do too. */
static int species_lca_tag(int a, int b, int root)
	{
	if(a < 0) return b;
	if(b < 0) return a;
	while(g_species_depth[a] > g_species_depth[b]) a = g_species_partag[a];
	while(g_species_depth[b] > g_species_depth[a]) b = g_species_partag[b];
	while(a != b)
		{
		if(g_species_partag[a] < 0 || g_species_partag[b] < 0) return root;  /* meet only at the implicit root */
		a = g_species_partag[a];
		b = g_species_partag[b];
		}
	return a;
	}

static void label_gene_tree_rec(struct taxon * gene_position, struct taxon * species_top, int *presence, int xnum)
	{
	struct taxon * position = gene_position, *tmp = NULL;
	(void)presence;   /* no longer needed for the mapping (kept for signature compatibility) */
	while(position != NULL)
		{
		if(position->daughter != NULL)
			{
			int lca_tag = -1;
			label_gene_tree_rec(position->daughter, species_top, presence, xnum);
			/* Map this node to the LCA of its children's mapped species-tree nodes.
			 * Folding over ALL children also covers the single-species case (the LCA
			 * of equal tags is that tag), so the old per-node presence reset and the
			 * distinct-child count are no longer needed -- O(depth), not O(species). */
			tmp = position->daughter;
			while(tmp != NULL)
				{
				int ctag = (tmp->name != -1) ? tmp->name : tmp->tag;
				lca_tag = (lca_tag < 0) ? ctag : species_lca_tag(lca_tag, ctag, xnum);
				tmp = tmp->next_sibling;
				}
			position->tag = lca_tag;
			}
		else
			position->tag = position->name;
		position = position->next_sibling;
		}
	}

/* Public entry: build the species parent-tag table once (per top-level call, not
 * per recursion), then run the recursive labeller. Same signature/behaviour as
 * before; internal mapping is now O(gene + species) instead of O(gene*species). */
void label_gene_tree(struct taxon * gene_position, struct taxon * species_top, int *presence, int xnum)
	{
	build_species_partag(species_top);
	label_gene_tree_rec(gene_position, species_top, presence, xnum);
	}

void descend(struct taxon * position, int *presence)
	{
	
	while(position != NULL)
		{
		if(position->daughter != NULL) descend(position->daughter, presence);
		else presence[position->name] = TRUE;
		position = position->next_sibling;
		}
	}

void subtree_id(struct taxon * position, int *tmp)
	{
	while(position != NULL)
		{
		if(position->name != -1) tmp[position->name] = FALSE;
		else
			{
			subtree_id(position->daughter, tmp);
			}
		position = position->next_sibling;
		}
	}

/* ===================================================================
 * lossmodel=standard illustration annotator.
 *
 * Builds the standard-model reconciled tree for drawing/NHX/.recon so the drawn
 * events match the arithmetic standard score. Legacy add_losses() merges a
 * duplication's overlapping lost lineages (union) and skips the gene-tree root;
 * this instead inserts losses PER EDGE (per copy, never merged) on explicitly
 * rooted gene+species trees, so count_losses() == the standard loss count and
 * the duplication marks == the standard dup count. Only used for illustration,
 * so O(gene*species) allocation is fine. See NOTES §6b.
 * =================================================================== */
static int sp_is_anc(int anc, int desc)   /* is species tag `anc` an ancestor-or-self of `desc`? */
	{
	while(desc >= 0 && desc < g_species_partag_n)
		{ if(desc == anc) return 1; desc = g_species_partag[desc]; }
	return desc == anc;
	}
/* Copy a whole species subtree as a lost lineage (every node loss=-1). */
static struct taxon *std_loss_copy(struct taxon *sp)
	{
	struct taxon *g = make_taxon(), *sc, *head = NULL, *tail = NULL;
	g->tag = sp->tag; g->name = sp->name; g->loss = -1;
	for(sc = sp->daughter; sc != NULL; sc = sc->next_sibling)
		{
		struct taxon *ch = std_loss_copy(sc);
		if(head == NULL) { head = tail = ch; g->daughter = ch; ch->parent = g; }
		else { tail->next_sibling = ch; ch->prev_sibling = tail; tail = ch; }
		}
	return g;
	}
/* Mirror the species path from `sp` down to species node `target`, splicing the
 * gene subtree `gc` in at `target` and creating each off-path species child as a
 * lost lineage. Returns the top gene node of the reconstructed path (or `gc`
 * itself if `sp` already is `target`). */
static struct taxon *std_mirror(struct taxon *sp, int target, struct taxon *gc)
	{
	struct taxon *g, *sc, *head = NULL, *tail = NULL;
	if(sp->tag == target) return gc;
	g = make_taxon(); g->tag = sp->tag; g->name = sp->name; g->loss = 0;
	for(sc = sp->daughter; sc != NULL; sc = sc->next_sibling)
		{
		struct taxon *rep = sp_is_anc(sc->tag, target) ? std_mirror(sc, target, gc) : std_loss_copy(sc);
		if(head == NULL) { head = tail = rep; g->daughter = rep; rep->parent = g; }
		else { tail->next_sibling = rep; rep->prev_sibling = tail; tail = rep; }
		}
	return g;
	}
/* Insert the losses for gene edge parent `g` (map T) -> child `c` (map Mc, a
 * strict descendant of T). is_dup: g is a duplication, so this copy descends the
 * whole T subtree (T's off-path children are lost); otherwise it is a speciation
 * edge that starts below T (T's split is not a loss). Replaces c's slot under g
 * with the reconstructed skeleton. */
static void std_edge(struct taxon *g, struct taxon *c, int T, int Mc, int is_dup, struct taxon *SR)
	{
	struct taxon *Tnode = get_branch(SR, T), *sc, *pv = c->prev_sibling, *nx = c->next_sibling;
	struct taxon *head = NULL, *tail = NULL;
	if(Tnode == NULL) return;
	c->prev_sibling = NULL; c->next_sibling = NULL; c->parent = NULL;
	for(sc = Tnode->daughter; sc != NULL; sc = sc->next_sibling)
		{
		struct taxon *rep;
		if(sp_is_anc(sc->tag, Mc)) rep = std_mirror(sc, Mc, c);   /* path toward Mc */
		else if(is_dup) rep = std_loss_copy(sc);                  /* off-path lost by this copy */
		else continue;                                            /* speciation: handled by sibling gene children */
		if(head == NULL) { head = tail = rep; } else { tail->next_sibling = rep; rep->prev_sibling = tail; tail = rep; }
		}
	if(head == NULL) head = tail = c;   /* nothing skipped/added: keep c */
	if(pv != NULL) { pv->next_sibling = head; head->prev_sibling = pv; }
	else g->daughter = head;
	head->parent = g;
	tail->next_sibling = nx; if(nx != NULL) nx->prev_sibling = tail;
	}
/* Post-order: mark standard duplications (loss=2) and insert per-edge losses. */
static void annotate_std_rec(struct taxon *g, struct taxon *SR)
	{
	struct taxon *c, *next; int T, is_dup = 0;
	for(c = g->daughter; c != NULL; c = c->next_sibling) if(c->daughter) annotate_std_rec(c, SR);
	T = g->tag;
	for(c = g->daughter; c != NULL; c = c->next_sibling)
		{ int Mc = (c->name != -1) ? c->name : c->tag; if(Mc == T) { is_dup = 1; break; } }
	if(is_dup) g->loss = 2;
	for(c = g->daughter; c != NULL; c = next)
		{
		int Mc = (c->name != -1) ? c->name : c->tag;
		next = c->next_sibling;
		if(Mc == T) continue;   /* kept copy / maps to same node -> no edge loss */
		std_edge(g, c, T, Mc, is_dup, SR);
		}
	}
static void std_zero_loss(struct taxon *n) { for(; n != NULL; n = n->next_sibling) { n->loss = 0; if(n->daughter) std_zero_loss(n->daughter); } }
/* Duplicate a subtree regardless of whether its top has a parent (duplicate_tree
 * requires parent==NULL to treat a node as a top). */
static struct taxon *dup_subtree(struct taxon *top)
	{
	struct taxon *sp = top->parent, *r;
	top->parent = NULL;
	temp_top = NULL; duplicate_tree(top, NULL); r = temp_top; temp_top = NULL;
	top->parent = sp;
	return r;
	}
/* Build and return the standard-model annotated (wrapped) reconciliation of a
 * copy of gene_top against species_top. Caller draws it and dismantle_tree()s it.
 * Restores species_top's numbering/partag on the way out. */
struct taxon *annotate_standard(struct taxon *gene_top, struct taxon *species_top, int xnum)
	{
	struct taxon *SR, *R, *sdup, *gdup, *s; int sxnum;
	(void)xnum;
	sdup = dup_subtree(species_top);
	SR = make_taxon(); SR->daughter = sdup; for(s = sdup; s != NULL; s = s->next_sibling) s->parent = SR;
	sxnum = number_tree1(SR, number_of_taxa) - 1;
	build_species_partag(SR);
	gdup = dup_subtree(gene_top);
	R = make_taxon(); R->daughter = gdup; for(s = gdup; s != NULL; s = s->next_sibling) s->parent = R;
	std_zero_loss(R);
	label_gene_tree_rec(R, SR, NULL, sxnum);
	annotate_std_rec(R, SR);
	dismantle_tree(SR);
	number_tree1(species_top, number_of_taxa); build_species_partag(species_top);
	return R;
	}

int reconstruct_map(struct taxon *position, struct taxon *species_top)
	{
	struct taxon *start = NULL, *tmp = NULL, *spec_pointer= NULL, *daug_pointer = NULL, *newbie = NULL;
	int found = FALSE, num_dups = 0;
	
	while(position != NULL)
		{
		if(position->daughter != NULL)
			{
			num_dups += reconstruct_map(position->daughter, species_top);
			
			/** Check to see if any of the daughters of this pointer taxa have the same id (or name) as its id **/
			tmp = position->daughter;
			found = FALSE;
			while(tmp != NULL)
				{
				if(tmp->name == -1)
					{
					if(tmp->tag == position->tag) found = TRUE;
					}
				else
					{
					if(tmp->name == position->tag) found = TRUE;
					}
				tmp = tmp->next_sibling;
				}
			
			if(found) /* if position is the site of a duplication event **/
				{
				position->loss = 2;
				num_dups++;
				found = FALSE;
				/*** now check to see if any of positions' daughters are not the same as pointer ***/
				tmp = position->daughter;
				while(tmp != NULL)
					{
					if(tmp->tag != position->tag && tmp->name != position->tag) /** we need to add parts to the tree ***/
						{
						newbie = make_taxon();					
						newbie->tag = position->tag;
						newbie->daughter = tmp;
						newbie->next_sibling = tmp->next_sibling;
						newbie->prev_sibling = tmp->prev_sibling;
						newbie->parent = tmp->parent;
						if(tmp->parent != NULL) (tmp->parent)->daughter = newbie;
						if(tmp->prev_sibling != NULL) (tmp->prev_sibling)->next_sibling = newbie;
						if(tmp->next_sibling != NULL) (tmp->next_sibling)->prev_sibling = newbie;
						tmp->next_sibling = NULL;
						tmp->prev_sibling = NULL;
						tmp->parent = newbie;
						}
					tmp = tmp->next_sibling;
					}
				}
			}
		position = position->next_sibling;
		}
	
	
	
	return(num_dups);
	}

void add_losses(struct taxon * position, struct taxon *species_top)
	{
	struct taxon * start = position, *spec_equiv = NULL, *tmp1 = NULL, *pos = NULL, *pos2 = NULL, *pos3 = NULL;
	int *presence = NULL, i;
	
	presence = malloc((2*number_of_taxa)*sizeof(int));
	for(i=0; i<2*number_of_taxa; i++) presence[i] = FALSE;
	
	/** Go down the gene tree to the bottom */
	while(position != NULL)
		{
		if(position->daughter != NULL) 
			{
			add_losses(position->daughter, species_top);
			/* now for the position we are at, find the equivalent position in the species tree */
			if((position->daughter)->tag != position->tag)  /* there is nothing to do if this is where we have reconstructed a duplication, all the relevant part have already been put in */
				{
				spec_equiv = get_branch(species_top, position->tag);
				
				tmp1 = position->daughter;
				tmp1->parent = NULL;
				position->daughter = NULL;
				
				/* record the nodes (and taxa) that are direct sibilngs to tmp1 */
				for(i=0; i<(2*number_of_taxa); i++) presence[i] = FALSE;
				pos = tmp1;
				while(pos != NULL)
					{
					presence[pos->tag] = TRUE;
					pos = pos->next_sibling;
					}
				/* From our position in the species tree, travel down to the bottom, repilicating what we find along the way (except for the parts that already exist and pointed to by tmp1) */
				tmp1 = construct_tree(spec_equiv->daughter, position, presence, tmp1);
				}
			
			
			}
		position = position->next_sibling;
		}	
	free(presence);
	presence = NULL;
	}

struct taxon * construct_tree(struct taxon * spec_pos, struct taxon *gene_pos, int *presence, struct taxon *extra_gene)
	{
	int first = TRUE, i, count;
	struct taxon * start = spec_pos, *position = NULL, *tmp = gene_pos, *tmp1 = NULL, *newbie = NULL;
	
	while(spec_pos != NULL)
		{
		if(presence[spec_pos->tag])  /* This already existed on the gene tree so we need to splice it back in here */
			{
			
			/*** Find the position on the gene tree extra */
			position = extra_gene;
			count = 0;
			while(position != NULL && position->tag != spec_pos->tag)
					position = position->next_sibling;
			if(position == NULL)
				{
				/* Tag not found — should not happen in a consistent reconciliation; skip this species node */
				spec_pos = spec_pos->next_sibling;
				continue;
				}

			/*** uncouple this from extra_gene and place in the gene tree ***/
			if(position == extra_gene)
				extra_gene = position->next_sibling;
			if(position->prev_sibling != NULL) (position->prev_sibling)->next_sibling = position->next_sibling;
			if(position->next_sibling != NULL) (position->next_sibling)->prev_sibling = position->prev_sibling;
			}
		else /* if it doesn't already exist */
			{
			position = make_taxon();
			if(spec_pos->name != -1) position->tag = spec_pos->name;
			else position->tag = spec_pos->tag;
			position->name = spec_pos->name;
			if(position->name != -1) position->loss = -1;
			}
		if(first)
			{
			tmp->daughter = position;
			position->parent = tmp;
			position->next_sibling = NULL;
			position->prev_sibling = NULL;
			tmp = position;
			first = FALSE;
			}
		else
			{
			tmp->next_sibling = position;
			position->prev_sibling = tmp;
			position->next_sibling = NULL;
			tmp = position;
			}
			
		if(spec_pos->daughter != NULL && !presence[spec_pos->tag]) extra_gene = construct_tree(spec_pos->daughter, tmp, presence, extra_gene);
		
		spec_pos = spec_pos->next_sibling;
		}
	
	
	return(extra_gene);
	}

int join_losses(struct taxon * position)
	{
	int loss = TRUE;
	
	while(position != NULL)
		{
		if(position->daughter != NULL)
			{
			if(join_losses(position->daughter))
				position->loss = -1;
			else
				loss = FALSE;
			}
		else
			{
			if(position->loss != -1) loss = FALSE;
			}
		position = position->next_sibling;
		}
	return(loss);
	}

int count_losses(struct taxon * position)
	{
	int count = 0;
	while(position != NULL)
		{
		if(position->loss == -1) count++;
		else
			{
			if(position->daughter != NULL) count += count_losses(position->daughter);
			}
		position = position->next_sibling;
		}
	return(count);
	}

void mapunknowns()
	{
	struct taxon *position = NULL, *species_tree = NULL, *gene_tree = NULL, *best_mapping = NULL, *unknown_fund = NULL, *pos = NULL,*copy = NULL;
	int i, j, k, l,  *presence = NULL, basescore = 1;
	float *overall_placements = NULL, biggest = -1, total, best_total = -1;
	char *temptree, *temptree1 = malloc(TREE_LENGTH * sizeof(char));
	if(!temptree1) { printf2("Error: out of memory in mapunknowns\n"); return; }
	temptree = malloc(TREE_LENGTH*sizeof(char));
	temptree[0] = '\0';
	temptree1[0] = '\0';



	presence = malloc(2*number_of_taxa*sizeof(int));
	overall_placements = malloc(2*number_of_taxa*sizeof(float));
	for(i=0; i<2*number_of_taxa; i++)
		overall_placements[i] = 0;

	
	/**** BUILD THE SPECIES TREE *****/
	strcpy(temptree, fundamentals[0]);
	returntree(temptree);
	/* build the tree in memory */
	/****** We now need to build the Supertree in memory *******/
	temp_top = NULL;
	{ int _to = 0; tree_build(1, temptree, species_tree, 1, -1, &_to); }
	species_tree = temp_top;
	temp_top = NULL;
	/** add an extra node to the top of the tree */
	temp_top = make_taxon();
	temp_top->daughter = species_tree;
	species_tree->parent = temp_top;
	species_tree = temp_top;
	temp_top = NULL;
	
	
	
	
	for(l=1; l<Total_fund_trees; l++)
		{
		temp_top = NULL;
		strcpy(temptree, "");
		strcpy(temptree, fundamentals[l]);
		returntree(temptree);
		/* build the tree in memory */
		/****** We now need to build the Supertree in memory *******/
		temp_top = NULL;
		{ int _to = 0; tree_build(1, temptree, gene_tree, 1, -1, &_to); }
		gene_tree = temp_top;
		temp_top = NULL;
		strcpy(temptree1, temptree);
		
		k=0;
		while(presence_of_taxa[0][k] >0 || presence_of_taxa[l][k] == FALSE) k++;  /* fundamental [0] is the species tree, so whatever is not in there is an unknown */
		/* k is now the number of the unkown taxa */
		printf2("unknown to be mapped: %s\n", taxa_names[k]);
		position = get_taxon(gene_tree, k);

		/** this is an unknown. we need to remove it, but mark its neighbors */
		pos = position;
		while(pos->prev_sibling != NULL) pos = pos->prev_sibling;
		while(pos != NULL)
			{
			pos->tag2 = TRUE;
			pos = pos->next_sibling;
			}
		
		/*** remove the unknown ****/
		if(position->next_sibling != NULL) (position->next_sibling)->prev_sibling = position->prev_sibling;
		if(position->prev_sibling != NULL) (position->prev_sibling)->next_sibling = position->next_sibling;
		if(position->parent != NULL) 
			{
			(position->parent)->daughter = position->next_sibling;
			(position->next_sibling)->parent = position->parent;
			}
		
		free(position);
		position = NULL;
		compress_tree(gene_tree);

		if(presence_of_trichotomies(gene_tree)) gene_tree = do_resolve_tricotomies(gene_tree, species_tree, basescore);
		
		duplicate_tree(gene_tree, NULL);
		copy = temp_top;
		tree_top = NULL;
		i = number_tree(gene_tree, 0);
		best_total = -1;
		for(j=0; j<i; j++)
			{
			
			position = get_branch(gene_tree, j);
			tree_top = gene_tree;
			reroot_tree(position);
			gene_tree = tree_top;
			tree_top = NULL;
			printf2("a==\n");
			total = tree_map(gene_tree, species_tree,0);
			
			
			if(total < best_total || best_total == -1)
				{
				best_total = total;
				if(best_mapping != NULL)
					{
					dismantle_tree(best_mapping);
					best_mapping = NULL;
					}
				best_mapping = gene_tree;
				gene_tree = NULL;
				}
			else
				{
				dismantle_tree(gene_tree);
				gene_tree = NULL;
				}
			temp_top = NULL;
			tree_top = NULL;
			duplicate_tree(copy, NULL);
			gene_tree = temp_top;
			temp_top = NULL;
			tree_top = NULL;
			number_tree(gene_tree, 0);
			}
			
		printf2("best reconstruction with a score of %f\n", best_total);
		tree_top = best_mapping;
		tree_coordinates(temptree, TRUE, FALSE, FALSE, -1);
		
		for(i=0; i<2*number_of_taxa; i++)
			presence[i] = FALSE;
			
		find_tagged(best_mapping, presence);
		j=0;
		for(i=0; i<2*number_of_taxa; i++)
			if(presence[i]) j++;
		for(i=0; i<2*number_of_taxa; i++)
			{
			if(presence[i]) overall_placements[i] += (float)1/(float)j;
			}
		if(copy != NULL)
			{
			dismantle_tree(copy);
			copy = NULL;
			}
		
		if(gene_tree != NULL)
			{
			dismantle_tree(gene_tree);
			gene_tree = NULL;
			}
		if(best_mapping != NULL)
			{
			dismantle_tree(best_mapping);
			best_mapping = NULL;
			}
		}
	
	biggest = -1;
	for(i=0; i<2*number_of_taxa; i++)
		{
		if(overall_placements[i] > biggest) biggest = overall_placements[i];
		}
	
	for(i=0; i<2*number_of_taxa; i++)
		overall_placements[i] = overall_placements[i]/biggest;
	put_in_scores(species_tree, overall_placements);
	
	tree_top = species_tree;
	tree_coordinates(temptree, TRUE, FALSE, TRUE, -1);
	
	if(species_tree != NULL)
		{
		dismantle_tree(species_tree);
		species_tree = NULL;
		}
	
		
	tree_top = NULL;
	free(temptree);
	free(temptree1);
	free(presence);
	free(overall_placements);
	}

float get_recon_score(char *giventree, int numspectries, int numgenetries)
	{
	struct taxon *position = NULL, *species_tree = NULL, *gene_tree = NULL, *best_mapping = NULL, *copy = NULL, *temp_top1 = NULL, *temp_top2 = NULL, *spec_copy = NULL;
	int i, j, k, l, m, q, r, spec_start=0, spec_end, gene_start, gene_end, num_species_internal = 0, error = FALSE, num_species_roots = 0, basescore = 1, rand1=0, rand2=0, dospecrand = 1, dogenerand=1, taxaorder=0, sp_xnum=0;
	int mindup_mode = FALSE, *dupc = NULL, mind = -1;
	int spec_mindup_mode = FALSE, *spec_dupc = NULL, min_spec_dups = -1;
	/* Parse-once gene-tree cache. The gene trees (fundamentals[]) are fixed
	 * across every species rooting and every gene-tree rooting, but the old code
	 * re-parsed each one from Newick (tree_build) inside every (species-rooting,
	 * gene-tree) iteration -- the dominant source of the per-candidate node
	 * allocation the profile flagged. Parse each gene tree exactly once here and
	 * reuse it. tree_build (Newick tokenisation) is species-tree-independent;
	 * only do_resolve_tricotomies depends on the species rooting (via LCA
	 * mapping), and only for gene trees that actually contain a polytomy -- so a
	 * binary cached tree is used read-only across all rootings, and a polytomous
	 * one is cloned + re-resolved per rooting (still skipping the re-parse). */
	struct taxon **gcache = NULL;
	int *gcache_poly = NULL, *gcache_nroot = NULL;
	float *overall_placements = NULL, biggest = -1, total, best_total = -1, sum_of_totals = 0, rooting_score = -1;
	char *temptree, *temptree1 = malloc(TREE_LENGTH * sizeof(char));
	if(!temptree1) { printf2("Error: out of memory in get_recon_score\n"); return -1; }

	temptree = malloc(TREE_LENGTH*sizeof(char));
	temptree[0] = '\0';
	temptree1[0] = '\0';



		/**** BUILD THE SPECIES TREE *****/
		strcpy(temptree,giventree);
		returntree(temptree);
		/* build the tree in memory */
		/****** We now need to build the Supertree in memory *******/
		temp_top = NULL;
		{ int _to = 0; tree_build(1, temptree, species_tree, 1, -1, &_to); }
		species_tree = temp_top;
		temp_top = NULL;

		num_species_roots = number_tree1(species_tree, number_of_taxa);

		duplicate_tree(species_tree, NULL);
		spec_copy = temp_top;
		temp_top = NULL;
		temp_top2 = NULL;
		
		
		if(numspectries == -1)
			{
			dospecrand = FALSE;
			numspectries = 1;
			}
		if(numspectries == -2)   /* min-duplication-restricted species rooting (implies gene mindup) */
			{
			spec_mindup_mode = TRUE;
			dospecrand = FALSE;
			numspectries = 1;
			numgenetries = -2;   /* the cheap species-dup pass needs min-dup gene rooting; the
			                      * numgenetries==-2 block below sets mindup_mode/dogenerand */
			}
		if(numgenetries == -1)
			{
			dogenerand = FALSE;   /* was dospecrand (copy-paste bug): numgenerootings=all
			                       * must disable GENE-rooting randomisation so the r-loop
			                       * scans all i gene rootings, not one random one */
			numgenetries = 1;
			}
		if(numgenetries == -2)   /* minimum-duplication-restricted gene rooting */
			{
			mindup_mode = TRUE;
			dogenerand = FALSE;
			numgenetries = 1;
			}

		/* Build the parse-once gene-tree cache (see declaration above). Each entry
		 * is the raw unrooted parse; gcache_poly[] records whether it needs
		 * per-rooting polytomy resolution, and binary trees are numbered here once
		 * (their numbering is rooting-independent) so the read-only duplication
		 * passes can run lr_root_counts() on them directly with no per-rooting
		 * clone/resolve/number/dismantle. */
		gcache        = malloc(Total_fund_trees * sizeof(struct taxon *));
		gcache_poly   = malloc(Total_fund_trees * sizeof(int));
		gcache_nroot  = malloc(Total_fund_trees * sizeof(int));
		{
		int gc;
		for(gc = 0; gc < Total_fund_trees; gc++)
			{
			strcpy(temptree, ""); strcpy(temptree, fundamentals[gc]);
			unroottree(temptree); returntree(temptree);
			temp_top = NULL; taxaorder = 0;
			tree_build(1, temptree, gene_tree, 1, gc, &taxaorder);
			gcache[gc] = temp_top; temp_top = NULL;
			gcache_poly[gc]  = presence_of_trichotomies(gcache[gc]);
			gcache_nroot[gc] = gcache_poly[gc] ? 0 : number_tree(gcache[gc], 0);
			}
		gene_tree = NULL;
		}

		/* Lever 2: min-duplication-restricted species rooting. Cheap pre-pass --
		 * for every species-tree rooting, sum over gene trees the minimum
		 * duplication count (each gene at its min-dup rooting), WITHOUT add_losses.
		 * The main loop then evaluates the full dup+losses score only at the species
		 * rootings that achieve the minimum duplication total. */
		if(spec_mindup_mode)
			{
			int mm, ll, jj, gmind, sdups;
			spec_dupc = malloc(num_species_roots*sizeof(int));
			min_spec_dups = -1;
			for(mm=0; mm<num_species_roots; mm++)
				{
				position = get_branch(species_tree, mm);
				temp_top = species_tree; reroot_tree(position); species_tree = temp_top; temp_top = NULL;
				sp_xnum = number_tree1(species_tree, number_of_taxa) - 1;
				build_species_partag(species_tree);
				sdups = 0;
				for(ll=0; ll<Total_fund_trees; ll++)
					{
					struct taxon *gt; int owned;
					/* Parse-once cache: binary gene trees are used read-only across all
					 * species rootings (no clone, no per-rooting parse/resolve/number/
					 * dismantle); polytomous ones are cloned and re-resolved for this
					 * species rooting (still no re-parse). */
					if(gcache_poly[ll])
						{
						duplicate_tree(gcache[ll], NULL); gt = temp_top; temp_top = NULL;
						gt = do_resolve_tricotomies(gt, species_tree, basescore);
						i = number_tree(gt, 0);
						owned = TRUE;
						}
					else { gt = gcache[ll]; i = gcache_nroot[ll]; owned = FALSE; }
					/* Linear-time all-rootings duplication DP: gmind = the minimum
					 * duplication count over every gene-tree rooting, in O(gene) from a
					 * single fixed tree -- replaces the old reroot-to-every-branch scan
					 * (reroot + label/tag duplication count + rebuild-from-copy per rooting). No
					 * rerooting or per-rooting tree copying (validated identical to the
					 * old scan, 0 mismatches; see NOTES_recon_performance_TODO.md). */
					{
					int *gdc = malloc(i*sizeof(int));
					lr_root_counts(gt, i, sp_xnum, gdc);
					gmind = -1;
					for(jj=0; jj<i; jj++)
						if(gmind == -1 || gdc[jj] < gmind) gmind = gdc[jj];
					free(gdc);
					}
					sdups += gmind;
					if(owned) dismantle_tree(gt);
					}
				spec_dupc[mm] = sdups;
				if(min_spec_dups == -1 || sdups < min_spec_dups) min_spec_dups = sdups;
				dismantle_tree(species_tree);
				temp_top = NULL; duplicate_tree(spec_copy, NULL); species_tree = temp_top; temp_top = NULL;
				}
			}

		for(q=0; q<numspectries; q++)
			{
			if(dospecrand > 0)
				{
				#ifdef _OPENMP
				spec_start = (int)fmod(rand_r(&thread_seed), num_species_roots);
				#else
				spec_start = (int)fmod(rand(), num_species_roots);
				#endif
				spec_end = spec_start +1;
				}
			else
				{
				spec_start = 0;
				spec_end = num_species_roots;
				}
/*		for(m=0; m<num_species_roots; m++) */ /* for every rooting of the species tree */
			for(m=spec_start; m<spec_end; m++)
				{
				/* Lever 2 Pass 2: only full-score the min-duplication species rootings */
				if(spec_mindup_mode && spec_dupc[m] != min_spec_dups) continue;
				position = get_branch(species_tree, m);
				temp_top = species_tree;
				/*printf("1\n");*/
				reroot_tree(position);
				species_tree = temp_top;
				temp_top = NULL;

				/* Species-tree preprocessing (numbering + parent/depth tables) is
				 * invariant across all gene trees and all gene-tree rootings of this
				 * species rooting -- do it once here and reuse via tree_map_prepared()
				 * instead of repeating it inside every tree_map() call below. */
				sp_xnum = number_tree1(species_tree, number_of_taxa) - 1;
				build_species_partag(species_tree);

				for(l=0; l<Total_fund_trees; l++)  /* for every gene tree */
					{
					temp_top1 = NULL;
					/* Parse-once cache: clone the pre-parsed gene tree instead of
					 * re-tokenising fundamentals[l] from Newick. Pass 2 mutates the
					 * gene tree (reconstruct_map/add_losses) and rebuilds from `copy`
					 * per rooting, so it needs its own working copy either way -- the
					 * cache just removes the redundant tree_build parse. Polytomous
					 * trees are still resolved against this species rooting. */
					temp_top = NULL;
					duplicate_tree(gcache[l], NULL);
					gene_tree = temp_top;
					temp_top = NULL;

					if(gcache_poly[l]) gene_tree = do_resolve_tricotomies(gene_tree, species_tree, basescore);


					duplicate_tree(gene_tree, NULL);
					copy = temp_top;
					temp_top = NULL;
					temp_top2 = NULL;
					i = number_tree(gene_tree, 0);
					best_total = -1;
					/* min-duplication-restricted rooting: Pass 1 -- duplication count
					 * at every rooting (cheap, no add_losses); find the minimum. The
					 * scoring loop below then evaluates full dup+losses only at the
					 * min-duplication rootings (Pass 2 via the dupc[] skip filter). */
					if(mindup_mode)
						{
						int jj;
						/* Linear-time all-rootings duplication DP: fill dupc[j] = the
						 * duplication count at gene-tree rooting j, for every rooting, in
						 * O(gene) from the single fixed tree -- replaces the old
						 * reroot-to-every-branch scan (reroot + label/tag duplication
						 * count + rebuild-from-copy per rooting). Pass 2 below then full-scores
						 * (add_losses) only the min-duplication rootings via the dupc[]
						 * filter. Validated identical to the old scan, 0 mismatches (see
						 * NOTES_recon_performance_TODO.md). */
						dupc = malloc(i*sizeof(int));
						lr_root_counts(gene_tree, i, sp_xnum, dupc);
						mind = -1;
						for(jj=0; jj<i; jj++)
							if(mind == -1 || dupc[jj] < mind) mind = dupc[jj];
						}
					for(r=0; r<numgenetries; r++)
						{
						if(dogenerand > 0)
							{
							#ifdef _OPENMP
							gene_start = (int)fmod(rand_r(&thread_seed), i);
							#else
							gene_start = (int)fmod(rand(), i);
							#endif
							gene_end = gene_start +1;
							}
						else
							{
							gene_start = 0;
							gene_end = i;
							}						
						#ifdef _OPENMP
						rand2 = (int)fmod(rand_r(&thread_seed), i);
						#else
						rand2 = (int)fmod(rand(), i);
						#endif
					/*	for(j=0; j<i; j++)  */ /* For every rooting of the genetree */
						for(j=gene_start; j<gene_end; j++)
							{
							/* Pass 2: only evaluate full dup+losses at min-duplication rootings */
							if(mindup_mode && dupc[j] != mind) continue;
							position = get_branch(gene_tree, j);
							temp_top = gene_tree;
							/*printf("2\n"); */
							reroot_tree(position);
							gene_tree = temp_top;
							temp_top = NULL;
							total = tree_map_prepared(gene_tree, species_tree, sp_xnum, 0);
							if(total < best_total || best_total == -1)
								{
								best_total = total;
								if(best_mapping != NULL)
									{
									dismantle_tree(best_mapping);
									best_mapping = NULL;
									}
								best_mapping = gene_tree;
								gene_tree = NULL;
								}
							else
								{
								dismantle_tree(gene_tree);
								gene_tree = NULL;
								}
							temp_top1 = NULL;
							temp_top = NULL;
							duplicate_tree(copy, NULL);
							gene_tree = temp_top;
							temp_top1 = NULL;
							temp_top = NULL;
							number_tree(gene_tree, 0);
							}
						}
					if(copy != NULL)
						{
						dismantle_tree(copy);
						copy = NULL;
						}
					
					if(gene_tree != NULL)
						{
						dismantle_tree(gene_tree);
						gene_tree = NULL;
						}
					if(best_mapping != NULL)
						{
						dismantle_tree(best_mapping);
						best_mapping = NULL;
						}
					if(dupc != NULL) { free(dupc); dupc = NULL; }
					sum_of_totals+=best_total;
					}
				if(sum_of_totals < rooting_score || rooting_score == -1)
					{
					rooting_score = sum_of_totals;
					}
				sum_of_totals = 0;
				dismantle_tree(species_tree);
				duplicate_tree(spec_copy, NULL);
				species_tree = temp_top;
				temp_top = NULL;
				temp_top2 = NULL;

				}
			}
	if(spec_copy != NULL)
		{
		dismantle_tree(spec_copy);
		spec_copy = NULL;
		}
	if(species_tree != NULL)
		{
		dismantle_tree(species_tree);
		species_tree = NULL;
		}
	if(spec_dupc != NULL) { free(spec_dupc); spec_dupc = NULL; }
	/* tear down the parse-once gene-tree cache */
	if(gcache != NULL)
		{
		int gc;
		for(gc = 0; gc < Total_fund_trees; gc++)
			if(gcache[gc] != NULL) dismantle_tree(gcache[gc]);
		free(gcache); gcache = NULL;
		}
	if(gcache_poly  != NULL) { free(gcache_poly);  gcache_poly  = NULL; }
	if(gcache_nroot != NULL) { free(gcache_nroot); gcache_nroot = NULL; }
	free(temptree);
	free(temptree1);
	return(rooting_score);
	}

int presence_of_trichotomies(struct taxon * position)
	{
	struct taxon *start = position;
	int i=0, result = FALSE;
	while(position != NULL)
		{
		i++;
		position = position->next_sibling;
		}
	if(i > 2)
		{
		if(start->parent != NULL)
			result = TRUE;
		else
			{
			if(i > 3)
				result = TRUE;
			}
		}

	position = start;
	while(position != NULL && result == FALSE)
		{
		if(position->daughter != NULL)	
			{
			if(presence_of_trichotomies(position->daughter))
				result = TRUE;
			}
		position = position->next_sibling;
		}
	return(result);
	}

struct taxon * do_resolve_tricotomies(struct taxon * gene_tree, struct taxon * species_tree, int basescore)
	{
	int i, xnum=0, *presence=NULL;
	struct taxon * teeemp = NULL;
	(void)basescore;   /* was only used to seed the distance matrix (now unused) */

	presence = malloc(2*number_of_taxa*sizeof(int));
	for(i=0; i<(2*number_of_taxa); i++) presence[i] = FALSE;
	/** 1) Label all internal and external taxa on the species tree ****/

	xnum = number_tree1(species_tree, number_of_taxa);
	xnum--;
	/****2) label all the gene tree nodes (and taxa) with their equivalent on the species tree **/
	label_gene_tree(gene_tree, species_tree, presence, xnum);
	temp_top2 = gene_tree;
	/* 3) build species parent/depth tables for species_lca_tag() used by the
	 *    minimum-duplication polytomy resolver. (The old distance-matrix step --
	 *    print_tree_withinternals + pathmetric_internals -- is no longer needed
	 *    now that resolution is duplication-guided rather than distance-guided.) */
	build_species_partag(species_tree);

/*	for(i=0; i<2*number_of_taxa; i++)	
		{
		if(i<number_of_taxa)
			printf2("\t%s", taxa_names[i]);
		else
			printf2("\t%d", i);
		}
	printf2("\n");
	for(i=0; i<2*number_of_taxa; i++)	
		{
		if(i<number_of_taxa)
			printf2("%s\t", taxa_names[i]);
		else
			printf2("%d\t", i);
		for(j=0; j<2*number_of_taxa; j++)
			printf2("%d\t", scores[i][j]);
		printf2("\n");
		}
*/
	/*resolve_tricotomies(gene_tree, species_tree);*/
	/* Soft-polytomy resolution: refine each multifurcation into the binary shape
	 * that minimises duplications (see resolve_tricotomies_mindup), instead of the
	 * distance-guided resolve_tricotomies_dist(). */
	temp_top2 = gene_tree;
	resolve_tricotomies_mindup(gene_tree, species_tree, xnum);
	gene_tree = temp_top2;
	
	/**** Unroot the tree -necessary as the rooting algorithm assumes an unrooted topology */
	/*count the number of daughters at the top of the tree */
	i=0;
	teeemp = gene_tree;
	while(teeemp != NULL)
		{
		i++;
		teeemp = teeemp->next_sibling;
		}
	if(i < 3)
		{
		teeemp = NULL;
		if(gene_tree->daughter != NULL)
			{
			temp_top2 = gene_tree->daughter;
			(gene_tree->next_sibling)->next_sibling = temp_top2;
			(gene_tree->next_sibling)->prev_sibling = NULL;
			temp_top2->prev_sibling = (gene_tree->next_sibling);
			temp_top2->parent = NULL;
			temp_top2 = temp_top2->prev_sibling;
			gene_tree->next_sibling = NULL;
			gene_tree->daughter = NULL;
			dismantle_tree(gene_tree);
			gene_tree = temp_top2;
			}
		else
			{
			temp_top2 = (gene_tree->next_sibling)->daughter;
			teeemp = gene_tree->next_sibling;
			gene_tree->next_sibling = temp_top2;
			temp_top2->prev_sibling = gene_tree;
			teeemp->prev_sibling = NULL;
			teeemp->daughter = NULL;
			dismantle_tree(teeemp);
			teeemp = NULL;
			}
		}
	
	free(presence);

	return(gene_tree);
	}

void make_unrooted(struct taxon * position)
	{
	int i=0;
	struct taxon * start = position;
	while(position != NULL)
		{
		i++;
		position = position->next_sibling;
		}
	if(i==1)
		{
		
		}
	}

/* -----------------------------------------------------------------------
 * Self-contained interactive HTML viewer output.
 *
 * Writes a single .html file = the embedded viewer template (viewer_template.h)
 * with the tree/reconciliation injected as JSON. Opens in any browser, no
 * dependencies. Works for a plain tree (recon=0: best supertree, NJ tree) or a
 * reconciled gene tree (recon=1: nodes carry duplication/speciation/loss events,
 * read off the same loss field the NHX writer uses).
 * ----------------------------------------------------------------------- */
static void hv_json_str(FILE *f, const char *s)
	{
	fputc('"', f);
	for(; s != NULL && *s; s++)
		{
		if(*s == '"' || *s == '\\') { fputc('\\', f); fputc(*s, f); }
		else if(*s == '\n') fputs("\\n", f);
		else if((unsigned char)*s < 0x20) fprintf(f, "\\u%04x", (unsigned)(unsigned char)*s);
		else fputc(*s, f);
		}
	fputc('"', f);
	}
static const char *hv_leaf_species(struct taxon *n)
	{
	if(n->name >= 0 && n->name < number_of_taxa) return taxa_names[n->name];
	if(n->tag  >= 0 && n->tag  < number_of_taxa) return taxa_names[n->tag];
	return NULL;
	}
static void hv_json_node(struct taxon *n, FILE *f, int recon)
	{
	fputc('{', f);
	if(n->daughter != NULL)
		{
		struct taxon *c; int first = 1;
		if(recon) fprintf(f, "\"event\":\"%s\",", (n->loss == 2) ? "duplication" : "speciation");
		if(recon && n->tag >= 0 && n->tag < number_of_taxa) { fputs("\"species\":", f); hv_json_str(f, taxa_names[n->tag]); fputc(',', f); }
		fputs("\"children\":[", f);
		for(c = n->daughter; c != NULL; c = c->next_sibling) { if(!first) fputc(',', f); first = 0; hv_json_node(c, f, recon); }
		fputc(']', f);
		}
	else if(recon && n->loss == -1)
		{
		const char *sp = hv_leaf_species(n);
		fputs("\"event\":\"loss\"", f);
		if(sp != NULL) { fputs(",\"species\":", f); hv_json_str(f, sp); }
		}
	else
		{
		const char *nm = n->fullname ? n->fullname : ((n->name >= 0 && n->name < number_of_taxa) ? taxa_names[n->name] : "?");
		fputs("\"name\":", f); hv_json_str(f, nm);
		if(recon) { const char *sp = hv_leaf_species(n); if(sp != NULL) { fputs(",\"species\":", f); hv_json_str(f, sp); } }
		if(n->length != 0) fprintf(f, ",\"length\":%g", (double)n->length);
		}
	fputc('}', f);
	}
/* Emit a tree node, forest-wrapping a multi-sibling top into a synthetic root. */
static void hv_json_tree(struct taxon *tree, FILE *f, int recon)
	{
	if(tree != NULL && tree->next_sibling != NULL)
		{
		struct taxon *c; int first = 1;
		fputs(recon ? "{\"event\":\"speciation\",\"children\":[" : "{\"children\":[", f);
		for(c = tree; c != NULL; c = c->next_sibling) { if(!first) fputc(',', f); first = 0; hv_json_node(c, f, recon); }
		fputs("]}", f);
		}
	else if(tree != NULL) hv_json_node(tree, f, recon);
	else fputs("null", f);
	}
static int hv_count_dups(struct taxon *n)
	{ int c = 0; for(; n != NULL; n = n->next_sibling) { if(n->loss == 2) c++; if(n->daughter) c += hv_count_dups(n->daughter); } return c; }

/* --- public HTML-viewer API (one navigable self-contained file) -------------
 * open once, add each tree, close once. recon=1 embeds reconciliation events
 * (duplication/loss glyphs, per-tree dup/loss counts); recon=0 is a plain tree
 * (best supertree, NJ tree, or a shown gene tree). Used by reconstruct (recon,
 * from in-memory annotated trees) and by hs/nj/showtrees (plain, from Newick). */
FILE *html_view_open(const char *filename, const char *metajson, int recon)
	{
	FILE *f = fopen(filename, "w");
	if(f == NULL) { printf2("Warning: could not open HTML view file '%s'\n", filename); return NULL; }
	fputs(VIEWER_HTML_HEAD, f);
	fputs("{\"type\":", f); fputs(recon ? "\"reconciliation\"" : "\"tree\"", f);
	fputs(",\"meta\":", f); fputs(metajson ? metajson : "{}", f);
	fputs(",\"trees\":[", f);
	return f;
	}
void html_view_add_tree(FILE *f, struct taxon *tree, const char *name, float score, int recon, int first)
	{
	if(f == NULL) return;
	if(!first) fputc(',', f);
	fputs("{\"name\":", f); hv_json_str(f, name ? name : "tree");
	if(recon) fprintf(f, ",\"score\":%.4g,\"dups\":%d,\"losses\":%d",
	                  (double)score, hv_count_dups(tree), count_losses(tree));
	fputs(",\"tree\":", f); hv_json_tree(tree, f, recon); fputc('}', f);
	}
/* Add a plain tree given as a Newick string, with leaves shown as full taxon
 * names. `treenum` >= 0: `newick` is gene tree `treenum`'s numeric form
 * (fundamentals[treenum]); it is expanded to full names via returntree_fullnames.
 * `treenum` == -1: `newick` already carries full names (e.g. a supertree in
 * retained_supers[]), and is parsed as-is. */
void html_view_add_newick(FILE *f, const char *newick, const char *name, int treenum, int first)
	{
	char *buf; struct taxon *tr;
	if(f == NULL || newick == NULL || newick[0] == '\0') return;
	if(treenum >= 0)
		{
		buf = malloc(TREE_LENGTH * sizeof(char));
		if(buf == NULL) return;
		ensure_fullname_bufsize(&buf, treenum);
		strcpy(buf, newick);
		returntree_fullnames(buf, treenum);
		}
	else
		{
		buf = malloc((strlen(newick) + 10) * sizeof(char));
		if(buf == NULL) return;
		strcpy(buf, newick);
		}
	temp_top = NULL; basic_tree_build(1, buf, NULL, TRUE);
	tr = temp_top; temp_top = NULL;
	if(tr != NULL) { html_view_add_tree(f, tr, name, 0, 0, first); dismantle_tree(tr); }
	free(buf);
	}
void html_view_close(FILE *f, const char *filename)
	{
	if(f == NULL) return;
	fputs("]}", f);
	fputs(VIEWER_HTML_TAIL, f);
	fclose(f);
	printf2("Interactive HTML view written to: %s\n", filename);
	}

/* -----------------------------------------------------------------------
 * NHX (New Hampshire eXtended) output for reconciled gene trees
 *
 * Conventions:
 *   Real gene-copy leaf:  GeneCopyName[&&NHX:S=Species:D=N]
 *   Loss leaf:            Species*LOST[&&NHX:S=Species]
 *   Duplication node:     (...)[&&NHX:D=Y]   (+ :S= when mapped to a leaf species)
 *   Speciation node:      (...)[&&NHX:D=N]   (+ :S= when mapped to a leaf species)
 *
 * Readable natively by Archaeopteryx (http://www.phylosoft.org/archaeopteryx/)
 * which renders D=Y nodes as filled squares and D=N as circles.
 * Also readable by ETE3 (Python) and ggtree (R/Bioconductor).
 * ----------------------------------------------------------------------- */

static void nhx_sibling_chain(struct taxon *pos, char *buf)
	{
	int first = 1;
	while (pos != NULL)
		{
		if (!first) strcat(buf, ",");
		first = 0;

		if (pos->daughter != NULL)
			{
			/* Internal node */
			strcat(buf, "(");
			nhx_sibling_chain(pos->daughter, buf);
			strcat(buf, ")");
			/* Build NHX annotation */
			char nhx[512];
			if (pos->loss == 2)
				snprintf(nhx, sizeof(nhx), "[&&NHX:D=Y");
			else
				snprintf(nhx, sizeof(nhx), "[&&NHX:D=N");
			/* Add S= only when this node maps to a named leaf species */
			if (pos->tag >= 0 && pos->tag < number_of_taxa)
				{
				size_t used = strlen(nhx);
				snprintf(nhx + used, sizeof(nhx) - used, ":S=%s", taxa_names[pos->tag]);
				}
			strncat(nhx, "]", sizeof(nhx) - strlen(nhx) - 1);
			strcat(buf, nhx);
			}
		else
			{
			/* Leaf node */
			char nhx[512];
			nhx[0] = '\0';
			if (pos->loss == -1)
				{
				/* Loss leaf: Species*LOST[&&NHX:S=Species] */
				const char *sname = NULL;
				if (pos->name >= 0 && pos->name < number_of_taxa)
					sname = taxa_names[pos->name];
				else if (pos->tag >= 0 && pos->tag < number_of_taxa)
					sname = taxa_names[pos->tag];
				if (sname != NULL)
					{
					strcat(buf, sname);
					strcat(buf, "*LOST");
					snprintf(nhx, sizeof(nhx), "[&&NHX:S=%s]", sname);
					}
				}
			else
				{
				/* Real gene-copy leaf */
				if (pos->fullname != NULL)
					strcat(buf, pos->fullname);
				else if (pos->name >= 0 && pos->name < number_of_taxa)
					strcat(buf, taxa_names[pos->name]);
				if (pos->name >= 0 && pos->name < number_of_taxa)
					snprintf(nhx, sizeof(nhx), "[&&NHX:S=%s:D=N]", taxa_names[pos->name]);
				}
			strcat(buf, nhx);
			}
		pos = pos->next_sibling;
		}
	}

/* Entry point — called with tree_top (the extra wrapper node), mirrors print_fullnamed_tree */
void print_nhx_tree(struct taxon *position, char *buf)
	{
	int first = 1;
	strcat(buf, "(");
	while (position != NULL)
		{
		if (!first) strcat(buf, ",");
		first = 0;
		if (position->daughter != NULL)
			nhx_sibling_chain(position->daughter, buf);
		else
			{
			/* Bare leaf at top level (unusual but handle gracefully) */
			if (position->fullname != NULL)
				strcat(buf, position->fullname);
			else if (position->name >= 0 && position->name < number_of_taxa)
				strcat(buf, taxa_names[position->name]);
			}
		position = position->next_sibling;
		}
	strcat(buf, ")");
	}

void reconstruct(int print_settings)  /* Carry out gene-tree reconciliation of source trees against a species tree */
	{
	struct taxon *position = NULL, *species_tree = NULL, *gene_tree = NULL, *best_mapping = NULL, *unknown_fund = NULL, *pos = NULL, *copy = NULL, *newbie = NULL;
	int i, j, k, l, xnum=0, *presence = NULL, **label_results = NULL, num_species_internal = 0, error = FALSE, printfiles = FALSE, how_many = 0, diff_overall =0, dorecon = FALSE, basescore = 1, taxaorder=0, species_source = 0, gene_tree_start = 0;
	float *overall_placements = NULL, biggest = -1, total, best_total = -1;
	char *temptree;
	char *temptree1 = malloc(TREE_LENGTH * sizeof(char));
	char *temptree2 = malloc(TREE_LENGTH * sizeof(char));
	char reconfilename[100], otherfilename[100], *tmp1 = NULL, c = '\0', speciestree_file[1000];
	char nhxfilename[1000];
	char htmlfilename[1000];
	FILE *nhxfile = NULL, *htmlfile = NULL;
	int htmlfirst = 1;
	if(!temptree1 || !temptree2) { free(temptree1); free(temptree2); printf2("Error: out of memory in reconstruct\n"); return; }
	FILE *reconstructionfile = NULL, *descendentsfile = NULL, *genebirthfile = NULL;

	temptree = malloc(TREE_LENGTH*sizeof(char));
	temptree[0] = '\0';
	speciestree_file[0] = '\0';
	temptree1[0] = '\0';
	reconfilename[0] ='\0';
	otherfilename[0] = '\0';
	nhxfilename[0] = '\0';
	htmlfilename[0] = '\0';
	
	tmp1 = malloc(TREE_LENGTH*sizeof(char));
	tmp1[0] = '\0';
	count_now = TRUE;

    for(i=0; i<num_commands; i++)
        {
        if(strcmp(parsed_command[i], "printfiles") == 0)
			{
			if(strcmp(parsed_command[i+1], "yes") == 0)
				printfiles = TRUE;
			}
		if(strcmp(parsed_command[i], "speciestree") == 0)
			{
			/* values: "memory", "first", or a filename */
			if(strcmp(parsed_command[i+1], "memory") == 0)
				species_source = 1;
			else if(strcmp(parsed_command[i+1], "first") == 0)
				species_source = 2;
			else
				{
				species_source = 3;
				strcpy(speciestree_file, parsed_command[i+1]);
				}
			}
		if(strcmp(parsed_command[i], "showrecon") == 0)
			{
			if(strcmp(parsed_command[i+1], "yes") == 0)
				dorecon = TRUE;
			}
		if(strcmp(parsed_command[i], "basescore") == 0)
			{
			basescore = atoi(parsed_command[i+1]);
			}
		if(strcmp(parsed_command[i], "nhxfile") == 0)
			{
			if(strcmp(parsed_command[i+1], "yes") == 0)
				{
				/* Default filename: inputfilename.nhx */
				strncpy(nhxfilename, inputfilename, sizeof(nhxfilename)-10);
				strncat(nhxfilename, ".nhx", sizeof(nhxfilename)-strlen(nhxfilename)-1);
				}
			else
				strncpy(nhxfilename, parsed_command[i+1], sizeof(nhxfilename)-1);
			}
		if(strcmp(parsed_command[i], "htmlview") == 0)
			{
			if(strcmp(parsed_command[i+1], "yes") == 0)
				{ strncpy(htmlfilename, inputfilename, sizeof(htmlfilename)-12); strncat(htmlfilename, ".recon.html", sizeof(htmlfilename)-strlen(htmlfilename)-1); }
			else
				strncpy(htmlfilename, parsed_command[i+1], sizeof(htmlfilename)-1);
			}
		}



	/* If autoprunemono was active, swap original (unpruned) trees in for reconstruct */
	int autoprunemono_swapped = FALSE;
	if(autoprunemono_active && original_fundamentals != NULL)
		{
		for(i = 0; i < Total_fund_trees; i++)
			{
			if(original_fundamentals[i] != NULL)
				{
				char *tmp_swap = fundamentals[i];
				fundamentals[i] = original_fundamentals[i];
				original_fundamentals[i] = tmp_swap;
				}
			}
		autoprunemono_swapped = TRUE;
		printf2("(Using original unpruned trees for reconstruct)\n");
		}

	/* If autodecompose is active, temporarily reload the pristine
	 * pre-decomposition dataset (NOTES_gene_tree_decomposition.md sec 6.3)
	 * in place of the live (decomposed) pool for the duration of this call,
	 * then reload the decomposed fragments back afterward. Unlike
	 * autoprunemono's element-by-element swap above, the tree count itself
	 * differs (one family -> k fragments) and per-tree bookkeeping arrays
	 * like fulltaxanames[] are indexed/sized for whichever dataset is
	 * CURRENTLY loaded -- a raw pointer-swap of fundamentals[]/tree_names[]/
	 * Total_fund_trees (tried first, see the comment on
	 * decompose_use_pristine_for_reconstruct() in main.c) segfaulted for
	 * exactly that reason. Going through a real execute_command() reload
	 * (like sec 6.2 already does for the initial decomposition) sidesteps
	 * it: every derived structure is freshly, correctly rebuilt for
	 * whichever dataset is loaded at the time. Without this, reconstruct
	 * would reconcile already-decomposed single-copy fragments against the
	 * species tree and silently report near-zero duplications/losses --
	 * wrong, not just imprecise. */
	int decompose_swapped = FALSE;
	int saved_trees_in_memory = trees_in_memory;
	if(decompose_active)
		{
		decompose_swapped = decompose_use_pristine_for_reconstruct();
		if(decompose_swapped)
			{
			/* execute_command()'s per-load teardown unconditionally zeroes
			 * trees_in_memory (it does NOT free/reallocate retained_supers[]
			 * itself, so the candidate species tree text a prior 'nj'/'hs'
			 * left there is still intact) -- restore the counter so
			 * species_source resolution just below still sees whatever
			 * supertree the user already had in memory before this reload. */
			trees_in_memory = saved_trees_in_memory;
			printf2("(Using original pre-decomposition gene trees for reconstruct)\n");
			}
		else
			{
			printf2("Error: could not switch to the pristine pre-decomposition gene trees for reconstruct -- aborting to avoid reconciling against already-decomposed fragments\n");
			free(temptree1); free(temptree2); free(tmp1); free(temptree);
			count_now = FALSE;
			return;
			}
		}

	if(printfiles)
		{
		strcpy(reconfilename, inputfilename);
		strcat(reconfilename, ".recon");
		reconstructionfile =fopen(reconfilename, "w");
		strcat(reconfilename, ".dist");
		distributionreconfile = fopen(reconfilename, "w");
		
		strcpy(otherfilename, inputfilename);
		strcat(otherfilename, ".species.decendents");
		descendentsfile = fopen(otherfilename, "w");

		strcpy(otherfilename, inputfilename);
		strcat(otherfilename, ".gene.births");
		genebirthfile = fopen(otherfilename, "w");	
		
		strcpy(otherfilename, inputfilename);
		strcat(otherfilename, ".onetoone");
		onetoonefile =  fopen(otherfilename, "w");

		strcpy(otherfilename, inputfilename);
		strcat(otherfilename, ".strict.onetoone");
		strictonetoonefile =  fopen(otherfilename, "w");

		
		fprintf(distributionreconfile, "Duplications\tLosses\n");
		}
	presence = malloc(2*number_of_taxa*sizeof(int));
	overall_placements = malloc(2*number_of_taxa*sizeof(float));
	for(i=0; i<2*number_of_taxa; i++)
		overall_placements[i] = 0;
	
	label_results = malloc(6*sizeof(int*));
	for(i=0; i<6; i++)
		{
		label_results[i] = malloc((2*number_of_taxa)*sizeof(int));
		for(j=0; j<(2*number_of_taxa); j++)
			{
			label_results[i][j] = 0;
			}
		}
	
	 for(i=0; i<num_commands; i++)
		{

		if(strcmp(parsed_command[i], "dups") == 0)
			{
			dup_weight = tofloat(parsed_command[i+1]);
			if(dup_weight < 0)
				{
				printf2("Error: '%s' is an invalid value for dups\n", parsed_command[i+1]);
				error = TRUE;
				}
			}
		if(strcmp(parsed_command[i], "losses") == 0)
			{
			loss_weight = tofloat(parsed_command[i+1]);
			if(loss_weight < 0)
				{
				printf2("Error: '%s' is an invalid value for losses\n", parsed_command[i+1]);
				error = TRUE;
				}
			}
		if(strcmp(parsed_command[i], "lossmodel") == 0)
			{
			if(strcmp(parsed_command[i+1], "standard") == 0) loss_model = 1;
			else if(strcmp(parsed_command[i+1], "legacy") == 0) loss_model = 0;
			else { printf2("Error: lossmodel must be 'standard' or 'legacy'\n"); error = TRUE; }
			}
		}
	/* Resolve species tree source and validate */
	if(species_source == 0)
		{
		if(trees_in_memory > 0)
			species_source = 1;
		else
			{
			printf2("Error: No supertree is in memory. Run 'hs', 'nj' or 'alltrees' first,\n");
			printf2("       or specify a species tree with: reconstruct speciestree <file>\n");
			error = TRUE;
			}
		}
	if(species_source == 1 && trees_in_memory == 0)
		{
		printf2("Error: No supertree in memory for 'speciestree memory'.\n");
		error = TRUE;
		}
	if(species_source == 3 && !error)
		{
		FILE *stcheck = fopen(speciestree_file, "r");
		if(stcheck == NULL)
			{
			printf2("Error: Cannot open species tree file '%s'\n", speciestree_file);
			error = TRUE;
			}
		else fclose(stcheck);
		}
	gene_tree_start = (species_source == 2) ? 1 : 0;
	if(!error)
		{
		if(print_settings)printf2("\nReconstruction of Most likely Duplications and Losses\n\tDuplication weight = %f\n\tLosses weight = %f\n\nTree name:\tReconstruction score\n", dup_weight, loss_weight);
		/**** BUILD THE SPECIES TREE *****/
		if(species_source == 1)
			{
			/* Use supertree from memory -- already has actual taxon names */
			strcpy(temptree, retained_supers[0]);
			printf2("Using supertree in memory as species tree.\n");
			}
		else if(species_source == 3)
			{
			/* Read species tree from user-specified file */
			FILE *stfile = fopen(speciestree_file, "r");
			int si = 0; int sc;
			while((sc = getc(stfile)) != EOF && sc != ';') { temptree[si++] = (char)sc; }
			temptree[si++] = ';'; temptree[si] = '\0';
			fclose(stfile);
			printf2("Using species tree from file: %s\n", speciestree_file);
			}
		else
			{
			/* Legacy: use first source tree (fundamentals[0]) */
			strcpy(temptree, fundamentals[0]);
			returntree(temptree);
			printf2("Using first source tree as species tree (legacy mode).\n");
			}
		/* build the tree in memory */
		/****** We now need to build the Species tree in memory *******/
		temp_top = NULL;
		{ int _to = 0; tree_build(1, temptree, species_tree, TRUE, -1, &_to); }
		species_tree = temp_top;
		temp_top = NULL;
		/** add an extra node to the top of the tree */
		temp_top = make_taxon();
		temp_top->daughter = species_tree;
		species_tree->parent = temp_top;
		species_tree = temp_top;
		temp_top = NULL;
		number_tree1(species_tree, number_of_taxa);
		num_species_internal = count_internal_branches(species_tree, 0);

		/* ---------------------------------------------------------------
		 * Find the optimal rooting of the species tree by trying every
		 * possible rooting and summing the minimum dup+loss score across
		 * all gene trees.  Only the inner scoring is done here (no output,
		 * no file writing).  The actual output pass that follows uses the
		 * best-rooted species tree.
		 * --------------------------------------------------------------- */
		{
		struct taxon *spec_orig = NULL, *spec_work = NULL, *spec_rr = NULL;
		struct taxon *gt = NULL, *gt_copy = NULL, *gt_root = NULL, *gt_pos = NULL;
		int n_spec_br, sr, ll, jr, gt_n;
		int best_spec_br = -1;
		float best_spec_total = -1.0f, sr_total, best_gt, sc, sto;
		int gt_to;

		/* Duplicate the entire species_tree wrapper (parent==NULL → temp_top is set) */
		duplicate_tree(species_tree, NULL);
		spec_orig = temp_top;   /* spec_orig = wrapper copy; spec_orig->parent == NULL */
		temp_top = NULL;

		/* Count branches on the copy, not the live tree */
		n_spec_br = number_tree(spec_orig->daughter, 0);

		printf2("Searching for optimal species tree rooting (%d rootings × %d gene trees)...\n",
		        n_spec_br, Total_fund_trees - gene_tree_start);

		for(sr = 0; sr < n_spec_br; sr++)
			{
			/* Fresh copy of species tree for this rooting trial.
			 * Duplicate from wrapper (parent==NULL) so temp_top is set. */
			{
			struct taxon *sw_wrap = NULL;
			duplicate_tree(spec_orig, NULL);
			sw_wrap = temp_top;
			temp_top = NULL;
			spec_work = sw_wrap->daughter;
			/* Detach from wrapper so reroot_tree parent-walk stops here */
			spec_work->parent = NULL;
			sw_wrap->daughter  = NULL;
			free(sw_wrap);
			}

			/* Index branches, get branch sr, reroot */
			number_tree(spec_work, 0);
			{
			struct taxon *sr_br = get_branch(spec_work, sr);
			temp_top = spec_work;
			reroot_tree(sr_br);
			}
			/* use the rerooted species tree directly (no extra root wrapper), as
			 * get_recon_score()/hs does */
			spec_rr = temp_top;
			temp_top = NULL;
			/* Tag nodes for reconciliation */
			number_tree1(spec_rr, number_of_taxa);

			sr_total = 0.0f;
			for(ll = gene_tree_start; ll < Total_fund_trees; ll++)
				{
				strcpy(temptree, fundamentals[ll]);
				unroottree(temptree);
				temp_top = NULL;
				gt_to = 0;
				tree_build(1, temptree, NULL, FALSE, -1, &gt_to);
				gt = temp_top;
				temp_top = NULL;
				/* Resolve gene-tree polytomies exactly as get_recon_score()/hs and the
				 * illustration pass (Pass B) below do, so the species-rooting selection
				 * and the reported total match the hs score (soft min-dup resolution). */
				if(presence_of_trichotomies(gt))
					gt = do_resolve_tricotomies(gt, spec_rr, basescore);
				tree_top = gt;
				compress_tree1(gt);
				gt = tree_top;
				tree_top = NULL;

				gt_n = number_tree(gt, 0);
				best_gt = -1.0f;

				if(gt_n > 2)
					{
					duplicate_tree(gt, NULL);
					gt_copy = temp_top;
					temp_top = NULL;

					for(jr = 0; jr < gt_n; jr++)
						{
						gt_pos = get_branch(gt, jr);
						temp_top = gt;
						reroot_tree(gt_pos);
						gt = temp_top;
						temp_top = NULL;

						/* score the rerooted gene tree directly (no extra root wrapper),
						 * matching get_recon_score()/hs; wrapping adds a spurious top-level
						 * duplication */
						sc = tree_map(gt, spec_rr, FALSE);
						if(sc < best_gt || best_gt < 0.0f) best_gt = sc;

						dismantle_tree(gt);
						gt = NULL;

						duplicate_tree(gt_copy, NULL);
						gt = temp_top;
						temp_top = NULL;
						number_tree(gt, 0);
						}
					dismantle_tree(gt_copy);
					gt_copy = NULL;
					/* gt was last restored from copy — dismantle it */
					dismantle_tree(gt);
					gt = NULL;
					}
				else
					{
					gt_root = make_taxon();
					gt_root->daughter = gt;
					if(gt) gt->parent = gt_root;
					sc = tree_map(gt_root, spec_rr, FALSE);
					if(sc < best_gt || best_gt < 0.0f) best_gt = sc;
					dismantle_tree(gt_root);
					gt_root = NULL;
					gt = NULL;
					}

				if(best_gt >= 0.0f) sr_total += best_gt;
				}  /* end gene tree loop */

			if(sr_total < best_spec_total || best_spec_total < 0.0f)
				{
				best_spec_total = sr_total;
				best_spec_br = sr;
				}

			dismantle_tree(spec_rr);
			spec_rr = NULL;
			}  /* end species tree rooting loop */

		/* Apply the best rooting to the live species_tree.
		 * Duplicate from wrapper (parent==NULL) then detach daughter before rerooting. */
		{
		struct taxon *sw_wrap = NULL;
		duplicate_tree(spec_orig, NULL);
		sw_wrap = temp_top;
		temp_top = NULL;
		spec_work = sw_wrap->daughter;
		spec_work->parent = NULL;
		sw_wrap->daughter  = NULL;
		free(sw_wrap);
		}
		number_tree(spec_work, 0);
		{
		struct taxon *best_br = get_branch(spec_work, best_spec_br);
		temp_top = spec_work;
		reroot_tree(best_br);
		}
		/* Replace the wrapper's daughter with the best-rooted tree */
		dismantle_tree(species_tree->daughter);
		species_tree->daughter = temp_top;
		if(temp_top) temp_top->parent = species_tree;
		temp_top = NULL;
		number_tree1(species_tree, number_of_taxa);
		num_species_internal = count_internal_branches(species_tree, 0);

		dismantle_tree(spec_orig);
		spec_orig = NULL;

		printf2("Optimal species tree rooting found (total dup+loss score: %.4f)\n", best_spec_total);
		}
		/* ---------------------------------------------------------------
		 * End optimal rooting search
		 * --------------------------------------------------------------- */

		if(printfiles)
			{
			fprintf(reconstructionfile, "Tree name\t");
			for(j=0; j<number_of_taxa+num_species_internal; j++)
				{
				check_tree(species_tree, j, reconstructionfile);
				}
			fprintf(reconstructionfile, "\n");
			
			fprintf(genebirthfile, "tree name\tgene birth\n");
			/* print out the descendents of each internal branch of the tree */
			do_descendents(species_tree, descendentsfile);
			fclose(descendentsfile);
			}
		
		
		/* Open NHX output file if requested */
		if(nhxfilename[0] != '\0')
			{
			nhxfile = fopen(nhxfilename, "w");
			if(nhxfile == NULL)
				{
				printf2("Warning: could not open NHX output file '%s' — NHX output disabled.\n", nhxfilename);
				nhxfilename[0] = '\0';
				}
			else
				printf2("NHX reconciliation output will be written to: %s\n", nhxfilename);
			}

		/* Open the single combined interactive HTML viewer if requested */
		if(htmlfilename[0] != '\0')
			{
			char meta[700];
			snprintf(meta, sizeof(meta),
			         "{\"dataset\":\"%s\",\"criterion\":\"recon\",\"lossmodel\":\"%s\"}",
			         inputfilename, (loss_model==1 ? "standard" : "legacy"));
			htmlfile = html_view_open(htmlfilename, meta, 1);
			htmlfirst = 1;
			}

		for(l=gene_tree_start; l<Total_fund_trees; l++)
			{
			temp_top = NULL;
			strcpy(temptree, "");
			strcpy(temptree, fundamentals[l]);
			unroottree(temptree);
			
			/*returntree(temptree); */ /* add in the actual names in the tree */
			/* build the tree in memory */
			/****** We now need to build the gene tree in memory *******/
			temp_top = NULL;
			how_many++;
			taxaorder=0;
			tree_build(1, temptree, gene_tree, FALSE, l, &taxaorder);
			gene_tree = temp_top;

			temp_top = NULL;
			strcpy(temptree1, temptree);
			i = count_taxa(gene_tree, 0);
			if(presence_of_trichotomies(gene_tree))
				{
				printf("Doing resolving of trichotomies\n");
				gene_tree = do_resolve_tricotomies(gene_tree, species_tree, basescore);
				}
			tree_top = gene_tree;
			compress_tree1(gene_tree);
			gene_tree = tree_top;
			temp_top = NULL;
			duplicate_tree(gene_tree, NULL);
			how_many++;
			copy = temp_top;
			tree_top = NULL;
			i = number_tree(gene_tree, 0); /* number every branch of the gene tree so we can use this to reroot */
			best_total = -1;

			if(i>2) /* if there are more than two branches in this gene tree */
				{
				for(j=0; j<i; j++) /* For every possible rerooting of the gene tree */
					{

					malloc_check=0;
					position = get_branch(gene_tree, j); /* find branch numbered j on the gene tree */ 
					temp_top = gene_tree;
					reroot_tree(position); /* reroot at this position */

					/* use the rerooted tree directly (no extra root wrapper), matching
						 * get_recon_score()/hs so the illustrated counts equal the hs score */
						gene_tree = temp_top;
						temp_top = NULL;

						diff_overall -= malloc_check;
						total = tree_map(gene_tree, species_tree->daughter, printfiles);
						diff_overall += malloc_check;
					if(total < best_total || best_total == -1)
						{
						best_total = total;
						if(best_mapping != NULL)
							{
							dismantle_tree(best_mapping);
							how_many--;
							best_mapping = NULL;
							}
						best_mapping = gene_tree;
						gene_tree = NULL;
						}
					else
						{
						dismantle_tree(gene_tree);
						how_many--;
						gene_tree = NULL;
						}
					temp_top = NULL;
					tree_top = NULL;
					duplicate_tree(copy, NULL);
					how_many++;
					gene_tree = temp_top;
					temp_top = NULL;
					tree_top = NULL;
					number_tree(gene_tree, 0);
					another_check += malloc_check;
					}
				}
			else
				{
				total = tree_map(gene_tree, species_tree->daughter, printfiles);				
				best_total = total;
				best_mapping = gene_tree;
				gene_tree = NULL;
				}
			if(strcmp(tree_names[l], "") == 0)
				printf2("Tree number: %d\t%f\n", l, best_total);
			else
				printf2("%s\t%f\n", tree_names[l], best_total);
			tree_top = best_mapping;
			/* lossmodel=standard: best_mapping is unannotated (the arithmetic scorer
			 * does no reconstruction), so build the standard-model reconciled tree for
			 * the illustration -- its dup marks and loss lineages match the standard
			 * score exactly (validated). Legacy mode leaves best_mapping as-is. */
			if(loss_model == 1)
				tree_top = annotate_standard(best_mapping, species_tree->daughter, 0);
			/** ADD IN PRINTFULLNAMED TREE HERE **/
			temptree[0] = '\0';
			print_fullnamed_tree(tree_top, temptree, l);
			/*print_named_tree(tree_top, temptree); */

			/* NHX output: write reconciled tree with duplication/loss annotations */
			if(nhxfile != NULL)
				{
				temptree1[0] = '\0';
				print_nhx_tree(tree_top, temptree1);
				/* Write tree name as a comment, then the NHX Newick on the next line */
				if(strcmp(tree_names[l], "") != 0)
					fprintf(nhxfile, "# %s  (score=%.4f)\n", tree_names[l], best_total);
				else
					fprintf(nhxfile, "# tree_%d  (score=%.4f)\n", l, best_total);
				fprintf(nhxfile, "%s;\n\n", temptree1);
				}

			/* Interactive HTML view: append this gene tree to the single combined viewer
			 * file. tree_top is the reconciled/annotated tree in both modes, so the drawn
			 * duplications and losses match the reported per-tree score. */
			if(htmlfile != NULL)
				{
				char tname[200];
				if(tree_names[l][0]) snprintf(tname, sizeof(tname), "%s", tree_names[l]);
				else snprintf(tname, sizeof(tname), "tree_%d", l);
				html_view_add_tree(htmlfile, tree_top, tname, best_total, 1, htmlfirst);
				htmlfirst = 0;
				}

			if(dorecon == TRUE)
				{
				tree_coordinates(temptree, TRUE, FALSE, FALSE, l); /* char *tree, int bootstrap, int build, int mapping, int fundnum */
				}
			/******* PRINT_TREE_LABELS TEST *******/	
			/****	label_results[0] = number of copies of the gene at this internal branch
					label_results[1] = number of duplications at this internal branch
					label_results[2] = number of losses at this internal branch
					label_results[3] = number of 1:1 orthologs after this internal branch (allowing duplications and losses in external taxa)
					label_results[4] = number of strict 1:1 orthologs after this internal branch (not allowing any duplications and losses)
					label_results[5] = (added Aug 17) Number of Duplications at this internal branch supported by splits in the orignal gene tree (Does not include those inferred from polytomies)
			****/

			if(printfiles)
				{	
				for(i=0; i<6; i++)
					{
					for(j=0; j<(2*number_of_taxa); j++)
						{
						label_results[i][j] = 0;
						}
					}
				
				if(tree_top->tag == number_of_taxa+num_species_internal-1)
					fprintf(genebirthfile, "%s\troot\n", tree_names[l]);
				else
					{
					fprintf(genebirthfile, "%s\t", tree_names[l]);
					check_tree(species_tree, tree_top->tag, genebirthfile);
					fprintf(genebirthfile, "\n");
					}
					
				print_tree_labels(tree_top, label_results, l, species_tree);
				fprintf(reconstructionfile, "%10s\t", tree_names[l]);
				for(j=0; j<(2*number_of_taxa)-1; j++)
					{
					
					fprintf(reconstructionfile, "%5d (%5d+ %5d+S %5d- %5d 1:1 %5d 1:1S)\t", label_results[0][j], label_results[1][j], label_results[5][j], label_results[2][j], label_results[3][j], label_results[4][j]);
					}
				fprintf(reconstructionfile, "\n");
				}
			/******* END TEST ********/
			for(i=0; i<2*number_of_taxa; i++)
				presence[i] = FALSE;
				
			find_tagged(best_mapping, presence);
			j=0;
			for(i=0; i<2*number_of_taxa; i++)
				if(presence[i]) j++;
			for(i=0; i<2*number_of_taxa; i++)
				{
				if(presence[i]) overall_placements[i] += (float)1/(float)j;
				}
			if(copy != NULL)
				{
				dismantle_tree(copy);
				copy = NULL;
				}
				how_many--;
			
			if(gene_tree != NULL)
				{
				dismantle_tree(gene_tree);
				gene_tree = NULL;
				}
			how_many--;
			/* lossmodel=standard: dismantle the separate annotated illustration tree */
			if(loss_model == 1 && tree_top != NULL && tree_top != best_mapping)
				{ dismantle_tree(tree_top); tree_top = NULL; }
			if(best_mapping != NULL)
				{
				dismantle_tree(best_mapping);
				best_mapping = NULL;
				}
				
				
			
			how_many--;
			how_many--;
			}
		
		if(species_tree != NULL)
			{
			dismantle_tree(species_tree);
			species_tree = NULL;
			}
		}
	if(printfiles)
		{
		fclose(distributionreconfile);
		fclose(reconstructionfile);
		fclose(genebirthfile);
		fclose(onetoonefile);
		fclose(strictonetoonefile);
		}
	if(nhxfile != NULL)
		{
		fclose(nhxfile);
		printf2("NHX output written to: %s\n", nhxfilename);
		nhxfile = NULL;
		}
	if(htmlfile != NULL) { html_view_close(htmlfile, htmlfilename); htmlfile = NULL; }
	tree_top = NULL;
	free(temptree);
	free(temptree1);
	free(temptree2);
	free(tmp1);
	free(presence);
	free(overall_placements);
	if(label_results != NULL)
		{
		for(i=0; i<6; i++)
			{
			free(label_results[i]);
			label_results[i] = NULL;
			}
		free(label_results);
		label_results = NULL;
		}

	/* Reload the decomposed dataset back in (undoes the pristine reload
	 * above) so other commands (hs, etc.) continue to see the decomposed
	 * fragments, before the autoprunemono per-index swap-back below (which
	 * relies on Total_fund_trees matching the live, decomposed count). */
	if(decompose_swapped)
		{
		decompose_restore_after_reconstruct();
		trees_in_memory = saved_trees_in_memory;   /* see comment at the swap-in site above */
		}

	/* Swap original trees back so other commands (hs, etc.) continue to use pruned versions */
	if(autoprunemono_swapped)
		{
		for(i = 0; i < Total_fund_trees; i++)
			{
			if(original_fundamentals[i] != NULL)
				{
				char *tmp_swap = fundamentals[i];
				fundamentals[i] = original_fundamentals[i];
				original_fundamentals[i] = tmp_swap;
				}
			}
		}

	count_now=FALSE;

	}

void put_in_scores(struct taxon * position, float * total)
	{
	while(position != NULL)
		{
		if(position->daughter != NULL) put_in_scores(position->daughter, total);
		if(total[position->tag] == 0) position->loss = -1;
		else position->loss = total[position->tag];
		position = position->next_sibling;
		}
	}

void assign_hgtdonors(struct taxon * position, int num, int part_num)
	{
	int i;
	while(position != NULL)
		{
		if(position->daughter != NULL) assign_hgtdonors(position->daughter, num, part_num);
		if(position->donor == NULL)
			{
			position->donor = malloc(num*sizeof(int));
			for(i=0; i<num; i++) position->donor[i] = FALSE;
			}
		if(position->tag2 == TRUE)
			{
			position->donor[part_num] = TRUE;
			}
		position = position->next_sibling;
		}
	}

void hgt_reconstruction()
	{
	struct taxon *position = NULL, *species_tree = NULL, *gene_tree = NULL, *best_mapping = NULL, *best_mapping1 = NULL, *best_mapping2 = NULL, *unknown_fund = NULL, *posit = NULL,*copy = NULL, *copy1 = NULL, **parts = NULL, *test_part = NULL, *pos = NULL, *best_donor = NULL, *best_HGT = NULL, *attached = NULL;
	int i, j, k, l,  *presence = NULL,*presence1 = NULL, *presence2 = NULL, hgt_receipient1, hgt_receipient2, **overall_presence = NULL, *overall_reconstruction = NULL, *overall_receptor = NULL, receptor, **tmp_presence1 = NULL, **tmp_presence2 = NULL, *before1 = NULL, *before2 = NULL, *after1 = NULL, *after2 = NULL, *temporary = NULL;
	float *overall_placements = NULL, biggest = -1,  total, best_total = -1,  best_total1 = -1, best_total2 = -1, HGT1 = 0, HGT2 = 0, original = 0, best_HGT_recon = -1, best_reconstruction = -1, sum, HGT_score = -1, donor_score = -1, tmp_allow = FALSE;
	char *temptree = NULL;
	char *temptree1 = malloc(TREE_LENGTH * sizeof(char));
	int **species_allowed = NULL, **dependent_species_allowed = NULL, *previous = NULL, xnum =0, x, y, z, partA, partB, q, r, s, allow_HGT1 = TRUE, allow_HGT2 = TRUE, numparts = 1,  place_marker = 1, found_better = FALSE, error = FALSE, taxaorder=0;
	int basescore = 1; /** see reconstrution command **/

	if(!temptree1) { printf2("Error: out of memory in recon_hgt\n"); return; }
	temptree = malloc(TREE_LENGTH*sizeof(char));
	temptree[0] = '\0';
	temptree1[0] = '\0';
	
	for(i=0; i<num_commands; i++)
		{
		if(strcmp(parsed_command[i], "dupweight") == 0)
			{
			dup_weight = tofloat(parsed_command[i+1]);
			}
		if(strcmp(parsed_command[i], "lossweight") == 0)
			{
			loss_weight = tofloat(parsed_command[i+1]);
			}
		if(strcmp(parsed_command[i], "hgtweight") == 0)
			{
			hgt_weight = tofloat(parsed_command[i+1]);
			}
		
		}
	
	printf2("Calculating best reconstruction of Duplications, losses and Horizontal gene transfers (HGT)\n\nDuplication weight set to %f\nLoss weight set to %f\nHGT weight set to %f\n\n", dup_weight, loss_weight, hgt_weight);
	
	tree_top = NULL;
	
	/**** First, take the species tree and calculate those part of the tree that would constitute travelling in time if a HGT was to occur (ie in ancestral or descendent branches) ***/
	species_allowed = malloc((2*number_of_taxa)* sizeof(int*));
	if(species_allowed == NULL) memory_error(120);
	for(i=0; i<(2*number_of_taxa); i++)
		{
		species_allowed[i] = malloc((2*number_of_taxa)*sizeof(int));
		if(species_allowed[i] == NULL) memory_error(121);
		for(j=0; j<(2*number_of_taxa); j++) species_allowed[i][j] = TRUE;
		}
	
	dependent_species_allowed = malloc((2*number_of_taxa)* sizeof(int*));
	if(!dependent_species_allowed) memory_error(122);
	for(i=0; i<(2*number_of_taxa); i++)
		{
		dependent_species_allowed[i] = malloc((2*number_of_taxa)*sizeof(int));
		if(dependent_species_allowed[i] == NULL) memory_error(123);
		for(j=0; j<(2*number_of_taxa); j++) dependent_species_allowed[i][j] = TRUE;
		}
	
	temporary = malloc(2*number_of_taxa*sizeof(int));
	presence = malloc(2*number_of_taxa*sizeof(int));
	presence1 = malloc(2*number_of_taxa*sizeof(int));
	presence2 = malloc(2*number_of_taxa*sizeof(int));
	
	previous = malloc((2*number_of_taxa)*sizeof(int));
	for(i=0; i<2*number_of_taxa; i++) previous[i] = FALSE;
	
	before1 = malloc((2*number_of_taxa)*sizeof(int));
	before2 = malloc((2*number_of_taxa)*sizeof(int));
	after1 = malloc((2*number_of_taxa)*sizeof(int));
	after2 = malloc((2*number_of_taxa)*sizeof(int));
	for(i=0; i<2*number_of_taxa; i++)
		{
		before1[i] = before2[i] = after1[i] = after2[i] = FALSE;
		}
	
	/**** BUILD THE SPECIES TREE *****/
	strcpy(temptree, fundamentals[0]);
	returntree(temptree);
	/* build the tree in memory */
	/****** We now need to build the Species tree in memory *******/
	temp_top = NULL;
	{ int _to = 0; tree_build(1, temptree, species_tree, 1, -1, &_to); }
	species_tree = temp_top;
	temp_top = NULL;
	/** add an extra node to the top of the tree */
	temp_top = make_taxon();
	temp_top->daughter = species_tree;
	species_tree->parent = temp_top;
	species_tree = temp_top;
	temp_top = NULL;
			printf2("2\n");

	xnum = number_tree1(species_tree, number_of_taxa); /* label the internal and external branches of the species tree */
	assign_ances_desc(species_tree, species_allowed, previous);
	free(previous);
	
	
	for(l=1; l<Total_fund_trees; l++)
		{
		best_reconstruction = -1;
		best_donor = NULL;
		best_HGT = NULL;
		temp_top = NULL;
		strcpy(temptree, "");
		strcpy(temptree, fundamentals[l]);
		returntree(temptree);
		unroottree(temptree);  /* we need the gene tree to be unrooted when we are breaking it up */
		/* build the tree in memory */
		/****** We now need to build the genetree in memory *******/
		temp_top = NULL;
		taxaorder=0;
		tree_build(1, temptree, gene_tree, 1, l, &taxaorder);
		gene_tree = temp_top;
		temp_top = NULL;
			
		strcpy(temptree1, temptree);
	
		/***** NOW THAT WE HAVE BUILT THE GENE TREE, WE NEED TO BREAK IT AT EVERY BRANCH IN TURN ****/

		num_gene_nodes = number_tree(gene_tree, 0);
		
		overall_receptor = malloc(num_gene_nodes*sizeof(int));
		for(i=0; i<num_gene_nodes; i++) overall_receptor[i] = -1;
		
		parts = malloc(num_gene_nodes*sizeof(struct taxon *));
		for(i=0; i<num_gene_nodes; i++) parts[i] = NULL;
		
		overall_placements = malloc(num_gene_nodes*sizeof(float));
		for(i=0; i<num_gene_nodes; i++)
			overall_placements[i] = 0;
		
		tmp_presence1 = malloc(num_gene_nodes*sizeof(int *));
		for(i=0; i<num_gene_nodes; i++) 
			{
			tmp_presence1[i] = malloc((2*number_of_taxa)*sizeof(int *));
			for(j=0; j<2*number_of_taxa; j++) tmp_presence1[i][j] = FALSE;
			}
		
		tmp_presence2 = malloc(num_gene_nodes*sizeof(int *));
		for(i=0; i<num_gene_nodes; i++) 
			{
			tmp_presence2[i] = malloc((2*number_of_taxa)*sizeof(int *));
			for(j=0; j<2*number_of_taxa; j++) tmp_presence2[i][j] = FALSE;
			}
		
		overall_presence = malloc(num_gene_nodes*sizeof(int *));
		for(i=0; i<num_gene_nodes; i++) 
			{
			overall_presence[i] = malloc((2*number_of_taxa)*sizeof(int *));
			for(j=0; j<2*number_of_taxa; j++) overall_presence[i][j] = FALSE;
			}

		if(presence_of_trichotomies(gene_tree)) gene_tree = do_resolve_tricotomies(gene_tree, species_tree, basescore);

		parts[0] = gene_tree;
		gene_tree = NULL;
		strcpy(temptree, "");
		print_tree(parts[0], temptree);
		
		duplicate_tree(parts[0], NULL);
		copy = tree_top;			
		i = number_tree(copy, 0);
		best_total = -1;
		for(j=0; j<i; j++) /* for every node (internal and external) on the tree */
			{
			position = get_branch(copy, j);
			tree_top = copy;
			reroot_tree(position);
		
			copy = tree_top;
			tree_top = NULL;
			printf2("e\n");
			total = tree_map(copy, species_tree,0);
			if(total < best_total || best_total == -1) best_total = total;
			dismantle_tree(copy);
			copy = NULL;
			temp_top = NULL;
			tree_top = NULL;
			duplicate_tree(parts[0], NULL);
			copy = tree_top;
			temp_top = NULL;
			tree_top = NULL;
			number_tree(copy, 0);
			}  /* WE now know the cost for the best reconstruction of dups and losses on this tree WITHOUT any HGTs */
		dismantle_tree(copy);
		copy = NULL;
		original = best_total;
		overall_placements[0] = original;
		numparts = 1;
		for(k=0; k<numparts; k++)
			{
			printf2("ON part %d out of %d parts\n", k, numparts);
			tree_top = parts[k];
			compress_tree1(parts[k]);
			parts[k] = tree_top;
			best_reconstruction = overall_placements[k];
			found_better = FALSE;
			j=0;
			while(parts[j] != NULL) j++;
			place_marker = j;   /* THis is to tell us the next available spce in parts[], so that if there is a HGT found we know where to attach it */
			reset_tag2(parts[k]);
			/* calculate the cost for duplications and losses for the tree as it is given  (trying all rootings )*/
			if(overall_placements[k] != 0 && overall_placements[k] > 1*hgt_weight)  /* If the best reconstruction of the unbroken tree requires no dups or losses, or if the cost of 1 HGT is greater than this score, we don't test for HGTs */
				{
				i = number_tree(parts[k], 0);
				for(j=0; j<i; j++) /* go through all the parts of the tree and for each part break it into two at every point and look for new hgts */
					{
					partA = partB = TRUE;
					/* Duplicate this tree (or part of), that we are at, and use this to find the best place to break the tree ***/

					tree_top = NULL;
					duplicate_tree(parts[k], NULL);
					copy = tree_top;
					
					number_tree(copy, 0);
					position = get_branch(copy, j); /* this is the branch at which we are going to break the tree :-) */
					if(position != NULL) /* no point in seperating at the top node of the tree (which in reality doesn't exist ) */
						{
						
						test_part = position;
						if(position == copy) copy = position->next_sibling;
						/* mark the position where this was removed from */
						posit = position;
						if(position->next_sibling != NULL) attached = position->next_sibling;
						else attached = position->prev_sibling;
						while(posit->prev_sibling != NULL) posit = posit->prev_sibling;
						while(posit != NULL)
							{
							posit->tag2 = TRUE;
							posit = posit->next_sibling;
							}
						
						/* remove this node from the tree */
						if(position->parent != NULL)
							{
							(position->parent)->daughter = position->next_sibling;
							(position->next_sibling)->parent = position->parent;
							position->parent = NULL;
							}
						if(position->next_sibling != NULL) (position->next_sibling)->prev_sibling = position->prev_sibling;
						if(position->prev_sibling != NULL) (position->prev_sibling)->next_sibling = position->next_sibling;
						position->next_sibling = NULL;
						position->prev_sibling = NULL;
						/* now we need to compress the tree we just pruned to make sure that there are no unnecessary pointer taxa */
						tree_top = copy;
						compress_tree1(copy);
						copy = tree_top;
						tree_top = test_part;
						compress_tree1(test_part);
						test_part = tree_top;
						
						
						/* now reconstruct the duplications and losses of "copy" and or "test_part" seperately to see if this helped the reconsctruction */
						/* Either "copy" or "test part"  could be the HGT part. We must treat both as the HGT part in turn and pick which one has the best dup+loss reconstruction */
						/* Whichever is being treated as the HGT cannot be rerooted because if this is the HGT part then where we broke the tree is when the HGT joined the new lineage */
						/* The other half is treated as the "donor" and can be rerooted to find the best rooting for dup+loss reconstruction */
						for(q=0; q<2*number_of_taxa; q++)
							{
							presence1[q] = FALSE;
							presence2[q] = FALSE;
							}
						
						/* first calculate the dup+loss score for the two parts, this is the score for each as if both were the HGT part, "test_part" doesn't have to be rerooted, but "copy" does, at the point that "test_part" was attached */
						if(k == 0)  /* we only need to assume "copy" was the HGT if k == 0, otherwise it must be the donor */
							{
							tree_top = NULL;
							duplicate_tree(copy, NULL);
							copy1 = tree_top;
							
							tree_top = copy;
							if(copy != attached) reroot_tree(attached);
							else
								{
								posit = make_taxon();
								posit->daughter = tree_top;
								tree_top->parent = posit;
								tree_top = posit;
								posit = NULL;
								} 
							copy = tree_top;  
							printf2("f\n");
							HGT1 = tree_map(copy, species_tree,0); /* HGT1 no has the score for the tree "copy" */
							hgt_receipient1 = copy->tag;  /* if this is the HGT part, then this is the hypothesised node on the species tree that received the HGT */

							dismantle_tree(copy);
							copy = copy1;
							copy1 = NULL;
							}
						
						/* If k != 0 then we don't reroot the tree before figuring out where the HGT was attached, if k != 0 then test_part cannot be the donor! (because this would mean we would have to reroot parts[k], which is a HGT already and as such is rooted :-) */
						tree_top = NULL;
						duplicate_tree(copy, NULL);
						copy1 = tree_top;
						printf2("g\n");
						best_total1 = tree_map(copy1, species_tree,0);
						find_tagged(copy1, presence1);
						best_mapping1 = copy1; /* this will be used later for checking HGT compatibilities (( this version is only used if k != 0 ))*/
						copy1 = NULL;
						
						tree_top = NULL;
						duplicate_tree(test_part, NULL);
						copy1 = tree_top;
						printf2("h\n");
						HGT2 = tree_map(copy1, species_tree,0); /* HGT2 no has the score for the tree "test_part" */
						hgt_receipient2 = copy1->tag; /* if this is the HGT part, then this is the hypothesised node on the species tree that received the HGT */
						/*find_tagged(copy1, presence2); This is commented out because if k != 0 then test_part cannot be the donor */
						dismantle_tree(copy1);
						copy1 = NULL;
						tree_top = NULL;
						/*** Next calculate the best rerooting of each tree in terms of dups+losses ***/
						
						if(k == 0) /* we can only reroot the original part of the tree, any subsequent ones are rooted because they are the HGTs */
							{
							dismantle_tree(best_mapping1);
							best_mapping1 = NULL;
							for(z=0; z<2; z++)
								{
								/**** CHECK TO SEE IF BREAKING THE TREE AT THIS POINT MAKES THE HGT "TRAVEL IN TIME" also calculate the best dups+loss reconstruction for each tree  */
								tree_top = NULL;
								if(z == 0) duplicate_tree(copy, NULL);
								else duplicate_tree(test_part, NULL);
								copy1 = tree_top;
								tree_top = NULL;
								tree_top = copy1;
								compress_tree1(copy1);
								copy1 = tree_top;
								
								x = number_tree(copy1, 0);
								best_total = -1;
								if(x >3)
									{
									for(y=0; y<x; y++)
										{
										position = get_branch(copy1, y);
										tree_top = copy1;
										if(copy1 != position) reroot_tree(position);
										else
											{
											if(z == 0)
												{
												posit = tree_top;
												q=0;
												while(posit != NULL)
													{
													q++;
													posit = posit->next_sibling;
													}
												posit = NULL;
												if(q != 1)
													{
													reroot_tree(tree_top);
													}
												}
											} 
										copy1 = tree_top;
										tree_top = NULL;
										
										printf2("I\n");
										total = tree_map(copy1, species_tree,0);
										if(total < best_total || best_total == -1)
											{
											best_total = total;
											if(best_mapping != NULL)
												{
												dismantle_tree(best_mapping);
												best_mapping = NULL;
												}
											best_mapping = copy1;
											copy1 = NULL;
											}
										else
											{
											dismantle_tree(copy1);
											copy1 = NULL;
											}
										temp_top = NULL;
										tree_top = NULL;
										if(z == 0) duplicate_tree(copy, NULL);
										else duplicate_tree(test_part, NULL);
										copy1 = tree_top;
										temp_top = NULL;
										tree_top = copy1;
										compress_tree1(copy1);
										copy1 = tree_top;
										tree_top = NULL;
										number_tree(copy1, 0);
										}
									}
								else
									{
									printf2("j\n");
									best_total = tree_map(copy1, species_tree,0);
									best_mapping = copy1;
									copy1 = NULL;
									}
								if(z==0) best_total1 = best_total;
								else best_total2 = best_total;
								
								
								for(q=0; q<2*number_of_taxa; q++)
									{
									if(z == 0) presence1[q] = FALSE;
									else presence2[q] = FALSE;
									}

								if(z == 0) find_tagged(best_mapping, presence1);  /* so after this function "presence" will have all the possible positions that the other part of the tree (if it was a HGT) could have come from. */
								else find_tagged(best_mapping, presence2);
								
								if(z==0) best_mapping1 = best_mapping;
								else best_mapping2 = best_mapping;
								best_mapping = NULL;
								
								if(copy1 != NULL)
									{
									dismantle_tree(copy1);
									copy1 = NULL;
									}
								}
							}
						
						/**** WE NOW HAVE THE FOLLOWING VALUES:  */
						/**** HGT1 = reconstruction cost of "copy" as it if was the HGT **/
						/**** HGT2 = reconstruction cost of "test_part" as it if was the HGT **/
						/**** best_total1 = reconstruction cost of "copy" as if it was the donor part */
						/**** best_total2 = reconstruction cost of "test_part" as if it was the donor part */
						/**** presence1 = contains the possible donor nodes if "copy" was the donor */
						/**** presence2 = contains the possible donor nodes if "test_part" was the donor */
						/**** hgt_receipient1 = is the node number of the receipient if "copy" wat the HGT part */
						/**** hgt_receipient2 = is the node number of the receipient if "test_part" wat the HGT part */
						/****/
						/**** With these we need to now test 1) if hypothesising if either was a HGT causes them to travel beack in time **/
						/**** 2) Whichever passes the above test must then have a reconstruction cost that is better than the original tree without ANY HGT hypothesised */
						/**** If both pass the first test and both have the same reconstruction cost which is better than the original tree, we arbitrarily choose one over the other (not sure how else to deal with ths :-) */
						
						/* First see if the placements of either contradicts the species tree array we made earlier */ 
						/* If the donor could have been the same as the receipient, then we don't allow the HGT */
						allow_HGT1 = allow_HGT2 = TRUE;
						if(presence1[hgt_receipient2] == TRUE) allow_HGT1 = FALSE;
						if(presence2[hgt_receipient1] == TRUE) allow_HGT2 = FALSE;  /* If we are able to hypothesise that the donor and the receipient was the same then we disregard the HGT */	
						
						if(k!= 0) allow_HGT2 = FALSE;	/* WE cannot reroot this if it isnot the original part of the tree (ie k == 0) and we already know that this is a HGT part, so it can only be a donor for a hgt event, so copy HAS TO BE the donor */
																												
						/* Test to see if any of the hypothesised HGTs would have conflicted with the structure of the Species tree */
						if(allow_HGT1 == TRUE || allow_HGT2 == TRUE)
							{
							s = FALSE;
							r = FALSE;
							for(q=0; q<(2*number_of_taxa); q++)
								{
								if(presence1[q] == TRUE && species_allowed[q][hgt_receipient2] == TRUE) s = TRUE;
								if(presence2[q] == TRUE && species_allowed[q][hgt_receipient1] == TRUE) r = TRUE;
								}
							if(s == FALSE) allow_HGT1 = FALSE;
							if(r == FALSE) allow_HGT2 = FALSE;
							}
						/*** NOW SEE IF THIS CONFLICTS WITH ANY PREVIOUS HGTS!!!!!!!****/
						
						if(numparts > 1) /* no point checking this if it is the first hgt defined */
							{
							for(q=0; q<num_gene_nodes; q++)
								{
								for(r=0; r<2*number_of_taxa; r++)
									{
									tmp_presence1[q][r] = FALSE;
									tmp_presence2[q][r] = FALSE;
									}
								}
							
							/* See if either of these two possible HGTs conflict with previous HGTs */
							
							if(allow_HGT1)
								{
								/* check both copy and test_part for any previous HGTs donor sites and then check what the possible donors are now, given the new HGT found */
								/* HGT1 assumes that copy is the donor and test_part is the HGT */
								/* we will use the best rooting of copy to calculate the possible donor positions of any hgts defined earlier, these will be saved in tmp_presence1 */
								/* best_mapping1  holds the optimal mapping of Copy as the donor */
								for(q=1; q<numparts; q++)
									{
									reset_tag2(best_mapping1);
									if(assign_tag2(best_mapping1, q))
										{
										find_tagged(best_mapping1, tmp_presence1[q]);
										}
									else
										{
										for(r=0; r<2*number_of_taxa; r++)
											{
											tmp_presence1[q][r] = overall_presence[q][r];
											}
										}
									}
								
								for(s=0; s<(2*number_of_taxa); s++)  /* initialise this to include the timetravelling HGTs that occur from the shape of the tree */
									{
									for(r=0; r<(2*number_of_taxa); r++) dependent_species_allowed[s][r] = species_allowed[s][r];
									}	
								
								/* NOW check to see if the new HGT causes any conflicts with the previous HGTs */
								
								for(q=1; q<numparts; q++)
									{
									
									 /* check to see if this HGT conflicts with any previously assigned ((not needed for the first hgt)) */
									
									if(allow_HGT1)
										{
										/* identify the nodes above and below the donor and receptor sites so cross hgts can be defined as being not allowed */
										for(r=0; r<(2*number_of_taxa); r++) before1[r]  = after1[r] = before2[r] = after2[r]= FALSE;
										for(r=0; r<2*number_of_taxa; r++) previous[r] = FALSE;
										for(r=0; r<(2*number_of_taxa); r++)
											{
											if(tmp_presence1[q][r] == TRUE)
												assign_before_after(species_tree, previous, before1, after1, r, FALSE);
											}
										for(r=0; r<(2*number_of_taxa); r++)
											{
											if(tmp_presence1[q][r] == TRUE) before1[r] = after1[r] = FALSE;
											}
										for(r=0; r<2*number_of_taxa; r++) previous[r] = FALSE;
										assign_before_after(species_tree, previous, before2, after2, overall_receptor[q], FALSE);
										before2[overall_receptor[q]] = after2[overall_receptor[q]] = FALSE;
										
										/** assign the new cross-hgt rules to the array dependent_species_allowed */
										for(r=0; r<2*number_of_taxa; r++)
											{
											if(before1[r] == TRUE)
												{
												for(s=0; s<2*number_of_taxa; s++)
													{
													if(after2[s] == TRUE)
														{
														dependent_species_allowed[r][s] = FALSE;
														dependent_species_allowed[s][r] = FALSE;
														}
													}
												}
											}
										for(r=0; r<2*number_of_taxa; r++)
											{
											if(before2[r] == TRUE)
												{
												for(s=0; s<2*number_of_taxa; s++)
													{
													if(after1[s] == TRUE)
														{
														dependent_species_allowed[r][s] = FALSE;
														dependent_species_allowed[s][r] = FALSE;
														}
													}
												}
											}
										}
									
									/* check this new HGT against each of the older hgts as they are added  */
									tmp_allow = FALSE;
									
									for(r=0; r<(2*number_of_taxa); r++)
										{
										if(presence1[r] == TRUE && dependent_species_allowed[r][hgt_receipient2] == TRUE)
											tmp_allow = TRUE;
										if(presence1[r] == TRUE && dependent_species_allowed[r][hgt_receipient2] == FALSE)
											presence1[r] = FALSE;
										}
									if(tmp_allow == FALSE)
										{
										allow_HGT1 = FALSE;
										q = (2*number_of_taxa);
										}
									/* this previous section checks the  newhgt against the rules put in place by previous HGTs and asks if this new HGT violates these rules */
									/* We need one other check. There is a chance that this new HGT may not conflict with any of the rules put in place by previous HGTs, but by its addition it causes some previous HGTs to violate the rules */
									/* this is what needs to be checked next */
									
									if(q>1)
										{
										/* dependent_species_allowed has all the rules imposed by each of the HGTs examined so far */
										
										for(s=1; s<q; s++)
											{
											tmp_allow = FALSE;
											for(r=0; r<(2*number_of_taxa); r++)
												{
												if(tmp_presence1[s][r] == TRUE)
													{
													if(dependent_species_allowed[r][overall_receptor[s]] == TRUE) tmp_allow = TRUE;
													}
												}
											if(tmp_allow == FALSE)
												{
												allow_HGT1 = FALSE;
												s = q;
												}
											}
										}
									
									} 
										
								}
							if(allow_HGT2)
								{
								/* check both copy and test_part for any previous HGTs donor sites and then check what the possible donors are now, given the new HGT found */
								/* HGT2 assumes that test_part is the donor and copy is the HGT */
								/* we will use the best rooting of copy to calculate the possible donor positions of any hgts defined earlier, these will be saved in tmp_presence2 */
								/* best_mapping2  holds the optimal mapping of Copy as the donor */
								for(q=1; q<numparts; q++)
									{
									reset_tag2(best_mapping2);
									if(assign_tag2(best_mapping2, q))
										{
										find_tagged(best_mapping2, tmp_presence2[q]);
										}
									else
										{
										for(r=0; r<2*number_of_taxa; r++)
											{
											tmp_presence2[q][r] = overall_presence[q][r];
											}
										}
									}
								
								for(s=0; s<(2*number_of_taxa); s++)  /* initialise this to include the timetravelling HGTs that occur from the shape of the tree */
									{
									for(r=0; r<(2*number_of_taxa); r++) dependent_species_allowed[s][r] = species_allowed[s][r];
									}	
								
								/* NOW check to see if the new HGT causes any conflicts with the previous HGTs */
								
								for(q=1; q<numparts; q++)
									{
									
									 /* check to see if this HGT conflicts with any previously assigned ((not needed for the first hgt)) */
									
									if(allow_HGT2)
										{
										/* identify the nodes above and below the donor and receptor sites so cross hgts can be defined as being not allowed */
										for(r=0; r<(2*number_of_taxa); r++) before1[r]  = after1[r] = before2[r] = after2[r]= FALSE;
										for(r=0; r<2*number_of_taxa; r++) previous[r] = FALSE;
										for(r=0; r<(2*number_of_taxa); r++)
											{
											if(tmp_presence2[q][r] == TRUE)
												assign_before_after(species_tree, previous, before1, after1, r, FALSE);
											}
										for(r=0; r<(2*number_of_taxa); r++)
											{
											if(tmp_presence2[q][r] == TRUE) before1[r] = after1[r] = FALSE;
											}
										for(r=0; r<2*number_of_taxa; r++) previous[r] = FALSE;
										assign_before_after(species_tree, previous, before2, after2, overall_receptor[q], FALSE);
										before2[overall_receptor[q]] = after2[overall_receptor[q]] = FALSE;
										
										/** assign the new cross-hgt rules to the array dependent_species_allowed */
										for(r=0; r<2*number_of_taxa; r++)
											{
											if(before1[r] == TRUE)
												{
												for(s=0; s<2*number_of_taxa; s++)
													{
													if(after2[s] == TRUE)
														{
														dependent_species_allowed[r][s] = FALSE;
														dependent_species_allowed[s][r] = FALSE;
														}
													}
												}
											}
										for(r=0; r<2*number_of_taxa; r++)
											{
											if(before2[r] == TRUE)
												{
												for(s=0; s<2*number_of_taxa; s++)
													{
													if(after1[s] == TRUE)
														{
														dependent_species_allowed[r][s] = FALSE;
														dependent_species_allowed[s][r] = FALSE;
														}
													}
												}
											}
										}
									
									/* check this new HGT against each of the older hgts as they are added  */
									tmp_allow = FALSE;
									
									for(r=0; r<(2*number_of_taxa); r++)
										{
										if(presence2[r] == TRUE && dependent_species_allowed[r][hgt_receipient1] == TRUE)
											tmp_allow = TRUE;
										if(presence2[r] == TRUE && dependent_species_allowed[r][hgt_receipient1] == FALSE)
											presence2[r] = FALSE;
										}
									if(tmp_allow == FALSE)
										{
										allow_HGT2 = FALSE;
										q = (2*number_of_taxa);
										}
									/* this previous section checks the  newhgt against the rules put in place by previous HGTs and asks if this new HGT violates these rules */
									/* We need one other check. There is a chance that this new HGT may not conflict with any of the rules put in place by previous HGTs, but by its addition it causes some previous HGTs to violate the rules */
									/* this is what needs to be checked next */
									
									if(q>1)
										{
										/* dependent_species_allowed has all the rules imposed by each of the HGTs examined so far */
										
										for(s=1; s<q; s++)
											{
											tmp_allow = FALSE;
											for(r=0; r<(2*number_of_taxa); r++)
												{
												if(tmp_presence2[s][r] == TRUE)
													{
													if(dependent_species_allowed[r][overall_receptor[s]] == TRUE) tmp_allow = TRUE;
													}
												}
											if(tmp_allow == FALSE)
												{
												allow_HGT2 = FALSE;
												s = q;
												}
											}
										}
									
									} 
										
								}							
							
							}
						
						if(best_mapping1 != NULL) dismantle_tree(best_mapping1);
							best_mapping1 = NULL;
						if(best_mapping2 != NULL) dismantle_tree(best_mapping2);
							best_mapping2 = NULL;
						
						/***************************************************************/
						/* now see if the reconstruction cost of either of the HGTs that are still allowed would be less than (or equal to?) the best so far (which may be the original) */
						if(allow_HGT1 && allow_HGT2)
							{
							
							if(best_total1 + HGT2 + (1*hgt_weight) < best_total2 + HGT1 + (1*hgt_weight) && best_total1 + HGT2 + (1*hgt_weight) < best_reconstruction)
								{
								found_better = TRUE;
								best_reconstruction = best_total1 + HGT2 + (1*hgt_weight);
								donor_score = best_total1;
								receptor = hgt_receipient2;
								HGT_score = HGT2;
								if(best_donor != NULL) dismantle_tree(best_donor);
								if(best_HGT != NULL) dismantle_tree(best_HGT);
								best_donor = copy;
								copy = NULL;
								best_HGT = test_part;
								test_part = NULL;
								for(q=0; q<2*number_of_taxa; q++) presence[q] = presence1[q];
								for(q=1; q<numparts; q++)
									{
									for(r=0; r<2*number_of_taxa; r++) overall_presence[q][r] = tmp_presence1[q][r];
									}
								}
							else
								{
								if(best_total2 + HGT1 + (1*hgt_weight) <= best_total1 + HGT2 + (1*hgt_weight) && best_total2 + HGT1 + (1*hgt_weight) < best_reconstruction)
									{
									found_better = TRUE;
									best_reconstruction = best_total2 + HGT1 + (1*hgt_weight);
									donor_score = best_total2;
									receptor = hgt_receipient1;
									HGT_score = HGT1;
									if(best_donor != NULL) dismantle_tree(best_donor);
									if(best_HGT != NULL) dismantle_tree(best_HGT);
									best_donor = test_part;
									test_part = NULL;
									best_HGT = copy;
									copy = NULL;
									for(q=0; q<2*number_of_taxa; q++) presence[q] = presence2[q];
									for(q=1; q<numparts; q++)
										{
										for(r=0; r<2*number_of_taxa; r++) overall_presence[q][r] = tmp_presence2[q][r];
										}
									}
								
								}
							}
						else
							{
							if(allow_HGT1) /* this assumes that "copy" is the donor and "test_part" is the HGT */
								{
								if(best_total1 + HGT2 + (1*hgt_weight) < best_reconstruction)
									{ /*do something */
									found_better = TRUE;
									best_reconstruction = best_total1 + HGT2 + (1*hgt_weight) ;
									donor_score = best_total1;
									receptor = hgt_receipient2;
									HGT_score = HGT2;
									if(best_donor != NULL) dismantle_tree(best_donor);
									if(best_HGT != NULL) dismantle_tree(best_HGT);
									best_donor = copy;
									copy = NULL;
									best_HGT = test_part;
									test_part = NULL;
									for(q=0; q<2*number_of_taxa; q++) presence[q] = presence1[q];
									for(q=1; q<numparts; q++)
										{
										for(r=0; r<2*number_of_taxa; r++) overall_presence[q][r] = tmp_presence1[q][r];
										}
									}
								}
							if(allow_HGT2) /* this assumes that "test_part" is the donor and "copy" is the HGT */
								{
								if(best_total2 + HGT1 + (1*hgt_weight) < best_reconstruction)
									{ /*do something */
									found_better = TRUE;
									best_reconstruction = best_total2 + HGT1 + (1*hgt_weight);
									donor_score = best_total2;
									receptor = hgt_receipient1;
									HGT_score = HGT1;
									if(best_donor != NULL) dismantle_tree(best_donor);
									if(best_HGT != NULL) dismantle_tree(best_HGT);
									best_donor = test_part;
									test_part = NULL;
									best_HGT = copy;
									copy = NULL;
									for(q=0; q<2*number_of_taxa; q++) presence[q] = presence2[q];
									for(q=1; q<numparts; q++)
										{
										for(r=0; r<2*number_of_taxa; r++) overall_presence[q][r] = tmp_presence2[q][r];
										}
									}
								}
							}
						}
					if(copy != NULL) dismantle_tree(copy);
					copy = NULL;
					if(test_part != NULL) dismantle_tree(test_part);
					test_part = NULL;
					}
				}
			if(found_better == TRUE) /* if we have found a better reconstruction that includes a HGT for this part of the tree */
				{
				/* dismantle the tree at parts[k] */
				printf2("HGT%d\n", numparts);
				overall_receptor[place_marker] = receptor;
				overall_placements[k] = donor_score;
				overall_placements[place_marker] = HGT_score;
				if(parts[k] != NULL) dismantle_tree(parts[k]);
				if(parts[place_marker] != NULL) dismantle_tree(parts[place_marker]);
				parts[k] = best_donor;
				best_donor = NULL;
				strcpy(temptree, "");
				print_tree(parts[k], temptree);
				printf2("donor:\n%s\n", temptree);
				parts[place_marker] = best_HGT;
				strcpy(temptree, "");
				print_tree(parts[place_marker], temptree);
				printf2("hgt:\n%s\n", temptree);

				best_HGT = NULL;
				assign_hgtdonors(parts[k], num_gene_nodes, place_marker);
				for(q=0; q<2*number_of_taxa; q++) overall_presence[place_marker][q] = presence[q]; /* record the possible donor nodes for this HGT */
				
				numparts++;
				k--; /* we need to re-evaluate the donor again to see if removing this HGT means that another HGT becomes apparant */
				}
			
			}
		printf2("\n\nprinting results\n");
			
		/**** print out the results */
		i=0;
		while(parts[i] != NULL && i < num_gene_nodes)
			{
			printf2("part %d\nScore = %f\nReceptor Node:%d\nPossible Donors: ",i,overall_placements[i], overall_receptor[i] );
			for(j=0; j<2*number_of_taxa; j++) printf2("%d,", overall_presence[i][j]);
			printf2("\n");
			strcpy(temptree, "");
			strcpy(temptree, "");
			print_tree(parts[i], temptree);
			printf2("parts[%d] %s\n",i, temptree);
			printf2("\n");
			/*print_tree(parts[i], temptree);
			tree_coordinates(temptree, TRUE, FALSE, FALSE);
			*/
			i++;
			}
		
			
		if(copy != NULL)
			{
			dismantle_tree(copy);
			copy = NULL;
			}
		if(gene_tree != NULL)
			{
			dismantle_tree(gene_tree);
			gene_tree = NULL;
			}
		if(best_mapping != NULL)
			{
			dismantle_tree(best_mapping);
			best_mapping = NULL;
			}
		for(i=0; i<num_gene_nodes; i++)
			{
			if(parts[i] != NULL) dismantle_tree(parts[i]);
			parts[i] = NULL;
			}
		
		free(parts);
		parts = NULL;
		free(overall_receptor);
		overall_receptor = NULL;
		for(i=0; i<num_gene_nodes; i++)
			{
			free(tmp_presence1[i]);
			tmp_presence1[i] = NULL;
			free(tmp_presence2[i]);
			tmp_presence2[i] = NULL;
			free(overall_presence[i]);
			overall_presence[i] = NULL;
			}
		
		free(overall_presence);
		overall_presence = NULL;
		free(tmp_presence1);
		tmp_presence1 = NULL;
		free(tmp_presence2);
		tmp_presence2 = NULL;
		
		free(overall_placements);
		overall_placements = NULL;
		}
	
	for(i=0; i<(2*number_of_taxa); i++) free(species_allowed[i]);
	free(species_allowed);
	species_allowed = NULL;
	for(i=0; i<(2*number_of_taxa); i++) free(dependent_species_allowed[i]);
	free(dependent_species_allowed);
	dependent_species_allowed = NULL;
	
	free(before1); 
	free(before2); 
	free(after1); 
	free(after2);

	
	if(species_tree != NULL)
		{
		dismantle_tree(species_tree);
		species_tree = NULL;
		}
	
	if(gene_tree != NULL)
		{
		dismantle_tree(gene_tree);
		gene_tree = NULL;
		}
	
	tree_top = NULL;
	free(temptree);
	free(temptree1);
	free(presence);
	free(presence1);
	free(presence2);
	free(temporary);



	}

void assign_before_after(struct taxon *position, int *previous, int *before, int *after, int num, int found)
	{
	int i=0;
	
	while(position != NULL)
		{
		if(position->tag == num)
			{
			for(i=0; i<(2*number_of_taxa); i++)
				before[i] = previous[i];  
			found = TRUE;
			}
		else
			{
			if(found) after[position->tag] = TRUE;
			}
			
		if(position->daughter != NULL)
			{
			previous[position->tag] = TRUE;
			assign_before_after(position->daughter,previous, before, after, num, found);
			previous[position->tag] = FALSE;
			}
		if(position->tag == num) found = FALSE;
		position = position->next_sibling;
		}
	}

void assign_ances_desc(struct taxon *position, int ** allowed_species, int * previous)
	{
	int i=0;
	
	while(position != NULL)
		{
		if(position->daughter != NULL)
			{
			previous[position->tag] = TRUE;
			assign_ances_desc(position->daughter, allowed_species, previous);
			previous[position->tag] = FALSE;
			}
		allowed_species[position->tag][position->tag] = FALSE;
		for(i=0; i<(2*number_of_taxa); i++)
			{
			if(previous[i])
				{
				allowed_species[i][position->tag] = FALSE;
				allowed_species[position->tag][i] = FALSE;
				}
			}
		position = position->next_sibling;
		}
	}

void resolve_tricotomies_dist (struct taxon *gene_tree, struct taxon *species_tree, int ** scores)
	{
	struct taxon *position = gene_tree, *start = gene_tree, *position2 = NULL, *new = NULL, *first = NULL, *second = NULL;
	int i=0, j, minscore, *presence = NULL;
	
	presence = malloc(number_of_taxa*sizeof(int));
	for(i=0; i<number_of_taxa; i++) presence[i] = FALSE;
	i=0;
	while(position != NULL)
		{
		i++;
		/* go though looking for internal branches and follow them down */
		if(position->daughter != NULL)
			resolve_tricotomies_dist(position->daughter, species_tree, scores);
		position = position->next_sibling;
		}
	
	/* for trichotomies at this level, idenfiy the internal branch on the species tree that best matches */
	while((temp_top2 != start && i > 2) || (temp_top2 == start && i > 3))
		{
		position = start;
		minscore=-1; first = NULL; second = NULL; 
		/* identify the branches at this level with the minimum distance between them */
		while(position->next_sibling != NULL)
			{
			position2 = position->next_sibling;
			while(position2 != NULL)
				{
				if(minscore==-1 || scores[position->tag][position2->tag] < minscore)
					{
					minscore = scores[position->tag][position2->tag];
					first = position;
					second = position2;
					}
				else
					{
					if((first->tag >= number_of_taxa || second->tag >= number_of_taxa) && (position->tag < number_of_taxa && position2->tag < number_of_taxa) && (scores[position->tag][position2->tag] == minscore))
						{
						minscore = scores[position->tag][position2->tag];
						first = position;
						second = position2;						
						}
					}
				position2=position2->next_sibling;
				}
			position = position->next_sibling;
			}
		position= first;
		position2= second;
		/* having identified the parts with the minimum distance, put them together */
		/* create a new internal branch */
		new = make_taxon();
		/* make the two min branches daughters of the new node */
		if(position->next_sibling != NULL) (position->next_sibling)->prev_sibling = position->prev_sibling;
		if(position->prev_sibling != NULL) (position->prev_sibling)->next_sibling = position->next_sibling;
		if(position->prev_sibling == NULL)
			{
			start = position->next_sibling;
			if(position->parent == NULL) 
				temp_top2 = start;
			else 
				(position->parent)->daughter = position->next_sibling;
			(position->next_sibling)->parent = position->parent;
			
			}
		position->parent = new;
		new->daughter = position;
		position->next_sibling = position2;
		position->prev_sibling = NULL;
		
		if(position2->next_sibling != NULL) (position2->next_sibling)->prev_sibling = position2->prev_sibling;
		if(position2->prev_sibling != NULL) (position2->prev_sibling)->next_sibling = position2->next_sibling;
		if(position2->prev_sibling == NULL)
			{
			start = position2->next_sibling;
			if(position2->parent == NULL)
				temp_top2 = start;
			else
				{
				(position2->parent)->daughter = position2->next_sibling;
				(position2->next_sibling)->parent = position2->parent;
				position2->parent = NULL;
				}
			
			
			}
		position2->prev_sibling = position;
		position2->next_sibling = NULL;
		/* put "new" into the tree at the present position (after the first) */
		if(start->next_sibling != NULL) (start->next_sibling)->prev_sibling = new;
		new->prev_sibling = start;
		new->next_sibling = start->next_sibling;
		start->next_sibling = new;
		
		/* Find the label for the new internal brnach from the species tree */
		if(position->tag == position2->tag)
			new->tag = position->tag;
		else
			{
			for(j=0; j<number_of_taxa; j++) presence[j] = FALSE;
			get_taxa(new->daughter, presence);
			new->tag = get_best_node(species_tree, presence, -1);
			
			}
		i--;
		}
	free(presence);
	}

/* resolve_tricotomies_mindup(): soft-polytomy resolution. Same greedy pairwise
 * machinery and pointer surgery as resolve_tricotomies_dist(), but the pair of
 * sibling branches merged at each step is chosen to MINIMISE duplications rather
 * than by species-tree path distance: prefer a pair whose LCA-map is a
 * speciation (the LCA equals neither child's map), and among those the deepest
 * LCA. A polytomy is "soft" -- treated as uncertainty about the true binary
 * shape -- so this yields the resolution with the fewest duplications, avoiding
 * the spurious duplications an arbitrary (distance-guided) resolution can force.
 *
 * Correctness: deepest-speciation-first greedy provably attains the minimum
 * number of duplications over ALL binary resolutions of the polytomy (verified
 * by exhaustive enumeration + an independent bottom-up max-speciation DP over
 * 400k random species-tree/polytomy instances, 0 mismatches). On a binary node
 * it does nothing, so binary gene trees are unaffected.
 *
 * Requires the species tree to be number_tree1-numbered and build_species_partag()
 * to have run (for species_lca_tag()/g_species_depth); xnum is the species root
 * sentinel. New internal nodes are mapped to the true LCA of their two children. */
void resolve_tricotomies_mindup(struct taxon *gene_tree, struct taxon *species_tree, int xnum)
	{
	struct taxon *position = gene_tree, *start = gene_tree, *position2 = NULL, *new = NULL, *first = NULL, *second = NULL;
	int i=0, bestspec, bestdepth, l, spec, d;
	(void)species_tree;

	i=0;
	while(position != NULL)
		{
		i++;
		/* resolve any polytomies deeper in the tree first (child maps are invariant) */
		if(position->daughter != NULL)
			resolve_tricotomies_mindup(position->daughter, species_tree, xnum);
		position = position->next_sibling;
		}

	/* internal level: reduce to 2 children; unrooted top level: reduce to 3 */
	while((temp_top2 != start && i > 2) || (temp_top2 == start && i > 3))
		{
		position = start;
		bestspec = -1; bestdepth = -1; first = NULL; second = NULL;
		/* pick the sibling pair whose merge best reduces duplications:
		 * speciation first, then deepest LCA (matches the verified greedy) */
		while(position->next_sibling != NULL)
			{
			position2 = position->next_sibling;
			while(position2 != NULL)
				{
				l    = species_lca_tag(position->tag, position2->tag, xnum);
				spec = (l != position->tag && l != position2->tag) ? 1 : 0;
				d    = (l == xnum) ? -1 : g_species_depth[l];
				if(first == NULL || spec > bestspec || (spec == bestspec && d > bestdepth))
					{ bestspec = spec; bestdepth = d; first = position; second = position2; }
				position2 = position2->next_sibling;
				}
			position = position->next_sibling;
			}
		position = first;
		position2 = second;
		/* merge `first` and `second` under a new internal node -- pointer surgery
		 * identical to resolve_tricotomies_dist() */
		new = make_taxon();
		if(position->next_sibling != NULL) (position->next_sibling)->prev_sibling = position->prev_sibling;
		if(position->prev_sibling != NULL) (position->prev_sibling)->next_sibling = position->next_sibling;
		if(position->prev_sibling == NULL)
			{
			start = position->next_sibling;
			if(position->parent == NULL)
				temp_top2 = start;
			else
				(position->parent)->daughter = position->next_sibling;
			(position->next_sibling)->parent = position->parent;
			}
		position->parent = new;
		new->daughter = position;
		position->next_sibling = position2;
		position->prev_sibling = NULL;

		if(position2->next_sibling != NULL) (position2->next_sibling)->prev_sibling = position2->prev_sibling;
		if(position2->prev_sibling != NULL) (position2->prev_sibling)->next_sibling = position2->next_sibling;
		if(position2->prev_sibling == NULL)
			{
			start = position2->next_sibling;
			if(position2->parent == NULL)
				temp_top2 = start;
			else
				{
				(position2->parent)->daughter = position2->next_sibling;
				(position2->next_sibling)->parent = position2->parent;
				position2->parent = NULL;
				}
			}
		position2->prev_sibling = position;
		position2->next_sibling = NULL;
		/* splice `new` into the sibling list after `start` */
		if(start->next_sibling != NULL) (start->next_sibling)->prev_sibling = new;
		new->prev_sibling = start;
		new->next_sibling = start->next_sibling;
		start->next_sibling = new;
		/* the new node maps to the LCA of its two children's species-tree nodes */
		new->tag = species_lca_tag((new->daughter)->tag, ((new->daughter)->next_sibling)->tag, xnum);
		i--;
		}
	}

