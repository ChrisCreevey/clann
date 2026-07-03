/*
 *  reconcile.h — Clann v5.0.0
 *  Gene-tree / species-tree reconciliation, duplication/loss/HGT inference,
 *  neighbour-joining supertree, and auxiliary tree-map helpers.
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

#ifndef CLANN_RECONCILE_H
#define CLANN_RECONCILE_H

#include "globals.h"
#include "utils.h"
#include "tree_ops.h"
#include "scoring.h"
#include "consensus.h"
#include "tree_io.h"

/* -----------------------------------------------------------------------
 * treecompare2.c functions called by reconcile functions
 * ----------------------------------------------------------------------- */
int  *apply_singlecopy_filter(void);
void  restore_singlecopy_filter(int *saved);
void  tree_coordinates(char *tree, int bootstrap, int build, int mapping, int fundnum);

/* -----------------------------------------------------------------------
 * Neighbour-joining supertree command
 * ----------------------------------------------------------------------- */
void  nj(void);
int   resolve_guide_tree(char *speciestree_opt, char *out_tree, char *error_msg);
int   decompose_gene_tree_stage2(int treenum, char *guide_tree_str, float dupsupport, int minfragtaxa, int minfragspecies, FILE *fragfile, FILE *infofile);

/* -----------------------------------------------------------------------
 * decomposegenetrees Stage 3: weighting and pool-integration prep
 * (see NOTES_gene_tree_decomposition.md §5 Stage 3 / §9 step 4). Does NOT
 * itself touch fundamentals[]/presence_of_taxa[][]/Total_fund_trees -- it
 * only produces an in-memory fragment list (with weights already assigned)
 * plus the exact bracket-annotated text ("[weight](tree);[name]") that a
 * caller can write to disk and reload via clann_load_trees()/execute_command,
 * per §6.2's "let the existing loader do pool registration" decision. That
 * reload wiring itself is step 5/6's job, not this one's. */
struct decompose_fragment
	{
	char  *newick;          /* fragment Newick, trailing ';', NO weight bracket */
	char   name[NAME_LENGTH];
	float  weight;          /* 1.0 for Stage 1 passthroughs, 1.0/k for Stage 2 output */
	int    source_treenum;  /* index into fundamentals[]/tree_names[] this came from */
	int    fragindex;       /* 1-based index of this fragment within its source family */
	int    k;               /* total surviving fragments from source_treenum (this->weight == 1.0f/k, or 1 for a passthrough) */
	};

int   decompose_gene_trees_stage3(char *guide_tree_str, float dupsupport, int minfragtaxa, int minfragspecies, FILE *infofile, struct decompose_fragment **out_frags);
char *build_decompose_output_text(struct decompose_fragment *frags, int n);
void  free_decompose_fragments(struct decompose_fragment *frags, int n);

/* -----------------------------------------------------------------------
 * decomposegenetrees CLI command (see NOTES_gene_tree_decomposition.md §3,
 * §9 step 5). Non-destructive: writes <filename>/<filename>_info.txt only.
 * ----------------------------------------------------------------------- */
void  decompose_gene_trees_cmd(void);

/* -----------------------------------------------------------------------
 * Gene-tree / species-tree reconciliation
 * ----------------------------------------------------------------------- */
void  reconstruct(int print_settings);
void  print_nhx_tree(struct taxon *position, char *buf);
float tree_map(struct taxon *gene_top, struct taxon *species_top, int print);
void  label_gene_tree(struct taxon *gene_position, struct taxon *species_top, int *presence, int xnum);
int   reconstruct_map(struct taxon *position, struct taxon *species_top);
void  add_losses(struct taxon *position, struct taxon *species_top);
struct taxon *construct_tree(struct taxon *spec_pos, struct taxon *gene_pos, int *presence, struct taxon *extra_gene);
int   join_losses(struct taxon *position);
int   count_losses(struct taxon *position);
float get_recon_score(char *giventree, int numspectries, int numgenetries);
void  put_in_scores(struct taxon *position, float *total);
void  assign_ances_desc(struct taxon *position, int **allowed_species, int *previous);
void  isittagged(struct taxon *position);

/* -----------------------------------------------------------------------
 * HGT reconstruction
 * ----------------------------------------------------------------------- */
void  hgt_reconstruction(void);
void  assign_hgtdonors(struct taxon *position, int num, int part_num);
void  assign_before_after(struct taxon *position, int *previous, int *before, int *after, int num, int found);

/* -----------------------------------------------------------------------
 * Tree-map / topology helpers
 * ----------------------------------------------------------------------- */
void  resolve_tricotomies(struct taxon *position, struct taxon *species_tree);
struct taxon *do_resolve_tricotomies(struct taxon *gene_tree, struct taxon *species_tree, int basescore);
int   presence_of_trichotomies(struct taxon *position);
void  resolve_tricotomies_dist(struct taxon *gene_tree, struct taxon *species_tree, int **scores);
void  print_tree_labels(struct taxon *position, int **results, int treenum, struct taxon *species_tree);
void  subtree_id(struct taxon *position, int *tmp);
void  descend(struct taxon *position, int *presence);
void  gene_content_parsimony(struct taxon *position, int *array);

/* -----------------------------------------------------------------------
 * Mapping unknowns
 * ----------------------------------------------------------------------- */
void  mapunknowns(void);

#endif /* CLANN_RECONCILE_H */
