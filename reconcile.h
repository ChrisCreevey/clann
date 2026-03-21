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

/* -----------------------------------------------------------------------
 * Gene-tree / species-tree reconciliation
 * ----------------------------------------------------------------------- */
void  reconstruct(int print_settings);
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
