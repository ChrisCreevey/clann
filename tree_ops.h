/*
 *  tree_ops.h — Clann v5.0.0
 *  Tree building, traversal, and manipulation functions.
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

#ifndef CLANN_TREE_OPS_H
#define CLANN_TREE_OPS_H

#include "globals.h"
#include "utils.h"

/* -----------------------------------------------------------------------
 * Tree building / parsing
 * ----------------------------------------------------------------------- */
int   tree_build(int c, char *treestring, struct taxon *parent, int fromfile, int fund_num, int taxaorder);
int   basic_tree_build(int c, char *treestring, struct taxon *parent, int fullnames);
struct taxon *make_taxon(void);
int   unroottree(char *tree);
int   treeToInt(char *array);
void  intTotree(int tree_num, char *array, int num_taxa);
void  totext(int c, char *array);

/* -----------------------------------------------------------------------
 * Tree structure manipulation
 * ----------------------------------------------------------------------- */
void  prune_tree(struct taxon *super_pos, int fund_num);
int   shrink_tree(struct taxon *position);
void  reset_tree(struct taxon *position);
void  dismantle_tree(struct taxon *position);
int   compress_tree(struct taxon *position);
int   compress_tree1(struct taxon *position);
void  duplicate_tree(struct taxon *orig_pos, struct taxon *prev_dup_pos);
void  reroot_tree(struct taxon *outgroup);
void  clean_pointer_taxa(struct taxon *position);
void  reallocate_retained_supers(void);

/* -----------------------------------------------------------------------
 * Tree query / traversal
 * ----------------------------------------------------------------------- */
int   count_taxa(struct taxon *position, int count);
int   find_taxa(struct taxon *position, char *query);
int   number_tree(struct taxon *position, int num);
int   number_tree1(struct taxon *position, int num);
int   number_tree2(struct taxon *position, int num);
void  check_tree(struct taxon *position, int tag_id, FILE *reconstructionfile);
int   count_internal_branches(struct taxon *position, int count);
void  identify_taxa(struct taxon *position, int *name_array);
int   check_taxa(struct taxon *position);
struct taxon *get_branch(struct taxon *position, int name);
struct taxon *get_taxon(struct taxon *position, int name);
struct taxon *find_remaining(struct taxon *position);
void  find_tagged(struct taxon *position, int *presence);
void  up_tree(struct taxon *position, int *presence);
void  down_tree(struct taxon *position, struct taxon *prev, int *presence);
int   are_siblings(struct taxon *position, int first, int second);
int   isit_onetoone(struct taxon *position, int onetoone);
void  print_onetoone_names(struct taxon *position, int onetoone);
int   get_min_node(struct taxon *position, int *presence, int num);
int   get_best_node(struct taxon *position, int *presence, int num);
void  find(struct taxon *position);
int   assign_tag2(struct taxon *position, int num);
void  reset_tag2(struct taxon *position);
void  get_taxa(struct taxon *position, int *presence);
void  get_taxa_details(struct taxon *position);
void  get_taxa_names(struct taxon *position, char **taxa_fate_names);
void  check_treeisok(struct taxon *position);

/* -----------------------------------------------------------------------
 * Tree printing
 * ----------------------------------------------------------------------- */
int   print_pruned_tree(struct taxon *position, int count, char *pruned_tree, int fullname, int treenum);
void  print_named_tree(struct taxon *position, char *tree);
void  print_fullnamed_tree(struct taxon *position, char *tree, int fundtreenum);
void  print_tree(struct taxon *position, char *tree);
void  print_single_subtree(struct taxon *node, char *buf);
void  print_tree_withinternals(struct taxon *position, char *tree);

/* -----------------------------------------------------------------------
 * Additional helpers
 * ----------------------------------------------------------------------- */
struct taxon *find_same(struct taxon *position, int tofind);

#endif /* CLANN_TREE_OPS_H */
