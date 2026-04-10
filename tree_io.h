/*
 *  tree_io.h — Clann v5.0.0
 *  Source-tree file reading, parsing, display, and filter commands.
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

#ifndef CLANN_TREE_IO_H
#define CLANN_TREE_IO_H

#include "clann.h"
#include "utils.h"

/* -----------------------------------------------------------------------
 * Globals defined in treecompare2.c used by tree_io functions.
 * ----------------------------------------------------------------------- */
extern char ***fulltaxanames;
extern int   *numtaxaintrees;
extern char **fundamentals;
extern char **stored_funds;
extern char **tree_names;
extern float *tree_weights;
extern char **taxa_names;
extern char **parsed_command;
extern char **stored_commands;
extern char  *tempsuper;
extern int    Total_fund_trees;
extern int    number_of_taxa;
extern int    num_excluded_trees;
extern int    num_excluded_taxa;
extern int    num_commands;
extern int    largest_tree;
extern int    smallest_tree;
extern int    fundamental_assignments;
extern int    tree_length_assignments;
extern int    name_assignments;
extern int    max_name_length;
extern int    got_weights;
extern int    delimiter;
extern int    parts;
extern int    trees_in_memory;
extern int    calculated_fund_scores;
extern int    select_longest;
extern int   *same_tree;
extern int   *number_of_comparisons;
extern int   *sourcetreetag;
extern int   *taxa_incidence;
extern int  **Cooccurrance;
extern int  **presence_of_taxa;
extern int  **stored_presence_of_taxa;
extern float *scores_retained_supers;
extern float *sourcetree_scores;
extern double sup;
extern char   delimiter_char;
extern char   inputfilename[10000];

/* The following variables are threadprivate in treecompare2.c.
 * The threadprivate pragma must be repeated in every TU that uses them,
 * per OpenMP §2.15.2. */
extern struct taxon *tree_top;
extern struct taxon *temp_top;
#ifdef _OPENMP
#pragma omp threadprivate(tree_top, temp_top, \
                          fundamentals, presence_of_taxa, \
                          number_of_comparisons, sourcetree_scores)
#endif

/* -----------------------------------------------------------------------
 * Functions defined in treecompare2.c called by tree_io functions.
 * ----------------------------------------------------------------------- */
void  totext(int c, char *array);
void  memory_error(int error_num);
int   tree_build(int c, char *treestring, struct taxon *parent, int fromfile, int fund_num, int taxaorder);
int   basic_tree_build(int c, char *treestring, struct taxon *parent, int fullnames);
void  dismantle_tree(struct taxon *position);
int   count_taxa(struct taxon *position, int count);
void  reset_tree(struct taxon *position);
int   find_taxa(struct taxon *position, char *query);
void  prune_tree(struct taxon *super_pos, int fund_num);
int   shrink_tree(struct taxon *position);
int   print_pruned_tree(struct taxon *position, int count, char *pruned_tree, int fullname, int treenum);
int   number_tree(struct taxon *position, int num);
int   check_taxa(struct taxon *position);
void  print_named_tree(struct taxon *position, char *tree);
void  print_fullnamed_tree(struct taxon *position, char *tree, int fundtreenum);
void  reallocate_retained_supers(void);
void  execute_command(char *filename, int do_all);
void  cal_fund_scores(int printfundscores);

/* -----------------------------------------------------------------------
 * Tree I/O and filter function declarations
 * ----------------------------------------------------------------------- */

/* File reading / parsing */
void     input_fund_tree(char *intree, int fundnum);
int      comment(FILE *file);
int      nexusparser(FILE *nexusfile);
void     input_file_summary(int do_all);
int      assign_taxa_name(char *name, int fund);

/* Display / filter commands */
void     showtrees(int savet);
void     exclude(int do_all);
void     include(int do_all);
void     exclude_taxa(int do_all);
void     restoretaxa(int do_all);
void     returntree(char *temptree);
void     returntree_fullnames(char *temptree, int treenum);

/* restoretaxa availability flag (set by exclude_taxa, cleared by restoretaxa) */
extern int restoretaxa_available;

#endif /* CLANN_TREE_IO_H */
