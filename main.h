/*
 *  main.h — Clann v5.0.0
 *  Entry point, command dispatch, and CLI interface declarations.
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

#ifndef CLANN_MAIN_H
#define CLANN_MAIN_H

#include "clann.h"
#include "utils.h"

/* -----------------------------------------------------------------------
 * Globals defined in treecompare2.c used by main.c functions.
 * ----------------------------------------------------------------------- */
extern FILE  *infile, *BR_file, *commands_file, *psfile, *logfile;
extern FILE  *distributionreconfile, *onetoonefile, *strictonetoonefile, *tempoutfile;
extern char **taxa_names, *commands_filename, ***fulltaxanames;
extern char **parsed_command, **fundamentals, **stored_funds;
extern char **retained_supers, **stored_commands, *tempsuper;
extern char **best_topology, **tree_names;
extern char **original_fundamentals;
extern int   autoprunemono_active;
extern float *tree_weights;
extern int  *numtaxaintrees, fullnamesnum, fullnamesassignments;
extern int   fundamental_assignments, tree_length_assignments, parsed_command_assignments;
extern int   name_assignments;
extern int  *taxa_incidence, number_of_taxa, Total_fund_trees;
extern int  *same_tree, **Cooccurrance, NUMSWAPS;
extern int ***fund_scores, ***stored_fund_scores, **super_scores;
extern int  *number_of_comparisons, *stored_num_comparisons;
extern int  **presence_of_taxa, **stored_presence_of_taxa;
extern int  *presenceof_SPRtaxa;
extern int   seed, num_commands, number_retained_supers, number_of_steps;
extern int   largest_tree, smallest_tree, criterion, parts;
extern int **total_coding, *coding_from_tree, total_nodes;
extern int   quartet_normalising, splits_weight, dweight;
extern int  *from_tree, method, tried_regrafts, hsprint;
extern int   max_name_length, got_weights, num_excluded_trees, num_excluded_taxa;
extern int   calculated_fund_scores, select_longest;
extern struct taxon *tree_top, *temp_top, *temp_top2, *branchpointer, *longestseq;
extern float *scores_retained_supers, *partition_number;
extern float  num_partitions, total_partitions, sprscore;
extern float *best_topology_scores, **weighted_scores, *sourcetree_scores;
extern float *score_of_bootstraps, *yaptp_results;
extern float  largest_length, dup_weight, loss_weight, hgt_weight, BESTSCORE;
extern float  ml_beta;
extern int    ml_scale;
extern time_t interval1, interval2;
extern double sup;
extern char   saved_supertree[TREE_LENGTH];
extern char  *test_array;
extern char   inputfilename[10000];
extern char   delimiter_char;
extern char   logfile_name[10000];
extern char   system_call[100000];
extern int    trees_in_memory, *sourcetreetag, remainingtrees, GC;
extern int    user_break, delimiter, print_log, num_gene_nodes, testarraypos;
extern int    malloc_check, count_now, another_check;
extern unsigned int thread_seed;
extern uint64_t *taxon_hash_vals;
extern VisitedSet    *visited_set;
extern VisitedSet    *thread_visited_acc;
extern LandscapeMap  *landscape_map;
extern time_t  rep_start_time;
extern int     hs_do_print;
extern float   last_status_score;
extern float   par_progress_best, par_last_print_score;
extern time_t  par_search_start;
extern int     skip_streak, hs_maxskips;

/* The following variables are threadprivate in treecompare2.c.
 * The threadprivate pragma must be repeated in every TU that uses them,
 * per OpenMP §2.15.2. */
#ifdef _OPENMP
#pragma omp threadprivate( \
    tree_top, temp_top, temp_top2, branchpointer, \
    super_scores, sourcetree_scores, presenceof_SPRtaxa, \
    sprscore, tried_regrafts, \
    retained_supers, scores_retained_supers, \
    best_topology, best_topology_scores, number_retained_supers, \
    BESTSCORE, NUMSWAPS, \
    thread_seed, \
    visited_set, thread_visited_acc, landscape_map, \
    rep_start_time, hs_do_print, last_status_score, \
    skip_streak, \
    fundamentals, presence_of_taxa, fund_scores, number_of_comparisons \
)
#endif

/* -----------------------------------------------------------------------
 * Functions defined in treecompare2.c called by main.c functions.
 * ----------------------------------------------------------------------- */
void  heuristic_search(int user, int print, int sample, int nreps);
void  bootstrap_search(void);
void  alltrees_search(int user);
void  usertrees_search(void);
void  yaptp_search(void);
void  do_consensus(void);
void  nj(void);
void  sourcetree_dists(void);
void  spr_dist(void);
void  reconstruct(int print_settings);
void  hgt_reconstruction(void);
void  mapunknowns(void);
void  generatetrees(void);
void  randomise_tree(char *tree);
void  random_prune(char *fund_tree);
void  prune_monophylies(void);
void  tips(int num);
void  exhaustive_SPR(char *string);
int   coding(int nrep, int scoring, int ptpreps);
void  cal_fund_scores(int printfundscores);
void  rf_precompute_fund_biparts(void);
void  dismantle_tree(struct taxon *position);
int   basic_tree_build(int c, char *treestring, struct taxon *parent, int fullnames);
void  reset_tree(struct taxon *position);
int   number_tree2(struct taxon *position, int num);
int   shrink_tree(struct taxon *position);
int   print_pruned_tree(struct taxon *position, int count, char *pruned_tree, int fullname, int treenum);
int   identify_species_specific_clades(struct taxon *position, int numt, int *taxa_fate, int clannID);

/* -----------------------------------------------------------------------
 * main.c public function declarations
 * ----------------------------------------------------------------------- */
int   main(int argc, char *argv[]);
void  clean_exit(int error);
int   seperate_commands(char *command);
void  print_splash(void);
void  print_commands(int num);
int   parse_command(char *command);
void  execute_command(char *commandline, int do_all);
void  set_parameters(void);
void  do_log(void);
void  controlc1(int signal);
void  controlc2(int signal);
void  controlc3(int signal);
void  controlc4(int signal);
void  controlc5(int signal);

#endif /* CLANN_MAIN_H */
