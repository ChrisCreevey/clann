/*
 *  globals.h — Clann v5.0.0
 *  Extern declarations for all file-scope globals defined in treecompare2.c,
 *  plus the OpenMP threadprivate pragma repeated for every translation unit
 *  that accesses a threadprivate variable (OpenMP §2.15.2).
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

#ifndef CLANN_GLOBALS_H
#define CLANN_GLOBALS_H

#include "clann.h"

/* -----------------------------------------------------------------------
 * File-scope globals defined in treecompare2.c
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
extern double ml_eta;    /* [experimental] tree-size scaling exponent: 0=Steel 2008, 1=normalised, >1=downweight large trees */
extern time_t interval1, interval2;
extern double sup;
extern char   saved_supertree[TREE_LENGTH];
extern char  *test_array;
extern char   inputfilename[10000];
extern char   delimiter_char;
extern char   logfile_name[10000];
extern char   system_call[100000];
extern int    trees_in_memory, *sourcetreetag, remainingtrees, GC;
extern volatile int user_break;
extern int    delimiter, print_log, num_gene_nodes, testarraypos;
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
extern int     hs_strategy;
extern int     hs_progress_interval;
extern time_t  par_last_progress_time;
extern float   hs_droprep;
extern int     rep_abandon;
extern int     hs_par_rep;
extern int     hs_thread_report_interval;
extern time_t  thread_report_last;

/* -----------------------------------------------------------------------
 * OpenMP threadprivate declarations — repeated per OpenMP §2.15.2 in
 * every TU that accesses a threadprivate variable.
 * ----------------------------------------------------------------------- */
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
    skip_streak, rep_abandon, hs_par_rep, thread_report_last, \
    fundamentals, presence_of_taxa, fund_scores, number_of_comparisons \
)
#endif

#endif /* CLANN_GLOBALS_H */
