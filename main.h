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

#include "globals.h"
#include "utils.h"

/* -----------------------------------------------------------------------
 * Functions defined in treecompare2.c (or other modules) called by main.c.
 * ----------------------------------------------------------------------- */
void  heuristic_search(int user, int print, int sample, int nreps);
void  bootstrap_search(void);
void  mlscores(void);
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
