/*
 *  scoring.h — Clann v5.0.0
 *  Supertree scoring criteria: distance fit, splits, quartet, RF, ML.
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

#ifndef CLANN_SCORING_H
#define CLANN_SCORING_H

#include "globals.h"
#include "utils.h"
#include "tree_ops.h"
#include "topology.h"

/* fund_bipart_sets: defined in treecompare2.c, used by rf scoring */
extern BipartSet *fund_bipart_sets;
#ifdef _OPENMP
#pragma omp threadprivate(fund_bipart_sets)
#endif

/* -----------------------------------------------------------------------
 * Distance and path-metric scoring
 * ----------------------------------------------------------------------- */
void  cal_fund_scores(int printfundscores);
void  pathmetric(char *string, int **scores);
void  pathmetric_internals(char *string, struct taxon *species_tree, int **scores);
void  calculate_withins(struct taxon *position, int **within, int *presence);
void  weighted_pathmetric(char *string, float **scores, int fund_num);

/* -----------------------------------------------------------------------
 * Supertree comparison / optimality criteria
 * ----------------------------------------------------------------------- */
float compare_trees(int spr);
void  rf_precompute_fund_biparts(void);
float compare_trees_rf(int spr);
float compare_trees_ml(int spr);
float compare_trees_sfit(int spr);
float compare_trees_qfit(int spr);
float MRC(char *supertree);
float quartet_compatibility(char *supertree);

#endif /* CLANN_SCORING_H */
