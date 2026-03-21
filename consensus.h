/*
 *  consensus.h — Clann v5.0.0
 *  Consensus tree construction, MRP matrix coding, average consensus.
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

#ifndef CLANN_CONSENSUS_H
#define CLANN_CONSENSUS_H

#include "globals.h"
#include "utils.h"
#include "tree_ops.h"
#include "scoring.h"
#include "tree_io.h"

/* -----------------------------------------------------------------------
 * treecompare2.c / viz.c functions called by consensus functions
 * ----------------------------------------------------------------------- */
void  neighbor_joining(int brlens, char *tree, int names);
void  tree_coordinates(char *tree, int bootstrap, int build, int mapping, int fundnum);

/* -----------------------------------------------------------------------
 * Consensus tree construction
 * ----------------------------------------------------------------------- */
void  consensus(int num_trees, char **trees, int num_reps, float percentage, FILE *outfile, FILE *guidetreefile);
void  do_consensus(void);

/* -----------------------------------------------------------------------
 * MRP / matrix coding
 * ----------------------------------------------------------------------- */
int   coding(int nrep, int scoring, int ptpreps);
int   MRP_matrix(char **trees, int num_trees, int consensus);
void  condense_coding(void);

/* -----------------------------------------------------------------------
 * Average consensus / NJ-based distances
 * ----------------------------------------------------------------------- */
int   average_consensus(int nrep, int missing_method, char *useroutfile, FILE *paupfile);

/* -----------------------------------------------------------------------
 * Descendant printing
 * ----------------------------------------------------------------------- */
void  print_descendents(struct taxon *position, FILE *outfile);
void  do_descendents(struct taxon *position, FILE *outfile);

#endif /* CLANN_CONSENSUS_H */
