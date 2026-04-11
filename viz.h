/*
 *  viz.h — Clann v5.0.0
 *  Tree visualisation: coordinate assignment, PostScript output, histogram
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

#ifndef CLANN_VIZ_H
#define CLANN_VIZ_H

#include "clann.h"
#include "utils.h"

/* -----------------------------------------------------------------------
 * Globals defined in treecompare2.c that the viz functions access.
 * ----------------------------------------------------------------------- */
extern FILE   *psfile;
extern float   largest_length;
extern int     number_of_taxa;
extern char  **taxa_names;

/* tree_top / temp_top are threadprivate in treecompare2.c.
 * viz functions are only ever called from the serial (non-parallel) main
 * loop, so this is safe.  The threadprivate pragma is repeated here so
 * that viz.c is standards-compliant per OpenMP §2.15.2. */
extern struct taxon *tree_top;
extern struct taxon *temp_top;
#ifdef _OPENMP
#pragma omp threadprivate(tree_top, temp_top)
#endif

/* -----------------------------------------------------------------------
 * Forward declarations for tree-structure functions (defined in
 * treecompare2.c) called by tree_coordinates.
 * ----------------------------------------------------------------------- */
void dismantle_tree(struct taxon *position);
int  basic_tree_build(int c, char *treestring, struct taxon *parent, int fullnames);
int  count_taxa(struct taxon *position, int count);
void memory_error(int error_num);

/* -----------------------------------------------------------------------
 * Viz function declarations
 * ----------------------------------------------------------------------- */

/* X-coordinate assignment */
int   xposition1(struct taxon *position, int count);
float middle_number(struct taxon *position);
void  xposition2(struct taxon *position);

/* Y-coordinate assignment */
int   yposition0(struct taxon *position, int level, int deepest);
int   yposition1(struct taxon *position, int level);
void  yposition2(struct taxon *position, int deepest);

/* PostScript colour and coordinate output */
void  printcolour(float real, int branch);
void  print_coordinates(struct taxon *position, char **treearray, int taxa_count, int mapping);
void  tree_coordinates(char *tree, int bootstrap, int build, int mapping, int fundnum);

/* Histogram rendering */
void  draw_histogram(FILE *outfile, int bins, float *results, int num_results);

#endif /* CLANN_VIZ_H */
