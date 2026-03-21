/*
 *  prune.h — Clann v5.0.0
 *  Source-tree pruning: monophylies, random prune, clade collapse.
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

#ifndef CLANN_PRUNE_H
#define CLANN_PRUNE_H

#include "globals.h"
#include "utils.h"
#include "tree_ops.h"
#include "tree_io.h"
#include "main.h"

/* -----------------------------------------------------------------------
 * Pruning operations
 * ----------------------------------------------------------------------- */
void  prune_monophylies(void);
void  random_prune(char *fund_tree);
void  collapse_clades(struct taxon *position, float user_limit, int *to_delete, FILE *rp_outfile);

/* -----------------------------------------------------------------------
 * Branch-length helpers
 * ----------------------------------------------------------------------- */
int   get_brlens(struct taxon *position, float *total, int *count);
float return_length(char *string);
long  extract_length(char *fullname);

/* -----------------------------------------------------------------------
 * Species-specific clade identification
 * ----------------------------------------------------------------------- */
int   identify_species_specific_clades(struct taxon *position, int numt, int *taxa_fate, int clannID);
void  untag_nodes_below(struct taxon *position, int *taxa_fate, int clannID);
void  untag_nodes_above(struct taxon *position, int *taxa_fate, int clannID);
long  list_taxa_in_clade(struct taxon *position, int *foundtaxa, struct taxon *longest, long seqlength);
long  list_taxa_above(struct taxon *position, int *foundtaxa, struct taxon *longest, long seqlength);

/* -----------------------------------------------------------------------
 * Taxon tagging helpers
 * ----------------------------------------------------------------------- */
int   untag_taxa(struct taxon *position, int *to_delete, int keep, int count, FILE *rp_outfile);
int   print_keep(struct taxon *position, int keep, int count, FILE *rp_outfile);

#endif /* CLANN_PRUNE_H */
