/*
 *  topology.h — Clann v5.0.0
 *  Topology hash set (VisitedSet), landscape map (LandscapeMap),
 *  and topology hash function (tree_topo_hash).
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

#ifndef CLANN_TOPOLOGY_H
#define CLANN_TOPOLOGY_H

#include "clann.h"
#include "utils.h"

/* -----------------------------------------------------------------------
 * Globals defined in treecompare2.c used by tree_topo_hash.
 * ----------------------------------------------------------------------- */
extern uint64_t *taxon_hash_vals;
extern int       number_of_taxa;

/* memory_error is defined in treecompare2.c */
void memory_error(int error_num);

/* -----------------------------------------------------------------------
 * VisitedSet — open-addressed uint64_t hash set for topology fingerprints
 * ----------------------------------------------------------------------- */
VisitedSet *vs_create(size_t cap);
void        vs_free(VisitedSet *vs);
void        vs_clear(VisitedSet *vs);
int         vs_contains(VisitedSet *vs, uint64_t key);
void        vs_insert(VisitedSet *vs, uint64_t key);
void        vs_merge(VisitedSet *dst, VisitedSet *src);

/* -----------------------------------------------------------------------
 * LandscapeMap — open-addressed hash map for visited-topology recording
 * ----------------------------------------------------------------------- */
LandscapeMap *lm_create(size_t cap);
void          lm_free(LandscapeMap *lm);
void          lm_record(LandscapeMap *lm, uint64_t hash, float score, const char *newick);
void          lm_update_score(LandscapeMap *lm, uint64_t hash, float new_score);
void          lm_merge(LandscapeMap *dst, LandscapeMap *src);
void          lm_write(LandscapeMap *lm, const char *filename);

/* -----------------------------------------------------------------------
 * Topology hash
 * ----------------------------------------------------------------------- */
uint64_t tree_topo_hash(struct taxon *root);

#endif /* CLANN_TOPOLOGY_H */
