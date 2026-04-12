/*
 *  treecluster.h — Clann v5.0.0
 *  Greedy RF-distance tree clustering on the hs landscape map.
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

#ifndef CLANN_TREECLUSTER_H
#define CLANN_TREECLUSTER_H

#include "clann.h"
#include "topology.h"
#include "globals.h"

/*
 * lm_cluster — greedy CD-HIT-style RF clustering of the landscape map.
 *
 *   lm         : fully-populated LandscapeMap (must be non-NULL)
 *   out_file   : output TSV filename
 *   threshold  : max normalized RF distance to join a cluster (0.0–1.0).
 *                Mirrors compare_trees_rf(): normalized by 2*(n-3) where
 *                n is the number of taxa in the tree.
 *                e.g. 0.1 = within 10% of the maximum possible RF distance.
 *   orderby    : 0 = sort by score (best first — descending for ML, ascending for all others),
 *                1 = sort by visit_count descending (most-visited first)
 *
 * Trees are sorted by the chosen criterion, then swept greedily: each tree
 * is compared against existing cluster representatives only.  First match
 * within threshold wins; no match opens a new cluster.
 *
 * Output columns (tab-separated, one header row):
 *   cluster_id  rep_newick  member_count  total_visits  best_score  worst_score
 *   score_mean  score_sd  rep_score  member_indices
 * Rows are sorted by member_count descending.
 * member_indices is a comma-separated list of the 1-based row indices from
 * the landscape TSV (index column) for every topology in this cluster.
 */
void lm_cluster(LandscapeMap *lm,
                const char   *out_file,
                float         threshold,
                int           orderby);

#endif /* CLANN_TREECLUSTER_H */
