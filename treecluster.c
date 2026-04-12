/*
 *  treecluster.c — Clann v5.0.0
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

#include "treecluster.h"
#include "utils.h"
#include <limits.h>

/* -----------------------------------------------------------------------
 * Internal helpers
 * ----------------------------------------------------------------------- */

static int cmp_uint64_tc(const void *a, const void *b)
    {
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return (x > y) - (x < y);
    }

/*
 * name_to_hash — look up a taxon name and return its splitmix64 weight.
 * Returns 0 for unrecognised names (no contribution to XOR bipart hash).
 */
static uint64_t name_to_hash(const char *name)
    {
    int i;
    if(!name || !taxa_names || !taxon_hash_vals) return 0;
    for(i = 0; i < number_of_taxa; i++)
        if(taxa_names[i] && strcmp(taxa_names[i], name) == 0)
            return taxon_hash_vals[i];
    return 0;
    }

/*
 * collect_biparts_named — parse a named-taxon Newick string into a sorted
 * array of canonical bipartition hashes.  Returns the number of bipartitions
 * written to `out`.  `out` must be at least (number_of_taxa + 1) uint64_t
 * elements.  `total_hash` is XOR of all taxon_hash_vals[] for this tree's
 * taxa (pre-computed by the caller).
 *
 * Uses the same XOR-min canonical hashing as topology.c / scoring.c.
 */
static int collect_biparts_named(const char *nwk, uint64_t total_hash,
                                  uint64_t *out)
    {
    /* stack depth bounded by tree depth ≤ number_of_taxa; 2*NAME_LENGTH+4 is
     * the same generous bound used in collect_biparts_newick() in scoring.c */
    uint64_t stack[2 * NAME_LENGTH + 4];
    int depth = 0, cnt = 0, i = 0;
    stack[0] = 0;
    while(nwk[i] && nwk[i] != ';')
        {
        if(nwk[i] == '(')
            { stack[++depth] = 0; i++; }
        else if(nwk[i] == ')')
            {
            uint64_t child_sh = stack[depth--];
            uint64_t comp     = total_hash ^ child_sh;
            uint64_t bh       = (child_sh < comp) ? child_sh : comp;
            if(bh != 0) out[cnt++] = bh;
            stack[depth] ^= child_sh;
            i++;
            /* skip optional internal node label */
            while(nwk[i] && nwk[i] != '(' && nwk[i] != ')' &&
                  nwk[i] != ',' && nwk[i] != ':' && nwk[i] != ';') i++;
            }
        else if(nwk[i] == ',')
            { i++; }
        else if(nwk[i] == ':')
            { while(nwk[i] && nwk[i] != ',' && nwk[i] != ')' &&
                    nwk[i] != ';') i++; }
        else
            {   /* taxon name — read until delimiter */
            char name[NAME_LENGTH + 1]; int j = 0;
            while(nwk[i] && nwk[i] != '(' && nwk[i] != ')' &&
                  nwk[i] != ',' && nwk[i] != ':' && nwk[i] != ';')
                { if(j < NAME_LENGTH) name[j++] = nwk[i]; i++; }
            name[j] = '\0';
            stack[depth] ^= name_to_hash(name);
            }
        }
    qsort(out, cnt, sizeof(uint64_t), cmp_uint64_tc);
    return cnt;
    }

/*
 * rf_distance — symmetric RF distance between two sorted bipartition arrays.
 * rf = |a| + |b| - 2 * |a ∩ b|   (sorted-merge intersection).
 */
static int rf_distance(const uint64_t *a, int na,
                       const uint64_t *b, int nb)
    {
    int i = 0, j = 0, shared = 0;
    while(i < na && j < nb)
        {
        if(a[i] == b[j])       { shared++; i++; j++; }
        else if(a[i] < b[j])   i++;
        else                   j++;
        }
    return na + nb - 2 * shared;
    }

/* -----------------------------------------------------------------------
 * Cluster record
 * ----------------------------------------------------------------------- */

typedef struct {
    uint64_t *biparts;    /* sorted bipartition hashes of representative */
    int       nb;         /* number of bipartitions in representative */
    char     *rep_newick; /* newick of representative */
    float     rep_score;  /* score of representative */
    float     best_score; /* best (lowest) score in cluster */
    int       member_count;
    int       total_visits;
} Cluster;

/* Comparator for sorting entries by score ascending */
typedef struct { const LandscapeEntry *e; } EntryPtr;

static int cmp_entry_score(const void *a, const void *b)
    {
    float sa = ((const EntryPtr *)a)->e->score;
    float sb = ((const EntryPtr *)b)->e->score;
    return (sa > sb) - (sa < sb);
    }

/* Comparator for sorting entries by visit_count descending */
static int cmp_entry_visits(const void *a, const void *b)
    {
    int va = ((const EntryPtr *)a)->e->visit_count;
    int vb = ((const EntryPtr *)b)->e->visit_count;
    return (vb > va) - (vb < va);
    }

/* Comparator for sorting clusters by member_count descending */
static int cmp_cluster_size(const void *a, const void *b)
    {
    int ca = ((const Cluster *)a)->member_count;
    int cb = ((const Cluster *)b)->member_count;
    return (cb > ca) - (cb < ca);
    }

/* -----------------------------------------------------------------------
 * lm_cluster
 * ----------------------------------------------------------------------- */

void lm_cluster(LandscapeMap *lm,
                const char   *out_file,
                int           threshold,
                int           orderby)
    {
    size_t    i;
    int       ti;
    uint64_t  total_hash = 0;
    size_t    n_entries;
    EntryPtr *entries    = NULL;
    Cluster  *clusters   = NULL;
    size_t    n_clusters = 0;
    size_t    c_alloc    = 0;
    uint64_t *tmp_biparts = NULL;
    FILE     *fp;

    if(!lm || !out_file || !out_file[0]) return;
    if(lm->count == 0)
        {
        printf2("  Landscape clustering: no entries to cluster.\n");
        return;
        }

    /* Pre-compute total_hash = XOR of all taxon_hash_vals */
    if(taxon_hash_vals)
        for(ti = 0; ti < number_of_taxa; ti++)
            total_hash ^= taxon_hash_vals[ti];

    /* Allocate temp bipartition array (one tree at a time).
     * A fully-bifurcating tree with N taxa has N-3 internal bipartitions;
     * +2 gives a safe margin for degenerate topologies. */
    tmp_biparts = malloc((size_t)(number_of_taxa + 2) * sizeof(uint64_t));
    if(!tmp_biparts)
        { printf2("  Landscape clustering: out of memory (tmp_biparts).\n"); return; }

    /* Collect valid (non-empty newick) entries into flat array */
    n_entries = 0;
    for(i = 0; i < lm->capacity; i++)
        if(lm->slots[i].hash != 0 && lm->slots[i].newick &&
           lm->slots[i].newick[0] != '\0')
            n_entries++;

    if(n_entries == 0)
        {
        printf2("  Landscape clustering: no entries with newick strings.\n");
        free(tmp_biparts);
        return;
        }

    entries = malloc(n_entries * sizeof(EntryPtr));
    if(!entries)
        {
        printf2("  Landscape clustering: out of memory (entries).\n");
        free(tmp_biparts);
        return;
        }

    {
    size_t idx = 0;
    for(i = 0; i < lm->capacity; i++)
        if(lm->slots[i].hash != 0 && lm->slots[i].newick &&
           lm->slots[i].newick[0] != '\0')
            entries[idx++].e = &lm->slots[i];
    }

    /* Sort entries */
    if(orderby == 1)
        qsort(entries, n_entries, sizeof(EntryPtr), cmp_entry_visits);
    else
        qsort(entries, n_entries, sizeof(EntryPtr), cmp_entry_score);

    /* Initial cluster allocation */
    c_alloc  = 64;
    clusters = malloc(c_alloc * sizeof(Cluster));
    if(!clusters)
        {
        printf2("  Landscape clustering: out of memory (clusters).\n");
        free(entries);
        free(tmp_biparts);
        return;
        }

    /* Greedy sweep */
    for(i = 0; i < n_entries; i++)
        {
        const LandscapeEntry *e = entries[i].e;
        int  nb = collect_biparts_named(e->newick, total_hash, tmp_biparts);
        int  best_rf = INT_MAX;
        int  best_c  = -1;
        size_t c;

        /* Compare against existing cluster representatives */
        for(c = 0; c < n_clusters; c++)
            {
            int rf = rf_distance(tmp_biparts, nb,
                                 clusters[c].biparts, clusters[c].nb);
            if(rf < best_rf)
                {
                best_rf = rf;
                best_c  = (int)c;
                if(rf == 0) break; /* exact match — no need to keep looking */
                }
            }

        if(best_c >= 0 && best_rf <= threshold)
            {
            /* Assign to existing cluster */
            clusters[best_c].member_count++;
            clusters[best_c].total_visits += e->visit_count;
            if(e->score < clusters[best_c].best_score)
                clusters[best_c].best_score = e->score;
            }
        else
            {
            /* Open a new cluster with this tree as representative */
            if(n_clusters == c_alloc)
                {
                size_t new_alloc = c_alloc * 2;
                Cluster *tmp_cl  = realloc(clusters,
                                           new_alloc * sizeof(Cluster));
                if(!tmp_cl)
                    {
                    printf2("  Landscape clustering: out of memory (grow clusters).\n");
                    break;
                    }
                clusters = tmp_cl;
                c_alloc  = new_alloc;
                }
            {
            Cluster *cl = &clusters[n_clusters++];
            cl->nb     = nb;
            cl->biparts = malloc((size_t)(nb > 0 ? nb : 1) * sizeof(uint64_t));
            if(cl->biparts)
                memcpy(cl->biparts, tmp_biparts,
                       (size_t)nb * sizeof(uint64_t));
            cl->rep_newick   = e->newick; /* points into lm — valid until lm_free */
            cl->rep_score    = e->score;
            cl->best_score   = e->score;
            cl->member_count = 1;
            cl->total_visits = e->visit_count;
            }
            }
        }

    free(entries);
    free(tmp_biparts);

    /* Sort clusters by member_count descending */
    qsort(clusters, n_clusters, sizeof(Cluster), cmp_cluster_size);

    /* Write output TSV */
    fp = fopen(out_file, "w");
    if(!fp)
        {
        printf2("  Landscape clustering: could not open output file '%s'\n",
                out_file);
        }
    else
        {
        fprintf(fp, "cluster_id\trep_newick\tmember_count\ttotal_visits"
                    "\tbest_score\trep_score\n");
        for(i = 0; i < n_clusters; i++)
            {
            fprintf(fp, "%zu\t%s\t%d\t%d\t%.6f\t%.6f\n",
                    i + 1,
                    clusters[i].rep_newick ? clusters[i].rep_newick : "",
                    clusters[i].member_count,
                    clusters[i].total_visits,
                    (double)clusters[i].best_score,
                    (double)clusters[i].rep_score);
            }
        fclose(fp);
        printf2("  Landscape clusters written to:   %s\n", out_file);
        printf2("  Clusters found: %zu  (threshold RF=%d, orderby=%s)\n",
                n_clusters, threshold,
                orderby == 1 ? "visits" : "score");
        }

    /* Free cluster bipart arrays */
    for(i = 0; i < n_clusters; i++)
        free(clusters[i].biparts);
    free(clusters);
    }
