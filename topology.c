/*
 *  topology.c — Clann v5.0.0
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

#include "topology.h"

/*  VisitedSet — open-addressed uint64_t hash set for topology fingerprints -----
 *  Key 0 is reserved as "empty"; if the computed hash is 0 we remap it to 1
 *  (one ghost collision possible; negligible in practice).
 */
VisitedSet *vs_create(size_t cap)
    {
    VisitedSet *vs = malloc(sizeof(VisitedSet));
    if(!vs) memory_error(201);
    vs->keys = calloc(cap, sizeof(uint64_t));
    if(!vs->keys) memory_error(202);
    vs->cap  = cap;
    vs->count = 0;
    vs->zero_present = 0;
    return vs;
    }

void vs_free(VisitedSet *vs)
    {
    if(!vs) return;
    free(vs->keys);
    free(vs);
    }

void vs_clear(VisitedSet *vs)
    {
    if(!vs) return;
    memset(vs->keys, 0, vs->cap * sizeof(uint64_t));
    vs->count = 0;
    vs->zero_present = 0;
    }

int vs_contains(VisitedSet *vs, uint64_t key)
    {
    size_t idx;
    if(!vs) return 0;
    if(key == 0) return vs->zero_present;
    idx = (size_t)(key & (vs->cap - 1));
    while(vs->keys[idx] != 0)
        {
        if(vs->keys[idx] == key) return 1;
        idx = (idx + 1) & (vs->cap - 1);
        }
    return 0;
    }

/* Maximum table size: 128 MB = 16M uint64_t slots */
#define VS_MAX_CAP (1u << 24)

void vs_insert(VisitedSet *vs, uint64_t key)
    {
    size_t idx;
    if(!vs) return;
    if(key == 0) { vs->zero_present = 1; return; }
    /* Grow if >75% full and below the cap */
    if(vs->count > vs->cap * 3 / 4 && vs->cap < VS_MAX_CAP)
        {
        size_t newcap = vs->cap * 2;
        uint64_t *newkeys = calloc(newcap, sizeof(uint64_t));
        size_t i;
        if(newkeys)
            {
            for(i = 0; i < vs->cap; i++)
                if(vs->keys[i] != 0)
                    {
                    size_t ni = (size_t)(vs->keys[i] & (newcap - 1));
                    while(newkeys[ni] != 0) ni = (ni + 1) & (newcap - 1);
                    newkeys[ni] = vs->keys[i];
                    }
            free(vs->keys);
            vs->keys = newkeys;
            vs->cap  = newcap;
            }
        /* if calloc failed, continue with old table (some false negatives) */
        }
    /* Stop inserting if at max cap and >75% full (act as probabilistic filter) */
    if(vs->count > vs->cap * 3 / 4) return;
    idx = (size_t)(key & (vs->cap - 1));
    while(vs->keys[idx] != 0)
        {
        if(vs->keys[idx] == key) return;  /* already present */
        idx = (idx + 1) & (vs->cap - 1);
        }
    vs->keys[idx] = key;
    vs->count++;
    }

void vs_merge(VisitedSet *dst, VisitedSet *src)
    {
    size_t i;
    if(!dst || !src) return;
    if(src->zero_present) vs_insert(dst, 0);
    for(i = 0; i < src->cap; i++)
        if(src->keys[i] != 0)
            vs_insert(dst, src->keys[i]);
    }


/*  LandscapeMap — open-addressed hash map for tree-space landscape recording -
 *  Records every unique topology visited during hs, with its score and the
 *  total number of times it was encountered (across all reps and threads).
 *  Key 0 is reserved as empty; hash==0 is remapped to 1 (one ghost possible).
 *  Enabled only when g_landscape_file[] is non-empty.
 */
LandscapeMap *lm_create(size_t cap)
    {
    LandscapeMap *lm = malloc(sizeof(LandscapeMap));
    if(!lm) memory_error(210);
    lm->slots = calloc(cap, sizeof(LandscapeEntry));
    if(!lm->slots) memory_error(211);
    lm->capacity = cap;
    lm->count = 0;
    return lm;
    }

void lm_free(LandscapeMap *lm)
    {
    size_t i;
    if(!lm) return;
    for(i = 0; i < lm->capacity; i++)
        if(lm->slots[i].hash != 0 && lm->slots[i].newick)
            { free(lm->slots[i].newick); lm->slots[i].newick = NULL; }
    free(lm->slots);
    free(lm);
    }

/* Grow lm to double capacity and re-insert all existing entries. */
static void lm_grow(LandscapeMap *lm)
    {
    size_t i, newcap = lm->capacity * 2;
    LandscapeEntry *newslots = calloc(newcap, sizeof(LandscapeEntry));
    if(!newslots) return; /* keep old table if OOM */
    for(i = 0; i < lm->capacity; i++)
        {
        if(lm->slots[i].hash == 0) continue;
        size_t idx = (size_t)(lm->slots[i].hash & (newcap - 1));
        while(newslots[idx].hash != 0) idx = (idx + 1) & (newcap - 1);
        newslots[idx] = lm->slots[i]; /* shallow copy — newick pointer transferred */
        }
    free(lm->slots);
    lm->slots = newslots;
    lm->capacity = newcap;
    }

/*  lm_record: insert on first visit (stores score and newick); increment
 *  visit_count on revisit.  newick may be NULL on revisit.
 *  key==0 is remapped to 1 to keep slot 0 as the sentinel.
 */
void lm_record(LandscapeMap *lm, uint64_t hash, float score, const char *newick)
    {
    size_t idx;
    if(!lm) return;
    if(hash == 0) hash = 1; /* remap sentinel */
    /* Grow if >75% full */
    if(lm->count >= lm->capacity * 3 / 4)
        lm_grow(lm);
    idx = (size_t)(hash & (lm->capacity - 1));
    while(lm->slots[idx].hash != 0)
        {
        if(lm->slots[idx].hash == hash)
            { lm->slots[idx].visit_count++; return; }
        idx = (idx + 1) & (lm->capacity - 1);
        }
    /* New entry */
    lm->slots[idx].hash        = hash;
    lm->slots[idx].score       = score;
    lm->slots[idx].visit_count = 1;
    lm->slots[idx].newick      = newick ? strdup(newick) : NULL;
    lm->count++;
    }

/*  lm_update_score: overwrite the stored score for an existing entry.
 *  Used after post-search rescoring to correct stale-cache scores recorded
 *  during the search.  No-op if the hash is not present in the map.
 */
void lm_update_score(LandscapeMap *lm, uint64_t hash, float new_score)
    {
    size_t idx;
    if(!lm) return;
    if(hash == 0) hash = 1; /* remap sentinel, mirrors lm_record */
    idx = (size_t)(hash & (lm->capacity - 1));
    while(lm->slots[idx].hash != 0)
        {
        if(lm->slots[idx].hash == hash)
            { lm->slots[idx].score = new_score; return; }
        idx = (idx + 1) & (lm->capacity - 1);
        }
    /* Not found — nothing to update */
    }

/*  lm_merge: merge src into dst.  visit_counts are summed; dst score is kept
 *  (same topology → same score, so first-seen is correct).
 */
void lm_merge(LandscapeMap *dst, LandscapeMap *src)
    {
    size_t i;
    if(!dst || !src) return;
    for(i = 0; i < src->capacity; i++)
        {
        LandscapeEntry *e = &src->slots[i];
        if(e->hash == 0) continue;
        /* Find slot in dst */
        size_t idx = (size_t)(e->hash & (dst->capacity - 1));
        /* Grow dst if needed before probing */
        if(dst->count >= dst->capacity * 3 / 4) lm_grow(dst);
        idx = (size_t)(e->hash & (dst->capacity - 1));
        while(dst->slots[idx].hash != 0)
            {
            if(dst->slots[idx].hash == e->hash)
                { dst->slots[idx].visit_count += e->visit_count; goto next_entry; }
            idx = (idx + 1) & (dst->capacity - 1);
            }
        /* New entry in dst */
        dst->slots[idx].hash        = e->hash;
        dst->slots[idx].score       = e->score;
        dst->slots[idx].visit_count = e->visit_count;
        dst->slots[idx].newick      = e->newick ? strdup(e->newick) : NULL;
        dst->count++;
        next_entry:;
        }
    }

/*  lm_read: read a landscape TSV written by lm_write() back into a new
 *  LandscapeMap.  Supports two header formats:
 *    new: index<TAB>newick<TAB>score<TAB>visit_count
 *    old: newick<TAB>score<TAB>visit_count
 *  Non-header lines that begin with '#' are silently skipped EXCEPT for
 *  the '# criterion=N' metadata line written by lm_write(), which is parsed
 *  and stored in *criterion_out (if non-NULL).  *criterion_out is initialised
 *  to -1 (meaning "not present in file") before any parsing begins.
 *  For old-format files without an index column, sequential 1-based indices
 *  are assigned as rows are read.
 *  Returns a freshly-allocated map on success, or NULL on error.
 */
LandscapeMap *lm_read(const char *filename, int *criterion_out)
    {
    FILE *fp;
    char *line = NULL;
    size_t line_alloc = 0;
    ssize_t nread;
    int lineno = 0;
    int has_index_col = 0;   /* 1 if header starts with "index\t" */
    int seq_index = 0;       /* counter used when file has no index column */
    LandscapeMap *lm;

    if(criterion_out) *criterion_out = -1;

    if(!filename || !filename[0]) return NULL;
    fp = fopen(filename, "r");
    if(!fp)
        { printf2("Error: could not open landscape file '%s' for reading\n", filename); return NULL; }

    lm = lm_create(8192);
    if(!lm) { fclose(fp); return NULL; }

    while((nread = getline(&line, &line_alloc, fp)) != -1)
        {
        char *p, *tab1, *tab2, *tab3;
        float score;
        int visit_count;
        int tree_index;
        const char *newick_ptr;
        /* Strip trailing newline / carriage return */
        while(nread > 0 && (line[nread-1] == '\n' || line[nread-1] == '\r'))
            { line[--nread] = '\0'; }

        lineno++;
        if(lineno == 1)
            {
            /* Detect format: new header starts with "index\t" */
            has_index_col = (strncmp(line, "index\t", 6) == 0);
            continue;  /* skip column header */
            }
        if(line[0] == '#')
            {
            /* Parse criterion metadata comment written by lm_write() */
            if(criterion_out && strncmp(line, "# criterion=", 12) == 0)
                *criterion_out = atoi(line + 12);
            continue;
            }
        if(line[0] == '\0') continue;

        if(has_index_col)
            {
            /* Format: index<TAB>newick<TAB>score<TAB>visit_count */
            tab1 = strchr(line, '\t');
            if(!tab1) continue;
            tree_index = atoi(line);
            *tab1 = '\0';
            newick_ptr = tab1 + 1;      /* newick starts here */

            tab2 = strchr(newick_ptr, '\t');
            if(!tab2) continue;
            *tab2 = '\0';

            tab3 = strchr(tab2 + 1, '\t');
            if(!tab3) continue;
            *tab3 = '\0';

            score       = (float)atof(tab2 + 1);
            visit_count = atoi(tab3 + 1);
            }
        else
            {
            /* Format: newick<TAB>score<TAB>visit_count */
            tab1 = strchr(line, '\t');
            if(!tab1) continue;
            *tab1 = '\0';
            p = tab1 + 1;

            tab2 = strchr(p, '\t');
            if(!tab2) continue;
            *tab2 = '\0';

            score       = (float)atof(p);
            visit_count = atoi(tab2 + 1);
            tree_index  = ++seq_index;  /* assign sequential index */
            newick_ptr  = line;
            }

        /* Intern the newick, compute a simple hash from the string content,
         * then insert.  We use a djb2-style hash so that lm_read entries
         * are self-consistent (visit_count and index are preserved via direct
         * slot write after lm_record sets them to 1 and 0 respectively). */
        {
        uint64_t h = 5381;
        const char *s = newick_ptr;
        while(*s) { h = h * 33 ^ (unsigned char)*s++; }
        if(h == 0) h = 1;

        lm_record(lm, h, score, newick_ptr);

        /* lm_record initialises visit_count to 1 and index to 0;
         * overwrite with stored values */
        {
        size_t idx = (size_t)(h & (lm->capacity - 1));
        while(lm->slots[idx].hash != 0)
            {
            if(lm->slots[idx].hash == h)
                {
                lm->slots[idx].visit_count = visit_count;
                lm->slots[idx].index       = tree_index;
                break;
                }
            idx = (idx + 1) & (lm->capacity - 1);
            }
        }
        }
        }

    free(line);
    fclose(fp);
    return lm;
    }


/*  lm_write: write TSV to filename.  Columns: index, newick, score, visit_count.
 *  A '# criterion=N' metadata comment is written after the column header so
 *  that 'recluster' can restore the correct score direction automatically.
 *  Each entry is assigned a 1-based sequential index (stored back into the
 *  LandscapeEntry.index field) so that lm_cluster can reference entries by
 *  row number in the member_indices column of the cluster TSV.            */
void lm_write(LandscapeMap *lm, const char *filename, int crit)
    {
    size_t i;
    int idx_counter = 0;
    FILE *fp;
    if(!lm || !filename || !filename[0]) return;
    fp = fopen(filename, "w");
    if(!fp) { printf2("Error: could not open landscape file '%s' for writing\n", filename); return; }
    fprintf(fp, "index\tnewick\tscore\tvisit_count\n");
    fprintf(fp, "# criterion=%d\n", crit);
    for(i = 0; i < lm->capacity; i++)
        {
        if(lm->slots[i].hash == 0) continue;
        lm->slots[i].index = ++idx_counter;
        fprintf(fp, "%d\t%s\t%.6f\t%d\n",
                lm->slots[i].index,
                lm->slots[i].newick ? lm->slots[i].newick : "",
                (double)lm->slots[i].score,
                lm->slots[i].visit_count);
        }
    fclose(fp);
    }


/*  tree_topo_hash -------------------------------------------------------
 *  Returns a 64-bit canonical hash of the unrooted topology rooted at
 *  `root`.  Uses the random-XOR bipartition method: each taxon i has a
 *  fixed weight taxon_hash_vals[i] (splitmix64); the hash of a bipartition
 *  is min(XOR of one side, XOR of other side); the tree hash is the XOR
 *  of all non-trivial bipartition hashes.  Rooting-independent.
 */
uint64_t sth_aux(struct taxon *pos, uint64_t total, uint64_t *tree_h)
    {
    uint64_t sh = 0;
    while(pos != NULL)
        {
        uint64_t child_sh;
        if(pos->daughter != NULL)
            {
            child_sh = sth_aux(pos->daughter, total, tree_h);
            /* non-trivial bipartition contribution (canonical: take smaller half) */
            uint64_t comp = total ^ child_sh;
            *tree_h ^= (child_sh < comp ? child_sh : comp);
            }
        else
            {
            /* leaf: trivial bipartition — do not add to tree_h */
            if(pos->name >= 0 && pos->name < number_of_taxa)
                child_sh = taxon_hash_vals[pos->name];
            else
                child_sh = (uint64_t)(pos->name + 1);  /* fallback */
            }
        sh ^= child_sh;
        pos = pos->next_sibling;
        }
    return sh;
    }

uint64_t tree_topo_hash(struct taxon *root)
    {
    uint64_t total = 0, tree_h = 0;
    int ti;
    if(!root || !taxon_hash_vals) return 0;
    for(ti = 0; ti < number_of_taxa; ti++) total ^= taxon_hash_vals[ti];
    /* Traverse from root itself.  In Clann's tree representation, tree_top
       is one branch of the unrooted trifurcation (not a virtual super-root):
       the trifurcation branches are tree_top, tree_top->next_sibling, etc.
       Calling sth_aux(root) processes all trifurcation branches correctly
       and is rooting-independent. */
    sth_aux(root, total, &tree_h);
    /* Remap 0 → 1 so that 0 remains the "empty" sentinel in VisitedSet */
    return (tree_h == 0) ? 1 : tree_h;
    }
