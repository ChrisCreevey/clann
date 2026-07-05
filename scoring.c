/*
 *  scoring.c — Clann v5.0.0
 *  Tree scoring and comparison metrics
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
#include "scoring.h"

/* --- RF scoring internal helpers (moved from treecompare2.c) ----------- */

static int cmp_uint64(const void *a, const void *b)
    {
    uint64_t x = *(const uint64_t*)a, y = *(const uint64_t*)b;
    return (x > y) - (x < y);
    }

static int collect_biparts_newick(const char *nwk, uint64_t total_hash, uint64_t *out)
    {
    uint64_t stack[2 * NAME_LENGTH + 4];  /* depth bounded by number of taxa */
    int depth = 0, cnt = 0, i = 0;
    stack[0] = 0;
    while(nwk[i] && nwk[i] != ';')
        {
        if(nwk[i] == '(')
            { stack[++depth] = 0; i++; }
        else if(nwk[i] == ')')
            {
            uint64_t child_sh = stack[depth--];
            uint64_t comp = total_hash ^ child_sh;
            uint64_t bh   = (child_sh < comp) ? child_sh : comp;
            if(bh != 0) out[cnt++] = bh;  /* skip trivial bipartitions */
            stack[depth] ^= child_sh;
            i++;
            /* skip optional internal node label (e.g. bootstrap support) to
             * prevent it being misread as a taxon index in the else branch */
            while(nwk[i] && nwk[i] != '(' && nwk[i] != ')' &&
                  nwk[i] != ',' && nwk[i] != ':' && nwk[i] != ';') i++;
            }
        else if(nwk[i] == ',')
            { i++; }
        else if(nwk[i] == ':')
            { while(nwk[i] && nwk[i] != ',' && nwk[i] != ')' && nwk[i] != ';') i++; }
        else
            {  /* taxon integer index — only valid at leaf position (after '(' or ',') */
            char num[64]; int j = 0;
            while(nwk[i] && nwk[i] != '(' && nwk[i] != ')' && nwk[i] != ',' && nwk[i] != ':')
                num[j++] = nwk[i++];
            num[j] = '\0';
            int tidx = atoi(num);
            if(tidx >= 0 && tidx < number_of_taxa)
                stack[depth] ^= taxon_hash_vals[tidx];
            }
        }
    qsort(out, cnt, sizeof(uint64_t), cmp_uint64);
    return cnt;
    }

static int bipart_intersection_count(const uint64_t *a, int na,
                                     const uint64_t *b, int nb)
    {
    int i = 0, j = 0, shared = 0;
    while(i < na && j < nb)
        {
        if(a[i] == b[j])      { shared++; i++; j++; }
        else if(a[i] < b[j])  i++;
        else                  j++;
        }
    return shared;
    }

/* count_cherries_newick: count internal nodes in pruned_nwk whose subtree
 * contains exactly 2 taxa.  Each such node is the shared ancestor of a cherry
 * pair (a split with one 2-taxon side) in the unrooted tree.
 * Uses the same Newick-parsing idiom as collect_biparts_newick. */
static int count_cherries_newick(const char *nwk)
    {
    int cnt_stack[2 * NAME_LENGTH + 4];
    int depth = 0, cherries = 0, i = 0;
    cnt_stack[0] = 0;
    while(nwk[i] && nwk[i] != ';')
        {
        if(nwk[i] == '(')
            { cnt_stack[++depth] = 0; i++; }
        else if(nwk[i] == ')')
            {
            int child_c = cnt_stack[depth--];
            if(child_c == 2) cherries++;
            cnt_stack[depth] += child_c;
            i++;
            /* skip optional internal node label */
            while(nwk[i] && nwk[i] != '(' && nwk[i] != ')' &&
                  nwk[i] != ',' && nwk[i] != ':' && nwk[i] != ';') i++;
            }
        else if(nwk[i] == ',')
            { i++; }
        else if(nwk[i] == ':')
            { while(nwk[i] && nwk[i] != ',' && nwk[i] != ')' && nwk[i] != ';') i++; }
        else
            {  /* taxon index at leaf — skip token, count one leaf */
            while(nwk[i] && nwk[i] != '(' && nwk[i] != ')' &&
                  nwk[i] != ',' && nwk[i] != ':') i++;
            cnt_stack[depth]++;
            }
        }
    return cherries;
    }

/* ---- Exact Bryant & Steel normalising constant ----------------------------
 * Computes log Z_T(beta) exactly via subset enumeration over B(T):
 *
 *   Z_T = exp(-2*beta*(n-3)) * sum_{S ⊆ B(T)} (exp(2*beta)-1)^|S| * N(S)
 *
 * N(S) = product over connected components C of the contracted-edge subgraph
 *        (edges of B(T) NOT in S) of (2*|C|-1)!!
 * where |C| = number of internal nodes merged into the super-node for C.
 *
 * Component sizes come from Union-Find: contracting edge (u,v) merges the
 * two internal node groups.  A group of k merged nodes has degree k+2 in the
 * partial tree, contributing (2k-1)!! binary resolutions.
 *
 * Limited to pruned trees with at most NORMCORR_EXACT_MAX_N taxa (2^(n-3)
 * subsets; for n=20 that is 2^17 = 131 072 iterations).  Larger trees fall
 * back to the large-beta approximation automatically.
 * -------------------------------------------------------------------------- */
#define NORMCORR_EXACT_MAX_N 20

/* Globals used by the recursive Newick parser (not thread-safe, but the
 * normcorr path in usertrees is single-threaded). */
static int _ncp_nc, _ncp_ec;
static int _ncp_eu[NORMCORR_EXACT_MAX_N], _ncp_ev[NORMCORR_EXACT_MAX_N];

/* Recursive descent: returns the internal node ID for this subtree, or -1
 * for a leaf.  Populates _ncp_eu/_ncp_ev for each internal edge found.     */
static int _ncp_parse(const char *nwk, int *p)
    {
    if(nwk[*p] != '(')
        {
        /* leaf — skip token and optional branch length */
        while(nwk[*p] && nwk[*p] != ',' && nwk[*p] != ')' && nwk[*p] != ';')
            (*p)++;
        return -1;
        }
    (*p)++;                           /* skip '(' */
    int my_id = _ncp_nc++;
    while(nwk[*p] != ')' && nwk[*p] && nwk[*p] != ';')
        {
        if(nwk[*p] == ',') { (*p)++; continue; }
        int child = _ncp_parse(nwk, p);
        if(child >= 0)                /* internal child → record edge */
            { _ncp_eu[_ncp_ec] = my_id; _ncp_ev[_ncp_ec] = child; _ncp_ec++; }
        }
    (*p)++;                           /* skip ')' */
    /* skip optional label */
    while(nwk[*p] && nwk[*p] != ',' && nwk[*p] != ')' && nwk[*p] != ';' && nwk[*p] != ':')
        (*p)++;
    /* skip optional branch length */
    if(nwk[*p] == ':')
        { (*p)++; while(nwk[*p] && nwk[*p] != ',' && nwk[*p] != ')' && nwk[*p] != ';') (*p)++; }
    return my_id;
    }

static double log_normconst_exact(const char *pruned_nwk, int n, double beta)
    {
    int m  = n - 3;   /* internal bipartitions */
    int nn = n - 2;   /* internal nodes        */
    int e, i;

    /* Parse tree topology */
    _ncp_nc = 0; _ncp_ec = 0;
    int pos = 0;
    _ncp_parse(pruned_nwk, &pos);
    /* _ncp_ec should equal m; _ncp_nc should equal nn */

    double x     = exp(2.0 * beta) - 1.0;
    double ex    = exp(-2.0 * beta * (double)m);

    /* Precompute powers of x */
    double pow_x[NORMCORR_EXACT_MAX_N + 1];
    pow_x[0] = 1.0;
    for(i = 1; i <= m; i++) pow_x[i] = pow_x[i-1] * x;

    double sum = 0.0;
    unsigned int total = 1u << m;

    for(unsigned int S = 0; S < total; S++)
        {
        /* Union-Find over nn internal nodes; contract edges NOT in S */
        int uf[NORMCORR_EXACT_MAX_N], sz[NORMCORR_EXACT_MAX_N];
        for(i = 0; i < nn; i++) { uf[i] = i; sz[i] = 1; }

        for(e = 0; e < m; e++)
            {
            if(S & (1u << e)) continue;   /* edge in S = kept, not contracted */
            int u = _ncp_eu[e], v = _ncp_ev[e];
            while(uf[u] != u) u = uf[u];
            while(uf[v] != v) v = uf[v];
            if(u != v)
                {
                if(sz[u] < sz[v]) { int t = u; u = v; v = t; }
                uf[v] = u; sz[u] += sz[v];
                }
            }

        /* N(S) = prod over components with size k>=2 of (2k-1)!! */
        double N = 1.0;
        for(i = 0; i < nn; i++)
            {
            if(uf[i] == i && sz[i] >= 2)
                {
                double df = 1.0;
                int f;
                for(f = 1; f <= 2*sz[i]-1; f += 2) df *= (double)f;
                N *= df;
                }
            }

        /* popcount(S) */
        int nbits = 0;
        unsigned int tmp = S;
        while(tmp) { nbits += (int)(tmp & 1u); tmp >>= 1; }

        sum += pow_x[nbits] * N;
        }

    return log(ex * sum);
    }

/* log_normconst_approx: approximate log(Z_T) for the Bryant & Steel (2008)
 * normalising constant using the truncated large-beta expansion:
 *
 *   Z_T = sum_m b_m(T) * e^{-beta*m}
 *       ≈ 1  +  b_2 * eps  +  b_4(c_T) * eps^2
 *
 * where
 *   eps   = e^{-2*beta}
 *   b_2   = 2*(n-3)                            [same for all binary n-leaf trees]
 *   b_4   = 4*C(n-3,2) + 6*(n-6+c_T)          [depends on cherry count c_T]
 *   c_T   = number of cherry pairs in T|X_i
 *
 * Accuracy: good for beta > ~1.5 (eps < 0.05).  Returns 0 for n < 4.
 *
 * Reference: Bryant & Steel (2008), arXiv:0810.0868, Section 3.
 * b_2 from their §3.2; b_4 from their Theorem 2.26 / equation after b_4(T). */
static double log_normconst_approx(int n, int cherries, double beta)
    {
    if(n < 4) return 0.0;
    double eps  = exp(-2.0 * beta);
    double eps2 = eps * eps;
    double b2   = 2.0 * (n - 3);
    double b4   = 4.0 * (n - 3) * (n - 4) * 0.5 + 6.0 * (n - 6 + cherries);
    double Z    = 1.0 + b2 * eps + b4 * eps2;
    if(Z < 1.0) Z = 1.0;   /* safety: b_4 can be slightly negative for n=4,5 */
    return log(Z);
    }

/* collect_biparts_newick_sz: like collect_biparts_newick but also fills
 * sizes_out[k] = min(|A_k|, |B_k|) (smaller-side leaf count) for each
 * non-trivial bipartition k.  ntaxa is the total leaf count in this tree.
 * supports_out[k]: if non-NULL, filled with the internal node label parsed
 * as a BS weight in [0,1] (labels in (1,100] are divided by 100); set to
 * 1.0f when no label is present.
 * out[], sizes_out[], and supports_out[] are co-sorted by hash on return. */
static int collect_biparts_newick_sz(const char *nwk, uint64_t total_hash,
                                     int ntaxa, uint64_t *out, int *sizes_out,
                                     float *supports_out)
    {
    uint64_t hash_stack[2 * NAME_LENGTH + 4];
    int      cnt_stack [2 * NAME_LENGTH + 4];
    int depth = 0, cnt = 0, i = 0;
    hash_stack[0] = 0;
    cnt_stack[0]  = 0;
    while(nwk[i] && nwk[i] != ';')
        {
        if(nwk[i] == '(')
            { hash_stack[++depth] = 0; cnt_stack[depth] = 0; i++; }
        else if(nwk[i] == ')')
            {
            uint64_t child_h = hash_stack[depth--];
            int      child_c = cnt_stack[depth + 1];
            uint64_t comp    = total_hash ^ child_h;
            uint64_t bh      = (child_h < comp) ? child_h : comp;
            int stored = 0;
            int other_c    = ntaxa - child_c;
            int min_side   = (child_c < other_c) ? child_c : other_c;
            /* Skip trivial bipartitions: bh==0 is the empty/full split; min_side<2
             * is a pendant (leaf) edge — non-informative for RF. A pendant split
             * can appear as a spurious internal edge when a source tree is stored
             * with a redundant degree-2 root (extra outer parentheses), e.g.
             * '(((A,(B,C)),D));'. Counting it inflates the source-tree bipartition
             * count and makes an otherwise-displayed tree look like a conflict,
             * producing scores inconsistent with the exhaustive/usertrees path. */
            if(bh != 0 && min_side >= 2)
                {
                out[cnt]       = bh;
                sizes_out[cnt] = min_side;
                cnt++;
                stored = 1;
                }
            hash_stack[depth] ^= child_h;
            cnt_stack[depth]  += child_c;
            i++;
            /* skip/capture optional internal node label (e.g. bootstrap support) */
            {
            char lbl[64]; int lj = 0;
            while(nwk[i] && nwk[i] != '(' && nwk[i] != ')' &&
                  nwk[i] != ',' && nwk[i] != ':' && nwk[i] != ';')
                lbl[lj++] = nwk[i++];
            if(supports_out && stored)
                {
                if(lj > 0)
                    {
                    lbl[lj] = '\0';
                    float bs = (float)atof(lbl);
                    supports_out[cnt-1] = (bs > 1.0f) ? bs / 100.0f : bs;
                    }
                else
                    supports_out[cnt-1] = 1.0f;
                }
            }
            }
        else if(nwk[i] == ',')
            { i++; }
        else if(nwk[i] == ':')
            { while(nwk[i] && nwk[i] != ',' && nwk[i] != ')' && nwk[i] != ';') i++; }
        else
            {  /* taxon integer index — only valid at leaf position (after '(' or ',') */
            char num[64]; int j = 0;
            while(nwk[i] && nwk[i] != '(' && nwk[i] != ')' && nwk[i] != ',' && nwk[i] != ':')
                num[j++] = nwk[i++];
            num[j] = '\0';
            int tidx = atoi(num);
            if(tidx >= 0 && tidx < number_of_taxa)
                { hash_stack[depth] ^= taxon_hash_vals[tidx]; cnt_stack[depth]++; }
            }
        }
    /* Co-sort (out, sizes_out, supports_out) by hash — insertion sort, n ≤ number_of_taxa */
    for(i = 1; i < cnt; i++)
        {
        uint64_t kh = out[i]; int ks = sizes_out[i];
        float    kp = supports_out ? supports_out[i] : 1.0f;
        int j = i - 1;
        while(j >= 0 && out[j] > kh)
            {
            out[j+1]  = out[j];  sizes_out[j+1] = sizes_out[j];
            if(supports_out) supports_out[j+1] = supports_out[j];
            j--;
            }
        out[j+1] = kh; sizes_out[j+1] = ks;
        if(supports_out) supports_out[j+1] = kp;
        }
    return cnt;
    }

/* quartet_intersect_score: merge-walk over sorted supertree hashes (a) and
 * sorted source-tree hashes (b, with parallel sizes sb).
 * Returns the sum of C(s,2)*C(n-s,2) over shared splits (= agreeing quartets).
 * Sets *q_total_out to the same sum over all source-tree splits. */
static int64_t quartet_intersect_score(const uint64_t *a, int na,
                                        const uint64_t *b, const int *sb, int nb,
                                        int ntaxa_i, int64_t *q_total_out)
    {
    int64_t q_agree = 0, q_total = 0;
    int i, j;
    for(j = 0; j < nb; j++)
        {
        int64_t s = sb[j], n = ntaxa_i;
        q_total += s*(s-1)/2 * (n-s)*(n-s-1)/2;
        }
    i = 0; j = 0;
    while(i < na && j < nb)
        {
        if(a[i] == b[j])
            {
            int64_t s = sb[j], n = ntaxa_i;
            q_agree += s*(s-1)/2 * (n-s)*(n-s-1)/2;
            i++; j++;
            }
        else if(a[i] < b[j]) i++;
        else                  j++;
        }
    *q_total_out = q_total;
    return q_agree;
    }

/* quartet_intersect_score_w: BS-weighted version of quartet_intersect_score.
 * Each split's quartet term is multiplied by its bootstrap support weight.
 * Uses double accumulators to handle fractional weights correctly.
 * Sets *q_total_out to the weighted sum over all source-tree splits. */
static double quartet_intersect_score_w(const uint64_t *a, int na,
                                         const uint64_t *b, const int *sb,
                                         const float *pb, int nb,
                                         int ntaxa_i, double *q_total_out)
    {
    double q_agree = 0.0, q_total = 0.0;
    int i, j;
    for(j = 0; j < nb; j++)
        {
        double s = sb[j], n = ntaxa_i;
        q_total += pb[j] * (s*(s-1)/2.0) * ((n-s)*(n-s-1)/2.0);
        }
    i = 0; j = 0;
    while(i < na && j < nb)
        {
        if(a[i] == b[j])
            {
            double s = sb[j], n = ntaxa_i;
            q_agree += pb[j] * (s*(s-1)/2.0) * ((n-s)*(n-s-1)/2.0);
            i++; j++;
            }
        else if(a[i] < b[j]) i++;
        else                  j++;
        }
    *q_total_out = q_total;
    return q_agree;
    }

/* ensure_fund_scores_alloc: lazily allocate the fund_scores path-distance
 * matrices (int[Total_fund_trees][N][N], zeroed). This array is O(trees*taxa^2)
 * and is only ever read by the distance-fit criterion (compare_trees), so it is
 * no longer allocated at load time. It is allocated on demand by the operations
 * that need it (dfit scoring via cal_fund_scores, and the boot/yaptp/
 * generatetrees/spr_dist snapshot bookkeeping). For large tree sets run under
 * the ML/RF/etc. criteria this avoids the whole array (e.g. ~110 GB for 121k
 * trees x 477 taxa). Idempotent: a no-op if already allocated. */
void ensure_fund_scores_alloc(void)
    {
    int i, j, k;
    if(fund_scores != NULL) return;
    fund_scores = malloc(Total_fund_trees * sizeof(int **));
    if(fund_scores == NULL) memory_error(35);
    for(i = 0; i < Total_fund_trees; i++)
        {
        fund_scores[i] = malloc(number_of_taxa * sizeof(int *));
        if(fund_scores[i] == NULL) memory_error(36);
        for(j = 0; j < number_of_taxa; j++)
            {
            fund_scores[i][j] = malloc(number_of_taxa * sizeof(int));
            if(fund_scores[i][j] == NULL) memory_error(37);
            for(k = 0; k < number_of_taxa; k++)
                fund_scores[i][j][k] = 0;
            }
        }
    }

void cal_fund_scores(int printfundscores)
    {
    int i =0, j=0, k=0;
    FILE *dists = NULL;

	ensure_fund_scores_alloc();   /* fund_scores is not allocated at load; allocate on first dfit use */
	calculated_fund_scores = TRUE;
    if(printfundscores)
            {
            dists = fopen("source_matrices.txt", "w");
            for(i=0; i<number_of_taxa; i++)
                fprintf(dists, "%s\t", taxa_names[i]);
            fprintf(dists, "\n\n");
            }
    for(i=0; i<Total_fund_trees; i++)
        {
        pathmetric(fundamentals[i], fund_scores[i]);
        if(printfundscores)  /* print the distance matrices for each fundamental tree */
            {
            
            for(j=0; j<number_of_taxa; j++)
                {
                fprintf(dists, "%s\t", taxa_names[j]);
                if(presence_of_taxa[i][j] > 0)
                    {
                    for(k=0; k<number_of_taxa; k++)
                        {
                        if(presence_of_taxa[i][k] > 0)
                            {
                            /** print out the distances **/
                            fprintf(dists, "%d\t", fund_scores[i][j][k]);
                            }
                        else
                            fprintf(dists, "\t");
                        }
                    fprintf(dists, "\n");
                    }
                else
                    {
                    fprintf(dists, "\n");
                    }
                }
            fprintf(dists, "\n\n\n");
            }
        
        
        }
    if(dists != NULL) fclose(dists);
    }

void pathmetric(char *string, int **scores)
    {
    int i=0, j=0, charactercount = -1, **closeP = NULL, variable = 0; 
    char number[30];


     
    /* The array characters is used to keep track, for each taxa, the open and closed brackets that has followed each */
    closeP = malloc((number_of_taxa)*(sizeof(int*)));
    for(i=0; i<number_of_taxa; i++)
        {
        closeP[i] = malloc(3*sizeof(int));
        closeP[i][0] = 0;  /* the number of open parentheses */
        closeP[i][1] = 0;  /* the number of close parentheses */
        closeP[i][2] = FALSE;  /* whether this taxa has been found yet */
        }

    unroottree(string);
    i=0;
    while(string[i] != ';')  /* until the end of the string */
        {
        switch(string[i])
            {
            case '(':	
                        for(j=0; j<number_of_taxa; j++)
                            {
                            if(closeP[j][2] == TRUE)
                                {
                                closeP[j][0] ++;
                                }
                            }
                        i++;
                        break;
            case ')':
                        for(j=0; j<number_of_taxa; j++)
                            {
                            if(closeP[j][2] == TRUE)
                                {
                                if(closeP[j][0] > 0)
                                    closeP[j][0]--;  /* if this close parenthesis cancels out a previously counted openparentheis */
                                else
                                    closeP[j][1]++;  
                                }
                            }
						i++;
						while(string[i] != ')' && string[i] != '(' && string[i] != ',' && string[i] != ';' && string[i] != ':')  /* skip the labels (if they are there) */
							i++;
                        break;
            case ',':
                        i++;
                        break;
            case ':':
                        while(string[i] != ')' && string[i] != '(' && string[i] != ',' && string[i] != ';')  /* skip the weights (if they are there) */
                            i++;
                        break;
            default:
                        /* this has to be a taxa number */
                        for(j=0; j<30; j++) number[j] = '\0';  /* initialise the array to hold the number in text form */
                        j=0;
                        while(string[i] != '(' && string[i] != ')' && string[i] != ',' && string[i] != ':')
                            {
                            number[j] = string[i];
                            i++; j++;
                            }
                        /* now need to change the string to an integer number */
                        charactercount = j;
                        variable = 0;
                        for(j=charactercount; j>0; j--)
                            {
                            variable += ( pow(10, (charactercount-j)) * (texttoint(number[j-1])));
                            }
                        closeP[variable][2] = TRUE;  /* tell the program that this caracter has now been passed and is to be counted from now on */
                        /* now need to assign the distance from any taxa we passed previously to this taxa */
                        for(j=0; j<number_of_taxa; j++)
                            {
                            if(closeP[j][2] == TRUE && j != variable) /* we have passed this taxa earlier */
                                {
                                scores[j][variable] = closeP[j][0]+closeP[j][1]+1;  /* the score is equal to the number of open parentheses not cancelled out plus the number of close parentheses not cancelled out + 1 */
                                scores[variable][j] = scores[j][variable];
                                }
                            }
           
                        break;
            }
            
        }
    for(i=0; i<number_of_taxa; i++)
        free(closeP[i]);
    free(closeP);
    closeP = NULL;
    }

void pathmetric_internals(char *string, struct taxon * species_tree, int **scores)
    {
    int i=0, j=0, charactercount = -1, **closeP = NULL, variable = 0, **within = NULL, *presence = NULL; 
    char number[30];

	/* WE need a matrix telling us which taxa anf branches are within other branches */
	within = malloc((2*number_of_taxa)*sizeof(int *));
	presence = malloc((2*number_of_taxa)*sizeof(int));
	for(i=0; i<(2*number_of_taxa); i++)
		{
		within[i] = malloc((2*number_of_taxa)*sizeof(int));
		presence[i] = FALSE;
		for(j=0; j<(2*number_of_taxa); j++)
			within[i][j] = 0;
		}
     i=0;
	 
	 calculate_withins(species_tree, within, presence);
	
	 
    /* The array characters is used to keep track, for each taxa, the open and closed brackets that has followed each */
    closeP = malloc((2*number_of_taxa)*(sizeof(int*)));
    for(i=0; i<2*number_of_taxa; i++)
        {
        closeP[i] = malloc(3*sizeof(int));
        closeP[i][0] = 0;  /* the number of open parentheses */
        closeP[i][1] = 0;  /* the number of close parentheses */
        closeP[i][2] = FALSE;  /* whether this taxa has been found yet */
        }

    unroottree(string);
    i=0;
    while(string[i] != ';')  /* until the end of the string */
        {
        switch(string[i])
            {
            case '(':	
                        for(j=0; j<2*number_of_taxa; j++)
                            {
                            if(closeP[j][2] == TRUE)
                                {
                                closeP[j][0] ++;
                                }
                            }
                        i++;
                        break;
            case ')':
                        for(j=0; j<2*number_of_taxa; j++)
                            {
                            if(closeP[j][2] == TRUE)
                                {
                                if(closeP[j][0] > 0)
                                    closeP[j][0]--;  /* if this close parenthesis cancels out a previously counted openparentheis */
                                else
                                    closeP[j][1]++;  
                                }
                            }
						i++;
						
						/* read the label on the internal branch */
                        for(j=0; j<30; j++) number[j] = '\0';  /* initialise the array to hold the number in text form */
                        j=0;
                        while(string[i] != '(' && string[i] != ')' && string[i] != ',' && string[i] != ':' && string[i] != ';')
                            {
                            number[j] = string[i];
                            i++; j++;
                            }
                        /* now need to change the string to an integer number */
                        charactercount = j;
                        variable = 0;
                        for(j=charactercount; j>0; j--)
                            {
                            variable += ( pow(10, (charactercount-j)) * (texttoint(number[j-1])));
                            }
                        closeP[variable][2] = TRUE;  /* tell the program that this caracter has now been passed and is to be counted from now on */
                        /* now need to assign the distance from any taxa we passed previously to this taxa */
                        for(j=0; j<2*number_of_taxa; j++)
                            {
                            if(closeP[j][2] == TRUE && j != variable) /* we have passed this taxa earlier */
                                {
								if(within[j][variable])
									scores[j][variable] = closeP[j][0]+closeP[j][1];
								else
									scores[j][variable] = closeP[j][0]+closeP[j][1]+1;  /* the score is equal to the number of open parentheses not cancelled out plus the number of close parentheses not cancelled out + 1 */
                                scores[variable][j] = scores[j][variable];
                                }
                            }
						
                        break;
            case ',':
                        i++;
                        break;
            case ':':
                        while(string[i] != ')' && string[i] != '(' && string[i] != ',' && string[i] != ';')  /* skip the weights (if they are there) */
                            i++;
                        break;
            default:
                        /* this has to be a taxa number */
                        for(j=0; j<30; j++) number[j] = '\0';  /* initialise the array to hold the number in text form */
                        j=0;
                        while(string[i] != '(' && string[i] != ')' && string[i] != ',' && string[i] != ':')
                            {
                            number[j] = string[i];
                            i++; j++;
                            }
                        /* now need to change the string to an integer number */
                        charactercount = j;
                        variable = 0;
                        for(j=charactercount; j>0; j--)
                            {
                            variable += ( pow(10, (charactercount-j)) * (texttoint(number[j-1])));
                            }
                        closeP[variable][2] = TRUE;  /* tell the program that this caracter has now been passed and is to be counted from now on */
                        /* now need to assign the distance from any taxa we passed previously to this taxa */
                        for(j=0; j<2*number_of_taxa; j++)
                            {
                            if(closeP[j][2] == TRUE && j != variable) /* we have passed this taxa earlier */
                                {
								if(within[j][variable])
									scores[j][variable] = closeP[j][0]+closeP[j][1];
								else
									scores[j][variable] = closeP[j][0]+closeP[j][1]+1;  /* the score is equal to the number of open parentheses not cancelled out plus the number of close parentheses not cancelled out + 1 */
                                scores[variable][j] = scores[j][variable];
                                }
                            }
           
                        break;
            }
            
        }
    for(i=0; i<2*number_of_taxa; i++)
		{
        free(closeP[i]);
		free(within[i]);
		}
	free(within);
	free(presence);
	free(closeP);
    closeP = NULL;
    }

void calculate_withins(struct taxon *position, int **within, int *presence)
	{
	int i=0;
	while(position != NULL)
		{
		for(i=0; i<2*number_of_taxa; i++)
			{
			if(presence[i] == TRUE)
				{
				within[i][position->tag] = TRUE;
				within[position->tag][i] = TRUE;
				}
			}
		if(position->daughter != NULL)
			{
			presence[position->tag] = TRUE;
			calculate_withins(position->daughter, within, presence);
			presence[position->tag] = FALSE;
			}
		position = position->next_sibling;
		}
	}

void weighted_pathmetric(char *string, float **scores, int fund_num)
    {
    int i=0, j=0, k=0, charactercount = -1, variable = 0, node_number = -1, open = 0; 
    char number[30];
    float *taxa_weights = NULL, *node_weights = NULL,  **closeP = NULL;

     
    /* The array characters is used to keep track, for each taxa, the open and closed brackets that has followed each */
    
	closeP = malloc((number_of_taxa)*(sizeof(float*)));
    if(!closeP) memory_error(90);
    for(i=0; i<number_of_taxa; i++)
        {
        closeP[i] = malloc(3*sizeof(float));
        if(!closeP[i]) memory_error(91);
        closeP[i][0] = 0;  /* the number of open parentheses */
        closeP[i][1] = 0;  /* the number of close parentheses */
        closeP[i][2] = FALSE;  /* whether this taxa has been found yet */
        }

    taxa_weights = malloc(number_of_taxa * sizeof(float));
        if(!taxa_weights) memory_error(92);
    node_weights = malloc(number_of_taxa * sizeof(float));
        if(!node_weights) memory_error(93);
    for(i=0; i<number_of_taxa; i++)
        {
        taxa_weights[i] = 1;
        node_weights[i] = 1;
        }
    unroottree(string);
    i=0;
    while(string[i] != ';')  /* until the end of the string */
        {
        switch(string[i])
            {
            case '(':	
                        if(i != 0)
                            {
                            k=i;
                            open = 1;
                            do
                                {
                                k++;
                                if(string[k] == '(') open++;
                                if(string[k] == ')') open--;
                                
                                }while(string[k] != ')' || open != 0);
							 node_number++;
							if(string[k+1] == ':')
								{
								k = k+2;  /* skip over the ':' to the start of the number */
								for(j=0; j<30; j++) number[j] = '\0';  /* initialise the array to hold the number in text form */
								j=0;
								while(string[k] != ')' && string[k] != ',' && string[k] != '(' && string[k] != ';')
									{
									number[j] = string[k];
									k++; j++;
									}  /* read in the weight */
								 node_weights[node_number] = tofloat(number);
								} /* otherwise the software will assume branch lengths of 1 for all */
                           
                           
                            for(j=0; j<number_of_taxa; j++)
                                {
                                if(closeP[j][2] == TRUE)
                                    {
                                    closeP[j][0] += node_weights[node_number];
                                    }
                                }
                            }
                        i++;
                        break;
            case ')':	
                        if(string[i+1] != ';')
                            {
                            i++;
                            while(string[i] != ')' && string[i] != '(' && string[i] != ',' && string[i] != ';')
                                i++;
                            for(j=0; j<number_of_taxa; j++)
                                {
                                if(closeP[j][2] == TRUE)
                                    {
                                    if(closeP[j][0] > 0)
                                        closeP[j][0] -= node_weights[node_number];  /* if this close parenthesis cancels out a previously counted openparentheis */
                                    else
                                        closeP[j][1] += node_weights[node_number];  
                                    }
                                }
                            node_weights[node_number] = 0;
                            node_number = 0;
                            }
                        else
                            i++;
                        break;
            case ',':
                        i++;
                        break;
            default:
                        /* this has to be a taxa number */
                        for(j=0; j<30; j++) number[j] = '\0';  /* initialise the array to hold the number in text form */ 
                        j=0;
                        while(string[i] != ':' && string[i] != ')' && string[i] != '(' && string[i] != ',')
                            {
                            number[j] = string[i];
                            i++; j++;
                            }

                        /* now need to change the string to an integer number */
                        charactercount = j;
                        variable = 0;
                        for(j=charactercount; j>0; j--)
                            {
                            variable += ( pow(10, (charactercount-j)) * (texttoint(number[j-1])));
                            }
                        closeP[variable][2] = TRUE;  /* tell the program that this caracter has now been passed and is to be counted from now on */
                        /* Now read in the weight for this taxa */
						if(string[i] == ':')
							{
							for(j=0; j<30; j++) number[j] = '\0';  /* initialise the array to hold the number in text form */
							j=0; i++;
							while(string[i] != ')' && string[i] != ',' && string[i] != '(' && string[i] != ';')
								{
								number[j] = string[i];
								i++; j++;
								}
							taxa_weights[variable] = tofloat(number);
							}
						
                        /* now need to assign the distance from any taxa we passed previously to this taxa */
                        
                        for(j=0; j<number_of_taxa; j++)
                            {
                            if(closeP[j][2] == TRUE && j != variable) /* we have passed this taxa earlier */
                                {
                                scores[j][variable] += (closeP[j][0]+closeP[j][1]+taxa_weights[variable]+taxa_weights[j])*tree_weights[fund_num];  /* the score is equal to the number of open parentheses not cancelled out plus the number of close parentheses not cancelled out + 1 */
                                scores[variable][j] = scores[j][variable];
                                }
                            }
           
                        break;
            }
            
    
        }
		
    for(i=0; i<number_of_taxa; i++)
        free(closeP[i]);
    free(closeP);
    closeP = NULL;
	free(taxa_weights);
	free(node_weights);

    }

float compare_trees(int spr)
    {
    int i=0, j=0, k=0, found = FALSE, here1 = FALSE, here2 = FALSE;
    float total =0, temp = 0;
    char *pruned_tree, *tmp;
    /********** allocate dynamic arrays **************/
    pruned_tree = malloc(TREE_LENGTH*sizeof(char));
    if(!pruned_tree)  memory_error(29);
        
    pruned_tree[0] = '\0';
    
    tmp = malloc(TREE_LENGTH*sizeof(char));
    if(!tmp)  memory_error(30);
    
    tmp[0] = '\0';
    
    
    
    /******* Next score this supertree against every fundamental tree in memory *********/
    for(i=0; i< Total_fund_trees; i++)
        {
        /* Allow Ctrl+C to interrupt long scoring loops; flush ensures the flag
         * is visible to worker threads on all architectures (e.g. ARM64). */
#ifdef _OPENMP
        #pragma omp flush(user_break)
#endif
        if(user_break) break;
		if(sourcetreetag[i]) /* if this sourcetree is to be used in the analysis ---- defined by the exclude command */
			{
			
			found = FALSE; here1 = TRUE; here2 = TRUE;
			
			/*** Next check to see if the pruning used affected this fundamental tree ***/
			if(sourcetree_scores[i] != -1 && spr)
				{
				
				/** check to see if any of the taxa from this source tree were in the subtree that was pruned **/
				for(j=0; j<number_of_taxa; j++)
					{
					if(presence_of_taxa[i][j] > 0  && presenceof_SPRtaxa[j] == TRUE) here1 = FALSE;
					if(presence_of_taxa[i][j] > 0 && presenceof_SPRtaxa[j] == FALSE) here2 = FALSE;
					}
				if(here1 || here2) found = TRUE;   /* This means that all the taxa from this source tree are all either in the main part of the supertree or all in the pruned subtree */
				}
			if(!found || !spr)
				{
				prune_tree(tree_top, i);  /* Prune the supertree so that it has the same taxa as the fundamental tree i */
				shrink_tree(tree_top);    /* Shrink the pruned tree by switching off any internal nodes that are not needed */
				
				pruned_tree[0] = '\0'; /* initialise the string */
				if(print_pruned_tree(tree_top, 0, pruned_tree, FALSE, 0) >1)
					{
					tmp[0] = '\0';
					strcpy(tmp, "(");
					strcat(tmp, pruned_tree);
					strcat(tmp, ")");
					strcpy(pruned_tree, tmp);
					}
				strcat(pruned_tree, ";");
				while(unroottree(pruned_tree));  /* remove bifurcating root from SPR-at-root grafts */
			
				/* initialise the super_scores array */
				for(j=0; j<number_of_taxa; j++)
					{
					for(k=j; k<number_of_taxa; k++)
						{
						super_scores[j][k] = 0;
						super_scores[k][j] = 0;
						}
					}

				pathmetric(pruned_tree, super_scores);  /*calculate the pathmetric for the pruned supertree */
				reset_tree(tree_top);  /* reset the supertree for the next comparison */
			
				temp = 0;
				/* score the differences between the two trees */
				for(j=0; j<number_of_taxa; j++)
					{
					for(k=(j+1); k<number_of_taxa; k++)
						{
						if(super_scores[j][k] != fund_scores[i][j][k])
							{
							if(super_scores[j][k] > fund_scores[i][j][k]) 
								{
								temp += super_scores[j][k] - fund_scores[i][j][k];
								}
							else
								{
								temp += fund_scores[i][j][k] - super_scores[j][k];
								}
							}
						}
					}
				if(dweight == 1) temp = temp/number_of_comparisons[i];  /* normalise the scores */
				
				/*** multiply the weight the tree has by the score */
				temp = temp*tree_weights[i];
				sourcetree_scores[i] = temp;
				
				total += temp;
				}
			else
				{
				total += sourcetree_scores[i];
				}
			}
        }
    free(tmp);
    free(pruned_tree);	
    return(total);
    }

void rf_precompute_fund_biparts(void)
    {
    int i, j;
    if(fund_bipart_sets != NULL)
        {
        for(i = 0; i < Total_fund_trees; i++)
            {
            free(fund_bipart_sets[i].hashes);
            free(fund_bipart_sets[i].sizes);
            free(fund_bipart_sets[i].supports);
            }
        free(fund_bipart_sets);
        }
    fund_bipart_sets = malloc(Total_fund_trees * sizeof(BipartSet));
    for(i = 0; i < Total_fund_trees; i++)
        {
        fund_bipart_sets[i].hashes   = malloc(number_of_taxa * sizeof(uint64_t));
        fund_bipart_sets[i].sizes    = malloc(number_of_taxa * sizeof(int));
        fund_bipart_sets[i].supports = malloc(number_of_taxa * sizeof(float));
        fund_bipart_sets[i].count    = 0;
        if(!sourcetreetag[i]) continue;
        int ntaxa_i = 0;
        uint64_t total_hash = 0;
        for(j = 0; j < number_of_taxa; j++)
            if(presence_of_taxa[i][j]) { ntaxa_i++; total_hash ^= taxon_hash_vals[j]; }
        char *tmp = malloc(strlen(fundamentals[i]) + 10);
        strcpy(tmp, fundamentals[i]);
        unroottree(tmp);
        fund_bipart_sets[i].count =
            collect_biparts_newick_sz(tmp, total_hash, ntaxa_i,
                                      fund_bipart_sets[i].hashes,
                                      fund_bipart_sets[i].sizes,
                                      fund_bipart_sets[i].supports);
        free(tmp);
        }
    }

/* -----------------------------------------------------------------------
 * compute_raw_rf_dists: raw (unnormalized, unweighted) RF distance from
 * tree_top to each source tree, stored in dists_out[i].
 * Requires rf_precompute_fund_biparts() to have been called.
 * dists_out must be allocated by caller (Total_fund_trees floats).
 * ----------------------------------------------------------------------- */
void compute_raw_rf_dists(float *dists_out)
    {
    int i, j;
    char *pruned_nwk = malloc(TREE_LENGTH * sizeof(char));
    if(!pruned_nwk) memory_error(90);
    char *tmp        = malloc(TREE_LENGTH * sizeof(char));
    if(!tmp) { free(pruned_nwk); memory_error(91); }
    uint64_t *super_bp = malloc(number_of_taxa * sizeof(uint64_t));
    if(!super_bp) { free(pruned_nwk); free(tmp); memory_error(92); }

    for(i = 0; i < Total_fund_trees; i++)
        {
        dists_out[i] = 0.0f;
        if(!sourcetreetag[i]) continue;

        int ntaxa_i = 0;
        uint64_t total_hash = 0;
        for(j = 0; j < number_of_taxa; j++)
            if(presence_of_taxa[i][j]) { ntaxa_i++; total_hash ^= taxon_hash_vals[j]; }
        if(ntaxa_i < 4) continue;

        prune_tree(tree_top, i);
        shrink_tree(tree_top);

        pruned_nwk[0] = '\0';
        if(print_pruned_tree(tree_top, 0, pruned_nwk, FALSE, 0) > 1)
            {
            tmp[0] = '\0';
            strcpy(tmp, "("); strcat(tmp, pruned_nwk); strcat(tmp, ")");
            strcpy(pruned_nwk, tmp);
            }
        strcat(pruned_nwk, ";");

        while(unroottree(pruned_nwk));
        int super_cnt = collect_biparts_newick(pruned_nwk, total_hash, super_bp);
        reset_tree(tree_top);

        int gene_cnt = fund_bipart_sets[i].count;
        int shared   = bipart_intersection_count(super_bp, super_cnt,
                                                 fund_bipart_sets[i].hashes, gene_cnt);
        dists_out[i] = (float)(super_cnt + gene_cnt - 2 * shared);
        }

    free(super_bp);
    free(tmp);
    free(pruned_nwk);
    }

/* -----------------------------------------------------------------------
 * compute_taxon_conflict: per-taxon conflict weight from MISSING source-tree
 * bipartitions (bipartition-level, not whole-tree).
 *
 * For each source tree that the current supertree (tree_top) does not perfectly
 * display (sourcetree_scores[i] > 0), enumerate that source tree's splits and,
 * for every split ABSENT from the induced supertree, add the source tree's
 * residual score to the taxa on the SMALLER side of the missing split. Only the
 * taxa actually involved in unsatisfied splits accrue weight — this pinpoints
 * the few genuinely-misplaced taxa rather than flagging every taxon in a
 * conflicting source tree.
 *
 * Requires: rf_precompute_fund_biparts() done, sourcetree_scores[] populated for
 * the current tree_top (call a compare_trees_*() first), tree_top = supertree.
 * Source trees carry integer taxon-id leaf labels, so membership is read
 * directly. Fills conflict_w[0..number_of_taxa-1]; returns count of nonzero.
 * ----------------------------------------------------------------------- */
int compute_taxon_conflict(float *conflict_w)
    {
    int i, j, cnt = 0, d;
    int maxdepth = 2 * number_of_taxa + 8;
    char     *pruned_nwk = malloc(TREE_LENGTH * sizeof(char));
    char     *tmp        = malloc(TREE_LENGTH * sizeof(char));
    uint64_t *super_bp   = malloc((number_of_taxa + 1) * sizeof(uint64_t));
    int      *tree_taxa  = malloc((number_of_taxa + 1) * sizeof(int));
    int     **memstack   = malloc(maxdepth * sizeof(int *));
    int      *memcount   = malloc(maxdepth * sizeof(int));

    for(j = 0; j < number_of_taxa; j++) conflict_w[j] = 0.0f;
    if(!pruned_nwk || !tmp || !super_bp || !tree_taxa || !memstack || !memcount)
        goto ctc_cleanup;
    for(d = 0; d < maxdepth; d++) memstack[d] = malloc((number_of_taxa + 1) * sizeof(int));

    for(i = 0; i < Total_fund_trees; i++)
        {
        if(!sourcetreetag[i]) continue;
        float resid = sourcetree_scores[i];
        if(resid <= 1e-6f) continue;   /* source tree already displayed — no conflict */

        int ntaxa_i = 0; uint64_t total_hash = 0;
        for(j = 0; j < number_of_taxa; j++)
            if(presence_of_taxa[i][j]) { tree_taxa[ntaxa_i++] = j; total_hash ^= taxon_hash_vals[j]; }
        if(ntaxa_i < 4) continue;

        /* Supertree bipartitions induced on this source tree's taxa. */
        prune_tree(tree_top, i);
        shrink_tree(tree_top);
        pruned_nwk[0] = '\0';
        if(print_pruned_tree(tree_top, 0, pruned_nwk, FALSE, 0) > 1)
            { tmp[0] = '\0'; strcpy(tmp, "("); strcat(tmp, pruned_nwk); strcat(tmp, ")"); strcpy(pruned_nwk, tmp); }
        strcat(pruned_nwk, ";");
        while(unroottree(pruned_nwk));
        int super_cnt = collect_biparts_newick(pruned_nwk, total_hash, super_bp);
        reset_tree(tree_top);

        /* Walk the source tree (integer labels), enumerating splits with membership. */
        const char *nwk = fundamentals[i];
        int p = 0; d = 0; memcount[0] = 0;
        while(nwk[p] && nwk[p] != ';')
            {
            char c = nwk[p];
            if(c == '(')
                { d++; if(d >= maxdepth) break; memcount[d] = 0; p++; }
            else if(c == ')')
                {
                int mc = memcount[d], k;
                uint64_t child_h = 0;
                for(k = 0; k < mc; k++) child_h ^= taxon_hash_vals[ memstack[d][k] ];
                uint64_t comp = total_hash ^ child_h;
                uint64_t bh   = (child_h < comp) ? child_h : comp;
                int other = ntaxa_i - mc;
                int minside = (mc < other) ? mc : other;
                if(bh != 0 && minside >= 2)   /* nontrivial split */
                    {
                    int found = 0, s;
                    for(s = 0; s < super_cnt; s++) if(super_bp[s] == bh) { found = 1; break; }
                    if(!found)   /* split present in source tree but MISSING from supertree */
                        {
                        if(mc <= other)                    /* members are the smaller side */
                            for(k = 0; k < mc; k++) conflict_w[ memstack[d][k] ] += resid;
                        else                               /* complement is the smaller side */
                            for(k = 0; k < ntaxa_i; k++)
                                {
                                int id = tree_taxa[k], inmem = 0, q;
                                for(q = 0; q < mc; q++) if(memstack[d][q] == id) { inmem = 1; break; }
                                if(!inmem) conflict_w[id] += resid;
                                }
                        }
                    }
                /* Pop: fold this subtree's members into the parent level. */
                if(d > 0)
                    { int par = d - 1, k2;
                      for(k2 = 0; k2 < mc; k2++)
                          if(memcount[par] < number_of_taxa) memstack[par][memcount[par]++] = memstack[d][k2]; }
                d--; if(d < 0) d = 0;
                p++;
                while(nwk[p] && nwk[p] != '(' && nwk[p] != ')' && nwk[p] != ',' && nwk[p] != ';') p++;
                }
            else if(c == ',') p++;
            else if(c == ':') { while(nwk[p] && nwk[p] != ',' && nwk[p] != ')' && nwk[p] != ';') p++; }
            else
                {  /* integer taxon-id leaf */
                char num[64]; int nj = 0;
                while(nwk[p] && nwk[p] != '(' && nwk[p] != ')' && nwk[p] != ',' && nwk[p] != ':') num[nj++] = nwk[p++];
                num[nj] = '\0';
                int id = atoi(num);
                if(id >= 0 && id < number_of_taxa && memcount[d] < number_of_taxa)
                    memstack[d][memcount[d]++] = id;
                }
            }
        }

    for(j = 0; j < number_of_taxa; j++) if(conflict_w[j] > 0.0f) cnt++;
    for(d = 0; d < maxdepth; d++) free(memstack[d]);

ctc_cleanup:
    free(memcount); free(memstack); free(tree_taxa);
    free(super_bp); free(tmp); free(pruned_nwk);
    return cnt;
    }

/* -----------------------------------------------------------------------
 * Smoothed search surrogate: transfer-distance ML.
 *
 * The RF-based ML score counts each source-tree split as present (0) or absent
 * (1). That is flat near the optimum — a supertree three splits from the true
 * tree scores the same, per split, as one thirty splits away — so the hill-climb
 * has no gradient to follow. compare_trees_transfer replaces the 0/1 term with
 * the TRANSFER INDEX of each source split: the minimum number of taxa that must
 * be moved for that split to appear in the induced supertree (Lemoine et al.
 * 2018). A split one taxon short scores 1, three taxa short scores 3, so the
 * surrogate gives a graded signal that funnels the search toward the right
 * basin. It shares the true optimum (every source split displayed => transfer 0
 * => score 0), so the search can optimise the surrogate and be re-scored with
 * true ML for the reported result.
 * ----------------------------------------------------------------------- */

/* enum_biparts_local: enumerate the nontrivial bipartitions of the integer-label
 * Newick nwk as bitsets over LOCAL taxon indices (g2l maps a global taxon id to
 * 0..ntaxa_i-1, or -1 if absent). Writes one bitset per branch into
 * bits[branch*nwords ...] (the branch's subtree side; transfer distance is
 * complement-invariant so either side works). Returns the branch count. */
static int enum_biparts_local(const char *nwk, const int *g2l, int ntaxa_i,
                              int nwords, uint64_t *bits)
    {
    int maxdepth = 2 * ntaxa_i + 8;
    uint64_t *ms = calloc((size_t)(maxdepth + 1) * nwords, sizeof(uint64_t));
    int cnt = 0, d = 0, p = 0, w;
    if(!ms) return 0;
    while(nwk[p] && nwk[p] != ';')
        {
        char c = nwk[p];
        if(c == '(')
            { d++; if(d > maxdepth) break; for(w = 0; w < nwords; w++) ms[d*nwords+w] = 0; p++; }
        else if(c == ')')
            {
            int pc = 0;
            for(w = 0; w < nwords; w++) pc += __builtin_popcountll(ms[d*nwords+w]);
            int other = ntaxa_i - pc, minside = (pc < other) ? pc : other;
            if(minside >= 2)
                { for(w = 0; w < nwords; w++) bits[cnt*nwords+w] = ms[d*nwords+w]; cnt++; }
            if(d > 0) for(w = 0; w < nwords; w++) ms[(d-1)*nwords+w] |= ms[d*nwords+w];
            d--; if(d < 0) d = 0;
            p++;
            while(nwk[p] && nwk[p] != '(' && nwk[p] != ')' && nwk[p] != ',' && nwk[p] != ';') p++;
            }
        else if(c == ',') p++;
        else if(c == ':') { while(nwk[p] && nwk[p] != ',' && nwk[p] != ')' && nwk[p] != ';') p++; }
        else
            {
            char num[64]; int nj = 0;
            while(nwk[p] && nwk[p] != '(' && nwk[p] != ')' && nwk[p] != ',' && nwk[p] != ':') num[nj++] = nwk[p++];
            num[nj] = '\0';
            int gid = atoi(num), lid = (gid >= 0 && gid < number_of_taxa) ? g2l[gid] : -1;
            if(lid >= 0) ms[d*nwords + (lid >> 6)] |= (1ULL << (lid & 63));
            }
        }
    free(ms);
    return cnt;
    }

float compare_trees_transfer(int spr)
    {
    (void)spr;   /* surrogate always recomputes fully */
    int i, j;
    float total = 0.0f;
    int nwmax = (number_of_taxa + 63) / 64;
    char     *pruned_nwk = malloc(TREE_LENGTH * sizeof(char));
    char     *tmp        = malloc(TREE_LENGTH * sizeof(char));
    int      *g2l        = malloc(number_of_taxa * sizeof(int));
    uint64_t *superbits  = malloc((size_t)(number_of_taxa + 1) * nwmax * sizeof(uint64_t));
    uint64_t *srcbits    = malloc((size_t)(number_of_taxa + 1) * nwmax * sizeof(uint64_t));
    if(!pruned_nwk || !tmp || !g2l || !superbits || !srcbits)
        { free(pruned_nwk); free(tmp); free(g2l); free(superbits); free(srcbits); return 0.0f; }

    for(i = 0; i < Total_fund_trees; i++)
        {
        if(user_break) break;
        if(!sourcetreetag[i]) continue;

        int ntaxa_i = 0;
        for(j = 0; j < number_of_taxa; j++) g2l[j] = -1;
        for(j = 0; j < number_of_taxa; j++) if(presence_of_taxa[i][j]) g2l[j] = ntaxa_i++;
        if(ntaxa_i < 4) { sourcetree_scores[i] = 0.0f; continue; }
        int nwords = (ntaxa_i + 63) / 64;

        /* supertree branches induced on this source tree's taxa */
        prune_tree(tree_top, i);
        shrink_tree(tree_top);
        pruned_nwk[0] = '\0';
        if(print_pruned_tree(tree_top, 0, pruned_nwk, FALSE, 0) > 1)
            { tmp[0] = '\0'; strcpy(tmp, "("); strcat(tmp, pruned_nwk); strcat(tmp, ")"); strcpy(pruned_nwk, tmp); }
        strcat(pruned_nwk, ";");
        while(unroottree(pruned_nwk));
        int ns = enum_biparts_local(pruned_nwk, g2l, ntaxa_i, nwords, superbits);
        reset_tree(tree_top);

        /* source-tree branches */
        int ng = enum_biparts_local(fundamentals[i], g2l, ntaxa_i, nwords, srcbits);

        /* For each source branch, its transfer index (min taxa-moves to appear
         * among the supertree branches). Exact matches (index 0) are the shared
         * splits, so we recover the true RF distance in the same pass:
         *   d_rf = ns + ng - 2*shared      (identical to compare_trees_ml)
         *   d_tr = sum of transfer indices (graded near-miss signal)
         * The search score keeps true RF as the PRIMARY term and uses transfer
         * only to break ties on RF plateaus (lexicographic via a large BIG),
         * so the true optimum and ranking are preserved — transfer supplies a
         * gradient only where RF is flat, without rewarding near-misses. */
        int shared = 0, a, b;
        double td = 0.0;
        for(a = 0; a < ng; a++)
            {
            int best = ntaxa_i;   /* upper bound */
            for(b = 0; b < ns; b++)
                {
                int hd = 0, w;
                for(w = 0; w < nwords; w++)
                    hd += __builtin_popcountll(srcbits[a*nwords+w] ^ superbits[b*nwords+w]);
                int t = (hd < ntaxa_i - hd) ? hd : ntaxa_i - hd;
                if(t < best) best = t;
                if(best == 0) break;
                }
            if(best == 0) shared++;
            td += (double)best;
            }
        double d_rf = (double)(ns + ng - 2 * shared);
        double BIG  = (double)Total_fund_trees * number_of_taxa * number_of_taxa + 1.0;
        double combined = d_rf * BIG + td;   /* RF dominates; transfer breaks ties */

        float sc = (float)(ml_beta * combined * (ml_scale == 1 ? LOG10E : 1.0f)) * tree_weights[i];
        sourcetree_scores[i] = sc;
        total += sc;
        }

    free(pruned_nwk); free(tmp); free(g2l); free(superbits); free(srcbits);
    return total;
    }

/* Search-time ML dispatcher: smoothed transfer surrogate when ml_smooth_search
 * is set (for guiding the hill-climb), else the true RF-based ML. The final
 * re-scoring and reporting always use compare_trees_ml() directly. */
float compare_trees_ml_search(int spr)
    {
    return ml_smooth_search ? compare_trees_transfer(spr) : compare_trees_ml(spr);
    }

float compare_trees_rf(int spr)
    {
    int i, j;
    float total = 0.0f;
    char *pruned_nwk    = malloc(TREE_LENGTH * sizeof(char));
    if(!pruned_nwk) memory_error(80);
    char *tmp           = malloc(TREE_LENGTH * sizeof(char));
    if(!tmp) { free(pruned_nwk); memory_error(81); }
    uint64_t *super_bp  = malloc(number_of_taxa * sizeof(uint64_t));
    if(!super_bp) { free(pruned_nwk); free(tmp); memory_error(82); }

    for(i = 0; i < Total_fund_trees; i++)
        {
#ifdef _OPENMP
        #pragma omp flush(user_break)
#endif
        if(user_break) break;
        if(!sourcetreetag[i]) continue;

        int found = FALSE, here1 = TRUE, here2 = TRUE;
        if(sourcetree_scores[i] != -1 && spr)
            {
            for(j = 0; j < number_of_taxa; j++)
                {
                if(presence_of_taxa[i][j] > 0 && presenceof_SPRtaxa[j] == TRUE)  here1 = FALSE;
                if(presence_of_taxa[i][j] > 0 && presenceof_SPRtaxa[j] == FALSE) here2 = FALSE;
                }
            if(here1 || here2) found = TRUE;
            }

        if(!found || !spr)
            {
            int ntaxa_i = 0;
            uint64_t total_hash = 0;
            for(j = 0; j < number_of_taxa; j++)
                if(presence_of_taxa[i][j]) { ntaxa_i++; total_hash ^= taxon_hash_vals[j]; }
            if(ntaxa_i < 4) { sourcetree_scores[i] = 0.0f; continue; }

            prune_tree(tree_top, i);
            shrink_tree(tree_top);

            pruned_nwk[0] = '\0';
            if(print_pruned_tree(tree_top, 0, pruned_nwk, FALSE, 0) > 1)
                {
                tmp[0] = '\0';
                strcpy(tmp, "("); strcat(tmp, pruned_nwk); strcat(tmp, ")");
                strcpy(pruned_nwk, tmp);
                }
            strcat(pruned_nwk, ";");

            while(unroottree(pruned_nwk));  /* remove bifurcating root from SPR-at-root grafts */
            int super_cnt = collect_biparts_newick(pruned_nwk, total_hash, super_bp);
            reset_tree(tree_top);

            int gene_cnt = fund_bipart_sets[i].count;
            int shared   = bipart_intersection_count(super_bp, super_cnt,
                                                     fund_bipart_sets[i].hashes, gene_cnt);
            float rf = (float)(super_cnt + gene_cnt - 2 * shared);
            int max_rf = 2 * (ntaxa_i - 3);
            if(max_rf > 0) rf /= (float)max_rf;
            rf *= tree_weights[i];
            sourcetree_scores[i] = rf;
            total += rf;
            }
        else
            {
            total += sourcetree_scores[i];
            }
        }

    free(super_bp);
    free(tmp);
    free(pruned_nwk);
    return total;
    }

float compare_trees_ml(int spr)
    {
    int i, j;
    float total = 0.0f;
    char *pruned_nwk   = malloc(TREE_LENGTH * sizeof(char));
    if(!pruned_nwk) memory_error(83);
    char *tmp          = malloc(TREE_LENGTH * sizeof(char));
    if(!tmp) { free(pruned_nwk); memory_error(84); }
    uint64_t *super_bp = malloc(number_of_taxa * sizeof(uint64_t));
    if(!super_bp) { free(pruned_nwk); free(tmp); memory_error(85); }

    for(i = 0; i < Total_fund_trees; i++)
        {
#ifdef _OPENMP
        #pragma omp flush(user_break)
#endif
        if(user_break) break;
        if(!sourcetreetag[i]) continue;

        int found = FALSE, here1 = TRUE, here2 = TRUE;
        if(sourcetree_scores[i] != -1 && spr)
            {
            for(j = 0; j < number_of_taxa; j++)
                {
                if(presence_of_taxa[i][j] > 0 && presenceof_SPRtaxa[j] == TRUE)  here1 = FALSE;
                if(presence_of_taxa[i][j] > 0 && presenceof_SPRtaxa[j] == FALSE) here2 = FALSE;
                }
            if(here1 || here2) found = TRUE;
            }

        if(!found || !spr)
            {
            int ntaxa_i = 0;
            uint64_t total_hash = 0;
            for(j = 0; j < number_of_taxa; j++)
                if(presence_of_taxa[i][j]) { ntaxa_i++; total_hash ^= taxon_hash_vals[j]; }
            if(ntaxa_i < 4) { sourcetree_scores[i] = 0.0f; continue; }

            prune_tree(tree_top, i);
            shrink_tree(tree_top);

            pruned_nwk[0] = '\0';
            if(print_pruned_tree(tree_top, 0, pruned_nwk, FALSE, 0) > 1)
                {
                tmp[0] = '\0';
                strcpy(tmp, "("); strcat(tmp, pruned_nwk); strcat(tmp, ")");
                strcpy(pruned_nwk, tmp);
                }
            strcat(pruned_nwk, ";");

            while(unroottree(pruned_nwk));  /* remove bifurcating root from SPR-at-root grafts */
            int super_cnt = collect_biparts_newick(pruned_nwk, total_hash, super_bp);

            /* Bryant & Steel (2008) normalising constant correction (usertrees normcorrect).
             * Computed here while pruned_nwk is still valid for T|X_i.
             * ml_do_normcorr==1: large-beta approx; ml_do_normcorr==2: exact subset enum. */
            if(ml_do_normcorr && ml_norm_logZ != NULL)
                {
                if(ml_do_normcorr == 2 && ntaxa_i <= NORMCORR_EXACT_MAX_N)
                    ml_norm_logZ[i] = log_normconst_exact(pruned_nwk, ntaxa_i, (double)ml_beta);
                else
                    {
                    if(ml_do_normcorr == 2)
                        {
                        static int _exact_warned = 0;
                        if(!_exact_warned)
                            { printf2("  Note: normcorrect=exact: trees with >%d taxa use large-beta approx\n", NORMCORR_EXACT_MAX_N); _exact_warned = 1; }
                        }
                    int nc = count_cherries_newick(pruned_nwk);
                    ml_norm_logZ[i] = log_normconst_approx(ntaxa_i, nc, (double)ml_beta);
                    }
                }

            reset_tree(tree_top);

            int gene_cnt = fund_bipart_sets[i].count;
            int shared   = bipart_intersection_count(super_bp, super_cnt,
                                                     fund_bipart_sets[i].hashes, gene_cnt);
            /* Raw RF distance — no normalisation */
            float d = (float)(super_cnt + gene_cnt - 2 * shared);
            /* [experimental] tree-size scaling: divide by k_i^eta (eta=0: Steel 2008, eta=1: normalised) */
            if(ml_eta != 0.0 && gene_cnt > 0)
                d /= (float)pow((double)gene_cnt, ml_eta);
            /* Negated log-likelihood: -ln L = beta * d  (Steel & Rodrigo 2008 formula).
             * mlscale=lust multiplies by log10(e) to match Akanni et al. 2014 (L.U.st). */
            float ml_score = (float)(ml_beta * d * (ml_scale == 1 ? LOG10E : 1.0f)) * tree_weights[i];
            sourcetree_scores[i] = ml_score;
            total += ml_score;
            }
        else
            {
            total += sourcetree_scores[i];
            }
        }

    free(super_bp);
    free(tmp);
    free(pruned_nwk);
    return total;
    }

/* -----------------------------------------------------------------------
 * compare_trees_sfit: splits-fit criterion.
 *   loss_i = (gene_cnt_i - shared_i) / gene_cnt_i
 * Fraction of source-tree splits absent from the pruned supertree.
 * Returns weighted sum of per-source-tree loss values (0 = perfect fit).
 * Requires rf_precompute_fund_biparts() to have been called.
 * ----------------------------------------------------------------------- */
float compare_trees_sfit(int spr)
    {
    int i, j;
    float total = 0.0f;
    char *pruned_nwk    = malloc(TREE_LENGTH * sizeof(char));
    if(!pruned_nwk) memory_error(86);
    char *tmp           = malloc(TREE_LENGTH * sizeof(char));
    if(!tmp) { free(pruned_nwk); memory_error(87); }
    uint64_t *super_bp  = malloc(number_of_taxa * sizeof(uint64_t));
    if(!super_bp) { free(pruned_nwk); free(tmp); memory_error(88); }

    for(i = 0; i < Total_fund_trees; i++)
        {
#ifdef _OPENMP
        #pragma omp flush(user_break)
#endif
        if(user_break) break;
        if(!sourcetreetag[i]) continue;

        int found = FALSE, here1 = TRUE, here2 = TRUE;
        if(sourcetree_scores[i] != -1 && spr)
            {
            for(j = 0; j < number_of_taxa; j++)
                {
                if(presence_of_taxa[i][j] > 0 && presenceof_SPRtaxa[j] == TRUE)  here1 = FALSE;
                if(presence_of_taxa[i][j] > 0 && presenceof_SPRtaxa[j] == FALSE) here2 = FALSE;
                }
            if(here1 || here2) found = TRUE;
            }

        if(!found || !spr)
            {
            int ntaxa_i = 0;
            uint64_t total_hash = 0;
            for(j = 0; j < number_of_taxa; j++)
                if(presence_of_taxa[i][j]) { ntaxa_i++; total_hash ^= taxon_hash_vals[j]; }
            if(ntaxa_i < 4) { sourcetree_scores[i] = 0.0f; continue; }

            prune_tree(tree_top, i);
            shrink_tree(tree_top);

            pruned_nwk[0] = '\0';
            if(print_pruned_tree(tree_top, 0, pruned_nwk, FALSE, 0) > 1)
                {
                tmp[0] = '\0';
                strcpy(tmp, "("); strcat(tmp, pruned_nwk); strcat(tmp, ")");
                strcpy(pruned_nwk, tmp);
                }
            strcat(pruned_nwk, ";");

            while(unroottree(pruned_nwk));
            int super_cnt = collect_biparts_newick(pruned_nwk, total_hash, super_bp);
            reset_tree(tree_top);

            int gene_cnt = fund_bipart_sets[i].count;
            float loss;
            if(bsweight && fund_bipart_sets[i].supports)
                {
                /* BS-weighted sfit: loss = (W_total - W_shared) / W_total
                 * where W = sum of per-split bootstrap support weights */
                float w_total = 0.0f, w_shared = 0.0f;
                int si = 0, gi = 0;
                for(gi = 0; gi < gene_cnt; gi++)
                    w_total += fund_bipart_sets[i].supports[gi];
                si = 0; gi = 0;
                while(si < super_cnt && gi < gene_cnt)
                    {
                    if(super_bp[si] == fund_bipart_sets[i].hashes[gi])
                        { w_shared += fund_bipart_sets[i].supports[gi]; si++; gi++; }
                    else if(super_bp[si] < fund_bipart_sets[i].hashes[gi]) si++;
                    else gi++;
                    }
                loss = (w_total > 0.0f) ? (w_total - w_shared) / w_total : 0.0f;
                }
            else
                {
                int shared = bipart_intersection_count(super_bp, super_cnt,
                                                       fund_bipart_sets[i].hashes, gene_cnt);
                loss = (gene_cnt > 0) ? (float)(gene_cnt - shared) / (float)gene_cnt : 0.0f;
                }
            loss *= tree_weights[i];
            sourcetree_scores[i] = loss;
            total += loss;
            }
        else
            {
            total += sourcetree_scores[i];
            }
        }

    free(super_bp);
    free(tmp);
    free(pruned_nwk);
    return total;
    }

/* -----------------------------------------------------------------------
 * compare_trees_qfit: quartet-fit criterion.
 *   Q_agree_i = Σ_{shared splits e} C(s_e,2)*C(n_i-s_e,2)
 *   Q_total_i = Σ_{all source splits e} C(s_e,2)*C(n_i-s_e,2)
 *   loss_i    = 1 - Q_agree_i / Q_total_i
 * Returns weighted sum of per-source-tree loss values (0 = all quartets agree).
 * Requires rf_precompute_fund_biparts() to have been called (needs sizes).
 * ----------------------------------------------------------------------- */
float compare_trees_qfit(int spr)
    {
    int i, j;
    float total = 0.0f;
    char *pruned_nwk    = malloc(TREE_LENGTH * sizeof(char));
    if(!pruned_nwk) memory_error(89);
    char *tmp           = malloc(TREE_LENGTH * sizeof(char));
    if(!tmp) { free(pruned_nwk); memory_error(90); }
    uint64_t *super_bp  = malloc(number_of_taxa * sizeof(uint64_t));
    if(!super_bp) { free(pruned_nwk); free(tmp); memory_error(91); }

    for(i = 0; i < Total_fund_trees; i++)
        {
#ifdef _OPENMP
        #pragma omp flush(user_break)
#endif
        if(user_break) break;
        if(!sourcetreetag[i]) continue;

        int found = FALSE, here1 = TRUE, here2 = TRUE;
        if(sourcetree_scores[i] != -1 && spr)
            {
            for(j = 0; j < number_of_taxa; j++)
                {
                if(presence_of_taxa[i][j] > 0 && presenceof_SPRtaxa[j] == TRUE)  here1 = FALSE;
                if(presence_of_taxa[i][j] > 0 && presenceof_SPRtaxa[j] == FALSE) here2 = FALSE;
                }
            if(here1 || here2) found = TRUE;
            }

        if(!found || !spr)
            {
            int ntaxa_i = 0;
            uint64_t total_hash = 0;
            for(j = 0; j < number_of_taxa; j++)
                if(presence_of_taxa[i][j]) { ntaxa_i++; total_hash ^= taxon_hash_vals[j]; }
            if(ntaxa_i < 4) { sourcetree_scores[i] = 0.0f; continue; }

            prune_tree(tree_top, i);
            shrink_tree(tree_top);

            pruned_nwk[0] = '\0';
            if(print_pruned_tree(tree_top, 0, pruned_nwk, FALSE, 0) > 1)
                {
                tmp[0] = '\0';
                strcpy(tmp, "("); strcat(tmp, pruned_nwk); strcat(tmp, ")");
                strcpy(pruned_nwk, tmp);
                }
            strcat(pruned_nwk, ";");

            while(unroottree(pruned_nwk));
            int super_cnt = collect_biparts_newick(pruned_nwk, total_hash, super_bp);
            reset_tree(tree_top);

            int gene_cnt = fund_bipart_sets[i].count;
            float loss;
            if(bsweight && fund_bipart_sets[i].supports)
                {
                double q_total_w;
                double q_agree_w = quartet_intersect_score_w(super_bp, super_cnt,
                                                              fund_bipart_sets[i].hashes,
                                                              fund_bipart_sets[i].sizes,
                                                              fund_bipart_sets[i].supports,
                                                              gene_cnt, ntaxa_i, &q_total_w);
                loss = (q_total_w > 0.0) ? (float)((q_total_w - q_agree_w) / q_total_w) : 0.0f;
                }
            else
                {
                int64_t q_total;
                int64_t q_agree = quartet_intersect_score(super_bp, super_cnt,
                                                           fund_bipart_sets[i].hashes,
                                                           fund_bipart_sets[i].sizes,
                                                           gene_cnt, ntaxa_i, &q_total);
                loss = (q_total > 0) ? (float)(q_total - q_agree) / (float)q_total : 0.0f;
                }
            loss *= tree_weights[i];
            sourcetree_scores[i] = loss;
            total += loss;
            }
        else
            {
            total += sourcetree_scores[i];
            }
        }

    free(super_bp);
    free(tmp);
    free(pruned_nwk);
    return total;
    }

float MRC(char *supertree)
    {
    int i=0, j=0, k=0, **super_matrix = NULL, count = 0, *tracking = NULL, parenthesis = 0, total=0, allsame1 = TRUE, allsame2 = TRUE;
     char number[100];
    float score = 0;


    /* assign the array super_matrix */
    super_matrix = malloc(number_of_taxa*sizeof(int *));
    if(!super_matrix) memory_error(64);
    for(i=0; i<number_of_taxa; i++)
        {
        super_matrix[i] = malloc((number_of_taxa)*sizeof(int));
        if(!super_matrix[i]) memory_error(65);
        for(j=0; j<number_of_taxa; j++)
            super_matrix[i][j] = 0;
        }


    /* Allocate the array that records the clade we are in during calculation */
    tracking = malloc(number_of_taxa*sizeof(int));
    for(i=0; i<number_of_taxa; i++) tracking[i] = FALSE;
    
    
    /** unroot the supertree (if necessary )  **/
    while(unroottree(supertree));
    
    i=0;j=0;k=0;
    /*** first we need to calculate the coding scheme for the supertree ***/
    /* we scan the nested parenthesis string once from right to left calculating the baum-ragan coding scheme */
    while(supertree[i] != ';')
        {
        switch(supertree[i])
            {

            case '(' :
                tracking[count] = TRUE;
                count++;
                i++;
                break;

            case ')':
                j = count;
                while(tracking[j] == FALSE) j--;
                if(tracking[j] == TRUE) tracking[j] = FALSE;
                i++;
                break;

            case ',':
                i++;
                break;

            default :
                for(j=0; j<100; j++)
                    number[j] = '\0';
                j=0;
                while(supertree[i] != '(' && supertree[i] != ')' && supertree[i] != ',')
                    {
                    number[j] = supertree[i];
                    i++;
                    j++;
                    }
                total = 0;
                j--; k=j;
                while(j>=0)
                    {
                    total+=(((int)pow(10,k-j))*texttoint(number[j]));  /* Turn the character to an integer */
                    j--;
                    }
                for(j=0; j<number_of_taxa; j++)
                    {
                    if(tracking[j] == TRUE)
                        {
                        super_matrix[j][total] = 1;
                        }
                    }
            }
        }
    total_partitions = 0;
    k=0;
    /** Next we are going to check for the presence of each of the partitions contained in total_coding in the super_matrix **/
    /** For each partition that is present in the supermatrix the score is incremented by 1 **/
   
    
    for(i=0; i<num_partitions; i++)  /* for every partition in total_coding */
        {
        for(j=1; j<number_of_taxa-2; j++) /* for every position in super_matrix */
            {
            allsame1 = TRUE;
            allsame2 = TRUE;
            for(k=0; k<number_of_taxa; k++)  /* for every position in both arrays */
                {
                if(total_coding[i][k] != 3)  
                    {
                    if(total_coding[i][k] != super_matrix[j][k]) allsame1 = FALSE;
                    if(total_coding[i][k] == super_matrix[j][k]) allsame2 = FALSE;
                    }
                }
            if(allsame1 || allsame2)
                {
                score += partition_number[i];
                /** if we have found a matching partition then we need to stop looking for this partition **/
                k=number_of_taxa;
                j=number_of_taxa-2;
                }
           
            }
        total_partitions += partition_number[i];
        }

    for(i=0; i<number_of_taxa; i++)
        {
        free(super_matrix[i]);
        }
    free(super_matrix);
    super_matrix = NULL;
    free(tracking);
    tracking = NULL;
    return(total_partitions-score);
    }

float quartet_compatibility(char * supertree)
    {
    int i=0, j=0, k=0,w=0, x=0, y=0, z=0, **super_matrix = NULL, count = 0, *tracking = NULL, parenthesis = 0, total=0,  allsame1 = TRUE, allsame2 = TRUE;
    char number[100];
    float score = 0, num_quartets = 0;
    int zerocount = 0, onecount = 0, *done_tree = NULL;
     
    /* assign the array super_matrix */
    super_matrix = malloc(number_of_taxa*sizeof(int *));
    if(!super_matrix) memory_error(64);
    for(i=0; i<number_of_taxa; i++)
        {
        super_matrix[i] = malloc((number_of_taxa)*sizeof(int));
        if(!super_matrix[i]) memory_error(65);
        for(j=0; j<number_of_taxa; j++)
            super_matrix[i][j] = 0;
        }


    done_tree = malloc(Total_fund_trees*sizeof(int));
    for(i=0; i<Total_fund_trees; i++)
        done_tree[i] = FALSE;

    /* Allocate the array that records the clade were are in during calculation */
    tracking = malloc(number_of_taxa*sizeof(int));
    for(i=0; i<number_of_taxa; i++) tracking[i] = FALSE;
    
    
    /** unroot the supertree (if necessary )  **/
    unroottree(supertree);
    
    i=0;j=0;k=0;
    /*** first we need to calculate the coding scheme for the supertree ***/
    /* we scan the nested parenthesis string once from right to left calculating the baum-ragan coding scheme */
    while(supertree[i] != ';')
        {
        switch(supertree[i])
            {

            case '(' :
                tracking[count] = TRUE;
                count++;
                i++;
                break;

            case ')':
                j = count;
                while(tracking[j] == FALSE) j--;
                if(tracking[j] == TRUE) tracking[j] = FALSE;
                i++;
                break;

            case ',':
                i++;
                break;

            default :
                for(j=0; j<100; j++)
                    number[j] = '\0';
                j=0;
                while(supertree[i] != '(' && supertree[i] != ')' && supertree[i] != ',')
                    {
                    number[j] = supertree[i];
                    i++;
                    j++;
                    }
                total = 0;
                j--;k=j;
                while(j>=0)
                    {
                    total+=(((int)pow(10,k-j))*texttoint(number[j]));  /* Turn the character to an integer */
                    j--;
                    }
                for(j=0; j<number_of_taxa; j++)
                    {
                    if(tracking[j] == TRUE)
                        {
                        super_matrix[j][total] = 1;
                        }
                    }
            }
        }
       


    
    k=0;
    /** Next we are going to check for the presence of each of the partitions contained in total_coding in the super_matrix **/
    /** For each partition that is present in the supermatrix the score is incremented by 1 **/
    
        for(w=0; w<number_of_taxa-3; w++)   /* first position in the quartet */
            {
            for(x=w+1; x<number_of_taxa-2; x++)   /* second position in the quartet */
                {
                for(y=x+1; y<number_of_taxa-1; y++)    /* third position in the quartet */
                    {
                    for(z=y+1; z<number_of_taxa; z++)     /* fourth position in the quartet  */
                        {
                        for(i=0; i<Total_fund_trees; i++) done_tree[i] = FALSE;
                        for(i=0; i<num_partitions; i++)  /* for every partition in total_coding */
                            {
                            if(!done_tree[from_tree[i]] && (total_coding[i][w]+total_coding[i][x]+total_coding[i][y]+total_coding[i][z]) == 2)  
                                {
                                for(j=1; j<number_of_taxa-2; j++) /* for every position in super_matrix */
                                    {
                                    if((super_matrix[j][w] +super_matrix[j][x] +super_matrix[j][y] +super_matrix[j][z]) == 2)
                                        {
                                        num_quartets += partition_number[i];  
                                        done_tree[from_tree[i]] = TRUE;                                
                                        if(total_coding[i][w] != super_matrix[j][w] || total_coding[i][x] != super_matrix[j][x] || total_coding[i][y] != super_matrix[j][y] || total_coding[i][z] != super_matrix[j][z]) allsame1 = FALSE;
                                        if(total_coding[i][w] == super_matrix[j][w] || total_coding[i][x] == super_matrix[j][x] || total_coding[i][y] == super_matrix[j][y] || total_coding[i][z] == super_matrix[j][z]) allsame2 = FALSE;
                                        if(allsame1 || allsame2)
                                            {
                                            score += partition_number[i];
                                            }
                                        j = number_of_taxa;  /** if we have found a matching partition then we need to stop looking for this quartet **/  
                                        allsame1 = TRUE;
                                        allsame2 = TRUE;
                                        }
                                    }

                                }
                            }
                        }
                    }
                }
            }
        
    
 
    for(i=0; i<number_of_taxa; i++)
        {
        free(super_matrix[i]);
        }
    free(super_matrix);
    super_matrix = NULL;
    free(tracking);
    if(done_tree != NULL) free(done_tree);
    tracking = NULL;
    return(num_quartets-score);
    }

