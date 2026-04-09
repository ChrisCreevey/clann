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
            if(bh != 0)
                {
                int other_c    = ntaxa - child_c;
                out[cnt]       = bh;
                sizes_out[cnt] = (child_c < other_c) ? child_c : other_c;
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

void cal_fund_scores(int printfundscores)
    {
    int i =0, j=0, k=0;
    FILE *dists = NULL;
        
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
            reset_tree(tree_top);

            int gene_cnt = fund_bipart_sets[i].count;
            int shared   = bipart_intersection_count(super_bp, super_cnt,
                                                     fund_bipart_sets[i].hashes, gene_cnt);
            /* Raw RF distance — no normalisation */
            float d = (float)(super_cnt + gene_cnt - 2 * shared);
            /* [experimental] tree-size scaling: divide by k_i^alpha (alpha=0: Steel 2008, alpha=1: normalised) */
            if(ml_alpha != 0.0 && gene_cnt > 0)
                d /= (float)pow((double)gene_cnt, ml_alpha);
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

