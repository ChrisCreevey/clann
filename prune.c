/*
 *  prune.c — Clann v5.0.0
 *  Tree pruning and collapsing utilities
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
#include "prune.h"

void random_prune(char *fund_tree)
	{
	float user_limit;
	int i=0, j=0;
	int *to_delete = NULL;
	char *prunecommand = NULL;
	FILE *rp_outfile = NULL;
	rp_outfile = fopen("prunedtaxa.txt", "w");
	
	prunecommand = malloc(1000000*sizeof(char));
	prunecommand[0] = '\0';
	
	to_delete = malloc(number_of_taxa*sizeof(int));
	for(i=0; i<number_of_taxa; i++) to_delete[i] = FALSE;
	
	for(i=0; i<num_commands; i++)
		{
		if(strcmp(parsed_command[i], "brlen") == 0)
			user_limit = atof(parsed_command[i+1]);
		}
	if(tree_top != NULL) dismantle_tree(tree_top);
	tree_top = NULL;
	
	temp_top = NULL;
	{ int _to = 0; tree_build(1, fund_tree, tree_top, 0, -1, &_to); }
	tree_top = temp_top;
	temp_top = NULL;
	
	collapse_clades(tree_top, user_limit, to_delete, rp_outfile);
	strcpy(prunecommand, "deletetaxa ");
	for(i=0; i<number_of_taxa; i++)
		{
		if(to_delete[i] == TRUE)
			{
			j++;
			strcat(prunecommand, taxa_names[i]);
			strcat(prunecommand, " ");
			
			}
		}
	if(j > 0)
		{
		num_commands = parse_command(prunecommand);
		exclude_taxa(FALSE);
		strcpy(prunecommand, "showtrees savetrees=yes display=no filename=prunedtree.ph");
		num_commands = parse_command(prunecommand);
		showtrees(FALSE);
		printf2("\n\tSummary of pruned taxa writtin to file prunedtaxa.txt\n");
		}
	else
		{
		printf2("0 clades met the criteria specified\n");
		}
	
	fclose(rp_outfile);
	}

void collapse_clades(struct taxon * position, float user_limit, int * to_delete, FILE *rp_outfile)
	{
	float total = 0;
	int count = 0, keep = 0, taxa_count = 0;
	
	
	while(position != NULL)
		{
		if(position->daughter != NULL)
			{
			total =0; count = 0;
			taxa_count = get_brlens(position->daughter, &total, &count);
			if(total/count <= user_limit)
				{
				/* collect names of taxa in clade for collapsing */
				keep = (int)fmod(rand(), taxa_count);
				print_keep(position->daughter, keep, 0, rp_outfile);
				untag_taxa(position->daughter, to_delete, keep, 0, rp_outfile);
				fprintf(rp_outfile, "\n");
				}
			else
				collapse_clades(position->daughter, user_limit, to_delete, rp_outfile);
			}
		position = position->next_sibling;
		}
	
	}

int get_brlens(struct taxon * position, float *total, int *count)
	{
	int taxa_count = 0, tmpcount = 0;
	while(position != NULL)
		{
		if(position->daughter != NULL)
			{
			taxa_count += get_brlens(position->daughter, total, count);
			}
		else
			{
			taxa_count++;
			}
		*total += return_length(position->weight);
		tmpcount ++;
		position = position->next_sibling;
		}
	*count+=tmpcount;
	return(taxa_count);
	}

float return_length(char *string)
	{
	int i=0, j=0;
	float length = 0;
	char *flt_length = malloc(TREE_LENGTH * sizeof(char));
	if(!flt_length) { printf2("Error: out of memory in return_length\n"); return 0; }

	while(i < strlen(string) && string[i] != ':' ) i++;
	
	if(i < strlen(string))
		{
		i++;
		j=0;
		while(string[i] != '\0')
			{
			flt_length[j] = string[i];
			j++;
			i++;
			}
		flt_length[j] = '\0';
		length = atof(flt_length);
		}
	free(flt_length);
	return(length);
	}

int print_keep(struct taxon *position, int keep, int count, FILE *rp_outfile)
	{
	while(position != NULL)
		{
		if(position->daughter != NULL)	
			count = print_keep(position->daughter, keep, count, rp_outfile);
		else
			{
			if(count == keep)
				{
				fprintf(rp_outfile, "%s:\t", taxa_names[position->name]);
				}
			count++;
			}
		position = position->next_sibling;
		}
	return(count);	
	}

int untag_taxa(struct taxon *position, int * to_delete, int keep, int count, FILE *rp_outfile)
	{
	
	
	while(position != NULL)
		{
		if(position->daughter != NULL)
			count = untag_taxa(position->daughter, to_delete, keep, count, rp_outfile);
		else
			{
			if(count != keep)
				{
				to_delete[position->name] = TRUE;
				fprintf(rp_outfile, "%s, ", taxa_names[position->name]);
				}
			count++;
			}
		position = position->next_sibling;
		}
	return(count);
	}

/* -----------------------------------------------------------------------
 * collapse_monophyly_in_tree: collapse species-specific (monophyletic
 * in-paralog) clades in a single fundamental tree down to one representative
 * per clade.  This is the per-tree logic shared by prune_monophylies()
 * (the "prunemonophylies" command) and decomposegenetrees()'s Stage 1
 * (see NOTES_gene_tree_decomposition.md §5, §9 step 1).
 *
 * treenum   : index into fundamentals[]/tree_names[] to collapse. Assumes
 *             select_longest has already been set by the caller (controls
 *             the representative-choice rule used by identify_species_specific_clades).
 * infofile  : if non-NULL, appends a "Tree # N [ name ]\n\tRemoved:...\tKEPT:...\n"
 *             report block, byte-identical to what prune_monophylies has
 *             always written to prunedtrees_info.txt.
 * out_tree  : caller-allocated buffer (>= TREE_LENGTH) that receives the
 *             collapsed Newick: numeric taxon ids, parenthesis-wrapped if
 *             more than one node remains, terminated with ';'.
 * out_num_commas : set to the number of top-level commas in *out_tree
 *             (the existing "4 or more taxa" gate used elsewhere is
 *             *out_num_commas > 2).
 * out_still_multicopy : if non-NULL, set to TRUE if any taxon still has
 *             more than one surviving copy in the collapsed tree (i.e. this
 *             tree needs further, non-topological decomposition), FALSE if
 *             the tree is now fully single-copy. May be NULL if the caller
 *             doesn't need this (prune_monophylies doesn't).
 * out_presence : if non-NULL, caller-allocated array of length
 *             number_of_taxa that receives the per-taxon surviving-leaf
 *             counts for the collapsed tree (same counts used internally
 *             for out_still_multicopy, exposed so callers -- e.g.
 *             decomposegenetrees Stage 3 -- can apply their own
 *             minfragtaxa/minfragspecies floor check without re-parsing
 *             *out_tree or re-deriving counts via print_pruned_tree's
 *             unreliable return value). Only populated if out_still_multicopy
 *             is also non-NULL; may be NULL if the caller doesn't need it.
 *
 * Returns the number of leaves retained (as returned by print_pruned_tree).
 * ----------------------------------------------------------------------- */
int collapse_monophyly_in_tree(int treenum, FILE *infofile, char *out_tree, int *out_num_commas, int *out_still_multicopy, int *out_presence)
	{
	int k=0, l=0, numt=0, *taxa_fate = NULL, clannID=0, report=FALSE, taxaorder=0, leaf_count=0, i=0, num_nodes=0;
	char *tmp = NULL, **taxa_fate_names = NULL;
	/* Size buffers to THIS tree's full-name form: returntree_fullnames() expands
	 * numeric ids to (longer) full names, which for large gene-family trees
	 * exceeds a fixed TREE_LENGTH. */
	size_t cbuf = returntree_fullname_buflen(fundamentals[treenum], treenum);
	if(cbuf < strlen(fundamentals[treenum]) + 1) cbuf = strlen(fundamentals[treenum]) + 1;
	if(cbuf < (size_t)TREE_LENGTH) cbuf = TREE_LENGTH;
	char *temptree = malloc(cbuf * sizeof(char));
	if(!temptree) { printf2("Error: out of memory in collapse_monophyly_in_tree\n"); return 0; }
	tmp = malloc(cbuf * sizeof(char));
	if(!tmp) { free(temptree); printf2("Error: out of memory in collapse_monophyly_in_tree\n"); return 0; }

	if(tree_top != NULL) dismantle_tree(tree_top);  /* Dismantle any trees already in memory */
	tree_top = NULL;

	temp_top = NULL;
	taxaorder=0;
	strcpy(temptree, fundamentals[treenum]);
	returntree_fullnames(temptree, treenum);

	/***  build the sourcetree in memory ***/
	if(tree_top != NULL)
        {
        dismantle_tree(tree_top);
        tree_top = NULL;
        }
    temp_top = NULL;
    basic_tree_build(1, temptree, tree_top, TRUE);

    tree_top = temp_top;
    temp_top = NULL;
    reset_tree(tree_top);

    /* Number all taxa in the tree for tracing which taxa were deleted or retained and record the number of taxa in the variable numt */
    numt=0; clannID=0;
    numt = number_tree2(tree_top, numt);
    /* alloc an arrary to hold all the names of hte taxa in the order they are found on the tree */
    taxa_fate_names = malloc(numt*sizeof(char*));
    for(k=0; k<numt; k++)
    	{
    	taxa_fate_names[k] = malloc(NAME_LENGTH*sizeof(char));
    	taxa_fate_names[k][0] = '\0';
    	}
    get_taxa_names(tree_top, taxa_fate_names); /* copy names into the array taxa_fate_names */

    /* Define an array to record the species that are pruned (and kept) in this tree */
    taxa_fate = malloc(numt*sizeof(int));
    for(k=0; k<numt; k++)
		{
		taxa_fate[k]=0;
		}

    clannID = identify_species_specific_clades(tree_top, numt, taxa_fate, clannID);  /* Call recursive function to travel down the tree looking for species-specific clades */
    if(infofile != NULL)
    	{
    	fprintf(infofile, "\nTree # %d [ %s ]\n\t", treenum, tree_names[treenum]);
    	/* Print out the list of deleted and pruned taxa from the taxa_fate array */
    	for(k=1; k<=numt; k++)
    		{
    		report=FALSE;
    		/* find out if there are any taxa in clann k for reporting */
    		for(l=0; l<numt; l++)
    			{
    			if(abs(taxa_fate[l]) == k)
    				{
    				report=TRUE;
    				l=numt;
    				}
    			}
    		if(report == TRUE)
    			{
    			for(l=0; l<numt; l++)
    				{
    				if(taxa_fate[l] == k*-1)
    					fprintf(infofile, "Removed:%s\t", taxa_fate_names[l]);
    				}
    			for(l=0; l<numt; l++)
    				{
    				if(taxa_fate[l] == k)
    					fprintf(infofile, "KEPT:%s\n\t", taxa_fate_names[l]);
    				}
    			}
    		}
    	}

	free(taxa_fate);
    for(k=0; k<numt; k++) free(taxa_fate_names[k]);
    free(taxa_fate_names);

    shrink_tree(tree_top);    /* Shrink the pruned tree by switching off any internal nodes that are not needed */

    if(out_still_multicopy != NULL)
    	{
    	int *presence = out_presence != NULL ? out_presence : malloc(number_of_taxa*sizeof(int));
    	for(k=0; k<number_of_taxa; k++) presence[k] = 0;
    	count_surviving_leaves(tree_top, presence);
    	*out_still_multicopy = FALSE;
    	for(k=0; k<number_of_taxa; k++)
    		if(presence[k] > 1) { *out_still_multicopy = TRUE; break; }
    	if(out_presence == NULL) free(presence);
    	}

    out_tree[0] = '\0'; /* initialise the string */
    leaf_count = print_pruned_tree(tree_top, 0, out_tree, TRUE, treenum);
    if(leaf_count > 1)
        {
        tmp[0] = '\0';
        strcpy(tmp, "(");
        strcat(tmp, out_tree);
        strcat(tmp, ")");
        strcpy(out_tree, tmp);
        }
    strcat(out_tree, ";");

    num_nodes=0;
    for(i=0; i<(int)strlen(out_tree); i++)
    	{if(out_tree[i] == ',') num_nodes++;	}
    if(out_num_commas != NULL) *out_num_commas = num_nodes;

    if(tree_top != NULL) { dismantle_tree(tree_top); tree_top = NULL; }
    temp_top = NULL;

    free(temptree);
    free(tmp);
    return(leaf_count);
	}

/* Recursively count surviving (tag != FALSE) leaves per taxon id, for the
 * still-multicopy check in collapse_monophyly_in_tree.  Mirrors the
 * tag-based skip logic used by print_pruned_tree/autoprunemono's
 * recount_from_tree (main.c). */
void count_surviving_leaves(struct taxon *position, int *presence)
	{
	while(position != NULL)
		{
		if(position->daughter != NULL)
			count_surviving_leaves(position->daughter, presence);
		else
			{
			if(position->tag != FALSE && position->name >= 0 && position->name < number_of_taxa)
				presence[position->name]++;
			}
		position = position->next_sibling;
		}
	}

void prune_monophylies(void)
    {
    int j=0, num_nodes=0, trees_included=0, trees_excluded=0;
    char *pruned_tree = NULL, filename2[10000];
    char filename3[10000];
    FILE *pm_outfile = NULL;
    select_longest=FALSE;
    filename2[0]= filename3[0] = '\0';
    strcpy(filename2, "prunedtrees.txt");
    strcpy(filename3, "prunedtrees_info.txt");

    size_t pt_bufsz = TREE_LENGTH;   /* grown per tree below to fit large source trees */
    pruned_tree = malloc(pt_bufsz*sizeof(char));
    pruned_tree[0] = '\0';

    for(j=0; j<num_commands; j++)
        {
        if(strcmp(parsed_command[j], "selection") == 0)
            {
            if(strcmp(parsed_command[j+1], "length")==0) select_longest=TRUE;
            }
        if(strcmp(parsed_command[j], "filename") == 0)
            {
            strcpy(filename2, parsed_command[j+1]);
            strcpy(filename3, filename2);
            strcat(filename3, "_info.txt");
            }

        }

    pm_outfile = fopen(filename2, "w");
    tempoutfile = fopen(filename3, "w");


    printf2("input trees will be pruned where clans of single species exist\nOne remaining representative will be chosen by ");
    if(select_longest == TRUE) printf2 ("the length of the sequence (in the name, after the species (delimited with a \"%c\")\n", delimiter_char);
        else printf2("random\n");

    for(j=0; j<Total_fund_trees; j++)
        {
        /* Grow pruned_tree to hold this tree's collapsed (numeric) form, which for
         * large gene-family trees can exceed a fixed TREE_LENGTH. */
        if(strlen(fundamentals[j]) + 128 > pt_bufsz)
            {
            char *np = realloc(pruned_tree, strlen(fundamentals[j]) + 128);
            if(np) { pruned_tree = np; pt_bufsz = strlen(fundamentals[j]) + 128; }
            }
        pruned_tree[0] = '\0';
        collapse_monophyly_in_tree(j, tempoutfile, pruned_tree, &num_nodes, NULL, NULL);

    	if(num_nodes > 2) /* if the remaining tree has more then 3 taxa */
    		{
	        if(strcmp(tree_names[j], "") != 0)
	        	fprintf(pm_outfile, "%s[%s]\n", pruned_tree, tree_names[j]);
	        else
	        	fprintf(pm_outfile, "%s[%d]\n", pruned_tree, j);
	        trees_included++;
        	}
        else
        	{
        	trees_excluded++;
        	}
        }

    printf2("\nPruning finished. %d pruned trees with 4 or more taxa written to the file \"%s\"\n", trees_included, filename2);
   printf2("Information on pruned and retained taxa in each tree written to the file \"%s\"\n",  filename3);
    printf2("%d pruned trees with less than 4 taxa were excluded\n",trees_excluded );
    free(pruned_tree);
    fclose(pm_outfile);
    fprintf(tempoutfile, "\n");
    fclose(tempoutfile);
    }

int identify_species_specific_clades(struct taxon * position, int numt, int * taxa_fate, int clannID)
    {
    float total = 0, keepnew=FALSE;
    int k=0, count = 0, taxa_count = 0, all_same_taxon = -1, *foundtaxa = NULL, i, *tmp_taxa_fate = NULL;
    long seqlength = 0;
    struct taxon * starting = position, *longest=NULL;
    foundtaxa = malloc(number_of_taxa*sizeof(int));

    /* Define an array to record the species that are pruned (and kept) in this tree */
    tmp_taxa_fate = malloc(numt*sizeof(int));
    for(k=0; k<numt; k++) tmp_taxa_fate[k]=0;
    
    /* 1) identify if all descendant nodes are from the same species */
    
    for(i=0; i<number_of_taxa; i++) foundtaxa[i] = FALSE;  /* array to keep track of taxa to delete */
    longestseq=NULL;
    seqlength= list_taxa_in_clade(position, foundtaxa, longest, seqlength);
    count=0;
    for(i=0; i<number_of_taxa; i++) {
        if(foundtaxa[i] == TRUE) count++;  /* count the number of different taxa we found here */
        }

    if(count > 1)
        { /* no point going further down this clade if we already know they are all the same species */
        while(position != NULL)
            {
            if(position->daughter != NULL)
                clannID = identify_species_specific_clades(position->daughter, numt, taxa_fate, clannID); /* if this node has a daughter, then call new instance of the function on the daughter */
            position = position->next_sibling;
            }
        }
    else /* if this was a species-specific clade, then untag everything, except for the longest */
        {
        clannID+=1; /* increment the clann ID */
        for(k=0; k<numt; k++) tmp_taxa_fate[k]=0;
        untag_nodes_below(position, tmp_taxa_fate, clannID);

    	/* check through tmp_taxa_fate to see if it overlaps with any clann definitions already in taxa_fate */
    	/* If there is overlap, it the new clann definition is a superset of the previous definition, then overwrite the previous definition in taxa_fate */
    	/* Otherwise do nothing */
    	keepnew=FALSE;
    	for(k=0; k<numt; k++)
    		{
    		if(tmp_taxa_fate[k] != 0 && taxa_fate[k] == 0) keepnew=TRUE;
    		}
    	if(keepnew == TRUE)
    		{
    		for(k=0; k<numt; k++) 
    			{
    			if(tmp_taxa_fate[k] != 0) 
					{
					taxa_fate[k] = tmp_taxa_fate[k];
					}
    			}
    		}
       /* fprintf (tempoutfile, "\n\t"); */
        }

    /* 2) Identify if all preceeding nodes are from the same species. */
        position = starting;
        for(i=0; i<number_of_taxa; i++) foundtaxa[i] = FALSE;  /* array to keep track of taxa to delete */
        count=0; seqlength = 0;
        if(position->parent != NULL)
            {   
            seqlength = list_taxa_above(position->parent, foundtaxa, longest, seqlength);
            }
        for(i=0; i<number_of_taxa; i++) {
            if(foundtaxa[i] == TRUE) count++;  /* count the number of different taxa we found here */
            }
        if(count == 1)
            {
            clannID+=1; /* increment the clann ID */
            for(k=0; k<numt; k++) tmp_taxa_fate[k]=0;
            untag_nodes_above(position->parent, tmp_taxa_fate, clannID);

	    	/* check through tmp_taxa_fate to see if it overlaps with any clann definitions already in taxa_fate */
	    	/* If there is overlap, it the new clann definition is a superset of the previous definition, then overwrite the previous definition in taxa_fate */
	    	/* Otherwise do nothing */
	    	keepnew=FALSE;
	    	for(k=0; k<numt; k++)
	    		{
	    		if(tmp_taxa_fate[k] != 0 && taxa_fate[k] == 0) keepnew=TRUE;
	    		}
	    	if(keepnew == TRUE)
	    		{
	    		for(k=0; k<numt; k++) 
	    			{
	    			if(tmp_taxa_fate[k] != 0) taxa_fate[k] = tmp_taxa_fate[k];
	    			}
	    		}


          /*  fprintf (tempoutfile, "\n\t"); */
            }

    free(foundtaxa);
    free(tmp_taxa_fate);
    return(clannID);
    }

long list_taxa_above(struct taxon * position, int * foundtaxa, struct taxon * longest, long seqlength) /* this function will check the parts of the tree above this position, it fixes aproblem that the tree rooting may mask monophylys */
    {
    long newseqlength=0;
    struct taxon * origin = position;
    /* go all the way to the last sibling */
    while(position->prev_sibling != NULL) position = position->prev_sibling;

    while(position != NULL)
        {   
        if(position->daughter != NULL && position != origin)
            {  /* go down throug this clade */
            seqlength = list_taxa_in_clade(position->daughter, foundtaxa, longest, seqlength);
            }
        if(position->name != -1 && position->tag == TRUE) /* check this node */
            {   
            foundtaxa[position->name]=TRUE;
            if(select_longest==TRUE)
                {   
                if((newseqlength = extract_length(position->fullname)) > seqlength) /* if this taxa has a longer length than the previously found longest */
                    {
                    seqlength = newseqlength;
                    longestseq = position;
                    }
                }
            else
                {   
                seqlength = newseqlength;
                longestseq = position;
                }
             }
        if(position->prev_sibling == NULL && position->parent != NULL)
            {   
            seqlength = list_taxa_above(position->parent, foundtaxa, longest, seqlength);
            }
        position = position->next_sibling;
        }

    return(seqlength);
    }

long list_taxa_in_clade(struct taxon * position, int * foundtaxa, struct taxon * longest, long seqlength) /* descend through the tree finding what taxa are there (and putting result into an array) and also identifying the longest sequence (the first number in the <<full>> name of the sequence, after the first "." and before the first "|") */
    {
    long newseqlength=0;
    
    while(position != NULL)
        {
        if(position->daughter != NULL)
            {
            seqlength = list_taxa_in_clade(position->daughter, foundtaxa, longest, seqlength);
            }
        else
            {
            if(position->tag == TRUE)
                {   
                foundtaxa[position->name]=TRUE;
                 if(select_longest==TRUE)
                    {   
                    if((newseqlength = extract_length(position->fullname)) > seqlength) /* if this taxa has a longer length than the previously found longest */
                        {
                        seqlength = newseqlength;
                        longestseq = position;
                        }
                    }
                else
                    {   
                    seqlength = newseqlength;
                    longestseq = position;
                    }
                }               
            }
        position = position->next_sibling;
        }
    return(seqlength);
    } 

long extract_length(char * fullname)
    {
    /* this function assumes that the length of the sequence is embedded in the full name of the taxa in the form "SPECIESNAME.SEQLEN|BLAHBLAHBLAH"
        the seperators can be a "." or a "|"
    */

    long seqlength = 0;
    char name[100000];
    char *eptr;

    name[0] = '\0';
    strcpy(name, fullname);

    const char s[3] = ".|";
    char *token;

    token=strtok(name, s); /* get first token (in this case the species name) */
    token = strtok(NULL, s); /* Get second token ( inthi case the squence length) */

    seqlength = strtol(token, &eptr, 10); /* extract string version of number as long int */

    return(seqlength);
    }

void untag_nodes_below(struct  taxon * position, int * taxa_fate, int clannID)
    {
    
    while(position != NULL)
        {
        if(position != longestseq) 
        	{
        		position->tag = FALSE;
        	if(position->daughter == NULL) 
        		{	
        		taxa_fate[position->tag2] = clannID*-1; /* to indicate its removal, assign tag2 to the negative of the clannID we are on */
        		/*fprintf(tempoutfile, "Removed:%s\t", position->fullname); */
        		}
        	}
        else
        	{
        	if(position->daughter == NULL) 
        		{	
        		taxa_fate[position->tag2] = clannID; /* to indicate that we are keeping this taxon, assign it to the positive of the clannID we are on */
        		/*fprintf(tempoutfile, "KEPT:%s\t", position->fullname); */
        		}
        	}
        if(position->daughter != NULL) untag_nodes_below(position->daughter, taxa_fate, clannID);
        position = position->next_sibling;
        }
    
    }

void untag_nodes_above(struct  taxon * position, int * taxa_fate, int clannID)
    {
    struct taxon * origin = position;
    while(position->prev_sibling != NULL) position = position->prev_sibling;

    while(position != NULL)
        {
        if(position != longestseq) 
        	{
        	position->tag = FALSE;
        	if(position->daughter == NULL) 
        		{
        		taxa_fate[position->tag2] = clannID*-1; /* to indicate its removal, assign tag2 to the negative of the clannID we are on */
        		/*fprintf(tempoutfile, "Removed:%s\t", position->fullname); */
        		}
        	}
        else
        	{	
        	if(position->daughter == NULL) 
        		{	
        		taxa_fate[position->tag2] = clannID; /* to indicate that we are keeping this taxon, assign it to the positive of the clannID we are on */
        		/*fprintf(tempoutfile, "KEPT:%s\t", position->fullname); */
        		}
        	}
        if(position->daughter != NULL && position != origin) untag_nodes_below(position->daughter, taxa_fate, clannID);
        if(position->parent != NULL) untag_nodes_above(position->parent, taxa_fate, clannID);
        position = position->next_sibling;
        }
    
    }

