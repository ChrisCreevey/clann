/* function for Karen to automatically collapse clades that have all of the same taxa in them, keeping only the one with the longest sequence (found in the full name) */

void prune_monophylies(char *fund_tree)
	{
	int i=0, j=0;
	char *pruned_tree = NULL;
	FILE *pm_outfile = NULL;
	rp_outfile = fopen("prunedtrees.txt", "w");
	
	pruned_tree = malloc(10000000*sizeof(char));
	pruned_tree[0] = '\0';

	for(j=0; j<Total_fund_trees; j++)
		{
  
		if(tree_top != NULL) dismantle_tree(tree_top);  /* Dismantle any trees already in memory */
		tree_top = NULL;
		
		temp_top = NULL;
		tree_build(1, fundamentals[j], tree_top, 0, -1); /* build the tree passed to the function */
		tree_top = temp_top;
		temp_top = NULL;
		
		identify_species_specific_clades(tree_top);  /* Call recursive function to travel down the tree looking for species-specific clades */
		shrink_tree(tree_top);    /* Shrink the pruned tree by switching off any internal nodes that are not needed */
		
		pruned_tree[0] = '\0'; /* initialise the string */
		if(print_pruned_tree(tree_top, 0, pruned_tree, TRUE) >1)
			{
			tmp[0] = '\0';
			strcpy(tmp, "(");
			strcat(tmp, pruned_tree);
			strcat(tmp, ")");
			strcpy(pruned_tree, tmp);
			}
		strcat(pruned_tree, ";");

		fprintf(pm_outfile, "%s\n", pruned_tree);
		}
	
	fclose(pm_outfile);
	}


void identify_species_specific_clades(struct taxon * position)
	{
	float total = 0;
	int count = 0, taxa_count = 0, all_same_taxon = -1, *foundtaxa = NULL;
	long seqlength = 0;
	struct taxon * starting = position, *longest=NULL;
	foundtaxa = malloc(number_of_taxa*sizeof(int));

	/* identify if all descendant nodes are from the same species */
	
	for(i=0; i<number_of_taxa; i++) foundtaxa[i] = FALSE;  /* array to keep track of taxa to delete */
	seqlength= list_taxa_in_clade(position, foundtaxa, longest, seqlength);
	for(i=0; i<number_of_taxa; i++) {
		if(foundtaxa[i] == TRUE) count++;  /* count the number of different taxa we found here */
		}

	if(count > 1)
		{
		reset_tree(position->daughter); /* if there was more than one species in this clade, then undo all the tagging from the previous */
			/* no point going further down this clade if we already know they are all the same species */
		while(position != NULL)
			{
			if(position->daughter != NULL)
				identify_species_specific_clades(position->daughter); /* if this node has a daughter, then call new instance of the function on the daughter */
			position = position->next_sibling;
			}
		}

	free(foundtaxa);

	}



long list_taxa_in_clade(struct taxon * position, int * foundtaxa, struct taxon * longest, long seqlength) /* descend through the tree finding what taxa are there (and putting result into an array) and also identifying the longest sequence (the first number in the <<full>> name of the sequence, after the first "." and before the first "|") */
	{

	long newseqlength=0;
	
	while(position != NULL)
		{
		if(position->daughter != NULL)
	        {
	        seqlength = list_taxa_in_clade(position->daughter);
	        }
	    else
	        {
	        if(position->name != -1)
	        	{
	        	foundtaxa[position->name]=TRUE;
	        	if((newseqlength = extract_length(position->fullname)) > seqlength) /* if this taxa has a longer length than the previously found longest */
	        		{
	        		seqlength = newseqlength;
	        		longest->tag = FALSE; /* mark the previously found largest for pruning */
	        		longest = position;
	        		}
	        	else
	        		{
	        		position->tag=FALSE;	
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



	
