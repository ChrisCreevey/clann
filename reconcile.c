/*
 *  reconcile.c — Clann v5.0.0
 *  Gene tree reconciliation and HGT reconstruction
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
#include "reconcile.h"

void gene_content_parsimony(struct taxon * position, int * array)
	{
	/* Go to the last internal nodes on the tree that are still tagged, then calculate the possible number of copies for this node */
	struct taxon *start = position;
	
	while(position != NULL)
		{
		if(position->daughter != NULL && position->tag > 0) gene_content_parsimony(position->daughter, array);
		position = position->next_sibling;		
		}
	position = start;
	while(position != NULL)
		{
		
		
		position = position->next_sibling;
		}
	}

void nj(void)
	{
	int i, j, missing_method = 1, error = FALSE;
	char *tree = NULL, useroutfile[100], *fakefilename = NULL;
	FILE *outfile = NULL;
	int *saved_tags = NULL;  /* for single-copy auto-filter */
	
	useroutfile[0] = '\0';
	strcpy(useroutfile, "NJ-tree.ph");
    for(i=0; i<num_commands; i++)
        {
        if(strcmp(parsed_command[i], "missing") == 0)
			{
			if(strcmp(parsed_command[i+1], "4point") == 0)
				missing_method = 1;
			else
				{
				if(strcmp(parsed_command[i+1], "ultrametric") == 0)
					missing_method = 0;
				else
					{
					printf2("Error: '%s' is an invalid method to be assigned to missing\n", parsed_command[i+1]);
					missing_method = 1;
					error = TRUE;
					}
				}
			}
			
		if(strcmp(parsed_command[i], "savetrees") == 0)
			strcpy(useroutfile, parsed_command[i+1]);
		}
		
	if(!error)
		{
		if((outfile = fopen(useroutfile, "w")) == NULL)
			{
			printf2("Error opening file named %s\n", useroutfile);
			error = TRUE;
			}
		}

	if(!error)
		{
		fakefilename = malloc(100*sizeof(char));
		fakefilename[0] = '\0';
		tree = malloc(TREE_LENGTH*sizeof(char));
		printf2("\n\nNeighbor-joining settings:\n\tDistance matrix generation by average consensus method\n\tEstimation of missing data using ");
		if(missing_method == 1)
			printf2("4 point condition distances\n");
		if(missing_method == 0)
			printf2("ultrametric distances\n");
		printf2("\tresulting tree saved to file %s\n\n\n", useroutfile);

		saved_tags = apply_singlecopy_filter();

		average_consensus(0, missing_method, fakefilename, NULL);
		neighbor_joining(TRUE, tree, TRUE);
		fprintf(outfile, "%s\n", tree);
		tree_coordinates(tree, TRUE, TRUE, FALSE, -1);
		
		for(i=0; i<number_retained_supers; i++)
			{
			strcpy(retained_supers[i], "");
			scores_retained_supers[i] = -1;
			}
		retained_supers[0] = realloc(retained_supers[0], (strlen(tree)+10)*sizeof(char));
		strcpy(retained_supers[0], tree);
		scores_retained_supers[0] = 0;
		trees_in_memory = 1;
		restore_singlecopy_filter(saved_tags);
		fclose(outfile);
		free(tree);
		}

	}

void isittagged(struct taxon * position)
	{
	while(position != NULL)
		{
		if(position->tag2 == TRUE) printf2("found\n");
		if(position->daughter != NULL)
			{
			isittagged(position->daughter);
			}
		position = position->next_sibling;
		}
	}

float tree_map(struct taxon * gene_top, struct taxon * species_top, int print)
	{
	int xnum =0, *presence = NULL, i, j, num_dups = 0, num_losses = 0, best = 0;
	char *temptree = NULL, reconfilename[100], *treetmp = NULL;

	treetmp = malloc(TREE_LENGTH*sizeof(int));
	treetmp[0] = '\0';
	
	presence = malloc(2*number_of_taxa*sizeof(int));
	for(i=0; i<2*number_of_taxa; i++)
		presence[i] = FALSE;
	
	/** 1) Label all internal and external taxa on the species tree ****/

	xnum = number_tree1(species_top, number_of_taxa);
	xnum--;
	
	/****2) label all the gene tree nodes (and taxa) with their equivalent on the species tree **/
	label_gene_tree(gene_top, species_top, presence, xnum);
	/*** 3) From the bottom-up, Identify those duplications at positions where the id of a pointer taxon is that same as any of its daughters **/
	num_dups = reconstruct_map(gene_top, species_top);
	/**** 4) we need to add the bits of the trees that are missing ****/
	add_losses(gene_top, species_top);
	join_losses(gene_top);
	num_losses = count_losses(gene_top);
	free(presence);
	free(treetmp);

	if(print)fprintf(distributionreconfile, "%d\t%d\n", num_dups, num_losses);
	return((dup_weight*(float)num_dups)+(loss_weight*(float)num_losses));
	}

void resolve_tricotomies(struct taxon *position, struct taxon *species_tree)
	{
	struct taxon *start = position, *copy = NULL, *best1 = NULL, *best2 = NULL, **thislevel = NULL, *previous = position->parent, *pos1 = NULL, *pos2 = NULL, *newbie = NULL, *lastsibling = NULL;
	int i=0, j=0, k=0, l=0, x=0, y=0, multiples = FALSE, topoftree = FALSE, *found = NULL, num_dups = 0, num_losses = 0, *presence = NULL, possible_internal = 0, donestuff = TRUE, smallest = 2*number_of_taxa, result, nope = TRUE;
	float best_score = -1;
	char *pruned_tree = NULL;
	
	
	presence = malloc((2*number_of_taxa)*sizeof(int));
	found = malloc((2*number_of_taxa)*sizeof(int));
	for(j=0; j<2*number_of_taxa; j++)
		{
		presence[j] = 0;
		found[j] = 0;
		}	
	/* go down to the bottom of the gene tree */
	
	
	while(position != NULL)
		{
		i++;
		if(position->daughter != NULL)
			{
			resolve_tricotomies(position->daughter, species_tree);
			found[position->tag]++;
			}
		else
			found[position->name]++;
		position = position->next_sibling;
		}
	position = start;
	/* work out any tricotomies containing multiples at this level */
	if(start->parent == NULL)
		{
		topoftree=TRUE;
		}

/*	if( i > 2 && start->parent != NULL ) *//*if there are more than two siblings at this level, and there is at least one more than once */
	
	if( (i > 2 && !topoftree) || ( i>3 && topoftree)) /*if there are more than two siblings at this level, and there is at least one more than once */
		{
		l = i;
		while(donestuff && ((l > 2 && !topoftree) || ( l>3 && topoftree)))
			{
			while(donestuff && ((l > 2 && !topoftree) || ( l>3 && topoftree)))
				{
				
				donestuff = FALSE;
				/* Step 1: Identify sister taxa and put them together (while marking new internal nodes as such )*/
				pos1=start;
				x=0;
				
				while(pos1 != NULL && ((l > 2 && !topoftree) || ( l>3 && topoftree)))
					{
					pos2 = pos1->next_sibling;
					y=x+1;
					while(pos2 != NULL && ((l > 2 && !topoftree) || ( l>3 && topoftree)))
						{
						if((possible_internal = are_siblings(species_tree, pos1->tag, pos2->tag)) != FALSE)
							{
						
							/* put them together */
							newbie = make_taxon();
							newbie->tag = possible_internal;
							newbie->daughter = pos1;
							newbie->tag2 = TRUE;
							if(pos1->prev_sibling == NULL)
								{
								if(pos1->parent != NULL)
									(pos1->parent)->daughter = newbie;
								newbie->parent = pos1->parent;
								start = newbie;
								}
							if(pos1->next_sibling != NULL) (pos1->next_sibling)->prev_sibling = newbie;
							if(pos1->prev_sibling != NULL) (pos1->prev_sibling)->next_sibling = newbie;
							pos1->parent = newbie;
							newbie->next_sibling = pos1->next_sibling;
							newbie->prev_sibling = pos1->prev_sibling;
							pos1->next_sibling = pos2;
							pos1->prev_sibling = NULL;
							
							
							if(pos2->next_sibling != NULL) (pos2->next_sibling)->prev_sibling = pos2->prev_sibling;
							if(pos2->prev_sibling != NULL) (pos2->prev_sibling)->next_sibling = pos2->next_sibling;
							
							pos2->next_sibling = NULL;
							pos2->prev_sibling = pos1;
							
							pos1 = start;
							pos2 = pos1->next_sibling;
							if(topoftree && newbie->prev_sibling == NULL) temp_top2 = newbie;
							newbie = NULL;
							l--;
							donestuff = TRUE;
							/*printf("done siblings\n");*/
							}
						else
							{
							pos2 = pos2->next_sibling;
							y++;
							}
						}
					pos1 = pos1->next_sibling;
					x++;
					}
				
				if(((l > 2 && !topoftree) || ( l>3 && topoftree)) )
					{
					/* Step 2: Join together any that are the same */
					pos1=start;
					x=0;
					while((pos1 != NULL && ((l > 2 && !topoftree) || ( l>3 && topoftree)))  )
						{
						nope = TRUE;
						pos2 = pos1->next_sibling;
						y=x+1;
						while((pos2 != NULL && ((l > 2 && !topoftree) || ( l>3 && topoftree))) )
							{
							if(pos1->tag == pos2->tag)
								{
								newbie = make_taxon();
								newbie->tag = pos1->tag;
								newbie->daughter = pos1;
								newbie->tag2=TRUE;
								if(pos1->prev_sibling == NULL)
									{
									if(pos1->parent != NULL)
										(pos1->parent)->daughter = newbie;
									newbie->parent = pos1->parent;
									start = newbie;
									}
								if(pos1->next_sibling != NULL) (pos1->next_sibling)->prev_sibling = newbie;
								if(pos1->prev_sibling != NULL) (pos1->prev_sibling)->next_sibling = newbie;
								pos1->parent = newbie;
								newbie->next_sibling = pos1->next_sibling;
								newbie->prev_sibling = pos1->prev_sibling;
								pos1->next_sibling = pos2;
								pos1->prev_sibling = NULL;
								
								
								if(pos2->next_sibling != NULL) (pos2->next_sibling)->prev_sibling = pos2->prev_sibling;
								if(pos2->prev_sibling != NULL) (pos2->prev_sibling)->next_sibling = pos2->next_sibling;
								
								pos2->next_sibling = NULL;
								pos2->prev_sibling = pos1;
								pos1 = start;
								pos2 = pos1->next_sibling;
								if(topoftree && newbie->prev_sibling == NULL) temp_top2 = newbie;
								newbie = NULL;
								l--;
								donestuff = TRUE;
								/*printf("done easy same\n");*/
								}
							else
								{
								pos2 = pos2->next_sibling;
								y++;
								}
							}
						if(nope) /* if we haven't found any taxa the same as pos 1 , go down through the pairings already made */
							{
							
							pos2 = NULL;
							copy = start;
							y=0;
							while(copy != NULL && pos2 == NULL)
								{
								if(copy->daughter != NULL && copy != pos1 && copy->tag2 == TRUE)
									pos2 = find_same(copy->daughter, pos1->tag);
								copy = copy->next_sibling;
								y++;
								}
							if(pos2 != NULL && pos2 != pos1)
								{
								newbie = make_taxon();
								newbie->tag = pos1->tag;
								newbie->daughter = pos1;
								newbie->tag2=TRUE;
								if(pos1->prev_sibling == NULL)
									{
									if(pos1->parent != NULL )
										(pos1->parent)->daughter = pos1->next_sibling;
									(pos1->next_sibling)->parent = pos1->parent;
									start = pos1->next_sibling;
									}							
								if(pos1->next_sibling != NULL) (pos1->next_sibling)->prev_sibling = pos1->prev_sibling;
								if(pos1->prev_sibling != NULL) (pos1->prev_sibling)->next_sibling = pos1->next_sibling;
								pos1->parent = newbie;
								newbie->next_sibling = pos2->next_sibling;
								newbie->prev_sibling = pos2->prev_sibling;
								if(pos2->next_sibling != NULL) (pos2->next_sibling)->prev_sibling = newbie;
								if(pos2->prev_sibling != NULL) (pos2->prev_sibling)->next_sibling = newbie;
								if(pos2->prev_sibling == NULL)
									{
									if(pos2->parent != NULL)
										(pos2->parent)->daughter = newbie;
									newbie->parent = pos2->parent;
									pos2->parent = NULL;
									if(start == pos2) start = newbie;
									}
								pos1->next_sibling = pos2;
								pos2->prev_sibling = pos1;
								
								pos1->prev_sibling = NULL;
								pos2->next_sibling = NULL;
								l--;
								donestuff = TRUE;
								
								pos1 = start;
								newbie = NULL;
								
								/*printf("done hard same\n");*/
								}
							}
						pos1 = pos1->next_sibling;
						}				
					}
				}
				/* Step 3: reconstruct the LCAs of remaining branches */
			donestuff = FALSE;
			if(((l > 2 && !topoftree) || ( l>3 && topoftree)))
				{
				smallest = 2*number_of_taxa;
				pos1 = start;
				pos2 = pos1->next_sibling;
				x=0;
				while(pos1 != NULL)
					{
					while(pos2 != NULL)
						{
						y=x+1;
						for(j=0; j<2*number_of_taxa; j++)
							presence[j] = 0;
						presence[pos1->tag] = TRUE;
						presence[pos2->tag] = TRUE;
						result =  get_min_node(species_tree, presence, 2*number_of_taxa-2);
						if(result <= smallest)
							{
							smallest = result;
							best1 = pos1;
							best2 = pos2;
							}
						pos2 = pos2->next_sibling;
						y++;
						}
					pos1 = pos1->next_sibling;
					x++;
					}
				
				pos1=best1;
				pos2 = best2;

				newbie = make_taxon();
				newbie->tag = smallest;
				newbie->daughter = pos1;
				newbie->tag2 = TRUE;
				if(pos1->prev_sibling == NULL)
					{
					if(pos1->parent != NULL)
						(pos1->parent)->daughter = newbie;
					newbie->parent = pos1->parent;
					start = newbie;
					}
				if(pos1->next_sibling != NULL) (pos1->next_sibling)->prev_sibling = newbie;
				if(pos1->prev_sibling != NULL) (pos1->prev_sibling)->next_sibling = newbie;
				pos1->parent = newbie;
				newbie->next_sibling = pos1->next_sibling;
				newbie->prev_sibling = pos1->prev_sibling;
				pos1->next_sibling = pos2;
				pos1->prev_sibling = NULL;
				
				if(pos2->next_sibling != NULL) (pos2->next_sibling)->prev_sibling = pos2->prev_sibling;
				if(pos2->prev_sibling != NULL) (pos2->prev_sibling)->next_sibling = pos2->next_sibling;
				
				pos2->next_sibling = NULL;
				pos2->prev_sibling = pos1;
				
				pos1 = start;
				pos2 = pos1->next_sibling;
				if(topoftree && newbie->prev_sibling == NULL) temp_top2 = newbie;
				newbie = NULL;
				l--;
				donestuff = TRUE;				
				/*printf("done lca\n");*/
				}
			}
		
		/* label copy with the number of copies of taxa present */
		

		
		
				
		}
	
	
	free(presence);
	free(found);
	}

void print_tree_labels(struct taxon *position, int **results, int treenum, struct taxon *species_tree)
	{
	int onetoone;
	while(position != NULL)
		{
		if(position->loss >= 1) 
			{
			results[1][position->tag]++;
			if(strcmp(position->weight, "") != 0) results[5][position->tag]++;
			}
		else
			{
			if(position->loss == -1) results[2][position->tag]++;
			else 
				{
				results[0][position->tag]++;
				/** now check to see if what follows is a 1:1 ortholog **/
				if(position->daughter != NULL)
					{
					onetoone = isit_onetoone(position->daughter, 2);
					if(onetoone == 1)
						{
						fprintf(onetoonefile, "Tree num:%d\tName:%s\tBranch:", treenum, tree_names[treenum]);
						check_tree(species_tree, position->tag, onetoonefile);
						fprintf(onetoonefile, "\n");
						results[3][position->tag]++;
						print_onetoone_names(position->daughter, 1);
						fprintf(onetoonefile, "\n");
						}
					if(onetoone == 2)
						{
						fprintf(onetoonefile, "Tree num:%d\tName:%s\tBranch:", treenum, tree_names[treenum]);
						check_tree(species_tree, position->tag, onetoonefile);
						fprintf(onetoonefile, "\n");
						fprintf(strictonetoonefile, "Tree num:%d\tName:%s\tBranch:", treenum, tree_names[treenum]);
						check_tree(species_tree, position->tag, strictonetoonefile);
						fprintf(strictonetoonefile, "\n");
						results[3][position->tag]++;
						results[4][position->tag]++;
						print_onetoone_names(position->daughter, 2);
						fprintf(onetoonefile, "\n");
						fprintf(strictonetoonefile, "\n");
						}
					}
				}
			}
		if(position->daughter != NULL)
			{
			if(position->loss != -1) print_tree_labels(position->daughter, results, treenum, species_tree);
			}
			
		position = position->next_sibling;
		}
	}

void label_gene_tree(struct taxon * gene_position, struct taxon * species_top, int *presence, int xnum)
	{
	struct taxon * position = gene_position, *tmp = NULL;
	int i =0, j=0, latest = -1;
/*	printf2("in Label_gene_tree\n");*/
	while(position != NULL)
		{
		if(position->daughter != NULL) 
			{
	/*		printf2("daughter\n"); */
			label_gene_tree(position->daughter, species_top, presence, xnum);
			for(i=0; i<2*number_of_taxa; i++) presence[i] = FALSE;
			tmp = position->daughter;
			j=0;
			while(tmp != NULL)
				{
				if(tmp->name != -1)
					{
					if(presence[tmp->name] == FALSE) j++;
					presence[tmp->name] = TRUE;
					latest = tmp->name;
					}
				else
					{
					if(presence[tmp->tag] == FALSE) j++;
					presence[tmp->tag] = TRUE;
					latest = tmp->tag;
					}
				tmp = tmp->next_sibling;
				}
			if(j>1)
				{
				for(i=0; i<2*number_of_taxa; i++) presence[i] = FALSE;
				descend(position->daughter, presence);
				position->tag = get_min_node(species_top, presence, xnum);
				}
			else
				position->tag = latest;
				
			}
		else
			position->tag = position->name;
		position = position->next_sibling;
		}
	}

void descend(struct taxon * position, int *presence)
	{
	
	while(position != NULL)
		{
		if(position->daughter != NULL) descend(position->daughter, presence);
		else presence[position->name] = TRUE;
		position = position->next_sibling;
		}
	}

void subtree_id(struct taxon * position, int *tmp)
	{
	while(position != NULL)
		{
		if(position->name != -1) tmp[position->name] = FALSE;
		else
			{
			subtree_id(position->daughter, tmp);
			}
		position = position->next_sibling;
		}
	}

int reconstruct_map(struct taxon *position, struct taxon *species_top)
	{
	struct taxon *start = NULL, *tmp = NULL, *spec_pointer= NULL, *daug_pointer = NULL, *newbie = NULL;
	int found = FALSE, num_dups = 0;
	
	while(position != NULL)
		{
		if(position->daughter != NULL)
			{
			num_dups += reconstruct_map(position->daughter, species_top);
			
			/** Check to see if any of the daughters of this pointer taxa have the same id (or name) as its id **/
			tmp = position->daughter;
			found = FALSE;
			while(tmp != NULL)
				{
				if(tmp->name == -1)
					{
					if(tmp->tag == position->tag) found = TRUE;
					}
				else
					{
					if(tmp->name == position->tag) found = TRUE;
					}
				tmp = tmp->next_sibling;
				}
			
			if(found) /* if position is the site of a duplication event **/
				{
				position->loss = 2;
				num_dups++;
				found = FALSE;
				/*** now check to see if any of positions' daughters are not the same as pointer ***/
				tmp = position->daughter;
				while(tmp != NULL)
					{
					if(tmp->tag != position->tag && tmp->name != position->tag) /** we need to add parts to the tree ***/
						{
						newbie = make_taxon();					
						newbie->tag = position->tag;
						newbie->daughter = tmp;
						newbie->next_sibling = tmp->next_sibling;
						newbie->prev_sibling = tmp->prev_sibling;
						newbie->parent = tmp->parent;
						if(tmp->parent != NULL) (tmp->parent)->daughter = newbie;
						if(tmp->prev_sibling != NULL) (tmp->prev_sibling)->next_sibling = newbie;
						if(tmp->next_sibling != NULL) (tmp->next_sibling)->prev_sibling = newbie;
						tmp->next_sibling = NULL;
						tmp->prev_sibling = NULL;
						tmp->parent = newbie;
						}
					tmp = tmp->next_sibling;
					}
				}
			}
		position = position->next_sibling;
		}
	
	
	
	return(num_dups);
	}

void add_losses(struct taxon * position, struct taxon *species_top)
	{
	struct taxon * start = position, *spec_equiv = NULL, *tmp1 = NULL, *pos = NULL, *pos2 = NULL, *pos3 = NULL;
	int *presence = NULL, i;
	
	presence = malloc((2*number_of_taxa)*sizeof(int));
	for(i=0; i<2*number_of_taxa; i++) presence[i] = FALSE;
	
	/** Go down the gene tree to the bottom */
	while(position != NULL)
		{
		if(position->daughter != NULL) 
			{
			add_losses(position->daughter, species_top);
			/* now for the position we are at, find the equivalent position in the species tree */
			if((position->daughter)->tag != position->tag)  /* there is nothing to do if this is where we have reconstructed a duplication, all the relevant part have already been put in */
				{
				spec_equiv = get_branch(species_top, position->tag);
				
				tmp1 = position->daughter;
				tmp1->parent = NULL;
				position->daughter = NULL;
				
				/* record the nodes (and taxa) that are direct sibilngs to tmp1 */
				for(i=0; i<(2*number_of_taxa); i++) presence[i] = FALSE;
				pos = tmp1;
				while(pos != NULL)
					{
					presence[pos->tag] = TRUE;
					pos = pos->next_sibling;
					}
				/* From our position in the species tree, travel down to the bottom, repilicating what we find along the way (except for the parts that already exist and pointed to by tmp1) */
				tmp1 = construct_tree(spec_equiv->daughter, position, presence, tmp1);
				}
			
			
			}
		position = position->next_sibling;
		}	
	free(presence);
	presence = NULL;
	}

struct taxon * construct_tree(struct taxon * spec_pos, struct taxon *gene_pos, int *presence, struct taxon *extra_gene)
	{
	int first = TRUE, i, count;
	struct taxon * start = spec_pos, *position = NULL, *tmp = gene_pos, *tmp1 = NULL, *newbie = NULL;
	
	while(spec_pos != NULL)
		{
		if(presence[spec_pos->tag])  /* This already existed on the gene tree so we need to splice it back in here */
			{
			
			/*** Find the position on the gene tree extra */
			position = extra_gene;
			count = 0;
			while(position != NULL && position->tag != spec_pos->tag)
					position = position->next_sibling;
			if(position == NULL)
				{
				/* Tag not found — should not happen in a consistent reconciliation; skip this species node */
				spec_pos = spec_pos->next_sibling;
				continue;
				}

			/*** uncouple this from extra_gene and place in the gene tree ***/
			if(position == extra_gene)
				extra_gene = position->next_sibling;
			if(position->prev_sibling != NULL) (position->prev_sibling)->next_sibling = position->next_sibling;
			if(position->next_sibling != NULL) (position->next_sibling)->prev_sibling = position->prev_sibling;
			}
		else /* if it doesn't already exist */
			{
			position = make_taxon();
			if(spec_pos->name != -1) position->tag = spec_pos->name;
			else position->tag = spec_pos->tag;
			position->name = spec_pos->name;
			if(position->name != -1) position->loss = -1;
			}
		if(first)
			{
			tmp->daughter = position;
			position->parent = tmp;
			position->next_sibling = NULL;
			position->prev_sibling = NULL;
			tmp = position;
			first = FALSE;
			}
		else
			{
			tmp->next_sibling = position;
			position->prev_sibling = tmp;
			position->next_sibling = NULL;
			tmp = position;
			}
			
		if(spec_pos->daughter != NULL && !presence[spec_pos->tag]) extra_gene = construct_tree(spec_pos->daughter, tmp, presence, extra_gene);
		
		spec_pos = spec_pos->next_sibling;
		}
	
	
	return(extra_gene);
	}

int join_losses(struct taxon * position)
	{
	int loss = TRUE;
	
	while(position != NULL)
		{
		if(position->daughter != NULL)
			{
			if(join_losses(position->daughter))
				position->loss = -1;
			else
				loss = FALSE;
			}
		else
			{
			if(position->loss != -1) loss = FALSE;
			}
		position = position->next_sibling;
		}
	return(loss);
	}

int count_losses(struct taxon * position)
	{
	int count = 0;
	while(position != NULL)
		{
		if(position->loss == -1) count++;
		else
			{
			if(position->daughter != NULL) count += count_losses(position->daughter);
			}
		position = position->next_sibling;
		}
	return(count);
	}

void mapunknowns()
	{
	struct taxon *position = NULL, *species_tree = NULL, *gene_tree = NULL, *best_mapping = NULL, *unknown_fund = NULL, *pos = NULL,*copy = NULL;
	int i, j, k, l,  *presence = NULL, basescore = 1;
	float *overall_placements = NULL, biggest = -1, total, best_total = -1;
	char *temptree, *temptree1 = malloc(TREE_LENGTH * sizeof(char));
	if(!temptree1) { printf2("Error: out of memory in mapunknowns\n"); return; }
	temptree = malloc(TREE_LENGTH*sizeof(char));
	temptree[0] = '\0';
	temptree1[0] = '\0';



	presence = malloc(2*number_of_taxa*sizeof(int));
	overall_placements = malloc(2*number_of_taxa*sizeof(float));
	for(i=0; i<2*number_of_taxa; i++)
		overall_placements[i] = 0;

	
	/**** BUILD THE SPECIES TREE *****/
	strcpy(temptree, fundamentals[0]);
	returntree(temptree);
	/* build the tree in memory */
	/****** We now need to build the Supertree in memory *******/
	temp_top = NULL;
	tree_build(1, temptree, species_tree, 1, -1, 0);
	species_tree = temp_top;
	temp_top = NULL;
	/** add an extra node to the top of the tree */
	temp_top = make_taxon();
	temp_top->daughter = species_tree;
	species_tree->parent = temp_top;
	species_tree = temp_top;
	temp_top = NULL;
	
	
	
	
	for(l=1; l<Total_fund_trees; l++)
		{
		temp_top = NULL;
		strcpy(temptree, "");
		strcpy(temptree, fundamentals[l]);
		returntree(temptree);
		/* build the tree in memory */
		/****** We now need to build the Supertree in memory *******/
		temp_top = NULL;
		tree_build(1, temptree, gene_tree, 1, -1, 0);
		gene_tree = temp_top;
		temp_top = NULL;
		strcpy(temptree1, temptree);
		
		k=0;
		while(presence_of_taxa[0][k] >0 || presence_of_taxa[l][k] == FALSE) k++;  /* fundamental [0] is the species tree, so whatever is not in there is an unknown */
		/* k is now the number of the unkown taxa */
		printf2("unknown to be mapped: %s\n", taxa_names[k]);
		position = get_taxon(gene_tree, k);

		/** this is an unknown. we need to remove it, but mark its neighbors */
		pos = position;
		while(pos->prev_sibling != NULL) pos = pos->prev_sibling;
		while(pos != NULL)
			{
			pos->tag2 = TRUE;
			pos = pos->next_sibling;
			}
		
		/*** remove the unknown ****/
		if(position->next_sibling != NULL) (position->next_sibling)->prev_sibling = position->prev_sibling;
		if(position->prev_sibling != NULL) (position->prev_sibling)->next_sibling = position->next_sibling;
		if(position->parent != NULL) 
			{
			(position->parent)->daughter = position->next_sibling;
			(position->next_sibling)->parent = position->parent;
			}
		
		free(position);
		position = NULL;
		compress_tree(gene_tree);

		if(presence_of_trichotomies(gene_tree)) gene_tree = do_resolve_tricotomies(gene_tree, species_tree, basescore);
		
		duplicate_tree(gene_tree, NULL);
		copy = temp_top;
		tree_top = NULL;
		i = number_tree(gene_tree, 0);
		best_total = -1;
		for(j=0; j<i; j++)
			{
			
			position = get_branch(gene_tree, j);
			tree_top = gene_tree;
			reroot_tree(position);
			gene_tree = tree_top;
			tree_top = NULL;
			printf2("a==\n");
			total = tree_map(gene_tree, species_tree,0);
			
			
			if(total < best_total || best_total == -1)
				{
				best_total = total;
				if(best_mapping != NULL)
					{
					dismantle_tree(best_mapping);
					best_mapping = NULL;
					}
				best_mapping = gene_tree;
				gene_tree = NULL;
				}
			else
				{
				dismantle_tree(gene_tree);
				gene_tree = NULL;
				}
			temp_top = NULL;
			tree_top = NULL;
			duplicate_tree(copy, NULL);
			gene_tree = temp_top;
			temp_top = NULL;
			tree_top = NULL;
			number_tree(gene_tree, 0);
			}
			
		printf2("best reconstruction with a score of %f\n", best_total);
		tree_top = best_mapping;
		tree_coordinates(temptree, TRUE, FALSE, FALSE, -1);
		
		for(i=0; i<2*number_of_taxa; i++)
			presence[i] = FALSE;
			
		find_tagged(best_mapping, presence);
		j=0;
		for(i=0; i<2*number_of_taxa; i++)
			if(presence[i]) j++;
		for(i=0; i<2*number_of_taxa; i++)
			{
			if(presence[i]) overall_placements[i] += (float)1/(float)j;
			}
		if(copy != NULL)
			{
			dismantle_tree(copy);
			copy = NULL;
			}
		
		if(gene_tree != NULL)
			{
			dismantle_tree(gene_tree);
			gene_tree = NULL;
			}
		if(best_mapping != NULL)
			{
			dismantle_tree(best_mapping);
			best_mapping = NULL;
			}
		}
	
	biggest = -1;
	for(i=0; i<2*number_of_taxa; i++)
		{
		if(overall_placements[i] > biggest) biggest = overall_placements[i];
		}
	
	for(i=0; i<2*number_of_taxa; i++)
		overall_placements[i] = overall_placements[i]/biggest;
	put_in_scores(species_tree, overall_placements);
	
	tree_top = species_tree;
	tree_coordinates(temptree, TRUE, FALSE, TRUE, -1);
	
	if(species_tree != NULL)
		{
		dismantle_tree(species_tree);
		species_tree = NULL;
		}
	
		
	tree_top = NULL;
	free(temptree);
	free(temptree1);
	free(presence);
	free(overall_placements);
	}

float get_recon_score(char *giventree, int numspectries, int numgenetries)
	{
	struct taxon *position = NULL, *species_tree = NULL, *gene_tree = NULL, *best_mapping = NULL, *copy = NULL, *temp_top1 = NULL, *temp_top2 = NULL, *spec_copy = NULL;
	int i, j, k, l, m, q, r, spec_start=0, spec_end, gene_start, gene_end, num_species_internal = 0, error = FALSE, num_species_roots = 0, basescore = 1, rand1=0, rand2=0, dospecrand = 1, dogenerand=1, taxaorder=0;
	float *overall_placements = NULL, biggest = -1, total, best_total = -1, sum_of_totals = 0, rooting_score = -1;
	char *temptree, *temptree1 = malloc(TREE_LENGTH * sizeof(char));
	if(!temptree1) { printf2("Error: out of memory in get_recon_score\n"); return -1; }

	temptree = malloc(TREE_LENGTH*sizeof(char));
	temptree[0] = '\0';
	temptree1[0] = '\0';



		/**** BUILD THE SPECIES TREE *****/
		strcpy(temptree,giventree);
		returntree(temptree);
		/* build the tree in memory */
		/****** We now need to build the Supertree in memory *******/
		temp_top = NULL;
		tree_build(1, temptree, species_tree, 1, -1, 0);
		species_tree = temp_top;
		temp_top = NULL;

		num_species_roots = number_tree1(species_tree, number_of_taxa);

		duplicate_tree(species_tree, NULL);
		spec_copy = temp_top;
		temp_top2 = NULL;
		
		
		if(numspectries == -1)
			{
			dospecrand = FALSE;
			numspectries = 1;
			}
		if(numgenetries == -1)
			{
			dospecrand = FALSE;
			numgenetries = 1;
			}
			
			
		for(q=0; q<numspectries; q++)
			{
			if(dospecrand > 0)
				{
				spec_start = (int)fmod(rand(), num_species_roots);
				spec_end = spec_start +1;
				}
			else
				{
				spec_start = 0;
				spec_end = num_species_roots;
				}
/*		for(m=0; m<num_species_roots; m++) */ /* for every rooting of the species tree */
			for(m=spec_start; m<spec_end; m++) 
				{
			
				position = get_branch(species_tree, m);
				temp_top = species_tree;
				/*printf("1\n");*/
				reroot_tree(position);
				species_tree = temp_top;
				temp_top = NULL;

				for(l=0; l<Total_fund_trees; l++)  /* for every gene tree */
					{
					temp_top1 = NULL;
					strcpy(temptree, "");
					strcpy(temptree, fundamentals[l]);
					unroottree(temptree);
					
					returntree(temptree);
					/* build the tree in memory */
					/****** We now need to build the genetree in memory *******/
					temp_top = NULL;
					taxaorder=0;
					tree_build(1, temptree, gene_tree, 1, l, 0);
					gene_tree = temp_top;
					temp_top = NULL;
					strcpy(temptree1, temptree);

					if(presence_of_trichotomies(gene_tree)) gene_tree = do_resolve_tricotomies(gene_tree, species_tree, basescore);	
				/*	else printf2("no resolving needed\n"); */


					duplicate_tree(gene_tree, NULL);
					copy = temp_top;
					temp_top2 = NULL;
					i = number_tree(gene_tree, 0);
					best_total = -1;
					for(r=0; r<numgenetries; r++)
						{
						if(dogenerand > 0)
							{
							gene_start = (int)fmod(rand(), i);
							gene_end = gene_start +1;
							}
						else
							{
							gene_start = 0;
							gene_end = i;
							}						
						rand2 = (int)fmod(rand(), i);
					/*	for(j=0; j<i; j++)  */ /* For every rooting of the genetree */
						for(j=gene_start; j<gene_end; j++)   
							{
							position = get_branch(gene_tree, j);
							temp_top = gene_tree;
							/*printf("2\n"); */
							reroot_tree(position);
							gene_tree = temp_top;
							temp_top = NULL;
							total = tree_map(gene_tree, species_tree,0);
							if(total < best_total || best_total == -1)
								{
								best_total = total;
								if(best_mapping != NULL)
									{
									dismantle_tree(best_mapping);
									best_mapping = NULL;
									}
								best_mapping = gene_tree;
								gene_tree = NULL;
								}
							else
								{
								dismantle_tree(gene_tree);
								gene_tree = NULL;
								}
							temp_top1 = NULL;
							temp_top = NULL;
							duplicate_tree(copy, NULL);
							gene_tree = temp_top;
							temp_top1 = NULL;
							temp_top = NULL;
							number_tree(gene_tree, 0);
							}
						}
					if(copy != NULL)
						{
						dismantle_tree(copy);
						copy = NULL;
						}
					
					if(gene_tree != NULL)
						{
						dismantle_tree(gene_tree);
						gene_tree = NULL;
						}
					if(best_mapping != NULL)
						{
						dismantle_tree(best_mapping);
						best_mapping = NULL;
						}
					sum_of_totals+=best_total;
					}
				if(sum_of_totals < rooting_score || rooting_score == -1)
					{
					rooting_score = sum_of_totals;
					}
				sum_of_totals = 0;
				dismantle_tree(species_tree);
				duplicate_tree(spec_copy, NULL);
				species_tree = temp_top;
				temp_top2 = NULL;

				}
			}
	if(spec_copy != NULL)
		{
		dismantle_tree(spec_copy);
		spec_copy = NULL;
		}
	if(species_tree != NULL)
		{
		dismantle_tree(species_tree);
		species_tree = NULL;
		}
	free(temptree);
	free(temptree1);
	return(rooting_score);
	}

int presence_of_trichotomies(struct taxon * position)
	{
	struct taxon *start = position;
	int i=0, result = FALSE;
	while(position != NULL)
		{
		i++;
		position = position->next_sibling;
		}
	if(i > 2)
		{
		if(start->parent != NULL)
			result = TRUE;
		else
			{
			if(i > 3)
				result = TRUE;
			}
		}

	position = start;
	while(position != NULL && result == FALSE)
		{
		if(position->daughter != NULL)	
			{
			if(presence_of_trichotomies(position->daughter))
				result = TRUE;
			}
		position = position->next_sibling;
		}
	return(result);
	}

struct taxon * do_resolve_tricotomies(struct taxon * gene_tree, struct taxon * species_tree, int basescore)
	{
	int i, j, xnum=0, *presence=NULL, **scores = NULL;
	char *temper=NULL;
	struct taxon * teeemp = NULL;
	
	temper = malloc(TREE_LENGTH*sizeof(int));
	temper[0] = '\0';
	
	scores = malloc((2*number_of_taxa)*sizeof(int *));
	for(i=0; i<2*number_of_taxa; i++)
		{
		scores[i] = malloc((2*number_of_taxa)*sizeof(int));
		for(j=0; j<2*number_of_taxa; j++)
			{
			scores[i][j] = basescore;
			}
		}
	
	presence = malloc(2*number_of_taxa*sizeof(int));
	for(i=0; i<(2*number_of_taxa); i++) presence[i] = FALSE;
	/** 1) Label all internal and external taxa on the species tree ****/

	xnum = number_tree1(species_tree, number_of_taxa);
	xnum--;
	/****2) label all the gene tree nodes (and taxa) with their equivalent on the species tree **/
	label_gene_tree(gene_tree, species_tree, presence, xnum);
	temp_top2 = gene_tree;
	/** 3) calculate the distance between all the branches (internal nad external) on the species tree (using the number of nodes as the criterion) */
	/* print out the species tree with the numbered labels on the internal branches */
	print_tree_withinternals(species_tree, temper);
	strcat(temper, ";");
	/* calculate the pathdistance between the branches (internal nad external) on the species tree */
	pathmetric_internals(temper, species_tree, scores);

/*	for(i=0; i<2*number_of_taxa; i++)	
		{
		if(i<number_of_taxa)
			printf2("\t%s", taxa_names[i]);
		else
			printf2("\t%d", i);
		}
	printf2("\n");
	for(i=0; i<2*number_of_taxa; i++)	
		{
		if(i<number_of_taxa)
			printf2("%s\t", taxa_names[i]);
		else
			printf2("%d\t", i);
		for(j=0; j<2*number_of_taxa; j++)
			printf2("%d\t", scores[i][j]);
		printf2("\n");
		}
*/
	/*resolve_tricotomies(gene_tree, species_tree);*/
	temp_top2 = gene_tree;
	resolve_tricotomies_dist(gene_tree, species_tree, scores);
	gene_tree = temp_top2;
	
	/**** Unroot the tree -necessary as the rooting algorithm assumes an unrooted topology */
	/*count the number of daughters at the top of the tree */
	i=0;
	teeemp = gene_tree;
	while(teeemp != NULL)
		{
		i++;
		teeemp = teeemp->next_sibling;
		}
	if(i < 3)
		{
		teeemp = NULL;
		if(gene_tree->daughter != NULL)
			{
			temp_top2 = gene_tree->daughter;
			(gene_tree->next_sibling)->next_sibling = temp_top2;
			(gene_tree->next_sibling)->prev_sibling = NULL;
			temp_top2->prev_sibling = (gene_tree->next_sibling);
			temp_top2->parent = NULL;
			temp_top2 = temp_top2->prev_sibling;
			gene_tree->next_sibling = NULL;
			gene_tree->daughter = NULL;
			dismantle_tree(gene_tree);
			gene_tree = temp_top2;
			}
		else
			{
			temp_top2 = (gene_tree->next_sibling)->daughter;
			teeemp = gene_tree->next_sibling;
			gene_tree->next_sibling = temp_top2;
			temp_top2->prev_sibling = gene_tree;
			teeemp->prev_sibling = NULL;
			teeemp->daughter = NULL;
			dismantle_tree(teeemp);
			teeemp = NULL;
			}
		}
	
	free(presence);
	free(temper);
	for(i=0; i<2*number_of_taxa; i++)
		free(scores[i]);
	free(scores);
	
	return(gene_tree);
	}

void make_unrooted(struct taxon * position)
	{
	int i=0;
	struct taxon * start = position;
	while(position != NULL)
		{
		i++;
		position = position->next_sibling;
		}
	if(i==1)
		{
		
		}
	}

void reconstruct(int print_settings)  /* Carry out gene-tree reconciliation of source trees against a species tree */
	{
	struct taxon *position = NULL, *species_tree = NULL, *gene_tree = NULL, *best_mapping = NULL, *unknown_fund = NULL, *pos = NULL, *copy = NULL, *newbie = NULL;
	int i, j, k, l, xnum=0, *presence = NULL, **label_results = NULL, num_species_internal = 0, error = FALSE, printfiles = FALSE, how_many = 0, diff_overall =0, dorecon = FALSE, basescore = 1, taxaorder=0, species_source = 0, gene_tree_start = 0;
	float *overall_placements = NULL, biggest = -1, total, best_total = -1;
	char *temptree;
	char *temptree1 = malloc(TREE_LENGTH * sizeof(char));
	char *temptree2 = malloc(TREE_LENGTH * sizeof(char));
	char reconfilename[100], otherfilename[100], *tmp1 = NULL, c = '\0', speciestree_file[1000];
	if(!temptree1 || !temptree2) { free(temptree1); free(temptree2); printf2("Error: out of memory in reconstruct\n"); return; }
	FILE *reconstructionfile = NULL, *descendentsfile = NULL, *genebirthfile = NULL;
	
	temptree = malloc(TREE_LENGTH*sizeof(char));
	temptree[0] = '\0';
	speciestree_file[0] = '\0';
	temptree1[0] = '\0';
	reconfilename[0] ='\0';
	otherfilename[0] = '\0';
	
	tmp1 = malloc(TREE_LENGTH*sizeof(char));
	tmp1[0] = '\0';
	count_now = TRUE;

    for(i=0; i<num_commands; i++)
        {
        if(strcmp(parsed_command[i], "printfiles") == 0)
			{
			if(strcmp(parsed_command[i+1], "yes") == 0)
				printfiles = TRUE;
			}
		if(strcmp(parsed_command[i], "speciestree") == 0)
			{
			/* values: "memory", "first", or a filename */
			if(strcmp(parsed_command[i+1], "memory") == 0)
				species_source = 1;
			else if(strcmp(parsed_command[i+1], "first") == 0)
				species_source = 2;
			else
				{
				species_source = 3;
				strcpy(speciestree_file, parsed_command[i+1]);
				}
			}
		if(strcmp(parsed_command[i], "showrecon") == 0)
			{
			if(strcmp(parsed_command[i+1], "yes") == 0)
				dorecon = TRUE;
			}
		if(strcmp(parsed_command[i], "basescore") == 0)
			{
			basescore = atoi(parsed_command[i+1]);
			}
		}



	/* If autoprunemono was active, swap original (unpruned) trees in for reconstruct */
	int autoprunemono_swapped = FALSE;
	if(autoprunemono_active && original_fundamentals != NULL)
		{
		for(i = 0; i < Total_fund_trees; i++)
			{
			if(original_fundamentals[i] != NULL)
				{
				char *tmp_swap = fundamentals[i];
				fundamentals[i] = original_fundamentals[i];
				original_fundamentals[i] = tmp_swap;
				}
			}
		autoprunemono_swapped = TRUE;
		printf2("(Using original unpruned trees for reconstruct)\n");
		}

	if(printfiles)
		{
		strcpy(reconfilename, inputfilename);
		strcat(reconfilename, ".recon");
		reconstructionfile =fopen(reconfilename, "w");
		strcat(reconfilename, ".dist");
		distributionreconfile = fopen(reconfilename, "w");
		
		strcpy(otherfilename, inputfilename);
		strcat(otherfilename, ".species.decendents");
		descendentsfile = fopen(otherfilename, "w");

		strcpy(otherfilename, inputfilename);
		strcat(otherfilename, ".gene.births");
		genebirthfile = fopen(otherfilename, "w");	
		
		strcpy(otherfilename, inputfilename);
		strcat(otherfilename, ".onetoone");
		onetoonefile =  fopen(otherfilename, "w");

		strcpy(otherfilename, inputfilename);
		strcat(otherfilename, ".strict.onetoone");
		strictonetoonefile =  fopen(otherfilename, "w");

		
		fprintf(distributionreconfile, "Duplications\tLosses\n");
		}
	presence = malloc(2*number_of_taxa*sizeof(int));
	overall_placements = malloc(2*number_of_taxa*sizeof(float));
	for(i=0; i<2*number_of_taxa; i++)
		overall_placements[i] = 0;
	
	label_results = malloc(6*sizeof(int*));
	for(i=0; i<6; i++)
		{
		label_results[i] = malloc((2*number_of_taxa)*sizeof(int));
		for(j=0; j<(2*number_of_taxa); j++)
			{
			label_results[i][j] = 0;
			}
		}
	
	 for(i=0; i<num_commands; i++)
		{

		if(strcmp(parsed_command[i], "dups") == 0)
			{
			dup_weight = tofloat(parsed_command[i+1]);
			if(dup_weight < 0)
				{
				printf2("Error: '%s' is an invalid value for dups\n", parsed_command[i+1]);
				error = TRUE;
				}
			}
		if(strcmp(parsed_command[i], "losses") == 0)
			{
			loss_weight = tofloat(parsed_command[i+1]);
			if(loss_weight < 0)
				{
				printf2("Error: '%s' is an invalid value for losses\n", parsed_command[i+1]);
				error = TRUE;
				}
			}
		}
	/* Resolve species tree source and validate */
	if(species_source == 0)
		{
		if(trees_in_memory > 0)
			species_source = 1;
		else
			{
			printf2("Error: No supertree is in memory. Run 'hs', 'nj' or 'alltrees' first,\n");
			printf2("       or specify a species tree with: reconstruct speciestree <file>\n");
			error = TRUE;
			}
		}
	if(species_source == 1 && trees_in_memory == 0)
		{
		printf2("Error: No supertree in memory for 'speciestree memory'.\n");
		error = TRUE;
		}
	if(species_source == 3 && !error)
		{
		FILE *stcheck = fopen(speciestree_file, "r");
		if(stcheck == NULL)
			{
			printf2("Error: Cannot open species tree file '%s'\n", speciestree_file);
			error = TRUE;
			}
		else fclose(stcheck);
		}
	gene_tree_start = (species_source == 2) ? 1 : 0;
	if(!error)
		{
		if(print_settings)printf2("\nReconstruction of Most likely Duplications and Losses\n\tDuplication weight = %f\n\tLosses weight = %f\n\nTree name:\tReconstruction score\n", dup_weight, loss_weight);
		/**** BUILD THE SPECIES TREE *****/
		if(species_source == 1)
			{
			/* Use supertree from memory -- already has actual taxon names */
			strcpy(temptree, retained_supers[0]);
			printf2("Using supertree in memory as species tree.\n");
			}
		else if(species_source == 3)
			{
			/* Read species tree from user-specified file */
			FILE *stfile = fopen(speciestree_file, "r");
			int si = 0; int sc;
			while((sc = getc(stfile)) != EOF && sc != ';') { temptree[si++] = (char)sc; }
			temptree[si++] = ';'; temptree[si] = '\0';
			fclose(stfile);
			printf2("Using species tree from file: %s\n", speciestree_file);
			}
		else
			{
			/* Legacy: use first source tree (fundamentals[0]) */
			strcpy(temptree, fundamentals[0]);
			returntree(temptree);
			printf2("Using first source tree as species tree (legacy mode).\n");
			}
		/* build the tree in memory */
		/****** We now need to build the Species tree in memory *******/
		temp_top = NULL;
		tree_build(1, temptree, species_tree, TRUE, -1, 0); 
		species_tree = temp_top;
		temp_top = NULL;
		/** add an extra node to the top of the tree */
		temp_top = make_taxon();
		temp_top->daughter = species_tree;
		species_tree->parent = temp_top;
		species_tree = temp_top;
		temp_top = NULL;
		number_tree1(species_tree, number_of_taxa);
		num_species_internal = count_internal_branches(species_tree, 0);
		if(printfiles)
			{
			fprintf(reconstructionfile, "Tree name\t");
			for(j=0; j<number_of_taxa+num_species_internal; j++)
				{
				check_tree(species_tree, j, reconstructionfile);
				}
			fprintf(reconstructionfile, "\n");
			
			fprintf(genebirthfile, "tree name\tgene birth\n");
			/* print out the descendents of each internal branch of the tree */
			do_descendents(species_tree, descendentsfile);
			fclose(descendentsfile);
			}
		
		
		for(l=gene_tree_start; l<Total_fund_trees; l++)
			{
			temp_top = NULL;
			strcpy(temptree, "");
			strcpy(temptree, fundamentals[l]);
			unroottree(temptree);
			
			/*returntree(temptree); */ /* add in the actual names in the tree */
			/* build the tree in memory */
			/****** We now need to build the gene tree in memory *******/
			temp_top = NULL;
			how_many++;
			taxaorder=0;
			tree_build(1, temptree, gene_tree, FALSE, l, 0);
			gene_tree = temp_top;

			temp_top = NULL;
			strcpy(temptree1, temptree);
			i = count_taxa(gene_tree, 0);
			if(presence_of_trichotomies(gene_tree))
				{
				printf("Doing resolving of trichotomies\n");
				gene_tree = do_resolve_tricotomies(gene_tree, species_tree, basescore);
				}
			tree_top = gene_tree;
			compress_tree1(gene_tree);
			gene_tree = tree_top;
			temp_top = NULL;
			duplicate_tree(gene_tree, NULL);
			how_many++;
			copy = temp_top;
			tree_top = NULL;
			i = number_tree(gene_tree, 0); /* number every branch of the gene tree so we can use this to reroot */
			best_total = -1;

			if(i>2) /* if there are more than two branches in this gene tree */
				{
				for(j=0; j<i; j++) /* For every possible rerooting of the gene tree */
					{

					malloc_check=0;
					position = get_branch(gene_tree, j); /* find branch numbered j on the gene tree */ 
					temp_top = gene_tree;
					reroot_tree(position); /* reroot at this position */

					/* add a node at the top of the re-rooted tree */
					gene_tree = make_taxon();
					gene_tree->daughter = temp_top;

					
					diff_overall -= malloc_check;
					total = tree_map(gene_tree, species_tree, printfiles); 
					diff_overall += malloc_check;
					if(total < best_total || best_total == -1)
						{
						best_total = total;
						if(best_mapping != NULL)
							{
							dismantle_tree(best_mapping);
							how_many--;
							best_mapping = NULL;
							}
						best_mapping = gene_tree;
						gene_tree = NULL;
						}
					else
						{
						dismantle_tree(gene_tree);
						how_many--;
						gene_tree = NULL;
						}
					temp_top = NULL;
					tree_top = NULL;
					duplicate_tree(copy, NULL);
					how_many++;
					gene_tree = temp_top;
					temp_top = NULL;
					tree_top = NULL;
					number_tree(gene_tree, 0);
					another_check += malloc_check;
					}
				}
			else
				{
				newbie = make_taxon();
				newbie->daughter=gene_tree;
				gene_tree->parent = newbie;
				gene_tree = newbie;
				newbie = NULL;
				
				total = tree_map(gene_tree, species_tree,printfiles);				
				best_total = total;
				best_mapping = gene_tree;
				gene_tree = NULL;
				}
			if(strcmp(tree_names[l], "") == 0)
				printf2("Tree number: %d\t%f\n", l, best_total);
			else
				printf2("%s\t%f\n", tree_names[l], best_total);
			tree_top = best_mapping;
			/** ADD IN PRINTFULLNAMED TREE HERE **/
			temptree[0] = '\0';
			print_fullnamed_tree(tree_top, temptree, l); 
			/*print_named_tree(tree_top, temptree); */

			if(dorecon == TRUE)
				{
				tree_coordinates(temptree, TRUE, FALSE, FALSE, l); /* char *tree, int bootstrap, int build, int mapping, int fundnum */
				}
			/******* PRINT_TREE_LABELS TEST *******/	
			/****	label_results[0] = number of copies of the gene at this internal branch
					label_results[1] = number of duplications at this internal branch
					label_results[2] = number of losses at this internal branch
					label_results[3] = number of 1:1 orthologs after this internal branch (allowing duplications and losses in external taxa)
					label_results[4] = number of strict 1:1 orthologs after this internal branch (not allowing any duplications and losses)
					label_results[5] = (added Aug 17) Number of Duplications at this internal branch supported by splits in the orignal gene tree (Does not include those inferred from polytomies)
			****/

			if(printfiles)
				{	
				for(i=0; i<6; i++)
					{
					for(j=0; j<(2*number_of_taxa); j++)
						{
						label_results[i][j] = 0;
						}
					}
				
				if(tree_top->tag == number_of_taxa+num_species_internal-1)
					fprintf(genebirthfile, "%s\troot\n", tree_names[l]);
				else
					{
					fprintf(genebirthfile, "%s\t", tree_names[l]);
					check_tree(species_tree, tree_top->tag, genebirthfile);
					fprintf(genebirthfile, "\n");
					}
					
				print_tree_labels(tree_top, label_results, l, species_tree);
				fprintf(reconstructionfile, "%10s\t", tree_names[l]);
				for(j=0; j<(2*number_of_taxa)-1; j++)
					{
					
					fprintf(reconstructionfile, "%5d (%5d+ %5d+S %5d- %5d 1:1 %5d 1:1S)\t", label_results[0][j], label_results[1][j], label_results[5][j], label_results[2][j], label_results[3][j], label_results[4][j]);
					}
				fprintf(reconstructionfile, "\n");
				}
			/******* END TEST ********/
			for(i=0; i<2*number_of_taxa; i++)
				presence[i] = FALSE;
				
			find_tagged(best_mapping, presence);
			j=0;
			for(i=0; i<2*number_of_taxa; i++)
				if(presence[i]) j++;
			for(i=0; i<2*number_of_taxa; i++)
				{
				if(presence[i]) overall_placements[i] += (float)1/(float)j;
				}
			if(copy != NULL)
				{
				dismantle_tree(copy);
				copy = NULL;
				}
				how_many--;
			
			if(gene_tree != NULL)
				{
				dismantle_tree(gene_tree);
				gene_tree = NULL;
				}
			how_many--;
			if(best_mapping != NULL)
				{
				dismantle_tree(best_mapping);
				best_mapping = NULL;
				}
				
				
			
			how_many--;
			how_many--;
			}
		
		if(species_tree != NULL)
			{
			dismantle_tree(species_tree);
			species_tree = NULL;
			}
		}
	if(printfiles)
		{
		fclose(distributionreconfile);
		fclose(reconstructionfile);	
		fclose(genebirthfile);
		fclose(onetoonefile);
		fclose(strictonetoonefile);
		}
	tree_top = NULL;
	free(temptree);
	free(temptree1);
	free(temptree2);
	free(tmp1);
	free(presence);
	free(overall_placements);
	if(label_results != NULL)
		{
		for(i=0; i<6; i++)
			{
			free(label_results[i]);
			label_results[i] = NULL;
			}
		free(label_results);
		label_results = NULL;
		}

	/* Swap original trees back so other commands (hs, etc.) continue to use pruned versions */
	if(autoprunemono_swapped)
		{
		for(i = 0; i < Total_fund_trees; i++)
			{
			if(original_fundamentals[i] != NULL)
				{
				char *tmp_swap = fundamentals[i];
				fundamentals[i] = original_fundamentals[i];
				original_fundamentals[i] = tmp_swap;
				}
			}
		}

	count_now=FALSE;

	}

void put_in_scores(struct taxon * position, float * total)
	{
	while(position != NULL)
		{
		if(position->daughter != NULL) put_in_scores(position->daughter, total);
		if(total[position->tag] == 0) position->loss = -1;
		else position->loss = total[position->tag];
		position = position->next_sibling;
		}
	}

void assign_hgtdonors(struct taxon * position, int num, int part_num)
	{
	int i;
	while(position != NULL)
		{
		if(position->daughter != NULL) assign_hgtdonors(position->daughter, num, part_num);
		if(position->donor == NULL)
			{
			position->donor = malloc(num*sizeof(int));
			for(i=0; i<num; i++) position->donor[i] = FALSE;
			}
		if(position->tag2 == TRUE)
			{
			position->donor[part_num] = TRUE;
			}
		position = position->next_sibling;
		}
	}

void hgt_reconstruction()
	{
	struct taxon *position = NULL, *species_tree = NULL, *gene_tree = NULL, *best_mapping = NULL, *best_mapping1 = NULL, *best_mapping2 = NULL, *unknown_fund = NULL, *posit = NULL,*copy = NULL, *copy1 = NULL, **parts = NULL, *test_part = NULL, *pos = NULL, *best_donor = NULL, *best_HGT = NULL, *attached = NULL;
	int i, j, k, l,  *presence = NULL,*presence1 = NULL, *presence2 = NULL, hgt_receipient1, hgt_receipient2, **overall_presence = NULL, *overall_reconstruction = NULL, *overall_receptor = NULL, receptor, **tmp_presence1 = NULL, **tmp_presence2 = NULL, *before1 = NULL, *before2 = NULL, *after1 = NULL, *after2 = NULL, *temporary = NULL;
	float *overall_placements = NULL, biggest = -1,  total, best_total = -1,  best_total1 = -1, best_total2 = -1, HGT1 = 0, HGT2 = 0, original = 0, best_HGT_recon = -1, best_reconstruction = -1, sum, HGT_score = -1, donor_score = -1, tmp_allow = FALSE;
	char *temptree = NULL;
	char *temptree1 = malloc(TREE_LENGTH * sizeof(char));
	int **species_allowed = NULL, **dependent_species_allowed = NULL, *previous = NULL, xnum =0, x, y, z, partA, partB, q, r, s, allow_HGT1 = TRUE, allow_HGT2 = TRUE, numparts = 1,  place_marker = 1, found_better = FALSE, error = FALSE, taxaorder=0;
	int basescore = 1; /** see reconstrution command **/

	if(!temptree1) { printf2("Error: out of memory in recon_hgt\n"); return; }
	temptree = malloc(TREE_LENGTH*sizeof(char));
	temptree[0] = '\0';
	temptree1[0] = '\0';
	
	for(i=0; i<num_commands; i++)
		{
		if(strcmp(parsed_command[i], "dupweight") == 0)
			{
			dup_weight = tofloat(parsed_command[i+1]);
			}
		if(strcmp(parsed_command[i], "lossweight") == 0)
			{
			loss_weight = tofloat(parsed_command[i+1]);
			}
		if(strcmp(parsed_command[i], "hgtweight") == 0)
			{
			hgt_weight = tofloat(parsed_command[i+1]);
			}
		
		}
	
	printf2("Calculating best reconstruction of Duplications, losses and Horizontal gene transfers (HGT)\n\nDuplication weight set to %f\nLoss weight set to %f\nHGT weight set to %f\n\n", dup_weight, loss_weight, hgt_weight);
	
	tree_top = NULL;
	
	/**** First, take the species tree and calculate those part of the tree that would constitute travelling in time if a HGT was to occur (ie in ancestral or descendent branches) ***/
	species_allowed = malloc((2*number_of_taxa)* sizeof(int*));
	if(species_allowed == NULL) memory_error(120);
	for(i=0; i<(2*number_of_taxa); i++)
		{
		species_allowed[i] = malloc((2*number_of_taxa)*sizeof(int));
		if(species_allowed[i] == NULL) memory_error(121);
		for(j=0; j<(2*number_of_taxa); j++) species_allowed[i][j] = TRUE;
		}
	
	dependent_species_allowed = malloc((2*number_of_taxa)* sizeof(int*));
	if(!dependent_species_allowed) memory_error(122);
	for(i=0; i<(2*number_of_taxa); i++)
		{
		dependent_species_allowed[i] = malloc((2*number_of_taxa)*sizeof(int));
		if(dependent_species_allowed[i] == NULL) memory_error(123);
		for(j=0; j<(2*number_of_taxa); j++) dependent_species_allowed[i][j] = TRUE;
		}
	
	temporary = malloc(2*number_of_taxa*sizeof(int));
	presence = malloc(2*number_of_taxa*sizeof(int));
	presence1 = malloc(2*number_of_taxa*sizeof(int));
	presence2 = malloc(2*number_of_taxa*sizeof(int));
	
	previous = malloc((2*number_of_taxa)*sizeof(int));
	for(i=0; i<2*number_of_taxa; i++) previous[i] = FALSE;
	
	before1 = malloc((2*number_of_taxa)*sizeof(int));
	before2 = malloc((2*number_of_taxa)*sizeof(int));
	after1 = malloc((2*number_of_taxa)*sizeof(int));
	after2 = malloc((2*number_of_taxa)*sizeof(int));
	for(i=0; i<2*number_of_taxa; i++)
		{
		before1[i] = before2[i] = after1[i] = after2[i] = FALSE;
		}
	
	/**** BUILD THE SPECIES TREE *****/
	strcpy(temptree, fundamentals[0]);
	returntree(temptree);
	/* build the tree in memory */
	/****** We now need to build the Species tree in memory *******/
	temp_top = NULL;
	tree_build(1, temptree, species_tree, 1, -1, 0);
	species_tree = temp_top;
	temp_top = NULL;
	/** add an extra node to the top of the tree */
	temp_top = make_taxon();
	temp_top->daughter = species_tree;
	species_tree->parent = temp_top;
	species_tree = temp_top;
	temp_top = NULL;
			printf2("2\n");

	xnum = number_tree1(species_tree, number_of_taxa); /* label the internal and external branches of the species tree */
	assign_ances_desc(species_tree, species_allowed, previous);
	free(previous);
	
	
	for(l=1; l<Total_fund_trees; l++)
		{
		best_reconstruction = -1;
		best_donor = NULL;
		best_HGT = NULL;
		temp_top = NULL;
		strcpy(temptree, "");
		strcpy(temptree, fundamentals[l]);
		returntree(temptree);
		unroottree(temptree);  /* we need the gene tree to be unrooted when we are breaking it up */
		/* build the tree in memory */
		/****** We now need to build the genetree in memory *******/
		temp_top = NULL;
		taxaorder=0;
		tree_build(1, temptree, gene_tree, 1, l, 0);
		gene_tree = temp_top;
		temp_top = NULL;
			
		strcpy(temptree1, temptree);
	
		/***** NOW THAT WE HAVE BUILT THE GENE TREE, WE NEED TO BREAK IT AT EVERY BRANCH IN TURN ****/

		num_gene_nodes = number_tree(gene_tree, 0);
		
		overall_receptor = malloc(num_gene_nodes*sizeof(int));
		for(i=0; i<num_gene_nodes; i++) overall_receptor[i] = -1;
		
		parts = malloc(num_gene_nodes*sizeof(struct taxon *));
		for(i=0; i<num_gene_nodes; i++) parts[i] = NULL;
		
		overall_placements = malloc(num_gene_nodes*sizeof(float));
		for(i=0; i<num_gene_nodes; i++)
			overall_placements[i] = 0;
		
		tmp_presence1 = malloc(num_gene_nodes*sizeof(int *));
		for(i=0; i<num_gene_nodes; i++) 
			{
			tmp_presence1[i] = malloc((2*number_of_taxa)*sizeof(int *));
			for(j=0; j<2*number_of_taxa; j++) tmp_presence1[i][j] = FALSE;
			}
		
		tmp_presence2 = malloc(num_gene_nodes*sizeof(int *));
		for(i=0; i<num_gene_nodes; i++) 
			{
			tmp_presence2[i] = malloc((2*number_of_taxa)*sizeof(int *));
			for(j=0; j<2*number_of_taxa; j++) tmp_presence2[i][j] = FALSE;
			}
		
		overall_presence = malloc(num_gene_nodes*sizeof(int *));
		for(i=0; i<num_gene_nodes; i++) 
			{
			overall_presence[i] = malloc((2*number_of_taxa)*sizeof(int *));
			for(j=0; j<2*number_of_taxa; j++) overall_presence[i][j] = FALSE;
			}

		if(presence_of_trichotomies(gene_tree)) gene_tree = do_resolve_tricotomies(gene_tree, species_tree, basescore);

		parts[0] = gene_tree;
		gene_tree = NULL;
		strcpy(temptree, "");
		print_tree(parts[0], temptree);
		
		duplicate_tree(parts[0], NULL);
		copy = tree_top;			
		i = number_tree(copy, 0);
		best_total = -1;
		for(j=0; j<i; j++) /* for every node (internal and external) on the tree */
			{
			position = get_branch(copy, j);
			tree_top = copy;
			reroot_tree(position);
		
			copy = tree_top;
			tree_top = NULL;
			printf2("e\n");
			total = tree_map(copy, species_tree,0);
			if(total < best_total || best_total == -1) best_total = total;
			dismantle_tree(copy);
			copy = NULL;
			temp_top = NULL;
			tree_top = NULL;
			duplicate_tree(parts[0], NULL);
			copy = tree_top;
			temp_top = NULL;
			tree_top = NULL;
			number_tree(copy, 0);
			}  /* WE now know the cost for the best reconstruction of dups and losses on this tree WITHOUT any HGTs */
		dismantle_tree(copy);
		copy = NULL;
		original = best_total;
		overall_placements[0] = original;
		numparts = 1;
		for(k=0; k<numparts; k++)
			{
			printf2("ON part %d out of %d parts\n", k, numparts);
			tree_top = parts[k];
			compress_tree1(parts[k]);
			parts[k] = tree_top;
			best_reconstruction = overall_placements[k];
			found_better = FALSE;
			j=0;
			while(parts[j] != NULL) j++;
			place_marker = j;   /* THis is to tell us the next available spce in parts[], so that if there is a HGT found we know where to attach it */
			reset_tag2(parts[k]);
			/* calculate the cost for duplications and losses for the tree as it is given  (trying all rootings )*/
			if(overall_placements[k] != 0 && overall_placements[k] > 1*hgt_weight)  /* If the best reconstruction of the unbroken tree requires no dups or losses, or if the cost of 1 HGT is greater than this score, we don't test for HGTs */
				{
				i = number_tree(parts[k], 0);
				for(j=0; j<i; j++) /* go through all the parts of the tree and for each part break it into two at every point and look for new hgts */
					{
					partA = partB = TRUE;
					/* Duplicate this tree (or part of), that we are at, and use this to find the best place to break the tree ***/

					tree_top = NULL;
					duplicate_tree(parts[k], NULL);
					copy = tree_top;
					
					number_tree(copy, 0);
					position = get_branch(copy, j); /* this is the branch at which we are going to break the tree :-) */
					if(position != NULL) /* no point in seperating at the top node of the tree (which in reality doesn't exist ) */
						{
						
						test_part = position;
						if(position == copy) copy = position->next_sibling;
						/* mark the position where this was removed from */
						posit = position;
						if(position->next_sibling != NULL) attached = position->next_sibling;
						else attached = position->prev_sibling;
						while(posit->prev_sibling != NULL) posit = posit->prev_sibling;
						while(posit != NULL)
							{
							posit->tag2 = TRUE;
							posit = posit->next_sibling;
							}
						
						/* remove this node from the tree */
						if(position->parent != NULL)
							{
							(position->parent)->daughter = position->next_sibling;
							(position->next_sibling)->parent = position->parent;
							position->parent = NULL;
							}
						if(position->next_sibling != NULL) (position->next_sibling)->prev_sibling = position->prev_sibling;
						if(position->prev_sibling != NULL) (position->prev_sibling)->next_sibling = position->next_sibling;
						position->next_sibling = NULL;
						position->prev_sibling = NULL;
						/* now we need to compress the tree we just pruned to make sure that there are no unnecessary pointer taxa */
						tree_top = copy;
						compress_tree1(copy);
						copy = tree_top;
						tree_top = test_part;
						compress_tree1(test_part);
						test_part = tree_top;
						
						
						/* now reconstruct the duplications and losses of "copy" and or "test_part" seperately to see if this helped the reconsctruction */
						/* Either "copy" or "test part"  could be the HGT part. We must treat both as the HGT part in turn and pick which one has the best dup+loss reconstruction */
						/* Whichever is being treated as the HGT cannot be rerooted because if this is the HGT part then where we broke the tree is when the HGT joined the new lineage */
						/* The other half is treated as the "donor" and can be rerooted to find the best rooting for dup+loss reconstruction */
						for(q=0; q<2*number_of_taxa; q++)
							{
							presence1[q] = FALSE;
							presence2[q] = FALSE;
							}
						
						/* first calculate the dup+loss score for the two parts, this is the score for each as if both were the HGT part, "test_part" doesn't have to be rerooted, but "copy" does, at the point that "test_part" was attached */
						if(k == 0)  /* we only need to assume "copy" was the HGT if k == 0, otherwise it must be the donor */
							{
							tree_top = NULL;
							duplicate_tree(copy, NULL);
							copy1 = tree_top;
							
							tree_top = copy;
							if(copy != attached) reroot_tree(attached);
							else
								{
								posit = make_taxon();
								posit->daughter = tree_top;
								tree_top->parent = posit;
								tree_top = posit;
								posit = NULL;
								} 
							copy = tree_top;  
							printf2("f\n");
							HGT1 = tree_map(copy, species_tree,0); /* HGT1 no has the score for the tree "copy" */
							hgt_receipient1 = copy->tag;  /* if this is the HGT part, then this is the hypothesised node on the species tree that received the HGT */

							dismantle_tree(copy);
							copy = copy1;
							copy1 = NULL;
							}
						
						/* If k != 0 then we don't reroot the tree before figuring out where the HGT was attached, if k != 0 then test_part cannot be the donor! (because this would mean we would have to reroot parts[k], which is a HGT already and as such is rooted :-) */
						tree_top = NULL;
						duplicate_tree(copy, NULL);
						copy1 = tree_top;
						printf2("g\n");
						best_total1 = tree_map(copy1, species_tree,0);
						find_tagged(copy1, presence1);
						best_mapping1 = copy1; /* this will be used later for checking HGT compatibilities (( this version is only used if k != 0 ))*/
						copy1 = NULL;
						
						tree_top = NULL;
						duplicate_tree(test_part, NULL);
						copy1 = tree_top;
						printf2("h\n");
						HGT2 = tree_map(copy1, species_tree,0); /* HGT2 no has the score for the tree "test_part" */
						hgt_receipient2 = copy1->tag; /* if this is the HGT part, then this is the hypothesised node on the species tree that received the HGT */
						/*find_tagged(copy1, presence2); This is commented out because if k != 0 then test_part cannot be the donor */
						dismantle_tree(copy1);
						copy1 = NULL;
						tree_top = NULL;
						/*** Next calculate the best rerooting of each tree in terms of dups+losses ***/
						
						if(k == 0) /* we can only reroot the original part of the tree, any subsequent ones are rooted because they are the HGTs */
							{
							dismantle_tree(best_mapping1);
							best_mapping1 = NULL;
							for(z=0; z<2; z++)
								{
								/**** CHECK TO SEE IF BREAKING THE TREE AT THIS POINT MAKES THE HGT "TRAVEL IN TIME" also calculate the best dups+loss reconstruction for each tree  */
								tree_top = NULL;
								if(z == 0) duplicate_tree(copy, NULL);
								else duplicate_tree(test_part, NULL);
								copy1 = tree_top;
								tree_top = NULL;
								tree_top = copy1;
								compress_tree1(copy1);
								copy1 = tree_top;
								
								x = number_tree(copy1, 0);
								best_total = -1;
								if(x >3)
									{
									for(y=0; y<x; y++)
										{
										position = get_branch(copy1, y);
										tree_top = copy1;
										if(copy1 != position) reroot_tree(position);
										else
											{
											if(z == 0)
												{
												posit = tree_top;
												q=0;
												while(posit != NULL)
													{
													q++;
													posit = posit->next_sibling;
													}
												posit = NULL;
												if(q != 1)
													{
													reroot_tree(tree_top);
													}
												}
											} 
										copy1 = tree_top;
										tree_top = NULL;
										
										printf2("I\n");
										total = tree_map(copy1, species_tree,0);
										if(total < best_total || best_total == -1)
											{
											best_total = total;
											if(best_mapping != NULL)
												{
												dismantle_tree(best_mapping);
												best_mapping = NULL;
												}
											best_mapping = copy1;
											copy1 = NULL;
											}
										else
											{
											dismantle_tree(copy1);
											copy1 = NULL;
											}
										temp_top = NULL;
										tree_top = NULL;
										if(z == 0) duplicate_tree(copy, NULL);
										else duplicate_tree(test_part, NULL);
										copy1 = tree_top;
										temp_top = NULL;
										tree_top = copy1;
										compress_tree1(copy1);
										copy1 = tree_top;
										tree_top = NULL;
										number_tree(copy1, 0);
										}
									}
								else
									{
									printf2("j\n");
									best_total = tree_map(copy1, species_tree,0);
									best_mapping = copy1;
									copy1 = NULL;
									}
								if(z==0) best_total1 = best_total;
								else best_total2 = best_total;
								
								
								for(q=0; q<2*number_of_taxa; q++)
									{
									if(z == 0) presence1[q] = FALSE;
									else presence2[q] = FALSE;
									}

								if(z == 0) find_tagged(best_mapping, presence1);  /* so after this function "presence" will have all the possible positions that the other part of the tree (if it was a HGT) could have come from. */
								else find_tagged(best_mapping, presence2);
								
								if(z==0) best_mapping1 = best_mapping;
								else best_mapping2 = best_mapping;
								best_mapping = NULL;
								
								if(copy1 != NULL)
									{
									dismantle_tree(copy1);
									copy1 = NULL;
									}
								}
							}
						
						/**** WE NOW HAVE THE FOLLOWING VALUES:  */
						/**** HGT1 = reconstruction cost of "copy" as it if was the HGT **/
						/**** HGT2 = reconstruction cost of "test_part" as it if was the HGT **/
						/**** best_total1 = reconstruction cost of "copy" as if it was the donor part */
						/**** best_total2 = reconstruction cost of "test_part" as if it was the donor part */
						/**** presence1 = contains the possible donor nodes if "copy" was the donor */
						/**** presence2 = contains the possible donor nodes if "test_part" was the donor */
						/**** hgt_receipient1 = is the node number of the receipient if "copy" wat the HGT part */
						/**** hgt_receipient2 = is the node number of the receipient if "test_part" wat the HGT part */
						/****/
						/**** With these we need to now test 1) if hypothesising if either was a HGT causes them to travel beack in time **/
						/**** 2) Whichever passes the above test must then have a reconstruction cost that is better than the original tree without ANY HGT hypothesised */
						/**** If both pass the first test and both have the same reconstruction cost which is better than the original tree, we arbitrarily choose one over the other (not sure how else to deal with ths :-) */
						
						/* First see if the placements of either contradicts the species tree array we made earlier */ 
						/* If the donor could have been the same as the receipient, then we don't allow the HGT */
						allow_HGT1 = allow_HGT2 = TRUE;
						if(presence1[hgt_receipient2] == TRUE) allow_HGT1 = FALSE;
						if(presence2[hgt_receipient1] == TRUE) allow_HGT2 = FALSE;  /* If we are able to hypothesise that the donor and the receipient was the same then we disregard the HGT */	
						
						if(k!= 0) allow_HGT2 = FALSE;	/* WE cannot reroot this if it isnot the original part of the tree (ie k == 0) and we already know that this is a HGT part, so it can only be a donor for a hgt event, so copy HAS TO BE the donor */
																												
						/* Test to see if any of the hypothesised HGTs would have conflicted with the structure of the Species tree */
						if(allow_HGT1 == TRUE || allow_HGT2 == TRUE)
							{
							s = FALSE;
							r = FALSE;
							for(q=0; q<(2*number_of_taxa); q++)
								{
								if(presence1[q] == TRUE && species_allowed[q][hgt_receipient2] == TRUE) s = TRUE;
								if(presence2[q] == TRUE && species_allowed[q][hgt_receipient1] == TRUE) r = TRUE;
								}
							if(s == FALSE) allow_HGT1 = FALSE;
							if(r == FALSE) allow_HGT2 = FALSE;
							}
						/*** NOW SEE IF THIS CONFLICTS WITH ANY PREVIOUS HGTS!!!!!!!****/
						
						if(numparts > 1) /* no point checking this if it is the first hgt defined */
							{
							for(q=0; q<num_gene_nodes; q++)
								{
								for(r=0; r<2*number_of_taxa; r++)
									{
									tmp_presence1[q][r] = FALSE;
									tmp_presence2[q][r] = FALSE;
									}
								}
							
							/* See if either of these two possible HGTs conflict with previous HGTs */
							
							if(allow_HGT1)
								{
								/* check both copy and test_part for any previous HGTs donor sites and then check what the possible donors are now, given the new HGT found */
								/* HGT1 assumes that copy is the donor and test_part is the HGT */
								/* we will use the best rooting of copy to calculate the possible donor positions of any hgts defined earlier, these will be saved in tmp_presence1 */
								/* best_mapping1  holds the optimal mapping of Copy as the donor */
								for(q=1; q<numparts; q++)
									{
									reset_tag2(best_mapping1);
									if(assign_tag2(best_mapping1, q))
										{
										find_tagged(best_mapping1, tmp_presence1[q]);
										}
									else
										{
										for(r=0; r<2*number_of_taxa; r++)
											{
											tmp_presence1[q][r] = overall_presence[q][r];
											}
										}
									}
								
								for(s=0; s<(2*number_of_taxa); s++)  /* initialise this to include the timetravelling HGTs that occur from the shape of the tree */
									{
									for(r=0; r<(2*number_of_taxa); r++) dependent_species_allowed[s][r] = species_allowed[s][r];
									}	
								
								/* NOW check to see if the new HGT causes any conflicts with the previous HGTs */
								
								for(q=1; q<numparts; q++)
									{
									
									 /* check to see if this HGT conflicts with any previously assigned ((not needed for the first hgt)) */
									
									if(allow_HGT1)
										{
										/* identify the nodes above and below the donor and receptor sites so cross hgts can be defined as being not allowed */
										for(r=0; r<(2*number_of_taxa); r++) before1[r]  = after1[r] = before2[r] = after2[r]= FALSE;
										for(r=0; r<2*number_of_taxa; r++) previous[r] = FALSE;
										for(r=0; r<(2*number_of_taxa); r++)
											{
											if(tmp_presence1[q][r] == TRUE)
												assign_before_after(species_tree, previous, before1, after1, r, FALSE);
											}
										for(r=0; r<(2*number_of_taxa); r++)
											{
											if(tmp_presence1[q][r] == TRUE) before1[r] = after1[r] = FALSE;
											}
										for(r=0; r<2*number_of_taxa; r++) previous[r] = FALSE;
										assign_before_after(species_tree, previous, before2, after2, overall_receptor[q], FALSE);
										before2[overall_receptor[q]] = after2[overall_receptor[q]] = FALSE;
										
										/** assign the new cross-hgt rules to the array dependent_species_allowed */
										for(r=0; r<2*number_of_taxa; r++)
											{
											if(before1[r] == TRUE)
												{
												for(s=0; s<2*number_of_taxa; s++)
													{
													if(after2[s] == TRUE)
														{
														dependent_species_allowed[r][s] = FALSE;
														dependent_species_allowed[s][r] = FALSE;
														}
													}
												}
											}
										for(r=0; r<2*number_of_taxa; r++)
											{
											if(before2[r] == TRUE)
												{
												for(s=0; s<2*number_of_taxa; s++)
													{
													if(after1[s] == TRUE)
														{
														dependent_species_allowed[r][s] = FALSE;
														dependent_species_allowed[s][r] = FALSE;
														}
													}
												}
											}
										}
									
									/* check this new HGT against each of the older hgts as they are added  */
									tmp_allow = FALSE;
									
									for(r=0; r<(2*number_of_taxa); r++)
										{
										if(presence1[r] == TRUE && dependent_species_allowed[r][hgt_receipient2] == TRUE)
											tmp_allow = TRUE;
										if(presence1[r] == TRUE && dependent_species_allowed[r][hgt_receipient2] == FALSE)
											presence1[r] = FALSE;
										}
									if(tmp_allow == FALSE)
										{
										allow_HGT1 = FALSE;
										q = (2*number_of_taxa);
										}
									/* this previous section checks the  newhgt against the rules put in place by previous HGTs and asks if this new HGT violates these rules */
									/* We need one other check. There is a chance that this new HGT may not conflict with any of the rules put in place by previous HGTs, but by its addition it causes some previous HGTs to violate the rules */
									/* this is what needs to be checked next */
									
									if(q>1)
										{
										/* dependent_species_allowed has all the rules imposed by each of the HGTs examined so far */
										
										for(s=1; s<q; s++)
											{
											tmp_allow = FALSE;
											for(r=0; r<(2*number_of_taxa); r++)
												{
												if(tmp_presence1[s][r] == TRUE)
													{
													if(dependent_species_allowed[r][overall_receptor[s]] == TRUE) tmp_allow = TRUE;
													}
												}
											if(tmp_allow == FALSE)
												{
												allow_HGT1 = FALSE;
												s = q;
												}
											}
										}
									
									} 
										
								}
							if(allow_HGT2)
								{
								/* check both copy and test_part for any previous HGTs donor sites and then check what the possible donors are now, given the new HGT found */
								/* HGT2 assumes that test_part is the donor and copy is the HGT */
								/* we will use the best rooting of copy to calculate the possible donor positions of any hgts defined earlier, these will be saved in tmp_presence2 */
								/* best_mapping2  holds the optimal mapping of Copy as the donor */
								for(q=1; q<numparts; q++)
									{
									reset_tag2(best_mapping2);
									if(assign_tag2(best_mapping2, q))
										{
										find_tagged(best_mapping2, tmp_presence2[q]);
										}
									else
										{
										for(r=0; r<2*number_of_taxa; r++)
											{
											tmp_presence2[q][r] = overall_presence[q][r];
											}
										}
									}
								
								for(s=0; s<(2*number_of_taxa); s++)  /* initialise this to include the timetravelling HGTs that occur from the shape of the tree */
									{
									for(r=0; r<(2*number_of_taxa); r++) dependent_species_allowed[s][r] = species_allowed[s][r];
									}	
								
								/* NOW check to see if the new HGT causes any conflicts with the previous HGTs */
								
								for(q=1; q<numparts; q++)
									{
									
									 /* check to see if this HGT conflicts with any previously assigned ((not needed for the first hgt)) */
									
									if(allow_HGT2)
										{
										/* identify the nodes above and below the donor and receptor sites so cross hgts can be defined as being not allowed */
										for(r=0; r<(2*number_of_taxa); r++) before1[r]  = after1[r] = before2[r] = after2[r]= FALSE;
										for(r=0; r<2*number_of_taxa; r++) previous[r] = FALSE;
										for(r=0; r<(2*number_of_taxa); r++)
											{
											if(tmp_presence2[q][r] == TRUE)
												assign_before_after(species_tree, previous, before1, after1, r, FALSE);
											}
										for(r=0; r<(2*number_of_taxa); r++)
											{
											if(tmp_presence2[q][r] == TRUE) before1[r] = after1[r] = FALSE;
											}
										for(r=0; r<2*number_of_taxa; r++) previous[r] = FALSE;
										assign_before_after(species_tree, previous, before2, after2, overall_receptor[q], FALSE);
										before2[overall_receptor[q]] = after2[overall_receptor[q]] = FALSE;
										
										/** assign the new cross-hgt rules to the array dependent_species_allowed */
										for(r=0; r<2*number_of_taxa; r++)
											{
											if(before1[r] == TRUE)
												{
												for(s=0; s<2*number_of_taxa; s++)
													{
													if(after2[s] == TRUE)
														{
														dependent_species_allowed[r][s] = FALSE;
														dependent_species_allowed[s][r] = FALSE;
														}
													}
												}
											}
										for(r=0; r<2*number_of_taxa; r++)
											{
											if(before2[r] == TRUE)
												{
												for(s=0; s<2*number_of_taxa; s++)
													{
													if(after1[s] == TRUE)
														{
														dependent_species_allowed[r][s] = FALSE;
														dependent_species_allowed[s][r] = FALSE;
														}
													}
												}
											}
										}
									
									/* check this new HGT against each of the older hgts as they are added  */
									tmp_allow = FALSE;
									
									for(r=0; r<(2*number_of_taxa); r++)
										{
										if(presence2[r] == TRUE && dependent_species_allowed[r][hgt_receipient1] == TRUE)
											tmp_allow = TRUE;
										if(presence2[r] == TRUE && dependent_species_allowed[r][hgt_receipient1] == FALSE)
											presence2[r] = FALSE;
										}
									if(tmp_allow == FALSE)
										{
										allow_HGT2 = FALSE;
										q = (2*number_of_taxa);
										}
									/* this previous section checks the  newhgt against the rules put in place by previous HGTs and asks if this new HGT violates these rules */
									/* We need one other check. There is a chance that this new HGT may not conflict with any of the rules put in place by previous HGTs, but by its addition it causes some previous HGTs to violate the rules */
									/* this is what needs to be checked next */
									
									if(q>1)
										{
										/* dependent_species_allowed has all the rules imposed by each of the HGTs examined so far */
										
										for(s=1; s<q; s++)
											{
											tmp_allow = FALSE;
											for(r=0; r<(2*number_of_taxa); r++)
												{
												if(tmp_presence2[s][r] == TRUE)
													{
													if(dependent_species_allowed[r][overall_receptor[s]] == TRUE) tmp_allow = TRUE;
													}
												}
											if(tmp_allow == FALSE)
												{
												allow_HGT2 = FALSE;
												s = q;
												}
											}
										}
									
									} 
										
								}							
							
							}
						
						if(best_mapping1 != NULL) dismantle_tree(best_mapping1);
							best_mapping1 = NULL;
						if(best_mapping2 != NULL) dismantle_tree(best_mapping2);
							best_mapping2 = NULL;
						
						/***************************************************************/
						/* now see if the reconstruction cost of either of the HGTs that are still allowed would be less than (or equal to?) the best so far (which may be the original) */
						if(allow_HGT1 && allow_HGT2)
							{
							
							if(best_total1 + HGT2 + (1*hgt_weight) < best_total2 + HGT1 + (1*hgt_weight) && best_total1 + HGT2 + (1*hgt_weight) < best_reconstruction)
								{
								found_better = TRUE;
								best_reconstruction = best_total1 + HGT2 + (1*hgt_weight);
								donor_score = best_total1;
								receptor = hgt_receipient2;
								HGT_score = HGT2;
								if(best_donor != NULL) dismantle_tree(best_donor);
								if(best_HGT != NULL) dismantle_tree(best_HGT);
								best_donor = copy;
								copy = NULL;
								best_HGT = test_part;
								test_part = NULL;
								for(q=0; q<2*number_of_taxa; q++) presence[q] = presence1[q];
								for(q=1; q<numparts; q++)
									{
									for(r=0; r<2*number_of_taxa; r++) overall_presence[q][r] = tmp_presence1[q][r];
									}
								}
							else
								{
								if(best_total2 + HGT1 + (1*hgt_weight) <= best_total1 + HGT2 + (1*hgt_weight) && best_total2 + HGT1 + (1*hgt_weight) < best_reconstruction)
									{
									found_better = TRUE;
									best_reconstruction = best_total2 + HGT1 + (1*hgt_weight);
									donor_score = best_total2;
									receptor = hgt_receipient1;
									HGT_score = HGT1;
									if(best_donor != NULL) dismantle_tree(best_donor);
									if(best_HGT != NULL) dismantle_tree(best_HGT);
									best_donor = test_part;
									test_part = NULL;
									best_HGT = copy;
									copy = NULL;
									for(q=0; q<2*number_of_taxa; q++) presence[q] = presence2[q];
									for(q=1; q<numparts; q++)
										{
										for(r=0; r<2*number_of_taxa; r++) overall_presence[q][r] = tmp_presence2[q][r];
										}
									}
								
								}
							}
						else
							{
							if(allow_HGT1) /* this assumes that "copy" is the donor and "test_part" is the HGT */
								{
								if(best_total1 + HGT2 + (1*hgt_weight) < best_reconstruction)
									{ /*do something */
									found_better = TRUE;
									best_reconstruction = best_total1 + HGT2 + (1*hgt_weight) ;
									donor_score = best_total1;
									receptor = hgt_receipient2;
									HGT_score = HGT2;
									if(best_donor != NULL) dismantle_tree(best_donor);
									if(best_HGT != NULL) dismantle_tree(best_HGT);
									best_donor = copy;
									copy = NULL;
									best_HGT = test_part;
									test_part = NULL;
									for(q=0; q<2*number_of_taxa; q++) presence[q] = presence1[q];
									for(q=1; q<numparts; q++)
										{
										for(r=0; r<2*number_of_taxa; r++) overall_presence[q][r] = tmp_presence1[q][r];
										}
									}
								}
							if(allow_HGT2) /* this assumes that "test_part" is the donor and "copy" is the HGT */
								{
								if(best_total2 + HGT1 + (1*hgt_weight) < best_reconstruction)
									{ /*do something */
									found_better = TRUE;
									best_reconstruction = best_total2 + HGT1 + (1*hgt_weight);
									donor_score = best_total2;
									receptor = hgt_receipient1;
									HGT_score = HGT1;
									if(best_donor != NULL) dismantle_tree(best_donor);
									if(best_HGT != NULL) dismantle_tree(best_HGT);
									best_donor = test_part;
									test_part = NULL;
									best_HGT = copy;
									copy = NULL;
									for(q=0; q<2*number_of_taxa; q++) presence[q] = presence2[q];
									for(q=1; q<numparts; q++)
										{
										for(r=0; r<2*number_of_taxa; r++) overall_presence[q][r] = tmp_presence2[q][r];
										}
									}
								}
							}
						}
					if(copy != NULL) dismantle_tree(copy);
					copy = NULL;
					if(test_part != NULL) dismantle_tree(test_part);
					test_part = NULL;
					}
				}
			if(found_better == TRUE) /* if we have found a better reconstruction that includes a HGT for this part of the tree */
				{
				/* dismantle the tree at parts[k] */
				printf2("HGT%d\n", numparts);
				overall_receptor[place_marker] = receptor;
				overall_placements[k] = donor_score;
				overall_placements[place_marker] = HGT_score;
				if(parts[k] != NULL) dismantle_tree(parts[k]);
				if(parts[place_marker] != NULL) dismantle_tree(parts[place_marker]);
				parts[k] = best_donor;
				best_donor = NULL;
				strcpy(temptree, "");
				print_tree(parts[k], temptree);
				printf2("donor:\n%s\n", temptree);
				parts[place_marker] = best_HGT;
				strcpy(temptree, "");
				print_tree(parts[place_marker], temptree);
				printf2("hgt:\n%s\n", temptree);

				best_HGT = NULL;
				assign_hgtdonors(parts[k], num_gene_nodes, place_marker);
				for(q=0; q<2*number_of_taxa; q++) overall_presence[place_marker][q] = presence[q]; /* record the possible donor nodes for this HGT */
				
				numparts++;
				k--; /* we need to re-evaluate the donor again to see if removing this HGT means that another HGT becomes apparant */
				}
			
			}
		printf2("\n\nprinting results\n");
			
		/**** print out the results */
		i=0;
		while(parts[i] != NULL && i < num_gene_nodes)
			{
			printf2("part %d\nScore = %f\nReceptor Node:%d\nPossible Donors: ",i,overall_placements[i], overall_receptor[i] );
			for(j=0; j<2*number_of_taxa; j++) printf2("%d,", overall_presence[i][j]);
			printf2("\n");
			strcpy(temptree, "");
			strcpy(temptree, "");
			print_tree(parts[i], temptree);
			printf2("parts[%d] %s\n",i, temptree);
			printf2("\n");
			/*print_tree(parts[i], temptree);
			tree_coordinates(temptree, TRUE, FALSE, FALSE);
			*/
			i++;
			}
		
			
		if(copy != NULL)
			{
			dismantle_tree(copy);
			copy = NULL;
			}
		if(gene_tree != NULL)
			{
			dismantle_tree(gene_tree);
			gene_tree = NULL;
			}
		if(best_mapping != NULL)
			{
			dismantle_tree(best_mapping);
			best_mapping = NULL;
			}
		for(i=0; i<num_gene_nodes; i++)
			{
			if(parts[i] != NULL) dismantle_tree(parts[i]);
			parts[i] = NULL;
			}
		
		free(parts);
		parts = NULL;
		free(overall_receptor);
		overall_receptor = NULL;
		for(i=0; i<num_gene_nodes; i++)
			{
			free(tmp_presence1[i]);
			tmp_presence1[i] = NULL;
			free(tmp_presence2[i]);
			tmp_presence2[i] = NULL;
			free(overall_presence[i]);
			overall_presence[i] = NULL;
			}
		
		free(overall_presence);
		overall_presence = NULL;
		free(tmp_presence1);
		tmp_presence1 = NULL;
		free(tmp_presence2);
		tmp_presence2 = NULL;
		
		free(overall_placements);
		overall_placements = NULL;
		}
	
	for(i=0; i<(2*number_of_taxa); i++) free(species_allowed[i]);
	free(species_allowed);
	species_allowed = NULL;
	for(i=0; i<(2*number_of_taxa); i++) free(dependent_species_allowed[i]);
	free(dependent_species_allowed);
	dependent_species_allowed = NULL;
	
	free(before1); 
	free(before2); 
	free(after1); 
	free(after2);

	
	if(species_tree != NULL)
		{
		dismantle_tree(species_tree);
		species_tree = NULL;
		}
	
	if(gene_tree != NULL)
		{
		dismantle_tree(gene_tree);
		gene_tree = NULL;
		}
	
	tree_top = NULL;
	free(temptree);
	free(temptree1);
	free(presence);
	free(presence1);
	free(presence2);
	free(temporary);



	}

void assign_before_after(struct taxon *position, int *previous, int *before, int *after, int num, int found)
	{
	int i=0;
	
	while(position != NULL)
		{
		if(position->tag == num)
			{
			for(i=0; i<(2*number_of_taxa); i++)
				before[i] = previous[i];  
			found = TRUE;
			}
		else
			{
			if(found) after[position->tag] = TRUE;
			}
			
		if(position->daughter != NULL)
			{
			previous[position->tag] = TRUE;
			assign_before_after(position->daughter,previous, before, after, num, found);
			previous[position->tag] = FALSE;
			}
		if(position->tag == num) found = FALSE;
		position = position->next_sibling;
		}
	}

void assign_ances_desc(struct taxon *position, int ** allowed_species, int * previous)
	{
	int i=0;
	
	while(position != NULL)
		{
		if(position->daughter != NULL)
			{
			previous[position->tag] = TRUE;
			assign_ances_desc(position->daughter, allowed_species, previous);
			previous[position->tag] = FALSE;
			}
		allowed_species[position->tag][position->tag] = FALSE;
		for(i=0; i<(2*number_of_taxa); i++)
			{
			if(previous[i])
				{
				allowed_species[i][position->tag] = FALSE;
				allowed_species[position->tag][i] = FALSE;
				}
			}
		position = position->next_sibling;
		}
	}

void resolve_tricotomies_dist (struct taxon *gene_tree, struct taxon *species_tree, int ** scores)
	{
	struct taxon *position = gene_tree, *start = gene_tree, *position2 = NULL, *new = NULL, *first = NULL, *second = NULL;
	int i=0, j, minscore, *presence = NULL;
	
	presence = malloc(number_of_taxa*sizeof(int));
	for(i=0; i<number_of_taxa; i++) presence[i] = FALSE;
	i=0;
	while(position != NULL)
		{
		i++;
		/* go though looking for internal branches and follow them down */
		if(position->daughter != NULL)
			resolve_tricotomies_dist(position->daughter, species_tree, scores);
		position = position->next_sibling;
		}
	
	/* for trichotomies at this level, idenfiy the internal branch on the species tree that best matches */
	while((temp_top2 != start && i > 2) || (temp_top2 == start && i > 3))
		{
		position = start;
		minscore=-1; first = NULL; second = NULL; 
		/* identify the branches at this level with the minimum distance between them */
		while(position->next_sibling != NULL)
			{
			position2 = position->next_sibling;
			while(position2 != NULL)
				{
				if(minscore==-1 || scores[position->tag][position2->tag] < minscore)
					{
					minscore = scores[position->tag][position2->tag];
					first = position;
					second = position2;
					}
				else
					{
					if((first->tag >= number_of_taxa || second->tag >= number_of_taxa) && (position->tag < number_of_taxa && position2->tag < number_of_taxa) && (scores[position->tag][position2->tag] == minscore))
						{
						minscore = scores[position->tag][position2->tag];
						first = position;
						second = position2;						
						}
					}
				position2=position2->next_sibling;
				}
			position = position->next_sibling;
			}
		position= first;
		position2= second;
		/* having identified the parts with the minimum distance, put them together */
		/* create a new internal branch */
		new = make_taxon();
		/* make the two min branches daughters of the new node */
		if(position->next_sibling != NULL) (position->next_sibling)->prev_sibling = position->prev_sibling;
		if(position->prev_sibling != NULL) (position->prev_sibling)->next_sibling = position->next_sibling;
		if(position->prev_sibling == NULL)
			{
			start = position->next_sibling;
			if(position->parent == NULL) 
				temp_top2 = start;
			else 
				(position->parent)->daughter = position->next_sibling;
			(position->next_sibling)->parent = position->parent;
			
			}
		position->parent = new;
		new->daughter = position;
		position->next_sibling = position2;
		position->prev_sibling = NULL;
		
		if(position2->next_sibling != NULL) (position2->next_sibling)->prev_sibling = position2->prev_sibling;
		if(position2->prev_sibling != NULL) (position2->prev_sibling)->next_sibling = position2->next_sibling;
		if(position2->prev_sibling == NULL)
			{
			start = position2->next_sibling;
			if(position2->parent == NULL)
				temp_top2 = start;
			else
				{
				(position2->parent)->daughter = position2->next_sibling;
				(position2->next_sibling)->parent = position2->parent;
				position2->parent = NULL;
				}
			
			
			}
		position2->prev_sibling = position;
		position2->next_sibling = NULL;
		/* put "new" into the tree at the present position (after the first) */
		if(start->next_sibling != NULL) (start->next_sibling)->prev_sibling = new;
		new->prev_sibling = start;
		new->next_sibling = start->next_sibling;
		start->next_sibling = new;
		
		/* Find the label for the new internal brnach from the species tree */
		if(position->tag == position2->tag)
			new->tag = position->tag;
		else
			{
			for(j=0; j<number_of_taxa; j++) presence[j] = FALSE;
			get_taxa(new->daughter, presence);
			new->tag = get_best_node(species_tree, presence, -1);
			
			}
		i--;
		}
	free(presence);
	}

