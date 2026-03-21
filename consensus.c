/*
 *  consensus.c — Clann v5.0.0
 *  Consensus tree construction methods
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
#include "consensus.h"

int average_consensus(int nrep, int missing_method, char * useroutfile, FILE *paupfile)
	{
	int **taxa_comp = NULL, i, j, k, l, found = FALSE, here = FALSE, error = FALSE;
	float used_weights = 0;
	char *temptree = NULL;
	
	
	temptree = malloc(TREE_LENGTH*sizeof(char));
	if(!temptree) printf2("out of memory'n");
	temptree[0] = '\0';
	taxa_comp = malloc(number_of_taxa*sizeof(int*));
	if(!taxa_comp) printf2("out of memory'n");
	for(i=0; i<number_of_taxa; i++)
		{
		taxa_comp[i] = malloc(number_of_taxa*sizeof(int));
		if(!taxa_comp[i]) printf2("out of memory'n");
		for(j=0; j<number_of_taxa; j++)
			taxa_comp[i][j] = FALSE;
		}
	for(i=0; i<Total_fund_trees; i++)
		{
		if(sourcetreetag[i])weighted_pathmetric(fundamentals[i], weighted_scores, i);
		}
	fflush(stdout);
	for(i=0; i<number_of_taxa; i++)
		{
		for(j=0; j<number_of_taxa; j++)
			{
			if(i != j)
				{
				used_weights = 0;
				for(k=0; k<Total_fund_trees; k++)
					{
					if(presence_of_taxa[k][i] > 0 && presence_of_taxa[k][j] > 0)
						used_weights = used_weights+tree_weights[k];
					}
				weighted_scores[i][j] = weighted_scores[i][j]/used_weights;
				}
			}
		}
	/** NEED something here to check that all the cells are filled out, and if not to calculate the ultrametric or the 4-point condition for the missing cell */
	for(k=0; k<Total_fund_trees; k++)
		{
		if(sourcetreetag[k])
			{
			for(i=0; i<number_of_taxa; i++)
				{
				for(j=0; j<number_of_taxa; j++) 
					{
					if(presence_of_taxa[k][i] > 0 && presence_of_taxa[k][j] > 0) taxa_comp[i][j] = TRUE; /* this records which taxa are never present together in the same tree */
					}
				}
			}
		}
	here = FALSE;
	for(i=0; i<number_of_taxa; i++)
		{
		for(j=0; j<number_of_taxa; j++)
			{
			if(!taxa_comp[i][j]) here = TRUE;
			}
		}
	found = TRUE;
	while(here && found)
		{
	
		for(i=0; i<number_of_taxa; i++)
			{
			for(j=0; j<number_of_taxa; j++)
				{
				if(!taxa_comp[i][j])
					{
					/** Calculate an approximation for the missing cell **/
					if(missing_method == 0) /*** Do the ultrametric calculation ***/
						{
						/** First find a cell shard by i and j that has a distace for to both **/
						found = FALSE;
						for(k=0; k<number_of_taxa; k++)
							{
							if(k != i && k != j && taxa_comp[i][k] && taxa_comp[j][k] && !found)
								{
								found = TRUE;
								if(weighted_scores[i][k] > weighted_scores[j][k])
									weighted_scores[i][j] = weighted_scores[j][i] = weighted_scores[i][k];
								else
									weighted_scores[i][j] = weighted_scores[j][i] = weighted_scores[j][k];
								taxa_comp[i][j] = taxa_comp[j][i] = TRUE;
								i = j = number_of_taxa;
								k = l = number_of_taxa;
								}
							}
						}
					if(missing_method == 1) /*** Do the 4 point condition calculation ***/
						{
						
						/** this method requires that i find two cells that have values for both i and j and for each other **/
						found = FALSE;
						for(k=0; k<number_of_taxa; k++)
							{
							for(l=0; l<number_of_taxa; l++)
								{
								if(k != i && k != j && l != i && l != j && l != k && taxa_comp[i][k] && taxa_comp[i][l] && taxa_comp[j][k] && taxa_comp[j][l] && taxa_comp[k][l] && !found)
									{
									found = TRUE;
									if((weighted_scores[i][l] + weighted_scores[j][k]) > (weighted_scores[i][k] + weighted_scores[j][l]))
										weighted_scores[i][j] = weighted_scores[j][i] = (weighted_scores[i][l] + weighted_scores[j][k]) - weighted_scores[k][l];
									else
										weighted_scores[i][j] = weighted_scores[j][i] = (weighted_scores[i][k] + weighted_scores[j][l]) - weighted_scores[k][l];
									
									taxa_comp[i][j] = taxa_comp[j][i] = TRUE;
									i = j = number_of_taxa;
									k = l = number_of_taxa;
									}
								}
							}
						}
					}
				}
			}
			
		here = FALSE;
		for(i=0; i<number_of_taxa; i++)
			{
			for(j=0; j<number_of_taxa; j++)
				{
				if(!taxa_comp[i][j]) here = TRUE;
				}
			}
		
		}
	
	if(!found)
		{
		printf2("\n\nERROR: the overlap in the data is too sparse to calculate missing cells using ");
		if(missing_method == 0) printf2("an ultrametric estimate\n");
		if(missing_method == 1) printf2("a 4 point condition estimate\n");
		error = TRUE;
		}
	else
		{
		if(paupfile != NULL)
			{
			if(nrep == 0 || nrep == 1) fprintf(paupfile,"#nexus\n");
			fprintf(paupfile, "\n\nbegin distances;\ndimensions ntax = %d;\nformat nodiagonal;\nmatrix\n", number_of_taxa);
			for(i=0; i<number_of_taxa; i++)
				{
				fprintf(paupfile, "\n\'%-20.20s\' ", taxa_names[i]);
				for(j=0; j<i; j++)
					{
					fprintf(paupfile, "%f ", weighted_scores[i][j]);
					}
				}
			fprintf(paupfile, "\n;\nend;\n");
			fprintf(paupfile, "begin paup;\nset increase=auto notifybeep=no errorbeep=no;\nset criterion = dist;\ndset objective = lsfit;\nhs;\nshowtrees;\n\nsavetrees file = %s brlens=yes format = nexus", useroutfile);
			if(nrep==0) fprintf(paupfile," Append=no replace=yes;\n\tquit;\n\nend;\n");  /* if there is one 1 rep */
			if(nrep==1) fprintf(paupfile," Append=no replace=yes;\n\nend;\n");     	/* if this is the forst of a few reps */
			if(nrep==2) fprintf(paupfile," Append=yes;\n\nend;\n");			/* if this is a middle rep */
			if(nrep==3) fprintf(paupfile," Append=yes;\n\tquit;\n\nend;\n");		/* if this is the last of a few reps */
			}
		
		}
	for(i=0; i<number_of_taxa; i++) free(taxa_comp[i]);
	free(taxa_comp);
	free(temptree);
	return(error);
	}

int coding(int nrep, int search, int ptpreps)
    {
    int i=0, j=0, k=0, nreps=10,  parenthesis = 0, count =0, *tracking = NULL, total = 0, **BR_coding = NULL, nodecount = 0, position = 0, split_count = 0, calculate_inhouse = FALSE;
    char number[100];
    char *string = malloc(TREE_LENGTH * sizeof(char));
    char filename[1000], **temptrees = NULL;
    int x=0, one_in_this = FALSE, zero_in_this = FALSE, swap=3, addseq=4, error=FALSE, *num_fund_taxa = NULL, njbuild = FALSE, weighted = FALSE;

    if(!string) { printf2("Error: out of memory in coding\n"); return FALSE; }
    filename[0] = '\0';
    for(i=0; i<num_commands; i++)
        {       
        if(strcmp(parsed_command[i], "nreps") == 0 && (strcmp(parsed_command[0], "boot") != 0 && strcmp(parsed_command[0], "bootstrap") != 0))
            {
            nreps = toint(parsed_command[i+1]);
            
            if(nreps == 0)
                {
                printf2("Error: '%s' is an invalid integer number to be assigned to nreps\n", parsed_command[i+1]);
                error = TRUE;
                }
            }

        if(strcmp(parsed_command[i], "hsreps") == 0 )
            {
            nreps = toint(parsed_command[i+1]);
            
            if(nreps == 0)
                {
                printf2("Error: '%s' is an invalid integer number to be assigned to nreps\n", parsed_command[i+1]);
                error = TRUE;
                }
            }
 
        if(strcmp(parsed_command[i], "analysis") == 0 )
            {
            if(strcmp(parsed_command[i+1], "parsimony") == 0)
                njbuild = FALSE;
            else
                {
                if(strcmp(parsed_command[i+1], "nj") == 0)
                    njbuild = TRUE;
                else
                    {
                    printf2("Error: option %s not valid analysis\n", parsed_command[i+1]);
                    error = TRUE;
                    }
                }
            }
            


        if(strcmp(parsed_command[i], "addseq") == 0)
            {
            if(strcmp(parsed_command[i+1], "simple") == 0)
                addseq = 1;
            if(strcmp(parsed_command[i+1], "closest") == 0)
                addseq = 2;
            if(strcmp(parsed_command[i+1], "asis") == 0)
                addseq = 3;
            if(strcmp(parsed_command[i+1], "random") == 0)
                addseq = 4;
            if(strcmp(parsed_command[i+1], "furthest") == 0)
                addseq = 5;
            if(strcmp(parsed_command[i+1], "simple") != 0 && strcmp(parsed_command[i+1], "closest") != 0 && strcmp(parsed_command[i+1], "asis") != 0 && strcmp(parsed_command[i+1], "random") != 0 && strcmp(parsed_command[i+1], "furthest") != 0)
                {
                printf2("Error addseq option %s not known\n", parsed_command[i+1]);
                error=TRUE;
                }
            }
        if(strcmp(parsed_command[i], "swap") == 0)
            {
            if(strcmp(parsed_command[i+1], "nni") == 0)
                swap = 1;
            if(strcmp(parsed_command[i+1], "spr") == 0)
                swap = 2;
            if(strcmp(parsed_command[i+1], "tbr") == 0)
                swap = 3;
            if(strcmp(parsed_command[i+1], "nni") != 0 && strcmp(parsed_command[i+1], "spr") != 0 && strcmp(parsed_command[i+1], "tbr") != 0)
                {
                printf2("Error swap option %s not known\n", parsed_command[i+1]);
                error=TRUE;
                }
            }
		
		if(strcmp(parsed_command[i], "weighted") == 0)
            {
            if(strcmp(parsed_command[i+1], "yes") == 0)
                weighted = TRUE;
            if(strcmp(parsed_command[i+1], "no") == 0)
                weighted = FALSE;
			}
		
		
        if(strcmp(parsed_command[i], "savetrees") == 0)
            {
            strcpy(filename, parsed_command[i+1]);
            }

		if(strcmp(parsed_command[i], "usepaup") == 0)
            {
			if(strcmp(parsed_command[i+1], "yes") == 0)
                calculate_inhouse = FALSE;
            if(strcmp(parsed_command[i+1], "no") == 0)
                calculate_inhouse = TRUE;
            }

    
        }




    if(!error)
        {
		temptrees = malloc(remainingtrees*sizeof(char *));
		for(i=0; i<remainingtrees; i++)
			{
			temptrees[i] = malloc(TREE_LENGTH*sizeof(char));
			temptrees[i][0] = '\0';
			}
		j=0;
		for(i=0; i<Total_fund_trees; i++)
			{
			if(sourcetreetag[i])
				{
				strcpy(temptrees[j], fundamentals[i]);
				j++;
				}
			}
				
		position = MRP_matrix(temptrees, remainingtrees, FALSE);

        num_partitions = (float)position;


        if(criterion==1)   /* if we are to use this coding scheme to do an MRP analysis in paup **/
            {
			if(calculate_inhouse)
				{
				fprintf(BR_file, " %d %d\n", number_of_taxa, position);
				for(i=0; i<number_of_taxa; i++)
					{
					fprintf(BR_file, "%-10d", i);
							
					for(j=0; j<position; j++)
						{
						if(total_coding[j][i] != 3)
							fprintf(BR_file, "%d ",total_coding[j][i]);
						else
							fprintf(BR_file, "? ");
						}
					fprintf(BR_file, "\n");
					}
				
				
				}
			else
				{
				if(nrep ==0 || nrep==1)fprintf(BR_file, "#NEXUS\n");     
				if(nrep ==0 || nrep==1)fprintf(BR_file, "[Baum-Ragan coding ]\n");
				fprintf(BR_file, "\n\nbegin data;\n");
				fprintf(BR_file, "\tdimensions ntax=%d nchar=%d;\n", number_of_taxa, position);
				fprintf(BR_file, "\tformat interleave missing=?;\n");
				fprintf(BR_file, "\tmatrix\n");
				
				for(i=0; i<number_of_taxa; i++)
					{
					fprintf(BR_file, "\'%-20.20s\' ", taxa_names[i]);
							
					for(j=0; j<position; j++)
						{
						if(total_coding[j][i] != 3)
							fprintf(BR_file, "%d ",total_coding[j][i]);
						else
							fprintf(BR_file, "? ");
						}
					fprintf(BR_file, "\n");
					}
				fprintf(BR_file, "\t;\nend;\n\nbegin paup;\n");
				
				if(weighted)
					{
					/*** Weight each of the characters depending on the size ofthe tree they come from ***/
					fprintf(BR_file, "\tweights ");
					j=0;
					for(i=0; i<Total_fund_trees; i++)
						{
						if(sourcetreetag[i])
							{
							k=0; nodecount = 0;
							while(fundamentals[i][k] != ';')
								{
								if(fundamentals[i][k] == '(') nodecount++;
								k++;
								}
							nodecount--;
							if(nodecount > 0)
								{
								if(i == Total_fund_trees-1)
									fprintf(BR_file, "%f:%d-%d", (float)(tree_weights[i]*1.0/(float)nodecount), j+1, j+nodecount);
								else
									fprintf(BR_file, "%f:%d-%d, ", (float)(tree_weights[i]*1.0/(float)nodecount), j+1, j+nodecount);
								j += nodecount;
								}
							}
						}
					fprintf(BR_file, ";\n");
					}
				
				if(search==2) /* if instead we want to do a ptp test on our data */
					{
					fprintf(BR_file,"\tset increase=auto notifybeep=no errorbeep=no;\n\tconstraints one=(%s,%s,(%s,%s));\n\tpermute NReps=%d;\n\tquit;\n\tend;\n", taxa_names[0], taxa_names[1], taxa_names[2], taxa_names[3], ptpreps);
					}
				else
					{
					if(search==1) 
						{
						fprintf(BR_file,"set increase=auto notifybeep=no errorbeep=no;\n");
						if(!njbuild)
							{
							fprintf(BR_file,"\ths nreps=%d swap=", nreps );
							if(swap==1) fprintf(BR_file,"nni "); if(swap==2) fprintf(BR_file, "spr "); if(swap==3) fprintf(BR_file, "tbr ");
							fprintf(BR_file, "addseq=");
							if(addseq==1) fprintf(BR_file, "simple"); if(addseq==2) fprintf(BR_file, "closest"); if(addseq==3) fprintf(BR_file, "asis"); if(addseq == 4) fprintf(BR_file, "random"); if(addseq==5) fprintf(BR_file, "furthest");
							fprintf(BR_file, ";\n\tshowtrees;\n");
							}
						else
							fprintf(BR_file, "nj;\n");
						fprintf(BR_file, "\tsavetrees FILE=");
						if(strcmp(filename, "") == 0) fprintf(BR_file,"MRP.tree Format=nexus treeWts=yes");    
						else fprintf(BR_file,"%s Format=nexus treeWts=yes", filename);
						}
					if(search==0) 
						{
						fprintf(BR_file,"set increase=auto notifybeep=no errorbeep=no;\n\talltrees;\n\tshowtrees;\n\tsavetrees FILE=");
						if(strcmp(filename, "") == 0) fprintf(BR_file,"MRP.tree Format=nexus treeWts=yes");    
						else fprintf(BR_file,"%s Format=nexus treeWts=yes", filename);
						}
					if(nrep==0) fprintf(BR_file," Append=no replace=yes;\n\tquit;\n\nend;\n");  /* if there is one 1 rep */
					if(nrep==1) fprintf(BR_file," Append=no replace=yes;\n\nend;\n");     	/* if this is the forst of a few reps */
					if(nrep==2) fprintf(BR_file," Append=yes;\n\nend;\n");			/* if this is a middle rep */
					if(nrep==3) fprintf(BR_file," Append=yes;\n\tquit;\n\nend;\n");		/* if this is the last of a few reps */
					}
				}
			}
        
        if(num_fund_taxa != NULL) free(num_fund_taxa);
        }
	fflush(BR_file);
	/*if(calculate_inhouse == FALSE) error = 3; */
	free(string);
    return(error);
    }

int MRP_matrix(char **trees, int num_trees, int consensus)
	{
	int i=0, j=0, k=0, count =0, x=0, *tracking = NULL, total = 0, **BR_coding = NULL, nodecount = 0, position = 0, split_count = 0;
    char number[100];
    char *string = malloc(TREE_LENGTH * sizeof(char));
	int *num_fund_taxa = NULL, one_in_this = FALSE, zero_in_this = FALSE;

	if(!string) { printf2("Error: out of memory in MRP_matrix\n"); return 0; }
	num_fund_taxa = malloc(num_trees*sizeof(int));
	if(!num_fund_taxa) memory_error(68);
	for(i=0; i<num_trees; i++) num_fund_taxa[i] = 0;
	
	if(coding_from_tree != NULL) 
		{
		free(coding_from_tree);
		coding_from_tree=NULL;
	}

	if(total_coding != NULL)
		{
		for(i=0; i<total_nodes; i++)
			free(total_coding[i]);
		free(total_coding);
		total_coding = NULL;
		}
	/* Calculated the number of nodes in all the fundamental trees */
	total_nodes= num_trees*(number_of_taxa - 2);
	/* dynamically allocate the array to hold the resulting the total Baum-ragan coding scheme */
	total_coding = malloc((total_nodes)*sizeof(int *));
	if(!total_coding) memory_error(80);
	for(i=0; i<total_nodes ; i++)
		{
		total_coding[i] = malloc((number_of_taxa)*sizeof(int));
		if(!total_coding[i]) memory_error(81);
		}
	for(i=0;i<total_nodes; i++)
			for(j=0; j<number_of_taxa; j++)
				total_coding[i][j] =0;

	/* Allcoate an aray that records which tree each partition comes from */
	coding_from_tree= malloc(total_nodes*sizeof(int));
	for(i=0; i<total_nodes ; i++) coding_from_tree=0;

	/* for this we need a new array that holds the number of times that any partition is repeated in the array */
	if(partition_number != NULL) free(partition_number);
	partition_number = malloc(total_nodes*sizeof(float));
	if(!partition_number) memory_error(67);
	for(i=0; i<total_nodes; i++) partition_number[i] = 1;   /* the score of 1 gives each partition the same weight, it means that larger trees will have more say, in the quartet sense it will be very biased, so will have to be modified for this case. **/
	/* for this we need a new array that records which source tree each partition comes from (so we don't include the same quartet from any tree more than once ) */
	if(from_tree != NULL) free(from_tree);
	from_tree = malloc(total_nodes*sizeof(int));
	if(!from_tree) memory_error(69);
	for(i=0; i<total_nodes; i++) from_tree[i] = 0;   
	for(x=0; x<num_trees; x++)
		{
		unroottree(trees[x]);
		string[0] = '\0';
		strcpy(string, trees[x]);

		nodecount = 0;
		i=0;
		while(string[i] != ';')
			{
			if(string[i] == '(') nodecount++;
			i++;
			}
		
		/* dynamically allocate the array to hold the resulting Baum-ragan coding scheme */
		BR_coding = NULL;
		BR_coding = malloc((nodecount)*sizeof(int *));
		if(!BR_coding) memory_error(77);
		for(i=0; i<nodecount; i++)
			{
			BR_coding[i] = malloc((number_of_taxa)*sizeof(int));
			if(!BR_coding[i]) memory_error(78);
			}
		for(i=0;i<nodecount; i++)
				for(j=0; j<number_of_taxa; j++)
					BR_coding[i][j] =0;
		/* Allocate the array that records the clade were are in during calculation */
		tracking = malloc((2*nodecount)*sizeof(int));
		if(!tracking) memory_error(79);
		for(i=0; i<(2*nodecount); i++) tracking[i] = FALSE;
		i=0;
		count = 0;
		/* we scan the nested parenthesis string once from left to right calculating the baum-ragan coding scheme */
		while(string[i] != ';')
			{
			
			switch(string[i])
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
					while(string[i] != '(' && string[i] != ')' && string[i] != ',' && string[i] != ';' &&  string[i] != ':')
						i++;
					break;
	
				case ',':
					i++;
					break;
				
				case ':':
					while(string[i] != '(' && string[i] != ')' && string[i] != ',' && string[i] != ';')
						i++;
					break;
	
				default :
					num_fund_taxa[x]++;
					for(j=0; j<100; j++)
						number[j] = '\0';
					j=0;
					while(string[i] != '(' && string[i] != ')' && string[i] != ',' && string[i] != ':')
						{
						number[j] = string[i];
						i++;
						j++;
						}
					total = 0;
					if(consensus)
						{
						total = assign_taxa_name(number, FALSE);
						}
					else
						{
						j--; k=j;
						while(j>=0)
							{
							total+=(((int)pow(10,k-j))*texttoint(number[j]));  /* Turn the character to an integer */
							j--;
							}
						}
					for(j=0; j<nodecount; j++)
						{
						if(tracking[j] == TRUE)
							{
							BR_coding[j][total] = 1;
							}
						}
				}
			}
		if(!consensus)
			{
			/*** Put ? at those taxa that are not present **/
			for(i=0; i<number_of_taxa; i++)
				{
				if(!presence_of_taxa[x][i])
					{
					for(j=0; j<nodecount; j++)
						{
						BR_coding[j][i] = 3;
						}
					}
				}
			}
		split_count = 0;
		for(i=0; i<nodecount; i++)
			{
			one_in_this = FALSE; zero_in_this = FALSE;
			for(j=0; j<number_of_taxa; j++)
				{
				if(BR_coding[i][j] == 1) one_in_this = TRUE;
				if(BR_coding[i][j] == 0) zero_in_this = TRUE;
				}
			if(one_in_this && zero_in_this) split_count++;
			
			}
		/*** Copy the coding from this tree into the total coding  **/
		for(i=0; i<nodecount; i++)
			{
			one_in_this = FALSE; zero_in_this = FALSE;
			for(j=0; j<number_of_taxa; j++)
				{
				if(BR_coding[i][j] == 1) one_in_this = TRUE;
				if(BR_coding[i][j] == 0) zero_in_this = TRUE;
				}
			
			if(one_in_this && zero_in_this)
				{
				from_tree[position]=x;
				for(j=0; j<number_of_taxa; j++)
					{
					total_coding[position][j] = BR_coding[i][j];
					}
				
				if(consensus)
					{
					partition_number[position] = score_of_bootstraps[x];   /** if this is being used for a consensus then we want the partition to have a weight equal to the inverse of the number of trees recorded during that replicate */
					}
				else
					{
					
					if(criterion == 3)
						{
						if(quartet_normalising == 1)
							partition_number[position] = 1*tree_weights[x];   /* equal weighting */
						if(quartet_normalising == 2)
							partition_number[position] = (1/((float)num_fund_taxa[x]-3))*tree_weights[x];    /* this is to try to normalise the score that any partition contributes ...... in the quartet method it will help to normalise for large tree bias, It may be ignored for the split compatibility analysis **/
						if(quartet_normalising == 3)
							{
							partition_number[position] = (1/((float)(num_fund_taxa[x]*(num_fund_taxa[x]-1)*(num_fund_taxa[x]-2)*(num_fund_taxa[x]-3))/(float)24))*tree_weights[x];
							}
						}
					if(criterion == 2)
						{
						if(splits_weight == 1)
							partition_number[position] = 1*tree_weights[x]; 	/* equal weighting */
						if(splits_weight == 2)
							partition_number[position] = (1/((float)split_count))*tree_weights[x];	/* weighted by number of splits in the source tree */
						}
					}
				position++;
				}
			}
		/** free up the memory space for BR_coding **/
		for(i=0; i<nodecount; i++)
			{
			free(BR_coding[i]);
			}
		free(BR_coding);
		BR_coding = NULL;
		free(tracking);
		tracking = NULL;

		}
	free(string);
	return(position);
	}

void condense_coding(void)
    {
    int i=0, j=0, k=0, l=0, same1 = TRUE, same2 = TRUE;
    
    
    /* what we need to do is check every defined partition in the total_coding array to see it it is repeated anywhere else in the array */
    for(i=0; i<num_partitions-1; i++)
        {
        if(total_coding[i][0] != 4)
            {
            for(j=i+1; j<num_partitions; j++)
                {
                if(total_coding[j][0] != 4)
                    {
                    same1 = TRUE; same2=TRUE;
                    for(k=0; k<number_of_taxa; k++)
                        {
                        if(total_coding[i][k] != total_coding[j][k])  /* if this not the same partition */
                            {
                            same1 = FALSE;
                            }
                        if(total_coding[i][k] == 3 || total_coding[j][k] == 3)  /* this section is to see if its the same partition but coded in the opposite manner (0s switched with 1s)  */
                            {
                            if(total_coding[j][k] != 3 || total_coding[i][k] != 3)
                                same2 = FALSE;
                            }
                        else
                            {
                            if(total_coding[i][k] == total_coding[j][k])
                                {
                                same2 = FALSE;
                                }
                            }
                        }
                    if(same1 || same2)
                        {
                        partition_number[i]++;
                        for(k=0; k<number_of_taxa; k++)
                            total_coding[j][k] = 4;
                        }
                    }
                }
            }
        }
 
 
    /* Next we need to bunch up the array, removing duplicates of the same partitions */
    k=0;l=1;
    for(i=0; i<num_partitions-1; i++)
        {
        if(total_coding[i][0] == 4)
            {
            l=i+1;
            while(l < num_partitions && total_coding[l][0] == 4) l++;
            if(l != num_partitions)
                {
                for(j=0; j<number_of_taxa; j++)
                    {
                    total_coding[i][j] = total_coding[l][j];
                    total_coding[l][j] = 4;
                    }
                partition_number[i] = partition_number[l];
                k++;
                }
            else
                i = num_partitions;
            }
        else
            k++;
        }

    k=0;
    for(i=0; i<num_partitions; i++)
        if(total_coding[i][0] != 4) k++;
        
     
    total_partitions= num_partitions;
    /* k isthe number of partitions that have been turned into 4's, this is the number of repeated partitions */            
    num_partitions = k;
    }

void do_consensus(void)
	{
	int tree_type = 0, i, j, k, l,m, numtrees = 0, present = TRUE, number, error = FALSE, useguide = FALSE;
	char **temptrees = NULL, c, tempname[1000], consensusfilename[1000], guidetreename[1000]; 
	float percentage = 0;
	FILE *consensusfile = NULL, *guidetreefile = NULL;
	
	guidetreename[0] = '\0';
	consensusfilename[0] = '\0';
	strcpy(consensusfilename, "consensus.ph");
	 for(i=0; i<num_commands; i++)
        {
		if(strcmp(parsed_command[i], "filename") == 0)
			{
			strcpy(consensusfilename, parsed_command[i+1]);
			}
		if(strcmp(parsed_command[i], "guidetree") == 0)
			{
			strcpy(guidetreename, parsed_command[i+1]);
			if(guidetreefile != NULL) fclose(guidetreefile);
			if((guidetreefile = fopen(guidetreename, "r")) == NULL)		/* check to see if the file is there */
				{								/* Open the source tree file */
				printf2("Cannot open file %s\n", guidetreename);
				error = TRUE;
				}
			useguide = TRUE;
			}
		if(strcmp(parsed_command[i], "method") == 0)
			{
			if(strcmp(parsed_command[i+1], "strict") == 0) percentage = 1.0;
			else
				{
				if(strcmp(parsed_command[i+1], "majrule") == 0) percentage = 0.5;
				else
					{
					if(strcmp(parsed_command[i+1], "minor") == 0) percentage = 0;
					else
						{
						percentage = tofloat(parsed_command[i+1]);
						if(percentage > 1)
							{
							printf2("Error: The cut off for a consensus tree must be between .5 and 1.0\n");
							error = TRUE;
							}
						}
					}
				}
			}


		if(strcmp(parsed_command[i], "data") == 0)
			{
			if(strcmp(parsed_command[i+1], "source") == 0)
				tree_type = 0;
			else
				{
				if(strcmp(parsed_command[i+1], "supertrees") == 0)
					tree_type = 1;
				else
					{
					if(strcmp(parsed_command[i+1], "bootstraps") == 0)
						tree_type = 2;
					else
						{
						printf2("Error: %s is an invalid option for data\n", parsed_command[i+1]);
						error = TRUE;
						}
					}
				}
			}
		}
	
	
	if(!error)
		{
		
					
			
			
		if(tree_type == 0)  /* do a consensus of the universally distributed source trees */
			{
			/*** count the nmber of universally distributed source trees */
			for(i=0; i<Total_fund_trees; i++)
				{
				if(sourcetreetag[i])
					{
					present = TRUE;
					for(k=0; k<number_of_taxa; k++)
						{
						if(presence_of_taxa[i][k] == FALSE)
							present = FALSE;
						}
					if(present) numtrees++;
					}
				}
				
			if(numtrees > 0)
				{
				
				if(score_of_bootstraps != NULL)
					free(score_of_bootstraps);
				score_of_bootstraps = malloc(numtrees*sizeof(float));
				tempname[0] = '\0';
				temptrees = malloc(numtrees*sizeof(char*));
				for(i=0; i<numtrees; i++)
					{
					score_of_bootstraps[i] = 1;
					temptrees[i] = malloc(TREE_LENGTH*sizeof(char));
					temptrees[i][0] = '\0';
					}
				j=0; k=0;
				for(i=0; i<Total_fund_trees; i++)
					{
					if(sourcetreetag[i])
						{
						present = TRUE;
						for(k=0; k<number_of_taxa; k++)
							{
							if(presence_of_taxa[i][k] == FALSE)
								present = FALSE;
							}
						if(present)
							{
							k=0;l=0;
							while(fundamentals[i][k] != ';')
								{
								switch(fundamentals[i][k])
									{
									case '(':
									case ',':
										temptrees[j][l] = fundamentals[i][k];
										k++; l++;
										break;
									case ')':
										temptrees[j][l] = fundamentals[i][k];
										k++; l++;
										while(fundamentals[i][k] != ')' && fundamentals[i][k] != '(' && fundamentals[i][k] != ':' && fundamentals[i][k] != ',' && fundamentals[i][k] != ';')
											{
											temptrees[j][l] = fundamentals[i][k];
											k++; l++;
											}
										break;
									case ':':
										temptrees[j][l] = fundamentals[i][k];
										k++; l++;
										while(fundamentals[i][k] != ')' && fundamentals[i][k] != '(' && fundamentals[i][k] != ':' && fundamentals[i][k] != ',' && fundamentals[i][k] != ';')
											{
											temptrees[j][l] = fundamentals[i][k];
											k++; l++;
											}
										
										break;
									default:
										m = 0;
										while(fundamentals[i][k] != ')' && fundamentals[i][k] != '(' && fundamentals[i][k] != ':' && fundamentals[i][k] != ',' && fundamentals[i][k] != ';')
											{
											tempname[m] = fundamentals[i][k];
											m++; k++;
											}
										tempname[m] = '\0';
										number = toint(tempname);
										m=0;
										while(taxa_names[number][m] != '\0')
											{
											temptrees[j][l] = taxa_names[number][m];
											m++; l++;
											}
										break;
									}
								}
							temptrees[j][l] = ';';
							temptrees[j][l+1] = '\0';
							j++;
						

							}
						}
					}
				printf2("\nConsensus settings:\n");
				printf2("\tConsensus of %d universally distributed source trees\n", numtrees);
				if(useguide)
					{
					printf2("\tOnly relationships as defined by the guidetree in %s ", guidetreename);
					}
				else
					{
					printf2("\tOnly relationships with ");
					if(percentage == 1) printf2("100%% support ");
					if(percentage == 0) printf2("50%% support or greater (including congruent minor components) ");
					if(percentage != 0 && percentage != 1) printf2("%0.0f%% support or greater ", percentage*100);
					}
				printf2("are included in the consensus\n");
				printf2("\tConsensus file = %s\n", consensusfilename);
				
				consensusfile = fopen(consensusfilename, "w");
				
				consensus(numtrees, temptrees, numtrees, percentage, consensusfile, guidetreefile);
				
				fclose(consensusfile);
				}
			else
				printf2("There are no tress that contain all the taxa, unable to construct a consensus tree\n");
					
			}
		}
	
	}

void consensus(int num_trees, char **trees, int num_reps, float percentage, FILE *outfile, FILE *guidetreefile)
    {
    int i, j, k, q, r, same1 = FALSE, same2 = FALSE, same3 = FALSE, same4 = FALSE,same5 = FALSE,same6 = FALSE,same7 = FALSE,same8 = FALSE, l, **sets, *in, *tmpcoding = NULL, end = FALSE, found = -1;
	/* The first thing needed is to create a Baum-Ragan coding scheme holding all the information from the bootstrapped trees */
	char **string, *tmp = NULL, name[100], *rest = malloc(TREE_LENGTH * sizeof(char)), value[100];
	if(!rest) { printf2("Error: out of memory in consensus\n"); return; }
	int count, first = -1, support = 0, **shorthand = NULL, subdivisions = ((int)(number_of_taxa/16))+1;
	float tmpnumber = 0;
	
	
    num_partitions = MRP_matrix(trees, num_trees, TRUE);
/*	shorthand = malloc(num_partitions*sizeof(int *));
	for(i=0; i<num_partitions; i++)
		{
		shorthand[i] = malloc(subdivisions*sizeof(int));
		for(j=0; j<subdivisions; j++) shorthand[i][j] = 0;
		}
*/		
	tmp = malloc(TREE_LENGTH*sizeof(char));
	tmp[0] = '\0';
	name[0] = '\0';
	rest[0] = '\0';
	value[0] = '\0';
	
	/** calculated the binary evivalent in blocks of 16 bits */
	/*  for(i=0; i<num_partitions-1; i++)
        {
        if(total_coding[i][0] != 4)
            {
			k=-1; q=0;
			for(j=0; j<number_of_taxa; j++)
				{
				if(j%16 == 0)
					{
					k++; shorthand[i][k]=0;
					if(k != 0) printf2("%d ", shorthand[k-1]);
					}
				shorthand[i][k] += (((int)(pow(2, j%16))) * total_coding[i][j]);
				}
			 printf2("%d ", shorthand[i][k]);
			printf2("\n");
			}
		}
	 */
	   for(i=0; i<num_partitions-1; i++)
        {
        if(total_coding[i][0] != 4)
            {
            for(j=i+1; j<num_partitions; j++)
                {
                if(total_coding[j][0] != 4)
                    {
                    same1 = TRUE; same2=TRUE;
                    for(k=0; k<number_of_taxa; k++)
                        {
                        if(total_coding[i][k] != total_coding[j][k])  /* if this not the same partition */
                            {
                            same1 = FALSE;
							if(same2 == FALSE) k = number_of_taxa;
                            }
                        if(total_coding[i][k] == 3 || total_coding[j][k] == 3)  /* this section is to see if its the same partition but coded in the opposite manner (0s switched with 1s)  */
                            {
                            if(total_coding[j][k] != 3 || total_coding[i][k] != 3)
								{
                                same2 = FALSE;
								if(same1 == FALSE) k = number_of_taxa;
								}
                            }
                        else
                            {
                            if(total_coding[i][k] == total_coding[j][k])
                                {
                                same2 = FALSE;
								if(same1 == FALSE) k = number_of_taxa;
                                }
                            }
                        }
                    if(same1 || same2)
                        {
                        partition_number[i] += partition_number[j];
                        for(k=0; k<number_of_taxa; k++)
                            total_coding[j][k] = 4;
                        }
                    }
                }
            }
        }
    /* Next we need to bunch up the array, removing duplicates of the same partitions */
    k=0;l=1;
    for(i=0; i<num_partitions-1; i++)
        {
        if(total_coding[i][0] == 4)
            {
            l=i+1;
            while(l < num_partitions && total_coding[l][0] == 4) l++;
            if(l != num_partitions)
                {
                for(j=0; j<number_of_taxa; j++)
                    {
                    total_coding[i][j] = total_coding[l][j];
                    total_coding[l][j] = 4;
                    }
                partition_number[i] = partition_number[l];
                k++;
                }
            else
                i = num_partitions;
            }
        else
            k++;
        }

    k=0;
    for(i=0; i<num_partitions; i++)
        if(total_coding[i][0] != 4) k++;
        
    total_partitions= num_partitions;
    /* k isthe number of partitions that have been turned into 4's, this is the number of repeated partitions */            
    num_partitions = k;

	/* sort the partitions according to their representation (score) */
	tmpcoding = malloc(number_of_taxa*sizeof(int));
	for(i=0; i<num_partitions; i++)
		{
		for(j=i+1; j<num_partitions; j++)
			{
			if(partition_number[i] < partition_number[j])
				{
				tmpnumber = partition_number[j];
				for(k=0; k<number_of_taxa; k++)
					tmpcoding[k] = total_coding[j][k];
				
				for(k=j; k>i; k--)
					{
					partition_number[k] = partition_number[k-1];
					for(q=0; q<number_of_taxa; q++)
						total_coding[k][q] = total_coding[k-1][q];
					}
				partition_number[i] = tmpnumber;
				for(k=0; k<number_of_taxa; k++)
						total_coding[i][k] = tmpcoding[k];
				j--;
				}
			
			}
		}

	free(tmpcoding);
/*** There are going to be two ways of summarising the results of a consensus, the user can ask for the best supported tree (default) or give a guidetree ("best tree") onto which the support tree is to be projected. */
	if(guidetreefile == NULL)	
		{
		/** From the remaining splits, only retain those that have at least x% representation (this would be 50% for maj-rule or 100% for strict) **/
		in = malloc(num_partitions*sizeof(int));  /* will record the partitions to be included */
		for(i=0; i<num_partitions; i++)
			in[i] = TRUE;
			
		if(percentage < 0.5)
			tmpnumber = 0.5;
		else
			tmpnumber = percentage;
		for(i=0; i<num_partitions; i++)
			{
			if((float)(partition_number[i]/num_reps) >= tmpnumber)
				{
				for(k=0; k<i; k++)
					{
					if(in[k] == TRUE)
						{
						same1 = same2 = same3 = same4 = same5 = same6 = same7 = same8 = TRUE;
						for(j=0; j<number_of_taxa; j++)
							{
							/* are all the 1s in this all 1s or all zeros in the earler one? */
							if(total_coding[i][j] == 1)
								{
								if(total_coding[i][j] != total_coding[k][j])
									same1 = FALSE;
								else
									same2 = FALSE;
								}
							if(total_coding[i][j] == 0)
								{
								if(total_coding[i][j] != total_coding[k][j])
									same3 = FALSE;
								else
									same4 = FALSE;
								}
							/* are all the 1s in this all 1s or all zeros in the earler one? */
							if(total_coding[k][j] == 1)
								{
								if(total_coding[i][j] != total_coding[k][j])
									same5 = FALSE;
								else
									same6 = FALSE;
								}
							if(total_coding[k][j] == 0)
								{
								if(total_coding[i][j] != total_coding[k][j])
									same7 = FALSE;
								else
									same8 = FALSE;
								}
							}
						if(same1 != TRUE && same2 != TRUE && same3 != TRUE && same4 != TRUE && same5 != TRUE && same6 != TRUE && same7 != TRUE && same8 != TRUE)
							in[i] = FALSE;
							
							
						}
					}
				}
			else
				in[i] = FALSE;
			}
			
	/**** Now do the same but identifying the minor components ***/
		if(percentage == 0)
			{
			for(i=0; i<num_partitions; i++)
				{
				if(in[i] == FALSE)
					{
					in[i] = TRUE;
					for(k=0; k<i; k++)
						{
						if(in[k] == TRUE)
							{
							same1 = same2 = same3 = same4 = same5 = same6 = same7 = same8 = TRUE;
							for(j=0; j<number_of_taxa; j++)
								{
								/* are all the 1s in this all 1s or all zeros in the earler one? */
								if(total_coding[i][j] == 1)
									{
									if(total_coding[i][j] != total_coding[k][j])
										same1 = FALSE;
									else
										same2 = FALSE;
									}
								if(total_coding[i][j] == 0)
									{
									if(total_coding[i][j] != total_coding[k][j])
										same3 = FALSE;
									else
										same4 = FALSE;
									}
								/* are all the 1s in this all 1s or all zeros in the earler one? */
								if(total_coding[k][j] == 1)
									{
									if(total_coding[i][j] != total_coding[k][j])
										same5 = FALSE;
									else
										same6 = FALSE;
									}
								if(total_coding[k][j] == 0)
									{
									if(total_coding[i][j] != total_coding[k][j])
										same7 = FALSE;
									else
										same8 = FALSE;
									}
								}
							if(same1 != TRUE && same2 != TRUE && same3 != TRUE && same4 != TRUE && same5 != TRUE && same6 != TRUE && same7 != TRUE && same8 != TRUE)
								in[i] = FALSE;
								
								
							}
						}
					}
				}
			}

		
		printf2("\n\n\nSets included in the consensus tree\n");
		for(i=0; i<num_partitions; i++)
			{
			if(in[i] && (float)(partition_number[i]/num_reps) >= percentage && (float)(partition_number[i]/num_reps) >= .5)
				{
				printf2("\t");
				for(j=0; j<number_of_taxa; j++)
					{
					if(total_coding[i][j] == 1)
						printf2("*");
					else
						printf2(".");
					}
				printf2("\t%.2f\n", (float)(partition_number[i]/num_reps));
				}
			}
		
		printf2("\n\n\nMinor Components included in the consensus tree\n");
		for(i=0; i<num_partitions; i++)
			{
			if(in[i] && (float)(partition_number[i]/num_reps) < .5)
				{
				printf2("\t");
				for(j=0; j<number_of_taxa; j++)
					{
					if(total_coding[i][j] == 1)
						printf2("*");
					else
						printf2(".");
					}
				printf2("\t%.2f\n", (float)(partition_number[i]/num_reps));
				}
			}
		
		printf2("\n\nSets not included in the consensus tree\n");
		for(i=0; i<num_partitions; i++)
			{
			if(!in[i])
				{
				printf2("\t");
				for(j=0; j<number_of_taxa; j++)
					{
					if(total_coding[i][j] == 1)
						printf2("*");
					else
						printf2(".");
					}
				printf2("\t%.2f\n", (float)(partition_number[i]/num_reps));
				}
			}
		/* build a nested parenthesis tree from the remaining splits and return. This is the result */
		
		/** next start to put together the tree from the matrix **/
		/** starting with those partitions with two 1s and working up to those partitions with n-2 partitions, put the tree together **/
		
		string = malloc(number_of_taxa*sizeof(char*));
		for(i=0; i<number_of_taxa; i++)
			{
			string[i] = malloc(TREE_LENGTH*sizeof(char));
			string[i][0] = '\0';
			strcpy(string[i], taxa_names[i]);
			}
		
		tmp[0] = '\0';
		for(i=2; i<number_of_taxa-2; i++)
			{
			same1 = FALSE;
			/** check each partition to see if it has i number of partitions */
			for(j=0; j<num_partitions; j++)
				{
				if(in[j])
					{
					count = 0;
					for(k=0; k<number_of_taxa; k++)
						{
						
						count += total_coding[j][k];
						}
					
					if(count == i || count == number_of_taxa-i)
						{
						/** if it is transposed we need to turn it back around **/
						if(count == number_of_taxa-i)
							{
							for(k=0; k<number_of_taxa; k++)
								{
								if(total_coding[j][k] == 1)
									total_coding[j][k] = 0;
								else
									total_coding[j][k] = 1;
								}
							}
						
						
						first = -1;
						for(k=0; k<number_of_taxa; k++)
							{
							if(total_coding[j][k] == 1 && first == -1)
								{
								first = k;
								strcpy(tmp, "(");
								strcat(tmp, string[first]);
								strcpy(string[first], tmp);
								}
							else
								{
								if(total_coding[j][k] == 1)
									{
									strcat(string[first], ",");
									strcat(string[first], string[k]);
									for(q=0; q<num_partitions; q++)
										total_coding[q][k] = 0;
									strcpy(string[k], "");
									
									}
								}
							}
						sprintf(tmp, ")%.2f", (float)(partition_number[j]/num_reps));
						strcat(string[first], tmp);
						same1 = TRUE;
						}
					}
				}
			if(same1) i =1;
			}
		strcpy(tmp, "(");
		first = FALSE;
		for(i=0; i<number_of_taxa; i++)
			{
			if(strcmp(string[i], "") != 0)
				{
				if(!first)
					{
					first = TRUE;
					strcat(tmp, string[i]);
					}
				else
					{
					strcat(tmp, ",");
					strcat(tmp, string[i]);
					}
				}
			}
		strcat(tmp, ");");	
		tree_coordinates(tmp, TRUE, TRUE, FALSE, -1);
		fprintf(outfile, "%s\n", tmp);
		free(in);
		for(i=0; i<number_of_taxa; i++)
			free(string[i]);
		free(string);
		}
	else  /* we are going to use a guide tree */
		{
		in = malloc(number_of_taxa*sizeof(int));
		for(i=0; i< number_of_taxa; i++) in[i] = FALSE;
		
		/** read in the guide tree **/
		i=0;
		tmp[i] = getc(guidetreefile);
		i++;
		while(!feof(guidetreefile))
			{
			tmp[i] = getc(guidetreefile);
			if(tmp[i-1] != ';')i++;
			}
		if(tmp[i] != ';' )tmp[i] = '\0';
		else tmp[i+1] = '\0';
		/** unroot the guide tree **/
		unroottree(tmp);
		/** identify the taxa inside each of the nested parentheses **/
		i=1; /* we don't do the first parenthesis */
		
		
		while(tmp[i] != ';' && tmp[i] != '\0')
			{
			switch(tmp[i])
				{
				case '(':
					/** find the corresponding closing parenthesis and remember all the taxa in between **/
					j=i+1;
					for(k=0; k<number_of_taxa; k++) in[k] = FALSE;
					k=1;
					end = FALSE;
					while(tmp[j] != '\0' && !end)
						{
						switch(tmp[j])
							{
							case '(':
								k++;
								j++;
								break;
							case ')':
								k--;
								j++;
								if(k == 0) end = TRUE;
								while(tmp[j] != '(' && tmp[j] != ')' && tmp[j] != ',' && tmp[j] != ';' && tmp[j]!= ':' && tmp[j] != '\0') j++;
								break;
							case ',':
								j++;
								break;
							case ':':
								while(tmp[j] != '(' && tmp[j] != ')' && tmp[j] != ',' && tmp[j] != ';' && tmp[j] != '\0') j++;
								break;
							default:
								q=0;
								while(tmp[j] != '(' && tmp[j] != ')' && tmp[j] != ',' && tmp[j] != ';' && tmp[j]!= ':' && tmp[j] != '\0')
									{
									name[q] = tmp[j];
									j++; q++;
									}
								name[q] = '\0';
								/** find this name **/
								for(q=0; q<number_of_taxa; q++)
									{
									if(strcmp(name, taxa_names[q]) == 0)
										{
										in[q] = TRUE;
										q = number_of_taxa;
										}
									}
								break;
							}
						}

					
					/** we now have a list of all the taxa inside this set of parentheses, j points to where we need to put the support found **/
					support = 0;
					/** find the partition that matches the set of taxa in "in" */
					
					found = -1;
					for(k=0; k<num_partitions; k++)
						{
						if(total_coding[k][0] != 4)
							{
							same1 = same2 = TRUE;
							for(q=0; q<number_of_taxa; q++)
								{
								/* are all the 1s in this all 1s or all zeros in the earler one? */
								if(in[q] == TRUE)
									{
									if(total_coding[k][q] == 1)
										same1 = FALSE;
									else
										same2 = FALSE;
									}
								if(in[q] == FALSE)
									{
									if(total_coding[k][q] == 0)
										same1 = FALSE;
									else
										same2 = FALSE;
									}
								}
							if(same1 == TRUE || same2 == TRUE)
								{
								/* this is the partition */
								found = k;
								k = num_partitions;
								}
							}
						}
					/* record everything from position j to the end into the array "rest" **/
					q=0; r = j;
					while(tmp[r] != '\0')
						{
						rest[q] = tmp[r];
						tmp[r] = '\0';
						q++; r++;
						}
					rest[q] = '\0';
					if(found != -1) sprintf(value, "%.2f", (float)(partition_number[found]/num_reps));
					else sprintf(value, "0.00");
					strcat(tmp, value);
					strcat(tmp, rest);
					i++;
					break;
				default:
					i++;
					break;
				}
			}
		tree_coordinates(tmp, TRUE, TRUE, FALSE, -1);
		fprintf(outfile, "%s\n", tmp);
		free(in);
		}
	free(rest);
	free(tmp);
    }

void print_descendents(struct taxon *position, FILE *outfile)
	{
	
	while(position != NULL)
		{
		if(position->daughter != NULL)
			{
			fprintf(outfile, "%s\t", position->weight);
			print_descendents(position->daughter, outfile);
			}
		else
			{
			fprintf(outfile, "%s\t", taxa_names[position->name]);
			}
		position = position->next_sibling;
		}
	}

void do_descendents(struct taxon *position, FILE *outfile)
	{
	struct taxon *start = position;
	
	if(start->parent != NULL)  /*If we are not at the top of the tree */
		{
		while(position != NULL)
			{
			if(position->daughter != NULL)
				{
				fprintf(outfile, "%s\t", position->weight);
				print_descendents(position->daughter, outfile);
				}
			else
				{
				fprintf(outfile, "%s\t", taxa_names[position->name]);
				}
			position = position->next_sibling;
			}
		}
	position = start;
	while(position != NULL)
		{
		if(position->daughter != NULL)
			{
			fprintf(outfile, "\n%s:\t", position->weight);
			do_descendents(position->daughter, outfile);
			}
		position = position->next_sibling;
		}
	
	
	}

