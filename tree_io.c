/*
 *  tree_io.c — Clann v5.0.0
 *  Source-tree file reading, parsing, display, and filter commands.
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

#include "tree_io.h"
#include "viz.h"

/* Private (static) forward declarations */
static void quick(float **items, int count);
static void qs(float **items, int left, int right);
static void prune_taxa_for_exclude(struct taxon *super_pos, int *tobeexcluded);

void input_fund_tree(char *intree, int fundnum)
	{
	int i=0,j=0, k=0, r=0, tree_length = 0, taxaposition = 0, tottaxaintree = 0;
	char temp[NAME_LENGTH], **temp_tree_names = NULL;
	float *temp_tree_weights = NULL;

	/* count the number of commas in the tree string to figure out the number of taxa in the tree */
	i=0;
	while(intree[i] != ';' )
		{
		if(intree[i] == ',') tottaxaintree++;
		i++;
		}
	tottaxaintree++;

	fulltaxanames[fundnum] = malloc(tottaxaintree*sizeof(char *));
	for(i=0; i<tottaxaintree; i++) fulltaxanames[fundnum][i] = NULL;
	numtaxaintrees[fundnum]=tottaxaintree;
	tree_length = strlen(intree);
	fundamentals[Total_fund_trees] = malloc((tree_length+100)*sizeof(char));
	if(!fundamentals[Total_fund_trees]) memory_error(119);
	fundamentals[Total_fund_trees][0] = '\0';
	i=0; r=0;
	while(intree[r] != ';' )
		{
		switch(intree[r])
			{
			case ')':
				fundamentals[Total_fund_trees][i] = intree[r];
				r++;
				i++;
				while(intree[r] != ')' && intree[r] != '(' && intree[r] != ',' && intree[r] != ';' && intree[r] != ':')
					{
					fundamentals[Total_fund_trees][i] = intree[r];
					r++; i++;
					}
				break;
			case '(':
			case ',':
			case ';':
				fundamentals[Total_fund_trees][i] = intree[r];
				r++;
				i++;
				break;
			case ' ':
			case '\n':
			case '\r':
			case '\t':
				r++;
				break;
			case ':':
				got_weights = TRUE;
				while(intree[r] != ')' && intree[r] != '(' && intree[r] != ',' && intree[r] != ';')
					{
					fundamentals[Total_fund_trees][i] = intree[r];
					r++;
					i++;
					}
				break;
			case '\'':
				r++;
				break;
			default :
				j=0;
				while(intree[r]  !=  ')' && intree[r] != '(' && intree[r] != ',' && intree[r] != ';' && intree[r] != ':' && intree[r] != '\'')
					{
					if(j<max_name_length) temp[j] = intree[r];
					j++;
					r++;
					}
				temp[j] = '\0';
				fulltaxanames[fundnum][taxaposition] = malloc((strlen(temp)+10)*sizeof(char));
				if(fulltaxanames[fundnum][taxaposition] == NULL) memory_error(444);
				fulltaxanames[fundnum][taxaposition][0] = '\0';
				strcpy(fulltaxanames[fundnum][taxaposition], temp);
				taxaposition++;
				k = assign_taxa_name(temp, TRUE);   /* this returns an int but we need it in a string format */
				totext(k, temp);		/* This will give us the string euivalent */
				j=0;
				while(temp[j] != '\0')
					{
					fundamentals[Total_fund_trees][i] = temp[j]; /* copy the string equivalent into the tree in memory */
					j++; i++;
					}

				break;
			}

		}

	fundamentals[Total_fund_trees][i] = ';';
	fundamentals[Total_fund_trees][i+1] = '\0';
	for(j=0; j<TAXA_NUM*name_assignments; j++)
		{
		if(same_tree[j] > 0)
			{
			number_of_comparisons[Total_fund_trees]++;  /* This will tell me how many taxa are in this tree  */
			for(k=j+1; k<TAXA_NUM*name_assignments; k++)
				{
				if(same_tree[k] > 0)
					{
					Cooccurrance[j][k]++;
					Cooccurrance[k][j]++;
					}
				}
			}
		same_tree[j] = 0;
		}
	k =number_of_comparisons[Total_fund_trees];
	if(k>largest_tree) largest_tree = k;
	if(k<smallest_tree) smallest_tree = k;
	number_of_comparisons[Total_fund_trees] = 0;
	for(j=0; j<k; j++)
		number_of_comparisons[Total_fund_trees] +=j;    /* using the number of taxa in the tree, calculate the total number of comparisons in this tree */


	Total_fund_trees++;
	/**************************/
	if(Total_fund_trees == (fundamental_assignments*FUNDAMENTAL_NUM)-2)  /* If we need to expand the size of the fundamentals array */
		{


		fundamental_assignments++;

		temp_tree_names = malloc(((fundamental_assignments-1)*FUNDAMENTAL_NUM)*sizeof(char *));
		for(i=0; i< ((fundamental_assignments-1)*FUNDAMENTAL_NUM); i++)
			{
			temp_tree_names[i] = malloc(1000*sizeof(char));
			temp_tree_names[i][0] = '\0';
			strcpy(temp_tree_names[i], tree_names[i]);
			free(tree_names[i]);
			}

		free(tree_names);
		tree_names = NULL;
		tree_names = malloc((fundamental_assignments*FUNDAMENTAL_NUM)*sizeof(char *));
		if(!tree_names) memory_error(102);
		temp_tree_weights = malloc(((fundamental_assignments-1)*FUNDAMENTAL_NUM)*sizeof(float));
		for(i=0; i<((fundamental_assignments-1)*FUNDAMENTAL_NUM); i++)
			temp_tree_weights[i] = tree_weights[i];
		free(tree_weights);
		tree_weights = NULL;

		tree_weights = malloc((fundamental_assignments*FUNDAMENTAL_NUM)*sizeof(float));
		if(!tree_weights) memory_error(103);
		for(i=0; i<(fundamental_assignments*FUNDAMENTAL_NUM); i++)
			{
			tree_names[i] = malloc(1000*sizeof(char));
			if(!tree_names[i]) memory_error(104);
			tree_names[i][0] = '\0';
			tree_weights[i] = 1;
			if(i<((fundamental_assignments-1)*FUNDAMENTAL_NUM))
				{
				strcpy(tree_names[i], temp_tree_names[i]);
				free(temp_tree_names[i]);
				tree_weights[i] = temp_tree_weights[i];
				}
			}
		free(temp_tree_names);
		temp_tree_names = NULL;
		free(temp_tree_weights);
		temp_tree_weights = NULL;
		number_of_comparisons = realloc(number_of_comparisons, (fundamental_assignments*FUNDAMENTAL_NUM)*sizeof(int));
		if(!number_of_comparisons) memory_error(105);
		for(i=((fundamental_assignments-1)*FUNDAMENTAL_NUM); i<(fundamental_assignments*FUNDAMENTAL_NUM); i++)
			number_of_comparisons[i] = 0;

		 /* store the contents of the fundamentals so far */
		stored_funds = NULL;
		stored_funds = malloc((((fundamental_assignments-1)*FUNDAMENTAL_NUM)-2)*sizeof(char *));
		if(!stored_funds) memory_error(43);
		for(i=0; i<((fundamental_assignments-1)*FUNDAMENTAL_NUM)-2; i++)
			{
			stored_funds[i] = NULL;
			j = (int)strlen(fundamentals[i]);
			stored_funds[i] = malloc(j+100*sizeof(char));
			if(!stored_funds[i]) memory_error(44);
			stored_funds[i][0] = '\0';
			strcpy(stored_funds[i], fundamentals[i]);
			}

		/***** free up the memory space taken by fundmaental ***/
		for(i=0; i<((fundamental_assignments-1)*FUNDAMENTAL_NUM); i++)
			{

			if(fundamentals[i] != NULL)free(fundamentals[i]);
			}
				free(fundamentals);
		fundamentals = NULL;

		/***** allocate the new size of fundamentals ********/
		fulltaxanames = realloc(fulltaxanames, (fundamental_assignments*FUNDAMENTAL_NUM)*sizeof(char **));
		for(i=(fundamental_assignments-1)*FUNDAMENTAL_NUM; i<fundamental_assignments*FUNDAMENTAL_NUM; i++)
			fulltaxanames[i] = NULL;

		numtaxaintrees =realloc(numtaxaintrees, (fundamental_assignments*FUNDAMENTAL_NUM)*sizeof(int));
		fundamentals = malloc((fundamental_assignments*FUNDAMENTAL_NUM)*sizeof(char *));
		if(!fundamentals) memory_error(12);

		for(i=0; i< fundamental_assignments*FUNDAMENTAL_NUM; i++)
			{
			fundamentals[i] = NULL;
			if(i<((fundamental_assignments-1)*FUNDAMENTAL_NUM)-2)
				fundamentals[i] = malloc((strlen(stored_funds[i])+100)*sizeof(char));
			else
				fundamentals[i] = malloc((tree_length_assignments*TREE_LENGTH)*sizeof(char));
			if(!fundamentals[i]) memory_error(13);
			fundamentals[i][0] = '\0';
			if(i <((fundamental_assignments-1)*FUNDAMENTAL_NUM)-2)
				{
				strcpy(fundamentals[i], stored_funds[i]);  /**** copying the values already read back into the fundamental array **/
				}
			}
		for(i=0; i<((fundamental_assignments-1)*FUNDAMENTAL_NUM)-2; i++)
			free(stored_funds[i]);
		free(stored_funds);
		stored_funds = NULL;
				/***** now do the same for presence_of_taxa *********/
		/*** store the present settings for presence_of_taxa ******/

		stored_presence_of_taxa = malloc(((fundamental_assignments-1)*FUNDAMENTAL_NUM)*sizeof(int *));
		for(i=0; i<((fundamental_assignments-1)*FUNDAMENTAL_NUM)-2; i++)
			{
			stored_presence_of_taxa[i] = NULL;
			stored_presence_of_taxa[i] = malloc((TAXA_NUM*name_assignments)*sizeof(int));
			for(j=0; j<TAXA_NUM*name_assignments; j++)
				stored_presence_of_taxa[i][j] = presence_of_taxa[i][j];
			}

		/****** free up presence of taxa ******/
		for(i=0; i<((fundamental_assignments-1)*FUNDAMENTAL_NUM); i++)
			free(presence_of_taxa[i]);
		free(presence_of_taxa);
		presence_of_taxa = NULL;
		/*** create the new sized array for presence_of_taxa ******/
		presence_of_taxa = malloc((fundamental_assignments*FUNDAMENTAL_NUM)*sizeof(int *));
	   if(!presence_of_taxa) memory_error(14);
		for(i=0; i<fundamental_assignments*FUNDAMENTAL_NUM; i++)
			{
			presence_of_taxa[i] = malloc((TAXA_NUM*name_assignments)*sizeof(int));
			if(!presence_of_taxa[i]) memory_error(15);

			for(j=0; j<TAXA_NUM*name_assignments; j++)
				{
				if(i<((fundamental_assignments-1)*FUNDAMENTAL_NUM)-2)
					presence_of_taxa[i][j] = stored_presence_of_taxa[i][j];
				else
					presence_of_taxa[i][j] = FALSE;
				}

			}
		for(i=0; i<((fundamental_assignments-1)*FUNDAMENTAL_NUM)-2; i++)
			{
			free(stored_presence_of_taxa[i]);
			}
		free(stored_presence_of_taxa);
		stored_presence_of_taxa = NULL;
				}
	}


int comment(FILE *file)
	{
	char c;
	c = getc(file);
	while(c != ']' && !feof(file))
		{
		if(c == '[') comment(file);
		c = getc(file);
		}
	if(!feof(file)) return(FALSE);
	else return(TRUE);
	}



int nexusparser(FILE *nexusfile)
	{
	char c, begin[6] = {'b', 'e', 'g', 'i', 'n', '\0'}, translate[10] = {'t','r','a','n','s','l','a','t','e','\0'}, tree[5] = {'t','r','e','e','\0'};
	int error = FALSE, i, j, k, l, translated = FALSE, num_taxa = 1000, num_trees = 0, found = FALSE, numtranslatedtaxa = 0;
	char *string = NULL, ***names = NULL, *newtree, single[1000];

	newtree = malloc(TREE_LENGTH*sizeof(char));
	newtree[0] = '\0';
	string = malloc(TREE_LENGTH*sizeof(char));
	string[0] = '\0';
	while(((c = getc(nexusfile)) == ' ' || c == '\t' || c == '\n' || c == '\r') && !feof(nexusfile));

	while(!feof(nexusfile) && !error)
		{
		while(((c = getc(nexusfile)) == ' ' || c == '\t' || c == '\n' || c == '\r') && !feof(nexusfile));
		switch(tolower(c))
			{
			case '[':

				comment(nexusfile);
				break;

			case 'b':   /* this is the beginning of a block **/
				for(i=1; i<5; i++)
					if((c = tolower(getc(nexusfile))) != begin[i]) error = TRUE;
				if(!error)
					{
					while(((c = getc(nexusfile)) == ' ' || c == '\t' || c == '\n' || c == '\r') && !feof(nexusfile));
					i = 0;
					do {	/* find the type of block it is **/
						string[i] = tolower(c);
						i++;
						c = getc(nexusfile);
						}while(c != ';' && c != ' ');
					string[i] = '\0';
					if(strcmp(string, "trees") == 0)
						{
						while(strcmp(string, "end;") != 0 && strcmp(string, "endblock;") != 0 && !feof(nexusfile) && !error)
							{
							/*** read in the trees ***/
							while(((c = getc(nexusfile)) == ' ' || c == '\t' || c == '\n' || c == '\r') && !feof(nexusfile));
							while(c == '[')
								{
								comment(nexusfile);
								while(((c = getc(nexusfile)) == ' ' || c == '\t' || c == '\n' || c == '\r') && !feof(nexusfile));
								}
							strcpy(string, "");
							i=1;
							string[0] = tolower(c);
							while(((c = getc(nexusfile)) != ' ' && c != '\n' && c != '\r') && !feof(nexusfile))
								{
								string[i] = tolower(c);
								i++;
								}
							string[i] = '\0';
							if(strcmp(string, tree) == 0)
								{
								/** read in trees **/
								/** read in the name of the tree **/
								while(((c = getc(nexusfile)) == ' ' || c == '\t' || c == '\n' || c == '\r') && !feof(nexusfile));
								i=1;
								string[0] = c;
								while(((c = getc(nexusfile)) != ' ' && c != '\n' && c != '\r') && !feof(nexusfile))
									{
									string[i] = c;
									i++;
									}
								string[i] = '\0';
								strcpy(tree_names[Total_fund_trees], string);

								/*** look out for any weights **/
								while(((c = getc(nexusfile)) == ' ' || c == '\t' || c == '=' || c == '\n' || c == '\r') && !feof(nexusfile));
								while(c == '[')
									{
									if((c = getc(nexusfile)) == '&')
										{
										if((c = tolower(getc(nexusfile))) == 'w')
											{
											while(((c = getc(nexusfile)) == ' ' || c == '\t' || c == '\n' || c == '\r') && !feof(nexusfile));
											i=1;
											string[0] = c;
											while((c = getc(nexusfile)) != ']' && c != '/' &&!feof(nexusfile))
												{
												string[i] = c;
												if(string[i] != ' ') i++;
												}
											string[i] = '\0';
											if(c == '/') /* if PAUP has saved tree weights as a fraction (i.e. 1/2 or 4/5 etc) */
												{
												tree_weights[Total_fund_trees] = tofloat(string);
												i=0;
												while((c = getc(nexusfile)) != ']' &&!feof(nexusfile))
													{
													string[i] = c;
													if(string[i] != ' ') i++;
													}
												string[i] = '\0';
												tree_weights[Total_fund_trees] = tree_weights[Total_fund_trees]/tofloat(string);
												}

											}
										else
											while((c = getc(nexusfile)) != ']' &&!feof(nexusfile));
										}
									else
										while((c = getc(nexusfile)) != ']' &&!feof(nexusfile));
									while(((c = getc(nexusfile)) == ' ' || c == '\t' || c == '\n' || c == '\r') && !feof(nexusfile));
									}


								/*** read in the tree **/
								if(!translated)
									{
									i=0;
									strcpy(string, "");
									while(c != ';')
										{
										if(string[i] != ' ' && string[i] != '\t' && string[i] != '\n' && string[i] !='\r') string[i] = c;
										i++;
										c = getc(nexusfile);
										}
									string[i] = ';';
									string[i+1] = '\0';
									input_fund_tree(string, Total_fund_trees);
									}
								else
									{
									/*** If there was a translation table then we need to put the names into the tree **/
									i=0;
									while(c != ';' && !feof(nexusfile) && !error)
										{
										switch(c)
											{
											case '(':
											case ',':
												newtree[i] = c;
												i++;
												c = getc(nexusfile);
												break;
											case ')':
											case ':':
												do
													{
													newtree[i] = c;
													i++;
													c = getc(nexusfile);
													}while(c != '(' && c != ')' && c != ',' && c != ';');
												break;
											default:
												j=0;
												while(c != '(' && c != ')' && c != ',' && c != ';' && c != ':' && !feof(nexusfile))
													{
													single[j] = c;
													j++;
													c = getc(nexusfile);
													}
												single[j] = '\0';
												found = -1;
												for(j=0; j<numtranslatedtaxa; j++)
													{
													if(strcmp(single, names[j][0])== 0)
														found = j;
													}
												if(found == -1)
													{
													printf2("Error: taxa %s is not in the translation table\n", single);
													error = TRUE;
													}
												else
													{
													j=0;
													while(names[found][1][j] != '\0')
														{
														newtree[i] = names[found][1][j];
														j++; i++;
														}
													}
										}
										}
									newtree[i] = ';';
									input_fund_tree(newtree, Total_fund_trees);
									}
								}
							if(strcmp(string, translate) == 0)
								{
								translated = TRUE;
								/** handle translation **/
								names = malloc(num_taxa*sizeof(char**));
								for(i=0; i<num_taxa; i++)
									{
									names[i] = malloc(2*sizeof(char*));
									for(j=0; j<2; j++)
										{
										names[i][j] = malloc(1000*sizeof(char));
										names[i][j][0] = '\0';
										}
									}
								i=0; j=0; c = getc(nexusfile);
								while(c != ';' && !feof(nexusfile))
									{
									while(c == ' ' || c == '\n' || c == '\r' || c == '\t' || c == ',') c = getc(nexusfile);
									if(j==2)
										{
										j = 0;
										i++;
										if(i+1 >= num_taxa)
											{
											names = realloc(names,(num_taxa+num_taxa)*sizeof(char**));
											for(k=num_taxa; k<num_taxa+num_taxa; k++)
												{
												names[k] = malloc(2*sizeof(char*));
												for(l=0; l<2; l++)
													{
													names[k][l] = malloc(1000*sizeof(char));
													names[k][l][0] = '\0';
													}
												}

											num_taxa += num_taxa;
											}
										}
									k=0;
									while(c != ' ' && c != '\n' && c != '\r' && c != '\t' && c != ';' && c != ',' && !feof(nexusfile))
										{
										names[i][j][k] = c;
										k++;
										c = getc(nexusfile);
										}
									names[i][j][k] = '\0';
									j++;
									}
								numtranslatedtaxa = i+1;
								}
							}
						}
					else
						{
						if(strcmp(string, "clann")==0)
							{
							strcpy(string, "");
							while(strcmp(string, "end;") != 0 && strcmp(string, "endblock;") != 0 && !feof(nexusfile) && !error)
								{
								i=0;
								while((c == ' ' || c == '\n' || c == '\r' || c == '\t' || c == ';') && !feof(nexusfile)) c = getc(nexusfile);
								while(c != ';' && !feof(nexusfile))
									{
									string[i] = c;
									i++;
									c = getc(nexusfile);
									}
								string[i] = ';';
								string[i+1] = '\0';
								strcpy(stored_commands[parts], string);
								parts++;


								}
							}
						else
							{
							/*** skip to the end of the block -- it must contain sequences or something **/
							strcpy(string, "");
							while(strcmp(string, "end;") != 0 && strcmp(string, "endblock;") != 0 && !feof(nexusfile) && !error)
								{
								i=0;
								while((c == ' ' || c == '\n' || c == '\r' || c == '\t' || c == ';') && !feof(nexusfile)) c = getc(nexusfile);
								while(c != ';' && c != ' ' && c != '\n' && c != '\r' && c != '\t' && !feof(nexusfile))
									{
									string[i] = tolower(c);
									if(i < 9997)i++;
									c = tolower(getc(nexusfile));
									}
								string[i] = ';';
								string[i+1] = '\0';
								}
							}
						}
					}
				else
					{
					/*** what to do if it doesn't say begin? ***/
					error = TRUE;
					}
				break;

			default:
				break;


			}
		}
	free(newtree);
/*	fclose(nexusfile); */
	free(string);
	if(translated)
		{
		for(i=0; i<num_taxa; i++)
			{
			for(j=0; j<2; j++)
				free(names[i][j]);
			free(names[i]);
			}
		free(names);
		}
	return(error);
	}


void input_file_summary(int do_all)
    {
    int i=0, j=0, k=0, l=0, limit, previous_limit, limit_reached, singlecopy=0;
    float *size_of_trees = NULL;

	size_of_trees = malloc((Total_fund_trees-num_excluded_trees)*sizeof(float));
	for(i=0; i<Total_fund_trees-num_excluded_trees; i++)
		size_of_trees[i] = 0;

    sup=1;
    printf2("\n\tSource tree summmary:\n\n");
    printf2("\t----------------------------------------------------\n");
    printf2("\t\tNumber of input trees: %d\n", Total_fund_trees-num_excluded_trees);
    printf2("\t\tNumber of unique taxa: %d\n", number_of_taxa-num_excluded_taxa);
    sup =1;
    for(i=4; i<=number_of_taxa-num_excluded_taxa; i++) sup*=((2*i)-5);
    printf2("\t\tTotal unrooted trees in Supertree space?\n\t\t\t%g\n\n", sup);


    printf2("\tOccurrence summary:\n\n");

    printf2("\t\tnumber\tTaxa name           \t\tOccurrence\n\n");
    for(i=0; i<number_of_taxa; i++)
        {
        printf2("\t\t%-4d\t%-30s\t%d\n", i, taxa_names[i], taxa_incidence[i]);
        }

	if(do_all)
		{
		printf2("\n\n\tCo-occurrence summary:\n\n");
		printf2("\t\tTaxa Number\n      ");


		limit_reached = FALSE;
		previous_limit = 0;
		limit = 12;
		while(!limit_reached)
			{
			if(limit >= number_of_taxa)
				{
				limit_reached = TRUE;
				limit = number_of_taxa;
				}
			if(previous_limit == 0)
				printf2("\n\t\t      ");
			else
				printf2("\n\tCo-occurence summary continued:\n\t\t      ");
			for(i=previous_limit; i<limit; i++)
				printf2("%-5d ", i);
			printf2("\n");
			for(i=previous_limit; i<number_of_taxa; i++)
				{
				printf2("\t\t%-5d ", i);
				for(j=previous_limit; j<=i; j++)
					{
					if(j<limit)
						{
						if(j== i)
							printf2("-     ");
						else
							printf2("%-5d ", Cooccurrance[i][j]);
						}
					}
				printf2("\n");
				}
			previous_limit = limit;
			limit = limit+12;
			if(limit > number_of_taxa) limit = number_of_taxa;
			}
		}

	printf2("\n\n\tSource tree size summary:\n\n  num leaves\n");
	l=0;
	smallest_tree = -1;
	largest_tree = -1;
	for(i=0; i<Total_fund_trees; i++)
		{
		if(sourcetreetag[i])
			{
			k=0;
			for(j=0; j<number_of_taxa; j++)
				{
				k+=presence_of_taxa[i][j];
				}
			if(k < smallest_tree || smallest_tree == -1) smallest_tree = k;
			if(k > largest_tree || largest_tree == -1) largest_tree = k;
			size_of_trees[l] = k;
			l++;
			}
		}

	draw_histogram(NULL, (largest_tree-smallest_tree)+1, size_of_trees, Total_fund_trees-num_excluded_trees);


	for(l=0; l<Total_fund_trees; l++)
		{
		i=0;
		for(k=0; k<number_of_taxa; k++)
			{
			/* count how many taxa are present more than once in each tree */
			if(presence_of_taxa[l][k] > 1) i++;
			}
		if(i==0) singlecopy++;
		}

	printf2("\tNumber of single copy trees:\t%d\n", singlecopy);
	printf2("\tnumber of multicopy trees:\t%d\n", Total_fund_trees-singlecopy);

    printf2("\t----------------------------------------------------\n");

	free(size_of_trees);
    }



int assign_taxa_name(char *inputname,int fund)
	{
	int i =0, j=0, k=0, answer = -1, taxa_on_tree = 0;
	char *delim = NULL, name[NAME_LENGTH];

	delim = malloc(10*sizeof(char));
	delim[0] = delimiter_char;
	delim[1] = '\0';
	name[0] = '\0';

        /* Reallocate memory if the number of taxa exceeds that originally defined */
        if(number_of_taxa+1 == TAXA_NUM*name_assignments)
            {
            name_assignments++;
            taxa_names = realloc(taxa_names, (TAXA_NUM*name_assignments)*sizeof(char *));
            if(!taxa_names) memory_error(16);

            for(j=TAXA_NUM*(name_assignments-1); j<TAXA_NUM*name_assignments; j++)
                    {
                    taxa_names[j] = malloc(NAME_LENGTH*sizeof(char));
                    if(!taxa_names[j]) memory_error(17);

                    taxa_names[j][0] = '\0';
                    }



            for(j=0; j<fundamental_assignments*FUNDAMENTAL_NUM; j++)
                {
                presence_of_taxa[j] = realloc(presence_of_taxa[j], (TAXA_NUM*name_assignments)*sizeof(int));
                if(!presence_of_taxa[j]) memory_error(19);

                for(k=(TAXA_NUM*(name_assignments-1)); k<(TAXA_NUM*name_assignments); k++)
                    presence_of_taxa[j][k] = FALSE;
                }

            Cooccurrance = realloc(Cooccurrance, (TAXA_NUM*name_assignments)*sizeof(int *));
            if(!Cooccurrance) memory_error(20);

            for(j=0; j<TAXA_NUM*(name_assignments-1); j++)
                    {
                    Cooccurrance[j] = realloc(Cooccurrance[j], (TAXA_NUM*name_assignments)*sizeof(int));
                    if(!Cooccurrance[j]) memory_error(21);

                    for(k=TAXA_NUM*(name_assignments-1); k<TAXA_NUM*name_assignments; k++)
                        Cooccurrance[j][k] = 0;
                    }

            for(j=TAXA_NUM*(name_assignments-1); j<TAXA_NUM*name_assignments; j++)
                    {
                    Cooccurrance[j] = malloc((TAXA_NUM*name_assignments)*sizeof(int));
                    if(!Cooccurrance[j]) memory_error(22);

                    for(k=0; k<TAXA_NUM*name_assignments; k++)
                        Cooccurrance[j][k] = 0;
                    }

            for(j=TAXA_NUM*(name_assignments-1); j<TAXA_NUM*name_assignments; j++)

            taxa_incidence = realloc(taxa_incidence, (TAXA_NUM*name_assignments)*sizeof(int));
            if(!taxa_incidence) memory_error(23);

            for(j=TAXA_NUM*(name_assignments-1); j<TAXA_NUM*name_assignments; j++)
                {
                taxa_incidence[j] = 0;
                }

            same_tree = realloc(same_tree, (TAXA_NUM*name_assignments)*sizeof(int));
            if(!same_tree)  memory_error(24);

            for(j=TAXA_NUM*(name_assignments-1); j<TAXA_NUM*name_assignments; j++)
                {
                same_tree[j] = 0;
                }
            }




        i=0;
	if(delimiter)
		{
		while(inputname[i] != delimiter_char && inputname[i] != '\0')
			{
			name[i]=inputname[i];
			i++;
			}
		name[i] = '\0';
		}
	else
		{
		strcpy(name, inputname);
		}
	i=0;
	while(answer == -1)
            {

            if(i != number_of_taxa)
                {
                if(strncmp(name, taxa_names[i], max_name_length) == 0)
                    {
                    answer = i;
                    taxa_incidence[i]++;
                    same_tree[i]++;
                    taxa_on_tree++;
                    if(fund == TRUE)
                        {
                        presence_of_taxa[Total_fund_trees][i]++;
                        }
                    }
                }
            else
                {
                strncpy(taxa_names[i], name, max_name_length);
                answer = i;
                taxa_incidence[i] = 1;
                same_tree[i]++;
                taxa_on_tree++;
                number_of_taxa++;
                if(fund == TRUE)
                    {
                    presence_of_taxa[Total_fund_trees][i]++;
                    }
                }
            i++;
            }
	free(delim);
	return(answer);
	}

void showtrees(int savet)
	{
	int worst = -2, best = -2,savetrees = FALSE, found = TRUE, taxachosen = 0, counter = 0, mode[5] = {TRUE, FALSE, FALSE, FALSE, FALSE}, start = 0, end = Total_fund_trees, error = FALSE, i=0, j=0, k=0, l=0, num=0, equalto = -1, greaterthan =0, lessthan = 1000000000, taxa_count = 0;
	char *temptree, string_num[10], namecontains[NAME_LENGTH], **containstaxa = NULL, savedfile[100];
	char *temptree1 = malloc(TREE_LENGTH * sizeof(char));
	char *tmp = malloc(TREE_LENGTH * sizeof(char));
	if(!temptree1 || !tmp) { free(temptree1); free(tmp); printf2("Error: out of memory in showtrees\n"); return; }
	FILE *showfile = NULL;
	float bestscore =10000000, worstscore = 0, **tempscores = NULL;
	int *tempsourcetreetag = NULL, display = TRUE, best_total = -1, total = 0, display_fullnames = FALSE, taxaorder=0;
	struct taxon *position = NULL, *species_tree = NULL, *gene_tree = NULL, *best_mapping = NULL, *unknown_fund = NULL, *pos = NULL,*copy = NULL;

	tempscores = malloc(Total_fund_trees*sizeof(float *));
	for(i=0; i<Total_fund_trees; i++)
		{
		tempscores[i] = malloc(2*sizeof(float));
		tempscores[i][0] = 0;
		tempscores[i][1] = i;
		}

	tempsourcetreetag = malloc(Total_fund_trees*sizeof(int));
	for(i=0; i<Total_fund_trees; i++) tempsourcetreetag[i] = sourcetreetag[i];
	savedfile[0] = '\0';
	if(savet == TRUE)
		{
		savetrees=TRUE;
		display=FALSE;
		strcpy(savedfile, "savedtrees.txt");
		}
	else
		strcpy(savedfile, "showtrees.txt");

	containstaxa = malloc(number_of_taxa*sizeof(char*));
	for(i=0; i<number_of_taxa; i++)
		{
		containstaxa[i] = malloc(NAME_LENGTH*sizeof(char));
		containstaxa[i][0] = '\0';
		}
	temptree = malloc(TREE_LENGTH*sizeof(char));
	temptree[0] = '\0';
	temptree1[0] = '\0';
	for(i=0; i<num_commands; i++)
		{
		if(strcmp(parsed_command[i], "fullnames") == 0)
			{
			if(strcmp(parsed_command[i+1], "no") == 0)
				display_fullnames = FALSE;
			else
				{
				if(strcmp(parsed_command[i+1], "yes") == 0)
					{
					display_fullnames = TRUE;
					}
				else
					{
					printf2("ERROR: %s not valid modifier of \"fullnames\"\n", parsed_command[i+1]);
					error = TRUE;
					}
				}
			}
		if(strcmp(parsed_command[i], "display") == 0)
			{
			if(strcmp(parsed_command[i+1], "no") == 0)
				display = FALSE;
			else
				{
				if(strcmp(parsed_command[i+1], "yes") == 0)
					{
					display = TRUE;
					}
				else
					{
					printf2("ERROR: %s not valid modifier of \"display\"\n", parsed_command[i+1]);
					error = TRUE;
					}
				}
			}
		if(savet==FALSE) /* If this has been calledf  from the showtrees command */
			{
			if(strcmp(parsed_command[i], "savetrees") == 0)
				{
				if(strcmp(parsed_command[i+1], "yes") == 0)
					savetrees = TRUE;
				else
					{
					if(strcmp(parsed_command[i+1], "no") == 0)
						savetrees = FALSE;
					else
						{
						printf2("Error: %s not a valid modifier of \"savetrees\"\n", parsed_command[i+1]);
						error = TRUE;
						}
					}
				}
			}
		if(strcmp(parsed_command[i], "filename") == 0)
			strcpy(savedfile, parsed_command[i+1]);

		if(strcmp(parsed_command[i], "range") == 0)
			{
			mode[0] = TRUE;
			start = toint(parsed_command[i+1])-1;
			end = toint(parsed_command[i+2])-1;
			if(start <0 || end > Total_fund_trees)
				{
				error = TRUE;
				printf2("Error: the start of the range must not be less than 1 and the end of the range must no be greater than the number of source trees\n");
				}
			}
		if(strcmp(parsed_command[i], "namecontains") == 0)
			{
			namecontains[0] = '\0';
			strcpy(namecontains, parsed_command[i+1]);
			mode[2] = TRUE;
			}
		if(strcmp(parsed_command[i], "containstaxa") == 0)
			{
			namecontains[0] = '\0';
			if(taxachosen < number_of_taxa)
				{
				strcpy(containstaxa[taxachosen], parsed_command[i+1]);
				taxachosen++;
				}
			else
				{
				printf2("Error: the number of taxa chosen os greater than the number of taxa in memory\n");
				error = TRUE;
				}
			mode[3] = TRUE;
			}

		if(strcmp(parsed_command[i], "score") == 0)
			{
			mode[4] = TRUE;
			bestscore = atof(parsed_command[i+1]);
			worstscore = atof(parsed_command[i+2]);
			}

		if(strcmp(parsed_command[i], "size") == 0)
			{
			mode[1] = TRUE;
			if(strcmp(parsed_command[i+1], "equalto") == 0)
				{
				equalto = toint(parsed_command[i+2]);

				if(equalto < 4 || equalto > number_of_taxa)
					{
					error = TRUE;
					printf2("Error in size \"equalto\"\n\n");
					}
				}
			else
				{
				if(strcmp(parsed_command[i+1], "greaterthan") == 0)
					{
					greaterthan = toint(parsed_command[i+2]);

					if(greaterthan >= number_of_taxa)
						{
						error = TRUE;
						printf2("Error in size \"greaterthan\"\n\n");
						}
					}
				else
					{
					if(strcmp(parsed_command[i+1], "lessthan") == 0)
						{
						lessthan = toint(parsed_command[i+2]);

						if(lessthan < 5)
							{
							error = TRUE;
							printf2("Error in size \"lessthan\"\n\n");
							}
						}
					else
						{
						printf2("Error: %s not valid option for \"size\"\n", parsed_command[i+1]);
						error = TRUE;
						}
					}
				}
			}



		}
	if(mode[4])
		{
		if(trees_in_memory == 0)
			{
			printf2("Error: There are no saved supertrees in memory from which to calculate scores\n");
			error = TRUE;
			}
		}

	if(savetrees)
		{
		if((showfile = fopen(savedfile, "w")) == NULL)		/* check to see if the file is there */
			{								/* Open the source tree file */
			printf2("Cannot open file %s\n", savedfile);
			error = TRUE;
			}
		}
	if(!error)
		{
		if(savet == TRUE)
			printf2("savetrees settings:\n\n");
		else
			printf2("showtrees settings:\n\n");
		if(savetrees)
			printf2("\tsaving selection of trees in phylip format to file: %s\n", savedfile);



		if(mode[0])  /* Specifies a particular range of trees - usually true by default */
			{
			for(i=0; i<start; i++)
				tempsourcetreetag[i] = FALSE;
			for(i=end+1; i<Total_fund_trees; i++)
				tempsourcetreetag[i] = FALSE;
			}

		if(mode[1]) /* Looks for trees of a particular size */
			{
			for(i=0; i<Total_fund_trees; i++)
				{
				strcpy(temptree, fundamentals[i]);
				returntree(temptree);
				/* build the tree in memory */
				/****** We now need to build the Supertree in memory *******/
				if(tree_top != NULL)
					{
					dismantle_tree(tree_top);
					tree_top = NULL;
					}
				temp_top = NULL;
				tree_build(1, temptree, tree_top, 1, -1, 0);

				tree_top = temp_top;
				temp_top = NULL;

				/* check how many taxa are on the tree */
				taxa_count = count_taxa(tree_top, 0);

				if(equalto != -1)
					{
					if(taxa_count != equalto)
						tempsourcetreetag[i] = FALSE;
					}
				else
					{
					if(taxa_count >= lessthan || taxa_count <= greaterthan)
						tempsourcetreetag[i] = FALSE;
					}
				}
			}
		if(mode[2])
			{
			for(i=0; i<Total_fund_trees; i++)
				{
				if(strstr(tree_names[i], namecontains) == NULL)
					{
					tempsourcetreetag[i] = FALSE;
					}
				}
			}
		if(mode[3]) /* Look for trees that contain a particular taxon */
			{
			for(i=0; i<Total_fund_trees; i++)
				{
				strcpy(temptree, fundamentals[i]);
				returntree(temptree);
				/* build the tree in memory */
				/****** We now need to build the Supertree in memory *******/
				if(tree_top != NULL)
					{
					dismantle_tree(tree_top);
					tree_top = NULL;
					}
				temp_top = NULL;
				tree_build(1, temptree, tree_top, 1, -1, 0);

				tree_top = temp_top;
				temp_top = NULL;


				found = TRUE;
				for(j=0; j<taxachosen; j++)
					{
					if(found == TRUE)
						found = find_taxa(tree_top, containstaxa[j]);
					}
				if(!found)
					tempsourcetreetag[i] = FALSE;
				}
			}

		if(mode[4])
			{
			for(i=0; i<Total_fund_trees; i++)
				{
				if(bestscore != worstscore)
					{
					if(sourcetree_scores[i] < bestscore || sourcetree_scores[i] > worstscore)
						tempsourcetreetag[i] = FALSE;
					}
				else
					{
					if(sourcetree_scores[i] < bestscore-0.000001 || sourcetree_scores[i] > worstscore+0.000001)
						tempsourcetreetag[i] = FALSE;
					}
				}
			}

		for(j=0; j<Total_fund_trees; j++)
	        {
	        if(tempsourcetreetag[j] && sourcetreetag[j])
				{
				if(savetrees)
					{
			        if(tree_top != NULL) dismantle_tree(tree_top);  /* Dismantle any trees already in memory */
			        tree_top = NULL;

			        temp_top = NULL;
			        taxaorder=0;
			        temptree[0] = '\0';
			        strcpy(temptree, fundamentals[j]);
					returntree_fullnames(temptree, j);
			        basic_tree_build(1, temptree, tree_top, TRUE);

			        tree_top = temp_top;
			        temp_top = NULL;
			        reset_tree(tree_top);

			        shrink_tree(tree_top);    /* Shrink the pruned tree by switching off any internal nodes that are not needed */
			        temptree[0] = '\0'; /* initialise the string */
			        if(print_pruned_tree(tree_top, 0, temptree, TRUE, j) >1)
			            {
			            tmp[0] = '\0';
			            strcpy(tmp, "(");
			            strcat(tmp, temptree);
			            strcat(tmp, ")");
			            strcpy(temptree, tmp);
			            }
			        strcat(temptree, ";");

			        fprintf(showfile, "%s[%f]", temptree, tree_weights[j]);
			        if(strcmp(tree_names[j], "") != 0)
			        	fprintf(showfile, "[%s]\n", tree_names[j]);
			        else
			        	fprintf(showfile, "[%d]\n", j+1);
					}
				if(display)
					{
					temptree[0] = '\0';
					strcpy(temptree, fundamentals[j]);
					returntree_fullnames(temptree, j);

					/*returntree(temptree); */
					printf2("\n\n\nTree number %d",j+1);
					if(strcmp(tree_names[j], "")!=0) printf2("\nTree name = %s", tree_names[j]);
					printf2("\nWeight = %f\n", tree_weights[j]);
					if(trees_in_memory > 0)printf2("Score = %f\n", sourcetree_scores[j]);
					tree_coordinates(temptree, TRUE, TRUE, FALSE, -1);
					}
				counter++;
		        }
			}


/*
		for(i=0; i<Total_fund_trees; i++)
			{
			if(tempsourcetreetag[i] && sourcetreetag[i])
				{
				strcpy(temptree, fundamentals[i]);
				returntree(temptree);
				if(savetrees)
					{
					fprintf(showfile, "%s ", temptree);
					if(strcmp(tree_names[i], "") != 0)
						fprintf(showfile, "[%s] ", tree_names[i]);
					if(trees_in_memory > 0)
						fprintf(showfile, " [%f]\n", sourcetree_scores[i]);
					else
						fprintf(showfile, "\n");
					}
				if(display)
					{
					printf2("\n\n\nTree number %d\nTree name = %s\n", i+1, tree_names[i]);
					printf2("Weight = %f\n", tree_weights[i]);
					if(trees_in_memory > 0)printf2("Score = %f\n", sourcetree_scores[i]);
					tree_coordinates(temptree, TRUE, TRUE, FALSE, -1);
					}
				counter++;
				}
			} */ /* Old save trees before names were delimited */
		if(savet == TRUE)
			printf2("\n%d source trees met with the criteria specified and saved to file\n", counter);
		else
			printf2("\n%d source trees met with the criteria specified\n", counter);
		}
	if(savetrees) fclose(showfile);
	free(temptree);
	free(temptree1);
	free(tmp);
	for(i=0; i<Total_fund_trees; i++)
		free(tempscores[i]);
	for(i=0; i<number_of_taxa; i++)
		free(containstaxa[i]);
	free(tempscores);
	free(containstaxa);
	free(tempsourcetreetag);
	}

/* quick sort */
static void quick(float ** items, int count)
	{
	qs(items, 0, count-1);
	}

/* the quick sort */
static void qs(float **items, int left, int right)
	{
	int i, j;
	float x, y, x1, y1;

	i= left; j=right;
	x = items[(left+right)/2][0];

	do {
		while((items[i][0] < x) && (i < right)) i++;
		while((x < items[j][0]) && (j > left)) j--;

		if(i <= j)
			{
			y = items[i][0];
			y1 = items[i][1];
			items[i][0] = items[j][0];
			items[i][1] = items[j][1];
			items[j][0] = y;
			items[j][1] = y1;
			i++; j--;
			}
		} while(i <= j);

	if(left < j) qs(items, left, j);
	if(i < right) qs(items, i, right);
	}


void exclude(int do_all)
	{
	int worst = -2, best = -2,savetrees = FALSE, found = TRUE, taxachosen = 0, counter = 0, mode[10] = {FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE}, start = 0, end = Total_fund_trees, error = FALSE, i=0, j=0, k=0, l=0, num=0, equalto = -1, greaterthan =1000000000, lessthan = 0, taxa_count = 0;
	char *temptree, string_num[10], namecontains[100], **containstaxa = NULL, savedfile[100], *command = NULL;
	char *tmp = malloc(TREE_LENGTH * sizeof(char));
	if(!tmp) { printf2("Error: out of memory in exclude\n"); return; }
	FILE *showfile = NULL, *tempfile = NULL;
	float bestscore =10000000, worstscore = 0, **tempscores = NULL;
	int *tempsourcetreetag = NULL, countedout =0, *temp_incidence = NULL, taxaorder=0;



	temp_incidence = malloc(number_of_taxa*sizeof(int));
	for(i=0; i<number_of_taxa; i++) temp_incidence[i] = 0;

	tempscores = malloc(Total_fund_trees*sizeof(float *));
	for(i=0; i<Total_fund_trees; i++)
		{
		tempscores[i] = malloc(2*sizeof(float));
		tempscores[i][0] = 0;
		tempscores[i][1] = i;
		}
	tempsourcetreetag = malloc(Total_fund_trees*sizeof(int));
	for(i=0; i<Total_fund_trees; i++) tempsourcetreetag[i] = sourcetreetag[i];
	for(i=0; i<Total_fund_trees; i++)
		if(tempsourcetreetag[i]) j++;
	containstaxa = malloc(number_of_taxa*sizeof(char*));
	for(i=0; i<number_of_taxa; i++)
		{
		containstaxa[i] = malloc(1000*sizeof(char));
		containstaxa[i][0] = '\0';
		}
	temptree = malloc(TREE_LENGTH*sizeof(char));
	temptree[0] = '\0';
	for(i=0; i<num_commands; i++)
		{
		if(strcmp(parsed_command[i], "singlecopy") == 0)
			{
			mode[5] = TRUE;
			}
		if(strcmp(parsed_command[i], "multicopy") == 0)
			{
			mode[6] = TRUE;
			}
		if(strcmp(parsed_command[i], "range") == 0)
			{
			mode[0] = TRUE;
			start = toint(parsed_command[i+1])-1;
			end = toint(parsed_command[i+2])-1;
			if(start <0 || end > Total_fund_trees)
				{
				error = TRUE;
				printf2("Error: the start of the range must not be less than 1 and the end of the range must no be greater than the number of source trees\n");
				}
			}
		if(strcmp(parsed_command[i], "namecontains") == 0)
			{
			namecontains[0] = '\0';
			strcpy(namecontains, parsed_command[i+1]);
			mode[2] = TRUE;
			}
		if(strcmp(parsed_command[i], "containstaxa") == 0)
			{
			namecontains[0] = '\0';
			if(taxachosen < number_of_taxa)
				{
				strcpy(containstaxa[taxachosen], parsed_command[i+1]);
				taxachosen++;
				}
			else
				{
				printf2("Error: the number of taxa chosen os greater than the number of taxa in memory\n");
				error = TRUE;
				}
			mode[3] = TRUE;
			}
		if(strcmp(parsed_command[i], "score") == 0)
			{
			mode[4] = TRUE;
			bestscore = atof(parsed_command[i+1]);
			worstscore = atof(parsed_command[i+2]);
			if(worstscore < bestscore)
				{
				error = TRUE;
				printf2("Error: the range must end with a larger or equal score to the start of the range\n");
				}
			}
		if(strcmp(parsed_command[i], "size") == 0)
			{
			mode[1] = TRUE;
			if(strcmp(parsed_command[i+1], "equalto") == 0)
				{
				equalto = toint(parsed_command[i+2]);

				if(equalto > number_of_taxa)
					{
					error = TRUE;
					printf2("Error in size \"equalto\"\n\n");
					}
				}
			else
				{
				if(strcmp(parsed_command[i+1], "greaterthan") == 0)
					{
					greaterthan = toint(parsed_command[i+2]);

					}
				else
					{
					if(strcmp(parsed_command[i+1], "lessthan") == 0)
						{
						lessthan = toint(parsed_command[i+2]);


						}
					else
						{
						printf2("Error: %s not valid option for \"size\"\n", parsed_command[i+1]);
						error = TRUE;
						}
					}
				}
			}



		}
	if(mode[4]) /* Look for trees that meet a particular score */
		{
		if(trees_in_memory == 0)
			{
			printf2("Error: There are no saved supertrees in memory from which to calculate scores\n");
			error = TRUE;
			}
		}
	if(mode[4])
		{
		if(trees_in_memory == 0)
			{
			printf2("Error: There are no saved supertrees in memory from which to calculate scores\n");
			error = TRUE;
			}
		}




	if(!error)
		{
		if(mode[0])
			{
			for(i=start-1; i<=end; i++)
				{
				tempsourcetreetag[i] = FALSE;
				countedout++;
				}
			}

		if(mode[1])
			{
			for(i=0; i<Total_fund_trees; i++)
				{
				strcpy(temptree, fundamentals[i]);
				returntree(temptree);
				/* build the tree in memory */
				/****** We now need to build the Supertree in memory *******/
				if(tree_top != NULL)
					{
					dismantle_tree(tree_top);
					tree_top = NULL;
					}
				temp_top = NULL;
				tree_build(1, temptree, tree_top, 1, -1, 0);

				tree_top = temp_top;
				temp_top = NULL;

				/* check how many taxa are on the tree */
				taxa_count = count_taxa(tree_top, 0);

				if(equalto != -1)
					{
					if(taxa_count == equalto)
						{
						tempsourcetreetag[i] = FALSE;
						countedout++;
						}
					}
				else
					{
					if(taxa_count < lessthan || taxa_count > greaterthan)
						{
						tempsourcetreetag[i] = FALSE;
						countedout++;
						}
					}
				}
			}
		if(mode[2])
			{
			for(i=0; i<Total_fund_trees; i++)
				{
				if(strstr(tree_names[i], namecontains) != NULL)
					{
					tempsourcetreetag[i] = FALSE;
					countedout++;
					}
				}
			}
		if(mode[3])
			{
			for(i=0; i<Total_fund_trees; i++)
				{
				strcpy(temptree, fundamentals[i]);
				returntree(temptree);
				/* build the tree in memory */
				/****** We now need to build the Supertree in memory *******/
				if(tree_top != NULL)
					{
					dismantle_tree(tree_top);
					tree_top = NULL;
					}
				temp_top = NULL;
				tree_build(1, temptree, tree_top, 1, -1, 0);

				tree_top = temp_top;
				temp_top = NULL;


				found = TRUE;
				for(j=0; j<taxachosen; j++)
					{
					if(found == TRUE)
						found = find_taxa(tree_top, containstaxa[j]);
					}
				if(found)
					{
					tempsourcetreetag[i] = FALSE;
					countedout++;
					}
				}
			}

		if(mode[4]) /* Look for trees that meet a particular score */
			{
			for(i=0; i<Total_fund_trees; i++)
				{
				if(bestscore != worstscore)
					{
					if(sourcetree_scores[i] >= bestscore || sourcetree_scores[i] <= worstscore)
						{
						tempsourcetreetag[i] = FALSE;
						countedout++;
						}
					}
				else
					{
					if(sourcetree_scores[i] > bestscore-0.000001 && sourcetree_scores[i] < worstscore+0.000001)
						{
						tempsourcetreetag[i] = FALSE;
						countedout++;
						}
					}
				}
			}
		if(mode[5] || mode[6]) /* Look for trees that are single copy of multicopy */
			{
			for(l=0; l<Total_fund_trees; l++)
				{
				i=0;
				for(k=0; k<number_of_taxa; k++)
					{
					/* count how many taxa are present more than once in each tree */
					if(presence_of_taxa[l][k] > 1) i++;
					}
				if(mode[5] && i==0) tempsourcetreetag[l] = FALSE; /*exlude the singlecopy*/
				if(mode[6] && i!=0) tempsourcetreetag[l] = FALSE; /*exlude the multicopy*/
				}
			}


/*
		if(mode[0] || mode[1] || mode[2] || mode[3] || mode[4] || mode[5] || mode[6])
			{
			for(i=0; i<number_of_taxa; i++) temp_incidence[i] = 0;
			for(i=0; i<Total_fund_trees; i++)
				{
				if(tempsourcetreetag[i])
					{
					strcpy(temptree, fundamentals[i]);
					returntree(temptree);
				*/	/* build the tree in memory */
					/****** We now need to build the Supertree in memory *******/
				/*	if(tree_top != NULL)
						{
						dismantle_tree(tree_top);
						tree_top = NULL;
						}
					temp_top = NULL;
					tree_build(1, temptree, tree_top, 1, -1, 0);

					tree_top = temp_top;
					temp_top = NULL;

					for(j=0; j<number_of_taxa; j++)
							temp_incidence[j] += find_taxa(tree_top, taxa_names[j]);
					}
				}
			l=0;
			for(i=0; i<number_of_taxa; i++)
				{
				if(temp_incidence[i] == 0) l++;
				}
				*/
			l=1;
			if(l> 0)
				{
				/*printf("\nWarning: %d Taxa are no longer represented in the included source trees\nThese taxa are as follows:\n", l);
				for(i=0; i<number_of_taxa; i++)
					{
					if(temp_incidence[i] == 0)
						printf2("\t%s\n", taxa_names[i]);
					}
				command = malloc(10000*sizeof(char));
				printf2("This will permenantly remove these trees from memory\nAre you sure you wish to continue: (yes/no) ");
				xgets(command);
				if(strcmp(command, "y") == 0 || strcmp(command, "Y") == 0 || strcmp(command, "yes") == 0 || strcmp(command, "Yes") == 0)
					{
				*/	tempfile = fopen("tempclannfile164.chr", "w");
					tmp[0] = '\0';
					countedout = 0;


					for(j=0; j<Total_fund_trees; j++)
					        {
					        if(tempsourcetreetag[j])
								{
						        if(tree_top != NULL) dismantle_tree(tree_top);  /* Dismantle any trees already in memory */
						        tree_top = NULL;

						        temp_top = NULL;
						        taxaorder=0;
						        strcpy(temptree, fundamentals[j]);
								returntree_fullnames(temptree, j);
			                    basic_tree_build(1, temptree, tree_top, TRUE);

						       /*  tree_build(1, fundamentals[j], tree_top, 0, j, 0); build the tree passed to the function */

						        tree_top = temp_top;
						        temp_top = NULL;
						        reset_tree(tree_top);

						        shrink_tree(tree_top);    /* Shrink the pruned tree by switching off any internal nodes that are not needed */
						        temptree[0] = '\0'; /* initialise the string */
						        if(print_pruned_tree(tree_top, 0, temptree, TRUE, j) >1)
						            {
						            tmp[0] = '\0';
						            strcpy(tmp, "(");
						            strcat(tmp, temptree);
						            strcat(tmp, ")");
						            strcpy(temptree, tmp);
						            }
						        strcat(temptree, ";");

						        if(strcmp(tree_names[j], "") != 0)
						        	fprintf(tempfile, "%s[%s]\n", temptree, tree_names[j]);
						        else
						        	fprintf(tempfile, "%s[%d]\n", temptree, j);
						        countedout++;
						        }
							else
								{
								if(sourcetreetag[i]) counter++;
								}
							}

		/*			for(i=0; i<Total_fund_trees; i++)
						{
						if(tempsourcetreetag[i])
							{
							temptree[0] = '\0';
							strcpy(temptree, fundamentals[i]);
							returntree(temptree);
							strcpy(tmp, "");
							strncpy(tmp, temptree, strlen(temptree)-1);
							tmp[strlen(temptree)-1] = '\0';
							strcpy(temptree, tmp);
							fprintf(tempfile, "%s", temptree);
							fprintf(tempfile, " [%f]; [ %s]\n", tree_weights[i], tree_names[i]);
							countedout++;
							}
						else
							{
							if(sourcetreetag[i]) counter++;
							}
						}  */ /* Old printing code, before implementation of delimited names input */
					fclose(tempfile);
					for(i=0; i<number_of_taxa; i++)
						{
						taxa_incidence[i] = 0;
						for(j=0; j<number_of_taxa; j++)
							{
							Cooccurrance[i][j] = 0;
							}
						}
					execute_command("tempclannfile164.chr", do_all);
					remove("tempclannfile164.chr");
				/*	}
				else
					{
					printf2("\nAction aborted\n");
					error = TRUE;
					} */
				}
			else
				{
				counter = 0;
				for(i=0; i<Total_fund_trees; i++)
					{
					if(tempsourcetreetag[i]== FALSE && sourcetreetag[i] == TRUE)
						{
						sourcetreetag[i] = FALSE;
						counter++;
						}
					}
				countedout = 0;
				for(i=0; i<Total_fund_trees; i++)
					{
					if(sourcetreetag[i])
						countedout++;
					}
				/*num_excluded_trees +=counter; */

				for(i=0; i<number_of_taxa; i++)
					{
					taxa_incidence[i] = 0;
					for(j=0; j<number_of_taxa; j++)
						{
						Cooccurrance[i][j] = 0;
						}
					for(k=0; k<Total_fund_trees; k++)
						{
						if(sourcetreetag[k])
							{
							if(presence_of_taxa[k][i] > 0)
								taxa_incidence[i]++;
							}
						}
					for(j=0; j<number_of_taxa; j++)
						{
						for(k=0; k<Total_fund_trees; k++)
							{
							if(sourcetreetag[k])
								{
								if(presence_of_taxa[k][i] > 0)
									{
									if(presence_of_taxa[k][j] > 0)
										Cooccurrance[i][j]++;
									}
								}
							}
						}
					}
				input_file_summary(do_all);
				}
			if(!error)
				{
				printf2("\n%d source trees were excluded, %d trees remain in memory\n", counter, countedout );
				}
			/*} */
		}

	free(temptree);
	free(tmp);
	for(i=0; i<Total_fund_trees; i++)
		free(tempscores[i]);
	for(i=0; i<number_of_taxa; i++)
		free(containstaxa[i]);
	free(tempscores);
	free(containstaxa);
	free(tempsourcetreetag);
	free(temp_incidence);

	free(command);

	}

void returntree(char *temptree) /* returns the tree with the names of the taxa included */
	{
	char string_num[10];
	char *string = malloc(TREE_LENGTH * sizeof(char));
	int i=0, j=0, k=0, l=0, num;

	if(!string) { printf2("Error: out of memory in returntree\n"); return; }
	string[0] = '\0';
	strcpy(string, temptree);
	strcpy(temptree, "");
	string_num[0] = '\0';
	j=0; k=0;
	while(string[j] != ';')
		{
		switch(string[j])
			{
			case ')':
				temptree[k] = ')';
				j++; k++;
				while(string[j] != '(' && string[j] != ')' && string[j] != ',' && string[j] != ';' && string[j] != ':')
					{
					temptree[k] = string[j];
					j++; k++;
					}
				break;
			case '(':
			case ',':
				temptree[k] = string[j];
				j++; k++;
				break;
			case ':':
				while(string[j] != '(' && string[j] != ')' && string[j] != ',' && string[j] != ';')
					{
					temptree[k] = string[j];
					j++; k++;
					}
				break;
			default:
				for(l=0; l<10; l++) string_num[l] = '\0';
				l=0;
				while(string[j] != '(' && string[j] != ')' && string[j] != ',' && string[j] != ';' && string[j] != ':')
					{
					string_num[l] = string[j];
					l++; j++;
					}
				string_num[l] = '\0';
				num = atoi(string_num);
				l=0;
				while(taxa_names[num][l] != '\0')
					{
					temptree[k] = taxa_names[num][l];
					k++;
					l++;
					}
				break;
			}
		}
	temptree[k] = ';';
	temptree[k+1] = '\0';
	free(string);

	}

void returntree_fullnames(char *temptree, int treenum) /* returns the tree with the names of the taxa included */
	{
	char string_num[10];
	char *string = malloc(TREE_LENGTH * sizeof(char));
	int i=0, j=0, k=0, l=0, num, taxaorder=0;

	if(!string) { printf2("Error: out of memory in returntree_fullnames\n"); return; }
	string[0] = '\0';
	strcpy(string, temptree);
	strcpy(temptree, "");
	string_num[0] = '\0';
	taxaorder=0;
	j=0; k=0;
	while(string[j] != ';')
		{
		switch(string[j])
			{
			case ')':
				temptree[k] = ')';
				j++; k++;
				while(string[j] != '(' && string[j] != ')' && string[j] != ',' && string[j] != ';' && string[j] != ':')
					{
					temptree[k] = string[j];
					j++; k++;
					}
				break;
			case '(':
			case ',':
				temptree[k] = string[j];
				j++; k++;
				break;
			case ':':
				while(string[j] != '(' && string[j] != ')' && string[j] != ',' && string[j] != ';')
					{
					temptree[k] = string[j];
					j++; k++;
					}
				break;
			default:
				for(l=0; l<10; l++) string_num[l] = '\0';
				l=0;
				while(string[j] != '(' && string[j] != ')' && string[j] != ',' && string[j] != ';' && string[j] != ':')
					{
					string_num[l] = string[j];
					l++; j++;
					}
				string_num[l] = '\0';
				num = atoi(string_num);
				l=0;

				while(fulltaxanames[treenum][taxaorder][l] != '\0')
					{
					temptree[k] = fulltaxanames[treenum][taxaorder][l];
					k++;
					l++;
					}
				taxaorder++;
				break;
			}
		}
	temptree[k] = ';';
	temptree[k+1] = '\0';
	free(string);

	}


void include(int do_all)
	{
	int worst = -2, best = -2,savetrees = FALSE, found = TRUE, taxachosen = 0, counter = 0, mode[5] = {FALSE, FALSE, FALSE, FALSE, FALSE}, start = 0, end = Total_fund_trees, error = FALSE, i=0, j=0, k=0, l=0, num=0, equalto = -1, greaterthan =number_of_taxa, lessthan = 3, taxa_count = 0;
	char *temptree, string_num[10], namecontains[100], **containstaxa = NULL, savedfile[100];
	FILE *showfile = NULL;
	float bestscore =10000000, worstscore = 0, **tempscores = NULL;
	int *tempsourcetreetag = NULL, countedout =0;


	tempscores = malloc(Total_fund_trees*sizeof(float *));
	for(i=0; i<Total_fund_trees; i++)
		{
		tempscores[i] = malloc(2*sizeof(float));
		tempscores[i][0] = 0;
		tempscores[i][1] = i;
		}
	tempsourcetreetag = malloc(Total_fund_trees*sizeof(int));
	for(i=0; i<Total_fund_trees; i++) tempsourcetreetag[i] = sourcetreetag[i];

	containstaxa = malloc(number_of_taxa*sizeof(char*));

	for(i=0; i<number_of_taxa; i++)
		{
		containstaxa[i] = malloc(1000*sizeof(char));
		containstaxa[i][0] = '\0';
		}
	temptree = malloc(TREE_LENGTH*sizeof(char));
	temptree[0] = '\0';
	for(i=0; i<num_commands; i++)
		{
		if(strcmp(parsed_command[i], "range") == 0)
			{
			mode[0] = TRUE;
			start = toint(parsed_command[i+1]);
			end = toint(parsed_command[i+2]);
			if(start <0 || end > Total_fund_trees)
				{
				error = TRUE;
				printf2("Error: the start of the range must not be less than 1 and the end of the range must no be greater than the number of source trees\n");
				}
			}
		if(strcmp(parsed_command[i], "namecontains") == 0)
			{
			namecontains[0] = '\0';
			strcpy(namecontains, parsed_command[i+1]);
			mode[2] = TRUE;
			}
		if(strcmp(parsed_command[i], "containstaxa") == 0)
			{
			namecontains[0] = '\0';
			if(taxachosen < number_of_taxa)
				{
				strcpy(containstaxa[taxachosen], parsed_command[i+1]);
				taxachosen++;
				}
			else
				{
				printf2("Error: the number of taxa chosen os greater than the number of taxa in memory\n");
				error = TRUE;
				}
			mode[3] = TRUE;
			}
		if(strcmp(parsed_command[i], "score") == 0)
			{
			mode[4] = TRUE;
			bestscore = atof(parsed_command[i+1]);
			worstscore = atof(parsed_command[i+2]);
			if(worstscore < bestscore)
				{
				error = TRUE;
				printf2("Error: the range must end with a larger score than the start of the range\n");
				}
			}
		if(strcmp(parsed_command[i], "size") == 0)
			{
			mode[1] = TRUE;
			if(strcmp(parsed_command[i+1], "equalto") == 0)
				{
				equalto = toint(parsed_command[i+2]);

				if(equalto < 4 || equalto > number_of_taxa)
					{
					error = TRUE;
					printf2("Error in size \"equalto\"\n\n");
					}
				}
			else
				{
				if(strcmp(parsed_command[i+1], "greaterthan") == 0)
					{
					greaterthan = toint(parsed_command[i+2]);

					if(greaterthan >= number_of_taxa)
						{
						error = TRUE;
						printf2("Error in size \"greaterthan\"\n\n");
						}
					}
				else
					{
					if(strcmp(parsed_command[i+1], "lessthan") == 0)
						{
						lessthan = toint(parsed_command[i+2]);

						if(lessthan < 5)
							{
							error = TRUE;
							printf2("Error in size \"lessthan\"\n\n");
							}
						}
					else
						{
						printf2("Error: %s not valid option for \"size\"\n", parsed_command[i+1]);
						error = TRUE;
						}
					}
				}
			}



		}
	if(mode[4]) /* Look for trees that meet a particular score */
		{
		if(trees_in_memory == 0)
			{
			printf2("Error: There are no saved supertrees in memory from which to calculate scores\n");
			error = TRUE;
			}
		}

	if(!error)
		{
		if(mode[0])
			{
			for(i=start-1; i<end; i++)
				{
				tempsourcetreetag[i] = TRUE;
				countedout++;
				}
			}

		if(mode[1])
			{
			for(i=0; i<Total_fund_trees; i++)
				{
				strcpy(temptree, fundamentals[i]);
				returntree(temptree);
				/* build the tree in memory */
				/****** We now need to build the Supertree in memory *******/
				if(tree_top != NULL)
					{
					dismantle_tree(tree_top);
					tree_top = NULL;
					}
				temp_top = NULL;
				tree_build(1, temptree, tree_top, 1, -1, 0);

				tree_top = temp_top;
				temp_top = NULL;

				/* check how many taxa are on the tree */
				taxa_count = count_taxa(tree_top, 0);

				if(equalto != -1)
					{
					if(taxa_count == equalto)
						{
						tempsourcetreetag[i] = TRUE;
						countedout++;
						}
					}
				else
					{
					if(taxa_count < lessthan || taxa_count > greaterthan)
						{
						tempsourcetreetag[i] = TRUE;
						countedout++;
						}
					}
				}
			}
		if(mode[2]) /* Looks for trees with a particular name, or part of name */
			{
			for(i=0; i<Total_fund_trees; i++)
				{
				if(strstr(tree_names[i], namecontains) != NULL)
					{
					tempsourcetreetag[i] = TRUE;
					countedout++;
					}
				}
			}
		if(mode[3])
			{
			for(i=0; i<Total_fund_trees; i++)
				{
				strcpy(temptree, fundamentals[i]);
				returntree(temptree);
				/* build the tree in memory */
				/****** We now need to build the Supertree in memory *******/
				if(tree_top != NULL)
					{
					dismantle_tree(tree_top);
					tree_top = NULL;
					}
				temp_top = NULL;
				tree_build(1, temptree, tree_top, 1, -1, 0);

				tree_top = temp_top;
				temp_top = NULL;


				found = TRUE;
				for(j=0; j<taxachosen; j++)
					{
					if(found == TRUE)
						found = find_taxa(tree_top, containstaxa[j]);
					}
				if(found)
					{
					tempsourcetreetag[i] = TRUE;
					countedout++;
					}
				}
			}
		if(mode[4]) /* Look for trees that meet a particular score */
			{
			for(i=0; i<Total_fund_trees; i++)
				{
				if(bestscore != worstscore)
					{
					if(sourcetree_scores[i] >= bestscore || sourcetree_scores[i] <= worstscore)
						{
						tempsourcetreetag[i] = TRUE;
						countedout++;
						}
					}
				else
					{
					if(sourcetree_scores[i] > bestscore-0.000001 && sourcetree_scores[i] < worstscore+0.000001)
						{
						tempsourcetreetag[i] = TRUE;
						countedout++;
						}
					}
				}
			}


		for(i=0; i<Total_fund_trees; i++)
			{
			if(tempsourcetreetag[i]== TRUE && sourcetreetag[i] == FALSE)
				{
				sourcetreetag[i] = TRUE;
				counter++;
				}
			}
		countedout = 0;
		for(i=0; i<Total_fund_trees; i++)
			{
			if(sourcetreetag[i])
				countedout++;
			}
		/*num_excluded_trees -=counter; */

		}
	l = 0;
	for(i=0; i<number_of_taxa; i++)
		{
		taxa_incidence[i] = 0;
		for(j=0; j<number_of_taxa; j++)
			{
			Cooccurrance[i][j] = 0;
			}
		for(k=0; k<Total_fund_trees; k++)
			{
			if(sourcetreetag[k])
				{
				if(presence_of_taxa[k][i] > 0)
					taxa_incidence[i]++;
				}
			}
		for(j=0; j<number_of_taxa; j++)
			{
			for(k=0; k<Total_fund_trees; k++)
				{
				if(sourcetreetag[k])
					{
					if(presence_of_taxa[k][i] > 0)
						{
						if(presence_of_taxa[k][j] > 0)
							Cooccurrance[i][j]++;
						}
					}
				}
			}
		}
	if(!error)
		{
		input_file_summary(do_all);
		printf2("\n%d source trees were re-included, %d trees remain in memory\n", counter, countedout );
		}

	free(temptree);
	for(i=0; i<Total_fund_trees; i++)
		free(tempscores[i]);
	for(i=0; i<number_of_taxa; i++)
		free(containstaxa[i]);
	free(tempscores);
	free(containstaxa);
	free(tempsourcetreetag);
	}

void exclude_taxa(int do_all)
	{
	char *temptree, *pruned_tree = NULL, *tmp = malloc(TREE_LENGTH * sizeof(char)), *command = NULL, tmpfilename[10000], previnputfilename[10000];
	if(!tmp) { printf2("Error: out of memory in exclude_taxa\n"); return; }
	int i=0, j=0, q=0, error = FALSE, taxachosen = 0, found = FALSE, *tobeexcluded = NULL, k=0, l=0, done = FALSE, num_left = 0, min_taxa = 4;
	FILE *tempfile = NULL;

	tmp[0] = '\0';
	tmpfilename[0] = '\0';
	previnputfilename[0] = '\0';
	pruned_tree = malloc(TREE_LENGTH*sizeof(char));

	tobeexcluded = malloc(number_of_taxa*sizeof(int));

	temptree = malloc(TREE_LENGTH*sizeof(char));
	temptree[0] = '\0';

	for(i=0; i<number_of_taxa; i++) tobeexcluded[i] = FALSE;
	for(i=0; i<num_commands; i++)
		{
		if(strcmp(parsed_command[i], "deletetaxa") == 0)
			{

			for(k=i+1; k<num_commands; k++)
				{
				if(strcmp(parsed_command[k], "mintaxa") == 0)
					{
					min_taxa = toint(parsed_command[k+1]);
					if(min_taxa < 1)
						{
						printf2("Error the minimum number of taxa cannot be set below 1\n");
						error = TRUE;
						}
					k = k+1;
					}
				else
					{
					found = FALSE;
					for(j=0; j<number_of_taxa; j++)
						{
						if(strcmp(taxa_names[j], parsed_command[k]) == 0)
							{
							found = TRUE;
							tobeexcluded[j] = TRUE;
							}
						}
					if(!found)
						{
						printf2("Error: Cannot find any taxa that match the string \"%s\"\n", parsed_command[k]);
						error = TRUE;
						}
					else
						q++;
					}
				}
			i=k;
			}

		}

	if(k == 0)
		{
		printf2("Error: you must specify the name of at least one taxa to delete\n");
		error = TRUE;
		}
	else
		printf2("\n%d taxa will be permanently deleted from the source trees\n", q);


		/*** go through the names given and identify the taxa ***/
	if(!error)
		{
		strcpy(previnputfilename, inputfilename);
		strcpy(tmpfilename, inputfilename);
		strcat(tmpfilename, ".tempclann.chr");
		tempfile = fopen(tmpfilename, "w");
		command = malloc(10000*sizeof(char));
		command[0] = '\0';
	/*	printf2("\nWarning: This command will permenantly delete any chosen taxa from each tree in memory\nThis will also delete any trees that contain less than 4 taxa after the pruning\n\nAre you sure you wish to continue? (yes/no): ");
		xgets(command);
		if(strcmp(command, "n") == 0 || strcmp(command, "no") == 0 || strcmp(command, "NO") == 0 || strcmp(command, "N") == 0 || strcmp(command, "No") == 0)
			error = TRUE;
	*/	free(command);
		}
	if(!error)
		{
		j=0;
		for(i=0; i<number_of_taxa; i++)
			{
			if(tobeexcluded[i]) j++;
			}
		for(i=0; i<Total_fund_trees; i++)
			{
			if(sourcetreetag[i])
				{
				l =0;
				/***   check that this tree when pruned will still have more than 3 taxa **/
				for(j=0; j<number_of_taxa; j++)
					if(presence_of_taxa[i][j] > 0 && !tobeexcluded[j]) l+=presence_of_taxa[i][j];

				/* if there will be more than 3 taxa **/
				if(l>= min_taxa)
					{
					num_left++;
					strcpy(temptree, fundamentals[i]);
					returntree_fullnames(temptree, i);
					/*printf("%s\n", temptree); */

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



					/***  prune the sourcetree of any of the taxa **/
					prune_taxa_for_exclude(tree_top, tobeexcluded);
					shrink_tree(tree_top);    /* Shrink the pruned tree by switching off any internal nodes that are not needed */

					for(j=0; j<TREE_LENGTH; j++)
						{
						pruned_tree[j] = '\0'; /* initialise the string */
						tmp[j] = '\0';
						}
					/*get_taxa_details(tree_top);
                   printf2("--\n"); */

					if(print_pruned_tree(tree_top, 0, pruned_tree, TRUE, i) >1)
						{
						tmp[0] = '\0';
						strcpy(tmp, "(");
						strcat(tmp, pruned_tree);
						strcat(tmp, ")");
						strcpy(pruned_tree, tmp);
						}

					fprintf(tempfile, "%s", pruned_tree);
					fprintf(tempfile, " [%f]; [ %s]\n", tree_weights[i], tree_names[i]);
					}
				else
					{
					if(!done)printf2("\n\nThe following trees have been deleted from memory because they now contain less than 4 taxa\nTree Number\tTree name\n");
					done = TRUE;
					printf2("%-10d\t%s\n", i+1, tree_names[i]);
					}
				}
			}
		}


	free(tobeexcluded);
	free(pruned_tree);
	free(temptree);
	free(tmp);

	if(!error)
		{
		fclose(tempfile);
		if(num_left == 0)
			printf2("\nError: This will delete all trees from memory..... aborting deletetaxa\n");
		else
			execute_command(tmpfilename, do_all);
		strcpy(inputfilename, previnputfilename);
		remove(tmpfilename);
		}

	}

	/* Prune tree: This is a recursive function that is called for every node position of the supertree
	it then checks to see if any of the siblings on this node are not contained in the fundamental tree, these siblings are then turned off.
	This only turns off taxa, pointer siblings will have to be turned off using a separate program */

static void prune_taxa_for_exclude(struct taxon * super_pos, int *tobeexcluded)
	{


	/* Traverse the supertree, visiting every taxa and checking if that taxa is on the fundamental tree */

	while(super_pos != NULL)
		{

		if(super_pos->daughter != NULL) prune_taxa_for_exclude(super_pos->daughter, tobeexcluded);  /* If this is pointer sibling, move down the tree */

		if(super_pos->name != -1) /* If there is an actual taxa on this sibling */
			{
			/* Check to see if that taxa is on the fundamental tree */

			if(tobeexcluded[super_pos->name] == TRUE)  /* presence of taxa is an array that is filled in as the fundamental trees are inputted */
				{  /* if its not there */
				super_pos->tag = FALSE;
				}
			}

		super_pos = super_pos->next_sibling;

		}

	}
