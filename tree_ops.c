/*
 *  tree_ops.c — Clann v5.0.0
 *  Tree construction, manipulation, and traversal operations
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
#include "tree_ops.h"

/* Forward declarations for cross-module functions called by tree_ops.c */
int   sort_tree(struct taxon *position);              /* treecompare2.c */
int   assign_taxa_name(char *name, int fund);         /* tree_io.c */
void  subtree_id(struct taxon *position, int *tmp);   /* reconcile.c */

int unroottree(char * tree)
    {
    int i=0, j=0, k=0, l=0, m=0, basecount = 0, parentheses=0, did_unrooting=FALSE;
    int foundopen = FALSE, foundclose = FALSE;
	float del_nodelen = 0;
	char length[100], *restof = malloc(TREE_LENGTH * sizeof(char));

	restof[0] = '\0';
	length[0] = '\0';
	/* scan through the tree counting the number of taxa/nodes at the base (for it to be unrooted there should be at least three) */
    while(tree[i] != ';')
        {
        switch(tree[i])
            {
            case '(':
                    parentheses++;
                    i++;
                    break;
            case ')':
                    parentheses--;
                    i++;
                    break;
            case ',':
                    if(parentheses == 1)
                        {
                        basecount++;
                        }
                    i++;
                    break;
            default:
                    i++;
                    break;
            }
        }
        
    if(basecount <2)  /* if the base of the tree is rooted */
        {
        did_unrooting=TRUE;
        i=0;
        parentheses = 0;
        while(tree[i] != ';')  /* delete the two parentheses to make the tree unrooted */
            {
            switch(tree[i])
                {
                case '(':
                        parentheses++;
                        if(parentheses == 2 && !foundopen)
                            {
                            tree[i] = '^';
                            foundopen = TRUE;
                            }
                        i++;
                        break;
                case ')':
                        if(parentheses == 2 && !foundclose)
                            {
                            tree[i] = '^';
							i++;
							while(tree[i] != ')' && tree[i] != '(' && tree[i] != ',' && tree[i] != ';' && tree[i] != ':')
								{
								tree[i] = '^';
								i++;
								}
                            if(tree[i] == ':')
                                {
								k=0;
								length[0] = '\0';
                                while(tree[i] != ')' && tree[i] != '(' && tree[i] != ',' && tree[i] != ';')
                                    {
									if(tree[i] != ':')
										{
										length[k] = tree[i];
										k++;
										}
									tree[i] = '^';
                                    i++; 
                                    }
								length[k] = '\0';
                                }
							if(length[0] != '\0') /* we have a branch length on the internal branch, so we need to add it to the branch length of the component that is to the direct right of this closing parenthesis */
								{
								del_nodelen = atof(length);
								k=i+1; /* This should be whatever is after the ',' which should be the next compoonent */
								if(tree[k] == '(') /* we need to find the end of this clade and add the value there */
									{
									l=1; k++;
									while((l != 0 || tree[k-1] != ')') && tree[k] != ';' )   /* CHANGED RECENTLY FROM while(l != 0 && tree[k-1] != ')' && tree[k] != ';' ) */
										{
										switch(tree[k])
											{
											case '(':
												l++;
												k++;
												break;
											case ')':
												l--;
												k++;
												break;
											default:
												k++;
												break;
											}
										}
									k--; /* k now points to the closing bracket */
									/* read in the length attached to this partenthsis */
									k++;
									while(tree[k] != ')' && tree[k] != '(' && tree[k] != ',' && tree[k] != ';' && tree[k] != ':') k++; /* k now points to the ":" at the start of the length (if there is a length) */
									}
								else
									{
									while(tree[k] != ')' && tree[k] != '(' && tree[k] != ',' && tree[k] != ';' && tree[k] != ':') k++; /* k now points to the ":" at the start of the length (if there is a length) */
									}
								if(tree[k] == ':') /* there is length attached to this */
									{
									m=k+1;
									length[0] = '\0'; l=0;
									while(tree[m] != ')' && tree[m] != '(' && tree[m] != ',' && tree[m] != ';')
										{
										length[l] = tree[m];
										l++; m++;
										}
									length[l] = '\0';
									del_nodelen += atof(length);
									
									}
								else
									m=k;
								/* now add del_nodelen to this point in the tree */
								 l=0;
								while(tree[m] != ';' && tree[m] != '\0')
									{
									restof[l] = tree[m];
									m++; l++;
									}
								restof[l] = ';';
								restof[l+1] = '\0';
								if(tree[k] == ':')
									tree[k] = '\0';
								else
									{
									tree[k] = '\0';
									}
								length[0] = '\0';
								sprintf(length, ":%f", del_nodelen);
								strcat(tree, length);
								strcat(tree, restof);
								}
                            foundclose = TRUE;
                            }
                        i++;
                        parentheses--;
                        break;
                default:
                        i++;
                        break;
                }
            }
                        
        /* scan through the string shifting up the characters to take into account those parentheses that have been deleted */
        i=0; j=0;
        while(tree[j] != ';')
            {
            if(tree[j] == '^')
                {
                while(tree[j] == '^')
                    j++;
                if(i!= j)tree[i] = tree[j];
                i++; j++;
                }
            else
                {
                if(i!=j)tree[i] = tree[j];
                i++;j++;
                }
            }
        tree[i] = tree[j];
        tree[i+1] = '\0';
        
        }
    free(restof);
    return(did_unrooting);
    }

int treeToInt(char *tree)
    {
    int i;
    struct taxon *sorted_tree = NULL;

    /* The first part is to order the taxa in the tree so that they are in the order that they would if the tree was produced by "intTOtree"
       This is done by starting at the the highest-level componenets of the tree and  find the lowest taxon in each, and reorder the compoentes so the taxa IDs are in order.
       Then go to each of the sub-components and continue the process.

       for example:
       (1,(3,0),(4,2));
       As usual, each component is eiather a taxon or a set of parired parentheses. 
       To re-order this tree, start with the first component, which is the outer-most pair of parentheses. 
       This has 3 components:
       		1
       		(3,0)
       		(4,2)
		These components ordered by the lowest taxon ID in each, would be:
			(3,2) - because it hosts taxon 0.
			1	  - because it has taxon 1.
			(4,2) - because it has taxon 2.
		This results in: ((3,0), 1, (4,2));

		Next take the next component (3,0). The two sub-components are 3 and 0 (obviously).
		These ordered result in (0,3) (obviuosly:)
		The tree now looks like: ((0,3),1,(4,2));

		Do this to the final component with more than one taxa (4,2) and you finally end up with the tree:

		((0,3),1,(2,4));

	*/


    /* The first part is to calculate the path of the tree, given the tree in nested parenthesis format */
    /*  this part of the code assumes that the tree being tested has been created using the intToTree function. The reason for this
        is that that function always puts the position of a new taxon after the parenthesis or taxon it is to be paired with.
        If a tree to be tested is made in any other way, it will lead to these assumptions to be violated and the function may crash or
        the wrong number may be calculated. cc 16/06/2003    */

    /* Step 1 Buid the tree in memory */

    temp_top = NULL;
    tree_build(1, tree, sorted_tree, FALSE, -1, 0);
    sorted_tree = temp_top;
    temp_top = NULL;


    /* Step 2: carry out a recursive depth-first search of the tree, finding and labelling each clade with the lowest taxon ID */

    sort_tree (sorted_tree);
    

    /* Step 3: reverse the tree build, noting the path taken */
    
    
    
    
    
    
    
    return(1);
    }

void intTotree(int tree_num, char *array, int num_taxa)
    {
    int  i = 0, j=0, k=0, l=0, exit = 0, bracket_count = 0;
    char *tmparray = malloc(TREE_LENGTH * sizeof(char));
    char *string= NULL;
	double supers = 1, max = 1, min = 0, oldmin = 0, *path = NULL;

	if(!tmparray) { printf2("Error: out of memory in intTotree\n"); return; }
	
	for(i=4; i<=num_taxa; i++) supers*=((2*i)-5);

    path = malloc((num_taxa-3)*sizeof(double));  /* this will store the path of the tree */
    string = malloc(30*sizeof(char));
    
    /************* calculate the path of the tree from a triplet to the full size of the tree ***************/
    
    min = 1;
    max = supers; /* max is now the number of possible trees in this tree space */
   /* printf2("path of tree number %d with %d taxa is:\n", tree_num, number_of_taxa);
   */
     for(i=4; i<=num_taxa; i++)  /* from a triplet to the full complement of taxa, we figure out the placement of each new taxa at each stage */
        {
        
        path[i-4] = div((tree_num - min), ((max - (min - 1))/((2 * i)- 5))).quot + 1;

   /*     printf2("min = %d, max = %d path = %d\n",min, max, path[i-4]);
     */  
        if(i != num_taxa)
            {
            oldmin = min;
            min = ((path[i-4]-1)*((max - (oldmin - 1))/((2 * i)- 5))) + oldmin;   /* we have to recalculate the min and max values depending on the path previously */
        
            max = ((path[i-4])*((max - (oldmin - 1))/((2 * i)- 5))) + (oldmin -1);
            }
        
        }
    
    /**************** next build the tree using the path calculated in the previous section ***************/
    
    /* we start with a triplet */
    strcpy(array, "(0,1,2);");
    
    for(i=3; i<num_taxa; i++)  /* for each taxa to be added */
        {
        /* the path array calculated previously, tells us which object the new taxa is to be paired with. an object can be a set of parenthesis or a taxa */
        j=0; k=0;
        while(j!= path[i-3])
            {
        /*    printf2("%s\n", array);  */
            exit = 0;
            while(exit == 0)
                {
                switch(array[k])
                    {
                    case '(':
                        break;
                    case ',':
                        break;
                    case ')':
                        j++;
                        exit = 1;
                        break;
                    default:
                        while(array[k+1] != ',' && array[k+1] != '(' && array[k+1] != ')' )
                            {
                            k++;
                            }
                        j++;
                        exit = 1;
                        break;
                    }
                k++;
                }
            } /* the element at position k-1 is the last character of the object the new taxa is to be paired with */
        
        for(j=0; j<k; j++)
            {
            tmparray[j] = array[j]; /* copy in all the elements up to this point into the tmp array */
            }
        
        /* insert the new taxa name preceeded by a comma */
        tmparray[j] = ',';
        tmparray[j+1] = '\0';
        
        string[0] = '\0'; /*initialise the string to hold the name of the taxa */
        totext(i, string);  /* turn the integer name of the taxa into its character equivalent */
        strcat(tmparray, string); /* append the name onto the tree */
        strcat(tmparray, ")" );  /* append on the closing bracket to hold the new pairing */
        
        if(array[j-1] == ')')  /* if we are appending the new taxa to an internal node */
            {
            /* find the end of the string in tmparray */
            k=j;
            while(tmparray[k] != '\0')
                {
                k++;
                }
                
            /* find the position we need to add an open bracket in tmparray */
            l=j;
            bracket_count = 1;
            do
                {
                l--;
                if(tmparray[l] == '(') bracket_count--;
                if(tmparray[l] == ')') bracket_count++;
                }while(bracket_count != 1);

            /* l now is the element in the tmparray we need to add the open bracket */
            /* k is now the last element in the tmparray */
            
            /* move every thing from the end of the string to where we need to add a open bracket forward one space */
            
            for(k=k; k>=l; k--)
                {
                tmparray[k+1] = tmparray[k];
                }
            /* now place the new open bracket in the position of l */
            tmparray[l] = '(';
            }
        else  /* else if we are appending the new taxa to another taxa */
            {
             /* find the end of the string in tmparray */
            k=j;
            while(tmparray[k] != '\0')
                {
                k++;
                }
                
            /* find the position we need to add an open bracket in tmparray */
            l=j;
            while(tmparray[l-1] != ',' && tmparray[l-1] != ')' && tmparray[l-1] != '(')
                {
                l--;
                }
            /* l now is the element in the tmparray we need to add the open bracket */
            /* k is now the last element in the tmparray */
            
            /* move every thing from the end of the string to where we need to add a open bracket forward one space */
            
            for(k=k; k>=l; k--)
                {
                tmparray[k+1] = tmparray[k];
                }
            /* now place the new open bracket in the position of l */
            tmparray[l] = '(';
            
            
            }
        /* now append the rest of the tree onto tmparray */
        /* travel to the end of the tmparray */
        k = j;
        while(tmparray[k] != '\0')
            {
            k++;
            }
        while(array[j] != '\0')
            {
            tmparray[k] = array[j];
            j++;k++;
            }
        tmparray[k] = '\0';
        
        /* copy the new tree onto the main array */
        strcpy(array, tmparray);
            
            
        }
    
    
    
    free(tmparray);
    free(path);
    free(string);
    }

int tree_build (int c, char *treestring, struct taxon *parent, int fromfile, int fundnum, int taxaorder)
	{
	
	char temp[NAME_LENGTH], tmptag[100];
	struct taxon *position = NULL, *extra = NULL;
	int i = 0, j =0, end = FALSE, out = FALSE, onlylength = FALSE;
	position = make_taxon();  /* make a new instance of the taxon for this level */

	tmptag[0] = '\0';
	for(i=0; i<NAME_LENGTH; i++)
		temp[i] = '\0';
	
	if(parent == NULL) temp_top = position;  /* assign the pointer in the array to the top of this tree */
	else
		{
		position->parent = parent;
		(position->parent)->daughter = position;
		}
		
	out = FALSE;
	while(!out && treestring[c] != ';')
		{
		switch (treestring[c])
			{
			case '(':
				/* If we assigned this taxon before */
				if(position->daughter != NULL || position->name != -1)
					{
					/* make a new taxon */
					extra = make_taxon();
					position->next_sibling = extra;
					extra->prev_sibling = position;
					position = extra;
					}
                                c++;
				c = tree_build(c, treestring, position, fromfile, fundnum, taxaorder);
				i=0;
				onlylength = FALSE;
				if(treestring[c] == ':') onlylength=TRUE;
				while(treestring[c+i] != ')' && treestring[c+i] != ',' && treestring[c+i] != ';')
					{
					position->weight[i] = treestring[c+i];
					i++;
					}
				c += i;
				position->weight[i] = '\0';
				break;
				
			case ',':
				c++;
				break;
				
			case ')':
				out = TRUE;
				break;
			
			case ':':
				i=0;
				while(treestring[c] != ')' && treestring[c] != '(' && treestring[c] != ',' && treestring[c] != ';')
					{
					position->weight[i] = treestring[c];
					i++;
					c++;
					}
				position->weight[i] = '\0';
				break;
			case ';':
				break;
				
			default:
					
				/* If we assigned this taxon before */
				if(position->daughter != NULL || position->name != -1)
					{
					/* make a new taxon */
					extra = make_taxon();
					position->next_sibling = extra;
					extra->prev_sibling = position;
					position = extra;
					}
                                        
				for(i=0; i<NAME_LENGTH; i++) temp[i] = '\0';
				i=0;
				end = FALSE;
				do{
					temp[i] = treestring[c];
						i++;
						c++;
					}while(treestring[c] != ')' && treestring[c] != ',' && treestring[c] != '(' && treestring[c] != ':');
				temp[i] = '\0';


				if(fromfile)  /* if the tree has come from a file, then the tree will contain the actual taxa names, unlike if it was created by all_trees, where the tree will contain taxa numbers */
					{
					/* assign the name to the taxon */
					j = number_of_taxa;
					position->name = assign_taxa_name(temp, FALSE);
					if(j != number_of_taxa)
						{
						printf2("Error, Taxa name %s is not contained in the source trees\n", temp);
						}
					}
				else  /* if the tree has been created internally by an alltrees search */
					{
					/* now need to change the string to an integer number */
					position->name = 0;
					for(j=i; j>0; j--)
						{
						position->name += ( pow(10, (i-j)) * (texttoint(temp[j-1])));
						}
					}
				if(fundnum != -1)
					{
					position->fullname = malloc((strlen(fulltaxanames[fundnum][taxaorder])+10)*sizeof(char));
					position->fullname[0] = '\0';
					strcpy(position->fullname, fulltaxanames[fundnum][taxaorder]);
					}
				taxaorder++;
				break;		
			}
		}
        c++;
        return(c);
	}

int basic_tree_build (int c, char *treestring, struct taxon *parent, int fullnames)
	{
	
	char temp[NAME_LENGTH], tmptag[100];
	struct taxon *position = NULL, *extra = NULL;
	int i = 0, j =0, end = FALSE, out = FALSE, onlylength = FALSE;
	position = make_taxon();  /* make a new instance of the taxon for this level */

	tmptag[0] = '\0';
	for(i=0; i<NAME_LENGTH; i++)
		temp[i] = '\0';
	
	if(parent == NULL) temp_top = position;  /* assign the pointer in the array to the top of this tree */
	else
		{
		position->parent = parent;
		(position->parent)->daughter = position;
		}
		
	out = FALSE;
	while(!out && treestring[c] != ';')
		{
		switch (treestring[c])
			{
			case '(':
				/* If we assigned this taxon before */
				if(position->daughter != NULL || position->name != -1)
					{
					/* make a new taxon */
					extra = make_taxon();
					position->next_sibling = extra;
					extra->prev_sibling = position;
					position = extra;
					}
                                c++;
				c = basic_tree_build(c, treestring, position, fullnames);
				i=0;
				onlylength = FALSE;
				if(treestring[c] == ':') onlylength=TRUE;
				while(treestring[c+i] != ')' && treestring[c+i] != ',' && treestring[c+i] != ';')
					{
					position->weight[i] = treestring[c+i];
					i++;
					}
				c += i;
				position->weight[i] = '\0';
				break;
				
			case ',':
				c++;
				break;
				
			case ')':
				out = TRUE;
				break;
			
			case ':':
				i=0;
				while(treestring[c] != ')' && treestring[c] != '(' && treestring[c] != ',' && treestring[c] != ';')
					{
					position->weight[i] = treestring[c];
					i++;
					c++;
					}
				position->weight[i] = '\0';
				break;
			case ';':
				break;
				
			default:
					
				/* If we assigned this taxon before */
				if(position->daughter != NULL || position->name != -1)
					{
					/* make a new taxon */
					extra = make_taxon();
					position->next_sibling = extra;
					extra->prev_sibling = position;
					position = extra;
					}
                                        
				for(i=0; i<NAME_LENGTH; i++) temp[i] = '\0';
				i=0;
				end = FALSE;
				do{
					temp[i] = treestring[c];
						i++;
						c++;
					}while(treestring[c] != ')' && treestring[c] != ',' && treestring[c] != '(' && treestring[c] != ':');
				temp[i] = '\0';


				if(fullnames == TRUE)  /* if the tree has come from a file, then the tree will contain the actual taxa names, unlike if it was created by all_trees, where the tree will contain taxa numbers */
					{
					/* assign the name to the taxon */
					position->name = assign_taxa_name(temp, FALSE);
					position->fullname = malloc((strlen(temp)+10)*sizeof(char));
					position->fullname[0] = '\0';
					strcpy(position->fullname,temp);
					}
				else  /* if the neame is a taxon ID (i.e. a number) */
					{
					/* now need to change the string to an integer number */
					position->name = 0;
					for(j=i; j>0; j--)
						{
						position->name += ( pow(10, (i-j)) * (texttoint(temp[j-1])));
						}
					}
				break;		
			}
		}
        c++;
        return(c);
	}

struct taxon * make_taxon(void)
	{
	struct taxon *position = NULL;
	
	if(count_now)malloc_check++;
	
	position = malloc(sizeof(struct taxon));
	if(!position)  memory_error(34);
		
	position->name = -1;
	position->fullname = NULL;
	position->daughter = NULL;
	position->parent = NULL;
	position->prev_sibling = NULL;
	position->next_sibling = NULL;
	position->tag = TRUE;
	position->tag2 = FALSE;
        position->xpos = 0;
        position->ypos = 0;
        position->spr = FALSE;
		position->weight[0] = '\0';
	position->loss = 0;
	position->length = 0;
	position->donor = NULL;
	return(position);
	}

void prune_tree(struct taxon * super_pos, int fund_num)
	{
	int _g = 0, _gmax = number_of_taxa * 4 + 16;

	/* Traverse the supertree, visiting every taxa and checking if that taxa is on the fundamental tree */

	while(super_pos != NULL && !user_break && _g++ < _gmax)
		{

		if(super_pos->daughter != NULL) prune_tree(super_pos->daughter, fund_num);  /* If this is pointer sibling, move down the tree */

		if(super_pos->name != -1) /* If there is an actual taxa on this sibling */
			{
			/* Check to see if that taxa is on the fundamental tree */

			if(presence_of_taxa[fund_num][super_pos->name] == FALSE)  /* presence of taxa is an array that is filled in as the fundamental trees are inputted */
				{  /* if its not there */
				super_pos->tag = FALSE;
				}
			}

		super_pos = super_pos->next_sibling;

		}
	if(_g >= _gmax) rep_abandon = 1;

	}

int shrink_tree (struct taxon * position)
	{
	int count = 0, tot = 0, i, j, k, l, havelabel = FALSE;
	float one, two;
	char tempstr[1000] , tmplabel[1000];
	struct taxon *tmppos = NULL;
	int _g = 0, _gmax = number_of_taxa * 4 + 16;

	tmplabel[0] = '\0';

	while(position != NULL && !user_break && _g++ < _gmax)
		{
		if(position->daughter != NULL) 
			{
			tot = shrink_tree(position->daughter);  /* tot will be equal to the number of daughters that this pointer has */
		
			if(tot < 2)
				{
				position->tag = FALSE; /* if the number of daughters is less than 2, then we need to switch off this pointer sibling */
				if(tot == 1 && strcmp(position->weight, "") != 0)
					{
					k=0;
					while( k < strlen(position->weight) && position->weight[k] != ':') k++;
					if(k < strlen(position->weight) && position->weight[k] == ':')  /* we only want to add the branch lengths if they are branchlengths */
						{
							
						k++;
						tempstr[0] = '\0'; i=1; j=0;
						while(position->weight[k] != '\0')
							{
							tempstr[j] = position->weight[k];
							j++; k++;
							}
						tempstr[j] = '\0';
						one = atof(tempstr);
						tempstr[0] = '\0'; i=1; j=0;
						tmppos = find_remaining(position->daughter);
						if(tmppos != NULL)  /* if we haven't deleted everything in this clade */
							{
							k=0;
							while( k < strlen(tmppos->weight) && tmppos->weight[k] != ':') k++;
							havelabel = FALSE;
							if(k != 0)
								{
								for(l=0; l<k; l++)
									tmplabel[l] = tmppos->weight[l];
								tmplabel[l] = '\0';
								havelabel = TRUE;
								}

							k++;
							while(tmppos->weight[k] != '\0')
								{
								tempstr[j] = tmppos->weight[k];
								j++; k++;
								}
							tempstr[j] = '\0';
							two = atof(tempstr) + one;
							if(!havelabel)
								sprintf(tmppos->weight, ":%f", two);
							else
								sprintf(tmppos->weight, "%s:%f", tmplabel, two);
							}
						}
					}
				count += tot;
				}
			else
				count++;
			}
		else  /* else this is a taxa node */
			{
			if(position->tag) count++;  /* if this taxa is not switched off then increment count */
			}
		position = position->next_sibling;
		}
	if(_g >= _gmax) rep_abandon = 1;
	return(count);

	}

int print_pruned_tree(struct taxon * position, int count, char *pruned_tree, int fullname, int treenum)
    {
    /* name only holds a leaf label: either a decimal integer (~12 chars) or a
     * fullname (bounded by NAME_LENGTH).  Use a stack buffer — the old 1 MB
     * malloc was called recursively for every internal node on every source
     * tree causing massive heap contention under parallel load. */
    char name[NAME_LENGTH + 16];
    int _g = 0, _gmax = number_of_taxa * 4 + 16;

    name[0] = '\0';

    while(position != NULL && !user_break && _g++ < _gmax)
        {
        if(position->daughter != NULL)
            {
            if(position->tag != FALSE)
                {
                if(count > 0)
                    {
                    strcat(pruned_tree, ",");
                    } 
                strcat(pruned_tree, "(");
                count++;
                print_pruned_tree(position->daughter, 0, pruned_tree, fullname, treenum);
                strcat(pruned_tree, ")");
				}
            else
                {
                count = print_pruned_tree(position->daughter, count, pruned_tree, fullname, treenum);
                }
            }
        else
            {
            if(position->tag != FALSE)
                {
                if(count > 0)
                    {
                    strcat(pruned_tree, ",");
                    }

                name[0] = '\0';
               /* totext(position->name, name); */ /* depreciated CC - June 2016*/ 

                if(fullname == FALSE)
                    {
                    sprintf(name, "%d", position->name);
                    }
                else
                    {
                    strcpy(name, position->fullname);
                    }  
                strcat(pruned_tree, name);
                count++;
                }
            }
		if(position->tag != FALSE && strcmp(position->weight, "") != 0)
			{
			strcat(pruned_tree, position->weight);
			}
        position = position->next_sibling;
        }
    if(_g >= _gmax) rep_abandon = 1;
    return(count);
    }

void totext(int c, char *array)
    {
    int count = 0, i =0;
    char tmp[30];
    
    if(c != 0)
        {
        while(pow(10, count) <= c)
            {
            tmp[count] = inttotext(fmod(c/pow(10, count), 10));
            count++;
            }
        tmp[count] = '\0';
        array[count] = '\0';
        for(i= count-1; i>=0; i--)
            array[(count-1)- i] = tmp[i];
	}
    else
        {
        array[0] = '0';
        array[1] = '\0';
        }
    }

void reset_tree(struct	taxon * position)
	{
	int _g = 0, _gmax = number_of_taxa * 4 + 16;

	while(position != NULL && !user_break && _g++ < _gmax)
		{
		position->tag = TRUE;
		if(position->daughter != NULL) reset_tree(position->daughter);
		position = position->next_sibling;
		}
	if(_g >= _gmax) rep_abandon = 1;

	}		

int count_taxa(struct taxon * position, int count)
	{
	struct taxon * start = position;
	int i =0;
	int _g = 0, _gmax = number_of_taxa * 4 + 16;

	while(position != NULL && !user_break && _g++ < _gmax)
		{
		if(position->daughter != NULL)
			{
			count = count_taxa(position->daughter, count);
			}
		else
			{
			if(position->name != -1) count++;
			}
		position = position->next_sibling;
		}
	return(count);
	}

int find_taxa(struct taxon * position, char *query)
	{
	struct taxon * start = position;
	int i =0, found = FALSE;
	
	while(position != NULL)
		{
		if(position->daughter != NULL)
                    {
					if(found == FALSE)
						found = find_taxa(position->daughter, query);
                    }
                else
                    {
                    if(position->name != -1)
						{
						if(strstr(taxa_names[position->name], query) != NULL)
							found = TRUE;
						}
                    }
                position = position->next_sibling;
		}
	return(found);
	}

int number_tree(struct taxon * position, int num)
	{
	struct taxon * start = position;
	int i =0;
	int _g, _gmax = number_of_taxa * 4 + 16;

	_g = 0;
	while(position != NULL && !user_break && _g++ < _gmax)
		{
		if(position->daughter != NULL)
			num = number_tree(position->daughter, num);
		position = position->next_sibling;
		}
	position = start;
	_g = 0;
	while(position != NULL && !user_break && _g++ < _gmax)
		{
		position->tag = num;
		num++;
		position = position->next_sibling;
		}
	return(num);
	}

void check_tree(struct taxon * position, int tag_id, FILE *reconstructionfile)
	{
	struct taxon * start = position;
	int i =0, found = FALSE;
	
	
	while(position != NULL && !found)
		{
		if(position->daughter != NULL)
			{
			if(!found)check_tree(position->daughter, tag_id, reconstructionfile);
			if(tag_id == position->tag)
				{
				found = TRUE;
				fprintf(reconstructionfile, "%d %s\t",position->tag, position->weight);
				}
			}
		else
			{
			if(position->name != -1) 
				{
				if(tag_id == position->tag)
					{
					found = TRUE;
					fprintf(reconstructionfile, "%d %s\t",position->tag, taxa_names[position->name]);
					}				
				}
			}
		position = position->next_sibling;
		}

	}

int count_internal_branches(struct taxon *position, int count)
	{
	
	while(position != NULL)
		{
		if(position->daughter != NULL)
			{
			count++;
			count = count_internal_branches(position->daughter, count);
			}
		position = position->next_sibling;
		}
	return(count);
	}

void identify_taxa(struct taxon * position, int *name_array)
	{
	int _g = 0, _gmax = number_of_taxa * 4 + 16;
	while(position != NULL && !user_break && _g++ < _gmax)
		{
		if(position->daughter != NULL)
			{
			identify_taxa(position->daughter, name_array);
			}
		else
			{
			name_array[position->name] = TRUE;
			}
		position = position->next_sibling;
		}
	}

int check_taxa(struct taxon * position)
	{
	struct taxon * start = position;
	int i =0, number = 0;
	int _g = 0, _gmax = number_of_taxa * 4 + 16;

	while(position != NULL && !user_break && _g++ < _gmax)
		{
		if(position->daughter != NULL)
                    {
                    number += check_taxa(position->daughter);
                    }
                else
                    {
                    if(position->name != -1) number++;
                    }
                position = position->next_sibling;
		}
        return(number);
	}

void dismantle_tree(struct taxon * position)
	{
	struct taxon * start = position;
	int _g, _gmax = number_of_taxa * 4 + 16;

	/* first scan through this level and go down any pointer there are here */
	_g = 0;
	while(position != NULL && !user_break && _g++ < _gmax)
		{
		if(position->daughter != NULL)
			dismantle_tree(position->daughter);

		position = position->next_sibling;
		}
	position = start;
	_g = 0;
	while(position->next_sibling != NULL && _g++ < _gmax) position= position->next_sibling;

	/* now remove every thing at this level */
	{ int _g3 = 0;
	while(position != NULL && _g3++ < _gmax)
		{
		if(count_now)malloc_check--;
		start = position;
		position = position->prev_sibling;
		if(start->donor != NULL) free(start->donor);
		if(start->fullname != NULL) free(start->fullname);
		free(start);
		}
	}
	}

void print_named_tree(struct taxon * position, char *tree)
	{
	struct taxon *place = position;
	int count = 0, j=0;
        char *name = NULL;
	int _g = 0, _gmax = number_of_taxa * 4 + 16;
	strcat(tree, "(");

	name = malloc(30*sizeof(char));
	name[0] = '\0';
	while(position != NULL && !user_break && _g++ < _gmax)
		{
		if(count >0) strcat(tree, ",");
		if(position->daughter != NULL)
			{
			print_named_tree(position->daughter, tree);
		/*	sprintf(name, "%d", position->tag);
			strcat(tree, name);  */ /* This is here for the development of the gene tree mapping */

			}
		else
			{
			strcat(tree, taxa_names[position->name]);
			}
		count++;
		position = position->next_sibling;

		}
	strcat(tree, ")");
	free(name);
	}

void print_fullnamed_tree(struct taxon * position, char *tree, int fundtreenum)
	{
	struct taxon *place = position;
	int count = 0, j=0;
        char *name = NULL;
	int _g = 0, _gmax = number_of_taxa * 4 + 16;
	strcat(tree, "(");

	name = malloc(30*sizeof(char));
	name[0] = '\0';
	while(position != NULL && !user_break && _g++ < _gmax)
		{
		if(count >0) strcat(tree, ",");
		if(position->daughter != NULL)
			{
			print_fullnamed_tree(position->daughter, tree, fundtreenum);
		/*	sprintf(name, "%d", position->tag);
			strcat(tree, name);  */ /* This is here for the development of the gene tree mapping */

			}
		else
			{
			if(position -> fullname != NULL)
			{	
				strcat(tree, position -> fullname);
			}
			else
			{
				strcat(tree, taxa_names[position->name]);
			}
			/*strcat(tree, fulltaxanames[fundtreenum][position->name]);*/
			}
		count++;
		position = position->next_sibling;

		}
	strcat(tree, ")");
	free(name);
	}

void print_tree(struct taxon * position, char *tree)
	{
	struct taxon *place = position;
	int count = 0, j=0;
        char *name = NULL;
	int _g = 0, _gmax = number_of_taxa * 4 + 16;
	strcat(tree, "(");

	name = malloc(30*sizeof(char));
	name[0] = '\0';
	while(position != NULL && !user_break && _g++ < _gmax)
		{
		if(count >0) strcat(tree, ",");
		if(position->daughter != NULL)
			{
			print_tree(position->daughter, tree);
			}
		else
			{
			strcpy(name, "");
			totext(position->name, name);
			strcat(tree, name);
			}
		count++;
		position = position->next_sibling;

		}
	strcat(tree, ")");
	free(name);
	}

/* Print exactly one node's subtree as a Newick fragment — does NOT traverse
 * next_sibling of `node` itself.  Unlike print_tree() which iterates the
 * whole sibling chain starting at `position`, this function treats `node`
 * as the sole root of the subtree being serialised.  Used by the TBR
 * string-based rerooting helpers in treecompare2.c. */
void print_single_subtree(struct taxon *node, char *buf)
	{
	char name[30];
	int _g = 0, _gmax = number_of_taxa * 4 + 16;

	if(node == NULL) return;

	if(node->daughter != NULL)
		{
		struct taxon *child = node->daughter;
		int first = 1;
		strcat(buf, "(");
		while(child != NULL && !user_break && _g++ < _gmax)
			{
			if(!first) strcat(buf, ",");
			first = 0;
			print_single_subtree(child, buf);   /* recurse on each child */
			child = child->next_sibling;
			}
		strcat(buf, ")");
		}
	else
		{
		name[0] = '\0';
		totext(node->name, name);
		strcat(buf, name);
		}
	}

void print_tree_withinternals(struct taxon * position, char *tree)
	{
	struct taxon *place = position;
	int count = 0, j=0;
        char *name = NULL;
	strcat(tree, "(");
	/*printf("(\n");*/
        
	name = malloc(30*sizeof(char));
	name[0] = '\0';
	while(position != NULL)
		{
		if(count >0) strcat(tree, ",");
		if(position->daughter != NULL)
			{
			print_tree_withinternals(position->daughter, tree);
			}
		else
			{
			strcpy(name, "");
			totext(position->name, name);
			strcat(tree, name);
		/*	printf2("%d\n", name);*/
			}
		count++;
		position = position->next_sibling;
		}
	strcat(tree, ")");	
	/*printf(")\n");*/
	if(place->parent != NULL)
		{
		strcpy(name, "");
		totext((place->parent)->tag, name);
		strcat(tree, name);
		}
	
	free(name);
	}

void reallocate_retained_supers(void)
    {

    int i=0;
    number_retained_supers = number_retained_supers+10;
    
    retained_supers = realloc(retained_supers, number_retained_supers*sizeof(char*));
    if(!retained_supers) memory_error(48);
    best_topology = realloc(best_topology, number_retained_supers*sizeof(char*));
    if(!best_topology) memory_error(87);
    scores_retained_supers = realloc(scores_retained_supers, number_retained_supers*sizeof(float));
    if(!scores_retained_supers) memory_error(50);
    best_topology_scores = realloc(best_topology_scores, number_retained_supers*sizeof(float));
    if(!best_topology_scores) memory_error(88);

    for(i=number_retained_supers - 10; i<number_retained_supers; i++)
        {
        retained_supers[i] = malloc(TREE_LENGTH*sizeof(char));
        if(!retained_supers[i]) memory_error(49);
        best_topology[i] = malloc(TREE_LENGTH*sizeof(char));
        if(!best_topology[i]) memory_error(89);
        best_topology[i][0] = '\0';
        best_topology_scores[i] = -1;
        retained_supers[i][0] = '\0';
        scores_retained_supers[i] = -1;
        }
    }

void reroot_tree(struct taxon *outgroup)
	{
	struct taxon *newbie = NULL, *temp_treetop = NULL, *position = NULL, *temp = NULL, *start = NULL, *parent_start = NULL, *next = NULL, *parent = NULL, *pointer = NULL;
	char *temptree = NULL;
	int _g, _gmax = number_of_taxa * 4 + 16;
	temptree = malloc(TREE_LENGTH*sizeof(char));
	temptree[0] = '\0';
	newbie = make_taxon();
	temp_treetop = newbie;
	/* every sibling must become a descendent and every descentdent as sibling of that level */
	/*** Assign newbie to point to any siblings of outgroup as a daughter **/
	position = outgroup;
	if(outgroup->prev_sibling != NULL)
		{
		_g = 0;
		while(position->prev_sibling != NULL && _g++ < _gmax) position = position->prev_sibling;
		if(_g >= _gmax) { rep_abandon = 1; free(temptree); return; }
		}
	else
		position = position->next_sibling;

	start = position;
	newbie->daughter = start;

	/***2: extract the ougroup ***/
	if(outgroup->next_sibling != NULL)(outgroup->next_sibling)->prev_sibling = outgroup->prev_sibling;
	if(outgroup->prev_sibling != NULL)(outgroup->prev_sibling)->next_sibling = outgroup->next_sibling;
	if(outgroup->prev_sibling == NULL && outgroup->parent != NULL)
		{
		(outgroup->parent)->daughter = outgroup->next_sibling;
		(outgroup->next_sibling)->parent = outgroup->parent;
		outgroup->parent = NULL;
		}
	outgroup->next_sibling = NULL;
	outgroup->prev_sibling = NULL;

	if(start!=NULL)parent = start->parent;
	else parent=NULL;

	/** create new pointer on first level */
	pointer = make_taxon();
	position = start;
	_g = 0;
	while(position->next_sibling != NULL && _g++ < _gmax) position = position->next_sibling;
	if(_g >= _gmax) { free(pointer); free(temptree); rep_abandon = 1; return; }
	position->next_sibling = pointer;
	pointer->prev_sibling = position;

	/* Walk up the parent chain inverting parent-child relationships.
	 * Guard against a cyclic parent chain (corrupt tree) with a depth counter. */
	_g = 0;
	while(parent != NULL && _g++ < _gmax)
		{
		position = parent;
		{ int _g2 = 0;
		  while(position->prev_sibling != NULL && _g2++ < _gmax) position = position->prev_sibling;
		  if(_g2 >= _gmax) { free(temptree); rep_abandon = 1; return; }
		}
		parent_start = position;
		pointer->daughter = parent_start;
		next = parent_start->parent;
		parent_start->parent = pointer;
		pointer = parent;
		start = parent_start;
		parent = next;
		}
	if(_g >= _gmax) { free(temptree); rep_abandon = 1; return; } /* cyclic parent chain */
	pointer->daughter = NULL;

	(newbie->daughter)->parent = newbie;

	/* now make outgroup a siblig to newbie */
	newbie->next_sibling = outgroup;
	outgroup->prev_sibling = newbie;

	temp_top = temp_treetop;
	/**** Now we need to remove any pointer taxa that are nolonger pointing to anything **/
	clean_pointer_taxa(temp_top);
	/*** Finally we need to put an extra node at the top */
/*	newbie = make_taxon();
	newbie->daughter = temp_top;
	temp_top->parent = newbie;
	temp_top = newbie; */
	free(temptree);
	}

void clean_pointer_taxa(struct taxon *position)
	{

	struct taxon *start = position, *tmp = NULL;
	int i=0;
	int _g, _gmax = number_of_taxa * 4 + 16;

	_g = 0;
	while(position != NULL && !user_break && _g++ < _gmax)
		{
		if(position->daughter != NULL) clean_pointer_taxa(position->daughter);
		position = position->next_sibling;
		}
	position = start;

	/*** Check to see if there are any pointer taxa not going anywhere ***/

	_g = 0;
	while(position != NULL && !user_break && _g++ < _gmax)
		{
		if(position->name == -1 && position->daughter == NULL)
			{
			if(start->parent != NULL) 
				{
				if(strcmp(position->weight, "") != 0) strcpy((start->parent)->weight, position->weight); /* copy weight to parent ponter node */
				}
			tmp = position->next_sibling;
			if(position->next_sibling != NULL) (position->next_sibling)->prev_sibling = position->prev_sibling;
			if(position->prev_sibling != NULL) (position->prev_sibling)->next_sibling = position->next_sibling;
			if(position->parent != NULL)
				{
				if(position->next_sibling != NULL)
					{
					start = position->next_sibling;
					(position->parent)->daughter = position->next_sibling;
					(position->next_sibling)->parent = position->parent;
					}
				else
					(position->parent)->daughter = NULL;
				position->parent = NULL;
				}
			if(position == temp_top)
				 temp_top = position->next_sibling;
			if(position == start) start = position->next_sibling;
			free(position);
			position = tmp;
			}
		else
			if(position != NULL) position = position->next_sibling;
		}
	position = start;
	
	/*** Count the number of children of each pointer taxa, if any only have 1 then delete that pointer taxa and replace it with its daughter ***/
	/** however if there are any weights (like branch lenghs or bootstraps) they need to be assigned to the parent of this level  (Aug 17)**/

	_g = 0;
	while(position != NULL && !user_break && _g++ < _gmax)
		{
		if(position->daughter != NULL)
			{
			int _g2 = 0;
			i=0;
			tmp = position->daughter;
			while(tmp != NULL && _g2++ < _gmax)
				{
				i++;
				tmp = tmp->next_sibling;
				}
			if(i == 1)
				{
				tmp = position->daughter;
				if(position->prev_sibling != NULL)(position->prev_sibling)->next_sibling = position->daughter;
				if(position->next_sibling != NULL)(position->next_sibling)->prev_sibling = position->daughter;
				(position->daughter)->next_sibling = position->next_sibling;
				(position->daughter)->prev_sibling = position->prev_sibling;
				(position->daughter)->parent = position->parent;
				if(position->parent != NULL) (position->parent)->daughter = position->daughter;
				if(position == temp_top) temp_top = position->daughter;
				if(position->fullname != NULL) free(position->fullname);
				free(position);
				position = tmp;
				}
			else
				position = position->next_sibling;
			}
		else
			position = position->next_sibling;
		}
		
		
	}

struct taxon * get_branch(struct taxon *position, int name)
	{
	struct taxon *start = position, *answer = NULL;
	int _g, _gmax = number_of_taxa * 4 + 16;

	_g = 0;
	while(position != NULL && answer == NULL && !user_break && _g++ < _gmax)
		{
		if(position->daughter != NULL)
			{
			answer = get_branch(position->daughter, name);
			}
		position = position->next_sibling;
		}
	if(answer == NULL)
		{
		position = start;
		_g = 0;
		while(position != NULL && answer == NULL && !user_break && _g++ < _gmax)
			{
			if(position->tag == name)
				{
				answer = position;

				}
			position = position->next_sibling;
			}
		}
	return(answer);
	}

void get_taxa_details(struct taxon *position)
	{
	struct taxon *start = position;
	
	
	while(position != NULL)
		{
		if(position->daughter != NULL)
			{
			get_taxa_details(position->daughter);
			}
		else
			{
			printf2("Taxon %d Fullname %s\n", position->name, position->fullname);
			}
		position = position->next_sibling;
		}	
	}

void get_taxa_names(struct taxon *position, char **taxa_fate_names)
	{
	struct taxon *start = position;
	
	
	while(position != NULL)
		{
		if(position->daughter != NULL)
			{
			get_taxa_names(position->daughter, taxa_fate_names);
			}
		else
			{
			strcpy(taxa_fate_names[position->tag2], position->fullname);
			}
		position = position->next_sibling;
		}	
	}

struct taxon * get_taxon(struct taxon *position, int name)
	{
	struct taxon *start = position, *answer = NULL;
	
	while(position != NULL && answer == NULL)
		{
		if(position->daughter != NULL)
			{
			answer = get_taxon(position->daughter, name);
			}
		position = position->next_sibling;
		}
	if(answer == NULL)
		{
		position = start;
		while(position != NULL && answer == NULL)
			{
			if(position->name == name && position->daughter == NULL)
				{
				answer = position;
				}
			position = position->next_sibling;
			}
		}
	return(answer);
	}

int compress_tree (struct taxon * position)
	{
	int count = 0, tot = 0, done = FALSE;
	struct taxon *pos = NULL;
	
	while(position != NULL)
		{
		done = FALSE;
		if(position->daughter != NULL) 
			{
			tot = compress_tree(position->daughter);  /* tot will be equal to the number of daughters that this pointer has */
			if(tot < 2)
				{
				pos = position->next_sibling;
				done = TRUE;
				if(position->next_sibling != NULL) (position->next_sibling)->prev_sibling = position->daughter;
				if(position->prev_sibling != NULL) (position->prev_sibling)->next_sibling = position->daughter;
				if(position->parent != NULL) 
					{
					(position->parent)->daughter = position->daughter;
					(position->daughter)->parent = position->parent;
					}
				else
					(position->daughter)->parent = NULL;
				(position->daughter)->next_sibling = position->next_sibling;
				(position->daughter)->prev_sibling = position->prev_sibling;
				if(position == tree_top) tree_top = position->daughter;
				if(position->donor != NULL) free(position->donor);
				free(position);
				position = pos;
				
				 /* if the number of daughters is less than 2, then we need to remove this pointer sibling */
				count += tot;
				}
			else
				count++;
			}
		else  /* else this is a taxa node */
			{
			if(position->tag) count++;  /* if this taxa is not switched off then increment count */
			}				
		if(!done) position = position->next_sibling;
		}
	return(count);
	}

int compress_tree1 (struct taxon * position)
	{
	int count = 0, tot = 0, done = FALSE, i;
	struct taxon *pos = NULL, *start = position;
	while(position != NULL)
		{
		done = FALSE;
		if(position->daughter != NULL) 
			{
			tot = compress_tree1(position->daughter);  /* tot will be equal to the number of daughters that this pointer has */
			if(tot < 2)
				{
				pos = position->next_sibling;
				done = TRUE;
				if(position->donor != NULL)
					{
					for(i=0; i<num_gene_nodes; i++)
						{
						if(position->donor[i] == TRUE) (position->daughter)->donor[i] = TRUE;
						}
					}
				if(position->next_sibling != NULL) (position->next_sibling)->prev_sibling = position->daughter;
				if(position->prev_sibling != NULL) (position->prev_sibling)->next_sibling = position->daughter;
				if(position->parent != NULL) 
					{
					(position->parent)->daughter = position->daughter;
					(position->daughter)->parent = position->parent;
					}
				else
					(position->daughter)->parent = NULL;
				(position->daughter)->next_sibling = position->next_sibling;
				(position->daughter)->prev_sibling = position->prev_sibling;
				if(position == tree_top) tree_top = position->daughter;
				if(position->donor != NULL) free(position->donor);
				free(position);
				position = pos;
				
				 /* if the number of daughters is less than 2, then we need to remove this pointer sibling */
				count += tot;
				}
			else
				count++;
			}
		else  /* else this is a taxa node */
			{
			count++;  /* if this taxa is not switched off then increment count */
			}				
		if(!done) position = position->next_sibling;
		}
	if(start == tree_top && count == 1 && start->daughter != NULL) /* then this is the top node and we need to delete it here and now */
		{
		tree_top = start->daughter;
		(start->daughter)->parent = NULL;
		if(start->donor != NULL)
			{
			for(i=0; i<num_gene_nodes; i++)
				{
				if(start->donor[i] == TRUE) (start->daughter)->donor[i] = TRUE;
				}
			free(start->donor);
			start->donor = NULL;
			}
		free(start);
		start = NULL;
		}
	return(count);
	}

void duplicate_tree(struct taxon * orig_pos, struct taxon * prev_dup_pos)
	{
	struct taxon *sibling = NULL, *dup_pos = NULL;
	int done = FALSE, i;
	
	if(orig_pos->parent == NULL) done = TRUE;
	
	while(orig_pos != NULL && !user_break)
		{
		dup_pos = make_taxon();
		if(done)
			{
			temp_top = dup_pos;
			done = FALSE;
			}
		dup_pos->name = orig_pos->name;
		if(orig_pos->fullname != NULL)
			{
			dup_pos->fullname = malloc((strlen(orig_pos->fullname)+10)*sizeof(char));
			dup_pos->fullname[0] = '\0';
			strcpy(dup_pos->fullname, orig_pos->fullname);
			}
		
		dup_pos->tag = orig_pos->tag;
		dup_pos->tag2 = orig_pos->tag2;
		dup_pos->xpos = orig_pos->xpos;
		dup_pos->ypos = orig_pos->ypos;
		dup_pos->spr = orig_pos->spr;
		strcpy(dup_pos->weight, orig_pos->weight);
		dup_pos->loss = orig_pos->loss;
		dup_pos->length = orig_pos->length;
		if(orig_pos->donor != NULL)
			{
			dup_pos->donor = malloc(num_gene_nodes*sizeof(int));
			for(i=0; i<num_gene_nodes; i++) dup_pos->donor[i] = orig_pos->donor[i];
			}
		if(orig_pos->prev_sibling != NULL)
			{
			dup_pos->prev_sibling = sibling;
			sibling->next_sibling = dup_pos;
			}
		if(orig_pos->parent != NULL)
			{
			dup_pos->parent = prev_dup_pos;
			prev_dup_pos->daughter = dup_pos;
			}
		sibling = dup_pos;
		if(orig_pos->daughter != NULL) duplicate_tree(orig_pos->daughter, dup_pos);
		orig_pos = orig_pos->next_sibling;
		}
	}

void printnamesandtags(struct taxon *position)
	{
	while(position != NULL)
		{
		if(position->daughter != NULL) 
			{
			printf2(" ->%d", position->tag);
			printnamesandtags(position->daughter);
			printf2(":");
			}
		else printf2(" %d", position->name);
		position=position->next_sibling;
		}
	}

struct taxon * find_same(struct taxon * position, int tofind)
	{
	struct taxon *answer = NULL;
	while(position != NULL && answer == NULL)
		{
		if(position->tag == tofind)
			{
			answer = position;
			}
		if(position->daughter != NULL && position->tag2 == TRUE && answer == NULL)
			answer = find_same(position->daughter, tofind);
		position = position->next_sibling;
		}	
	return(answer);
	}

int are_siblings(struct taxon *position, int first, int second)
	{
	struct taxon *start = position;
	int i=0, answer = FALSE;
	while(position != NULL)
		{
		if(position->tag == first || position->tag == second || position->name == first || position->name == second)
			{
			i++;
			}
		position=position->next_sibling;
		}
	if(i == 2)
		{
		 answer = start->parent->tag;
		}
	if(!answer)
		{
		position = start;
		while(position != NULL)
			{
			if(position->daughter != NULL) answer = are_siblings(position->daughter, first, second);
			position=position->next_sibling;
			}
		}
	return(answer);
	}

int number_tree1(struct taxon * position, int num)
	{
	struct taxon * start = position;
	int i =0;
	
	
	while(position != NULL)
		{
		if(position->daughter != NULL)
			num = number_tree1(position->daughter, num);
		position = position->next_sibling;
		}
	position = start;
	while(position != NULL)
		{
		if(position->daughter != NULL)
			{
			position->tag = num;
			num++;
			}
		else
			{
			position->tag = position->name;
			}
		position = position->next_sibling;
		}
	return(num);
	}

int number_tree2(struct taxon * position, int num)
	{
	struct taxon * start = position;
	int i =0;
	
	
	while(position != NULL)
		{
		if(position->daughter != NULL)
			num = number_tree2(position->daughter, num);
		position = position->next_sibling;
		}
	position = start;
	while(position != NULL)
		{
		if(position->daughter == NULL)
			{
			position->tag2 = num;
			num++;
			}
		position = position->next_sibling;
		}
	return(num);
	}

void print_onetoone_names(struct taxon *position, int onetoone)
	{
    while(position != NULL)
		{
		if(position->daughter != NULL)
			{
			print_onetoone_names(position->daughter, onetoone);
			}
		if(position->fullname != NULL)
			{
			fprintf(onetoonefile, "%s\t", position->fullname);
			if(onetoone == 2) fprintf(strictonetoonefile, "%s\t", position->fullname);
			}
		position=position->next_sibling;
		}
	}

int isit_onetoone(struct taxon *position, int onetoone)  /* This will return 0 if not a 1:1, 1 if it was a relaxed 1:1 and 2 if it was a strict 1:1 */
	{
	while(position != NULL && onetoone != 0)
		{
		if(position->loss >= 1)
			{
			if(position->daughter != NULL) onetoone = 0;
			else onetoone = 1;
			}
		else
			{
			if(position->loss == -1)
				{
				if(position->daughter != NULL) onetoone = 0;
				else onetoone = 1;
				}
			}
		if(position->daughter != NULL && onetoone != 0)
			{
			if(position->loss != -1) onetoone = isit_onetoone(position->daughter, onetoone);
			}
		position = position->next_sibling;
		}
	return(onetoone);
	}

int get_min_node(struct taxon * position, int *presence, int num)
	{
	struct taxon * start = position;
	int *tmp, i, ans;
	
	tmp = malloc(number_of_taxa*sizeof(int));
	
	while(position != NULL)
		{
		if(position->daughter != NULL)
			{
			num = get_min_node(position->daughter, presence, num);
			for(i=0; i<number_of_taxa; i++) tmp[i] = presence[i];
			subtree_id(position->daughter, tmp);
			ans = TRUE;
			for(i=0; i<number_of_taxa; i++)
				{
				if(tmp[i] == TRUE) ans = FALSE;
				}
			if(ans == TRUE && position->tag < num) num = position->tag;
			}
		position = position->next_sibling;
		}
	free(tmp);
	return(num);
	}

void find(struct taxon * position)
	{
	while(position != NULL)
		{
		if(position->daughter != NULL) find(position->daughter);
		if(position->tag2 == TRUE)
			printf2("!%d\n", position->tag);
		position = position->next_sibling;
		}
	}

struct taxon * find_remaining(struct taxon * position)  /* This looks for the next tagged branch in the tree down from here and returns a pointer to it, otherwise it returns null pointer */
	{
	struct taxon * start = position, *pos = NULL, *found = NULL;
	int i = 0;
	
	while(position != NULL && found == NULL)
		{
		if(position->tag == TRUE) 
			{
			found = position;			
			}
		position = position->next_sibling;
		}
	position = start;
	while(position != NULL && found == NULL)
		{
		if(position->daughter != NULL) found = find_remaining(position->daughter);
		position = position->next_sibling;
		}
	return(found);
	}

void find_tagged(struct taxon * position, int *thispresence)
	{
	struct taxon * start = position, *pos = NULL;
	int i = 0;
	
	while(position != NULL)
		{
		if(position->daughter != NULL) find_tagged(position->daughter, thispresence);
		
		if(position->tag2 == TRUE) 
			{
			/* Mark this as being present */
			thispresence[position->tag] = TRUE;
			if(position->daughter != NULL) down_tree(position->daughter, position, thispresence);
			/**** count how many siblings on this level are not "lost"   */
			pos = start;
			i=0;
			down_tree(start, position, thispresence);
						
			}
		position = position->next_sibling;
		}
	}

void up_tree(struct taxon * position, int * presence)
	{
	int i=0;
	struct taxon * prev = position , *start = NULL;
	presence[position->tag] = TRUE;
	
	while(position->next_sibling != NULL) position = position->next_sibling;
	while(position != NULL)
		{
		if(position->loss != -1) i++;
		start = position;
		position = position->prev_sibling;
		}
	position = start;	
	if(i<2)
		{
		if(position->parent != NULL) up_tree(position->parent, presence);
		while(position != NULL)
			{
			presence[position->tag] = TRUE;
			if(position != prev && position->daughter != NULL) down_tree(position->daughter, position, presence);
			position =position->next_sibling;
			}
		}
	}

void down_tree(struct taxon * position, struct taxon *prev, int * presence)
	{
	int i=0;
	struct taxon *start = position;
	while(position != NULL)
		{
		if(position->loss == -1)
			{
			presence[position->tag] = TRUE;
			if(position->daughter != NULL) down_tree(position->daughter, position, presence);
			}
		else
			i++;
		position = position->next_sibling;
		}
	if(i < 2)
		{
		position = start;
		if(start->parent != NULL && start->parent != prev) up_tree(start->parent, presence);
		while(position != NULL)
			{
			if(position->loss != -1 && position->tag2 != TRUE)
				{
				presence[position->tag] = TRUE;
				if(position->daughter != NULL) down_tree(position->daughter, position, presence);
				}
			position = position->next_sibling;
			}
		}
	}

void reset_tag2(struct taxon * position)
	{
	while(position != NULL)
		{
		position->tag2 = FALSE;
		if(position->daughter != NULL) reset_tag2(position->daughter);
		position = position->next_sibling;
		}
	}

int assign_tag2(struct taxon * position, int num)
	{
	int found = FALSE, i;
	while(position != NULL)
		{
		if(position->daughter != NULL)
			{
			if(assign_tag2(position->daughter, num)) found = TRUE;
			}
		if(position->donor == NULL)
			{
			position->donor = malloc(num_gene_nodes*sizeof(int));
			for(i=0; i<num_gene_nodes; i++) position->donor[i] = FALSE;
			}
		if(position->donor[num] == TRUE)
			{
			position->tag2 = TRUE;
			found = TRUE;
			}
		position = position->next_sibling;
		}
	return(found);
	}

void get_taxa(struct taxon *position, int *presence)
	{
	while(position != NULL)
		{
		if(position->daughter != NULL) 
			get_taxa(position->daughter, presence);
		else
			presence[position->name] = TRUE;
		
		position = position->next_sibling;
		}
	}

int get_best_node(struct taxon * position, int *presence, int num)
	{
	struct taxon * start = position;
	int *tmp, *tmp1, i, ans, ans1, num1;
	
	tmp = malloc(number_of_taxa*sizeof(int));
	tmp1 = malloc(number_of_taxa*sizeof(int));
	while(position != NULL)
		{
		if(position->daughter != NULL)
			{
			num = get_best_node(position->daughter, presence, num);
			for(i=0; i<number_of_taxa; i++)
				{
				tmp[i] = presence[i];
				if(presence[i] == TRUE ) tmp1[i] = FALSE;
				else tmp1[i] = TRUE;
				}
			subtree_id(position->daughter, tmp);
			subtree_id(position->daughter, tmp1);
			ans = TRUE; ans1=TRUE;
			for(i=0; i<number_of_taxa; i++)
				{
				if(tmp[i] == TRUE) ans = FALSE;
				if(tmp1[i] == TRUE) ans1 = FALSE;
				}
			if(ans == TRUE && (num == -1 || position->tag < num)) 
				num = position->tag;
			else
				{
				if(num == -1)
					{
					if(ans1 == TRUE && position->tag > num)
						num = position->tag;
					}
				}
			}
		position = position->next_sibling;
		}
	free(tmp);
	free(tmp1);
	return(num);
	}

void check_treeisok(struct taxon *position)
	{
	struct taxon *start = position;
	
	if(position->parent != NULL) printf2("^^ position->parent %d\n", position->parent->tag);
	while(position != NULL)
		{
		if(position->prev_sibling != NULL) printf2("[prev %d\t", position->prev_sibling->tag);
		else printf2("[no prev\t");
		if(position->parent != NULL) printf2("^^ parent = %d ^^\t ", position->parent->tag);
		if(position->daughter != NULL)
			printf2("pointer vv%d (daughter = %d)\t", position->tag, position->daughter->tag);
		else
			printf2("taxa %d\t", position->name);
		if(position->next_sibling != NULL) printf2("\tnext %d]\t", position->next_sibling->tag);
		else printf2("no next]\n");
		position = position->next_sibling;
		}
	position = start;
	while(position != NULL)
		{
		if(position->daughter != NULL) check_treeisok(position->daughter);
		position = position->next_sibling;
		}
	}

