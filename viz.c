/*
 *  viz.c — Clann v5.0.0
 *  Tree visualisation: coordinate assignment, PostScript output, histogram
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

#include "viz.h"

/** this labels all the taxa nodes */
int xposition1 (struct taxon *position, int count)
    {
    
    while(position != NULL)
        {
        if(position->daughter != NULL)
            count = xposition1(position->daughter, count);
        else
            {
            position->xpos = (float)count;
            count++;
            }
        position = position->next_sibling;
        }
    return(count);
    }
    
    
/* find the smallest xpos label on the next level */
float middle_number(struct taxon *position)
    {
    float start, end;
    
    start= position->xpos;
    while(position->next_sibling != NULL)
        {
        position = position->next_sibling;
        }
    end = position->xpos;
    return((start+end)/2);
    }


/* set all of the pointer taxa xpos */
void xposition2(struct taxon *position)
    {
    while(position != NULL)
        {
        if(position->daughter != NULL)
            {
            xposition2(position->daughter);
            position->xpos = middle_number(position->daughter);
            }
        position = position->next_sibling;
        }
    }

/** the next functions are  for determining the y position of each node of the tree **/
int yposition0(struct taxon *position, int level, int deepest)
    {
    int tmp, deephere = level;
	struct taxon *start = position;
	
    while(position != NULL)
        {
        if(position->daughter != NULL)
            {
            tmp = yposition0(position->daughter, level+1, deepest);
			if(tmp > deephere) deephere = tmp;
            }
        position = position->next_sibling;
        }
	position = start;
	while(position != NULL)
        {
		position->ypos = (deepest - (deephere-level)) -1;
		position = position->next_sibling;
		}
	return(deephere);	
    }



/** the next functions are  for determining the y position of each node of the tree **/
int yposition1(struct taxon *position, int level)
    {
    int tmp, deepest = level;
    while(position != NULL)
        {
        if(position->daughter != NULL)
            {
            tmp = yposition1(position->daughter, level+1);
            if(tmp > deepest) deepest = tmp;
            }
        position = position->next_sibling;
        }
    return(deepest);
    }


void yposition2(struct taxon *position, int deepest)
    {
    while(position != NULL)
        {
        if(position->daughter != NULL)
            yposition2(position->daughter, deepest);
            
        position->ypos = position->ypos/deepest;
        position = position->next_sibling;
        }
    }

void printcolour(float real, int branch)
        {
        float value;
		
		if(branch)value = real/largest_length;
		else value = real;
		
        if(value >= 1)
            fprintf(psfile, "\t0.933 0.000 0.000 setrgbcolor\n");
        else
            {
            if(value >= .89)
                fprintf(psfile, "\t1.000 0.100 0.000 setrgbcolor\n");
            else
                {
                if(value >= .79)
                    fprintf(psfile, "\t1.000 0.200 0.000 setrgbcolor \n");
                else
                    {
                    if(value >= .69)
                        fprintf(psfile, "\t1.000 0.300 0.000 setrgbcolor \n");
                    else
                        {
                        if(value >= .59)
                            fprintf(psfile, "\t1.000 0.400 0.000 setrgbcolor \n");
                        else
                            {
                            if(value >= .49)
                                fprintf(psfile, "\t1.000 0.500 0.000 setrgbcolor \n");
                            else
                                {
                                if(value >= .39)
                                    fprintf(psfile, "\t1.000 0.600 0.000 setrgbcolor \n");
                                else
                                    {
                                    if(value >=.29)
                                        fprintf(psfile, "\t1.000 0.700 0.000 setrgbcolor \n");
                                    else
                                        {
                                        if(value >=.19)
                                            fprintf(psfile, "\t1.000 0.800 0.000 setrgbcolor \n");
                                        else
                                            {
                                            if(value >=.09)
                                                fprintf(psfile, "\t1.000 9.000 0.000 setrgbcolor \n");
											else
												{
												if(value > 0)
													fprintf(psfile, "\t1.000 1.000 0.000 setrgbcolor \n");
												else
													{
													if(value > -1)
														fprintf(psfile, "\t0 0 0 setrgbcolor \n");
													else
														fprintf(psfile, "\t.8 .8 .8 setrgbcolor \n");
													}
												}
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }


void print_coordinates(struct taxon *position, char **treearray, int taxa_count, int mapping)
    {
    int i, j, x, y, v1, v2, y1;
    float x1, x2;
    struct taxon *start = position;
    
        
    while(position != NULL)
        {
        x = (int)((position->xpos/taxa_count)*(taxa_count*2));
        y = (int)(position->ypos*(60));

        if(position->daughter != NULL)
            {
            v1 = (int)(position->ypos*500); if(v1 == 0) v1 = 30;
			printcolour(position->loss, 0);
            fprintf(psfile, "%d %d moveto\n%d %d lineto stroke\n", v1, (int)((position->xpos/taxa_count)*750)+15, (int)((position->daughter)->ypos*500), (int)((position->xpos/taxa_count)*750)+15);
			fprintf(psfile, "\t0 0 0 setrgbcolor \n");
            if(position->weight[0] != ':') fprintf(psfile, "%d %d moveto\n (%s) show\n", (int)((position->daughter)->ypos*500)+5, (int)((position->xpos/taxa_count)*750)+12, position->weight);
			treearray[x][y] = '+';
            for(i=y+1; i< (int)((position->daughter)->ypos*60); i++)
                {
				if(!mapping)
					{
					if(position->loss == 0) treearray[x][i] = '-';
					if(position->loss == -1) treearray[x][i] = '.';
					if(position->loss >= 1) treearray[x][i] = '=';
					}
				else
					treearray[x][i] = '-';
				}

            if(strcmp(position->weight, "") != 0 && position->weight[0] != ':')  /* if there is a weight on the tree */
				{
				i = i+2;
				j=0;
				
				while(position->weight[j] != '\0' && position->weight[j] != ':')
					{
					treearray[x][i+j] = position->weight[j];
					j++;
					}
				}
            print_coordinates(position->daughter, treearray, taxa_count, mapping);
            }
        else
            {
            v1 = (int)(position->ypos*500);
            if(v1 == 0) v1 = 30;
			printcolour(position->loss, 0);
            fprintf(psfile, "%d %d moveto\n%d %d lineto stroke\n", v1, (int)((position->xpos/taxa_count)*750)+15 , 500, (int)((position->xpos/taxa_count)*750)+15);
			fprintf(psfile, "\t0 0 0 setrgbcolor \n");
            fprintf(psfile, "%d %d moveto\n (%s) show\n", 505, (int)((position->xpos/taxa_count)*750)+12, taxa_names[position->name]);
            treearray[x][y] = '+';
            for(i=y+1; i<60; i++)
				{
				if(!mapping)
					{
					if(position->loss == 0) treearray[x][i] = '-';
					if(position->loss == -1) treearray[x][i] = '.';
					if(position->loss >= 1) treearray[x][i] = '=';
					}
				else
					treearray[x][i] = '-';
				}
            j=0; i++;
			if(position->fullname != NULL)
				{
				while(i<80 && position->fullname[j] != '\0')
					{
					treearray[x][i] = position->fullname[j];
					i++; j++;
					}
				}
			else
				{
				while(i<80 && taxa_names[position->name][j] != '\0')
					{
					treearray[x][i] = taxa_names[position->name][j];
					i++; j++;
					}
				}
            }
        
        position = position->next_sibling;
        }
    
    /** Now do the vertical lines **/
    
    position = start;
    
    v1 = (int)((position->xpos/taxa_count)*(taxa_count*2));
    x1 = (position->xpos/taxa_count);
    if(x1 == 0) v1 =1;
    
    while(position->next_sibling != NULL)
        position = position->next_sibling;
    
    v2 = (int)((position->xpos/taxa_count)*(taxa_count*2));
    x2 = (position->xpos/taxa_count);
    if(x2 == 0) x2 = 30;
    y1 = (int)(position->ypos*500);
    if(y1 ==0) y1 = 30;

	printcolour(start->loss, 0);
    fprintf(psfile, "%d %d moveto\n%d %d lineto stroke\n", y1,(int)(x1*750)+15, y1,(((int)(x2*750)+15+(int)(x1*750)+15)/2)+1 );
	printcolour(position->loss, 0);
	fprintf(psfile, "%d %d moveto\n%d %d lineto stroke\n", y1,(int)(x2*750)+15, y1,(((int)(x2*750)+15+(int)(x1*750)+15)/2)+1 );

    for(i=v1+1; i<v2; i++)
        treearray[i][(int)(position->ypos*(60))] = '|';
    
    }

/* Print the tree graphically to screen */
void tree_coordinates(char *tree, int bootstrap, int build, int mapping, int fundnum)
    {
    char **treearray = NULL;
    int i=0, j=0, taxa_count = 0, deepest = 0, acrosssize = 20, upsize = 5, taxaorder=0;
    
	if(build)
		{
		 /* build the tree in memory */
		/****** We now need to build the Supertree in memory *******/
		if(tree_top != NULL)
			{
			dismantle_tree(tree_top);
			tree_top = NULL;
			}
		temp_top = NULL;
		largest_length = 0;
		/*tree_build(1, tree, tree_top, 1, fundnum, 0); */ /* removed 21/03/25 to allow fullnames to be displayed in showtrees */ 
		basic_tree_build(1, tree, tree_top, TRUE);
		tree_top = temp_top;
		temp_top = NULL;
		}
	/* check how many taxa are on the tree */
	taxa_count = count_taxa(tree_top, 0);
	
    /* create the array to hold the tree graphically */
    treearray = malloc((taxa_count*2+1)*sizeof(char*));
    if(!treearray) memory_error(70);
    for(i=0; i<taxa_count*2+1; i++)
        {
        treearray[i] = malloc(80*sizeof(char));
        if(!treearray[i]) memory_error(71);
        
        for(j=0; j<80; j++) treearray[i][j] = ' ';
        }


	if(bootstrap && !mapping)
		psfile = fopen("trees.ps", "w");
	else
		{
		if(mapping)
			psfile = fopen("trees.ps", "w");
		}
    
    /* now deal with the x coordinates */
    xposition1(tree_top, 1);   /* label all the taxa nodes with xpos. */
    xposition2(tree_top); /* label all the pointer nodes with xpos. */
    
	
	deepest = yposition1(tree_top, 0)+1;
	yposition0(tree_top, 0, deepest);
    /* now deal with the y coordinates */
    yposition2(tree_top, deepest);
  
	 /** open the post-script file and start the blurb **/
    
    fprintf(psfile, "%%!PS-Adobe-3 Created by Clann\n");
    fprintf(psfile, "\n\n2 setlinecap\n2 setlinejoin\n1.5 setlinewidth\n\n");
    fprintf(psfile, "\n\n%%Draw lines\n");
    if(number_of_taxa <=50)fprintf(psfile,"\n\n/Times findfont 12 scalefont setfont\n");
	else
		{
		if(number_of_taxa <=100)fprintf(psfile,"\n\n/Times findfont 8 scalefont setfont\n");
		else
			fprintf(psfile,"\n\n/Times findfont 5 scalefont setfont\n");
		}
    
   
	if(mapping)
		{
		  /*** Put in scale of colours **/
		fprintf(psfile, "\t.8 .8 .8 setrgbcolor\n");
		fprintf(psfile, "%d %d moveto\n%d %d lineto\n%d %d lineto\n%d %d lineto\n closepath\nfill\n", acrosssize+5, 3*upsize, 11*acrosssize+5, 3*upsize, 11*acrosssize+5, 3*upsize+10, acrosssize+5, 3*upsize);
		
		for(i=0; i<=10; i++)
			{
			if(i == 0)fprintf(psfile, "\t.8 .8 .8 setrgbcolor\n");
			else printcolour((float)i/(float)10, 0);
			fprintf(psfile, "%d %d moveto\n%d %d lineto\n%d %d lineto\n%d %d lineto\nclosepath\nfill\n\n",(acrosssize*i)+5, 5, ((i+1)*acrosssize)+5, 5, ((i+1)*acrosssize)+5, (upsize)+5, ((i)*acrosssize)+5, (upsize)+5 );
			}
		}
		
    /* now print the coordinates */
    print_coordinates(tree_top, treearray, taxa_count, mapping);
    fprintf(psfile, "\nshowpage\n");
    
    
    /** print out the tree graphically **/
    for(i=0; i<taxa_count*2+1; i++)
        {
        for(j=0; j<80; j++)
            {
            printf2("%c", treearray[i][j]);
            }
        printf2("\n");
        }
	
	if(bootstrap)
		fclose(psfile);
		
    for(i=0; i<taxa_count*2+1; i++)
        free(treearray[i]);
    free(treearray);
    }
    
void draw_histogram(FILE *outfile, int bins, float *results, int num_results)
	{
	int i, j, x, y, z, value;
	double max1, min = -1, max = -1, mean = 0, std_dev,  variance, skewness,  k;
	int *hist = NULL;
	for(i=0; i<num_results; i++)
		{
		if(min == -1)
			{
			min = results[i];
			max = results[i];
			}
		else
			{
			if(results[i] < min)
				min = results[i];
			if(results[i] > max)
				max = results[i];
			}
		/* calculate mean */
			
		mean += results[i];
		}

	mean = mean/num_results;

	if(outfile != NULL)fprintf(outfile, "range\tnumber of trees\n");
	if(outfile != NULL)
		{
		if(bins ==0)
			{
			bins = 20;
			k = 1;
			}
		else
			{
			k = ((max - min) / bins);
			}
		}
	else
		k =1;
	hist = malloc(bins*sizeof(int));
	for(i=0; i<bins; i++) hist[i] = 0;
	for(i=0; i<num_results; i++)
		{
		value = (int)((results[i]-min)/k);
		if(value >= bins) value = bins-1;
		if(value < 0) value = 0;

		hist[value]++;

		}
	max1 = -1;
	for(i=0; i<bins; i++)
		{
		if(hist[i] > max1) max1 = hist[i];
		}
		
	if(outfile != NULL)printf2("\nResults as follows:\n\n\n");
	for(i=0; i<bins; i++)
		{
		 x = (hist[i]/max1)*40;
		if(outfile != NULL)fprintf(outfile, "%.2f-%.2f\t%d\n", min+(i*k), (min+((i+1)*k))-0.01, hist[i]);
		if(outfile != NULL)printf2("%.2f - %.2f\t|", min+(i*k), (min+((i+1)*k))-0.01);
		else printf2("\t%-4.0f|", min+(i*k));
		for(j=0; j<x; j++)
			printf2("=");
		printf2(" (%d)\n",hist[i]);
		y=i;
		z =0;
		while(y < bins && hist[y] == 0)
			{
			y++;
			z++;
			}
		if(z > 15)
			{
			printf2("\t    .\n\t    .\n\t    .\n\t    \n");
			i = y-3;
			}
		
		}
		printf2("\n\n\n");

	/* calculate the rest of the statistics */
	variance = 0;
	for(i=0; i<num_results; i++)
		{
		variance += pow((results[i] - mean), 2);
		}
	variance = variance/(num_results -1);
	std_dev = sqrt(variance);
	skewness = 0;
	for(i=0; i<num_results; i++)
		{
		skewness += pow((results[i] - mean)/std_dev, 3);
		}
	skewness = skewness/num_results;
	
	
	if(outfile != NULL)printf2("Moments of the Distribution:\n\n\tMean = %f\n\tVariance = %f\n\tStandard Deviation = %f\n\tSkewness = %f\n\tStandard deviation of skewness = %f\n\n", mean, variance, std_dev, skewness, sqrt(((double)6)/(double)num_results));
	
	
	
	
	
	}
