/*
 *  treecompare2.c
 *  
 *
 *  Created by Chris Creevey on Sat Oct 04 2003.
 *  Copyright (c) 2003 - 2016 Chris Creevey. All rights reserved.
 *
 */

#include "config.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <time.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>

#ifdef HAVE_READLINE /* from config.h */
#include <readline/readline.h>
#include <readline/history.h>
#endif

/*#include "readline.h"
#include "history.h"
*/

/****** Define  ********/

#define FUNDAMENTAL_NUM 1000
#define TREE_LENGTH 1000000
#define TAXA_NUM 1000
#define NAME_LENGTH 1000



#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif
#ifndef TEMP
#define TEMP 2
#endif



/******** Structure definitions *************/

struct taxon {
	
	int name;					/* This holds number which refers to the name of the taxon at this node if it doesn't point to a daughter node */
	char *fullname;				/* This holds the position of the original name in the array fulltaxanames as it was given in the tree file, this allows us to record if there was parts of the name excluded earlier */
	float score;
	struct taxon *daughter;   	/* This points to a daughter node, but is only set if name is null */
	struct taxon *parent;		/* This points to the parent node, however its only set on the first sibling at each level */
	struct taxon *prev_sibling;	/* Points to the previous sibling at this node */
	struct taxon *next_sibling;	/* Points to the next sibling at this node */
	int tag;					/* A tag which is included as something which might be useful at some point */
	int tag2;					/* sometimes one of something just isn't enough ;o) */
	float loss;					/* A tag which will tell whether this node was lost.... only used for tree-mapping procedure */
        float xpos;   /* used to fix the position of the node for graphical representation */
        float ypos;	
        int spr;
		char weight[100];  /* this will be used when we have a weight to add to the node like with a Bootstrap proportion */
		float length;
		int *donor; /* THis is for the HGT analysis */
	} taxon_type;
	
	





/*********** Function definitions ********************/

void find(struct taxon * position);
int run_main(int argc, char *argv[]);
void input_file_summary(int do_all);
void clean_exit(int error);
char *xgets(char *s);
char inttotext(int c);
void totext(int c, char *array);
int assign_taxa_name(char *name, int fund);
void execute_command(char *filename, int do_all);
int seperate_commands(char *command);
int parse_command(char *command);
void print_commands(int num);
int texttoint(char c);
void cal_fund_scores(int printfundscores);
void pathmetric(char *string, int **scores);
void weighted_pathmetric(char *string, float **scores, int fund_num);
void unroottree(char * tree);
void alltrees_search(int user);
float compare_trees(int spr);
struct taxon * make_taxon(void);
void intTotree(int tree_num, char *array, int num_taxa);
int tree_build (int c, char *treestring, struct taxon *parent, int fromfile, int fund_num);
void prune_tree(struct taxon * super_pos, int fund_num);
int treeToInt(char *array);
int shrink_tree (struct taxon * position);
int print_pruned_tree(struct taxon * position, int count, char *pruned_tree, int fullname);
void reset_tree(struct	taxon * position);
int count_taxa(struct taxon * position, int count);
void check_tree(struct taxon * position, int tag_id, FILE *reconstructionfile);
int check_taxa(struct taxon * position);
int find_taxa(struct taxon * position, char *query);
int number_tree(struct taxon * position, int num);
void dismantle_tree(struct taxon * position);
void bootstrap_search(void);
void memory_error(int error_num);
void print_named_tree(struct taxon * position, char *tree);
void print_tree(struct taxon * position, char *tree);
int toint(char *number);
void reallocate_retained_supers(void);
void usertrees_search(void);
void heuristic_search(int user, int print, int sample, int nreps);
int average_consensus(int nrep, int missing_method, char * useroutfile, FILE *paupfile);
int do_search(char *tree, int user, int print, int maxswaps, FILE *outfile, int numspectries, int numgenetries);
int branchswap(int number_of_swaps, float score, int numspectries, int numgenetries);
int find_swaps(float * number, struct taxon * position, int number_of_swaps, int numspectries,int numgenetries);
void do_swap(struct taxon * first, struct taxon * second);
int swapper(struct taxon * position,struct taxon * prev_pos, int stepstaken, struct taxon * first_swap, struct taxon * second_swap, float * number, int * swaps, int number_of_swaps, int numspectries, int numgenetries);
void yaptp_search(void);
void randomise_tree(char *tree);
void randomise_taxa(char *tree);
void random_star_decom(char *tree);
int check_if_diff_tree(char *tree);
int coding(int nrep, int scoring, int ptpreps);
int MRP_matrix(char **trees, int num_trees, int consensus);
void set_parameters(void);
float MRC(char *supertree);
float quartet_compatibility(char *supertree);
void condense_coding(void);
void reset_spr (struct taxon *position);
int remaining_spr (struct taxon *position);
int spr(struct taxon * position, int maxswaps, int numspectries, int numgenetries);
int regraft(struct taxon * position, struct taxon * newbie, struct taxon * last, int steps, int maxswaps,int numspectries,int numgenetries);
void get_lengths(struct taxon *position);
int xposition1 (struct taxon *position, int count);
float middle_number(struct taxon *position);
void xposition2(struct taxon *position);
int yposition0(struct taxon *position, int level, int deepest);
int yposition1(struct taxon *position, int level);
void yposition2(struct taxon *position, int deepest);
void print_coordinates(struct taxon *position, char **treearray, int taxa_count, int mapping);
void tree_coordinates(char *tree, int bootstrap, int build, int mapping, int fundnum);
float tofloat(char *number);
void generatetrees(void);
void draw_histogram(FILE *outfile, int bins, float *results, int num_results);
void consensus(int num_trees, char **trees, int num_reps, float percentage, FILE *outfile, FILE *guidetreefile);
void input_fund_tree(char *intree, int fundnum);
int nexusparser(FILE *nexusfile);
void do_consensus(void);
int comment(FILE *file);
void showtrees(void);
void exclude(int do_all);
void returntree(char *temptree);
void quick(float ** items, int count);
void qs(float **items, int left, int right);
void include(int do_all);
void exclude_taxa(int do_all);
void sourcetree_dists();
void prune_taxa_for_exclude(struct taxon * super_pos, int *tobeexcluded);
void spr_dist(void);
int string_SPR(char * string); /* carries out random SPR operations on the tree */
void neighbor_joining( int brlens, char *tree, int names);
void nj(void);
void identify_taxa(struct taxon * position, int *name_array);
void reroot_tree(struct taxon *outgroup);
void clean_pointer_taxa(struct taxon *position);
struct taxon * get_branch(struct taxon *position, int name);
float tree_map(struct taxon * gene_top, struct taxon * species_top, int print);
int number_tree1(struct taxon * position, int num);
void label_gene_tree(struct taxon * gene_position, struct taxon * species_top, int *presence, int xnum);
int get_min_node(struct taxon * position, int *presence, int num);
void subtree_id(struct taxon * position, int *tmp);
void descend(struct taxon * position, int *presence);
int reconstruct_map(struct taxon *position, struct taxon *species_top);
void add_losses(struct taxon * position, struct taxon *species_top);
int join_losses(struct taxon * position);
int count_losses(struct taxon * position);
struct taxon * construct_tree(struct taxon * spec_pos, struct taxon *gene_pos, int *presence, struct taxon *extra_gene);
int compress_tree (struct taxon * position);
int compress_tree1 (struct taxon * position);
struct taxon * get_taxon(struct taxon *position, int name);
void duplicate_tree(struct taxon * orig_pos, struct taxon * prev_dup_pos);
void find_tagged(struct taxon * position, int *presence);
void up_tree(struct taxon * position, int *presence);
void down_tree(struct taxon * position, struct taxon *prev, int *presence);
void mapunknowns();
void reconstruct(int print_settings);
void put_in_scores(struct taxon * position, float * total);
void assign_ances_desc(struct taxon *position, int ** allowed_species, int * previous);
void isittagged(struct taxon * position);
void hgt_reconstruction();
void assign_hgtdonors(struct taxon * position, int num, int part_num);
void reset_tag2(struct taxon * position);
int assign_tag2(struct taxon * position, int num);
void assign_before_after(struct taxon *position, int *previous, int *before, int *after, int num, int found);
struct taxon * find_remaining(struct taxon * position);
void exhaustive_SPR(char * string);
void print_tree_labels(struct taxon *position, int **results, int treenum, struct taxon *species_tree);
int count_internal_branches(struct taxon *position, int count);
float get_recon_score(char *giventree, int numspectries, int numgenetries);
void print_descendents(struct taxon *position, FILE *outfile);
void do_descendents(struct taxon *position, FILE *outfile);
void resolve_tricotomies(struct taxon *position, struct taxon *species_tree);
void gene_content_parsimony(struct taxon * position, int * array);
struct taxon * do_resolve_tricotomies(struct taxon * gene_tree, struct taxon * species_tree, int basescore);
int presence_of_trichotomies(struct taxon * position);
int are_siblings(struct taxon *position, int first, int second);
int isit_onetoone(struct taxon *position, int onetoone);
void print_onetoone_names(struct taxon *position, int onetoone);
int get_best_node(struct taxon * position, int *presence, int num);

void random_prune(char *fund_tree);
void collapse_clades(struct taxon * position, float user_limit, int * to_delete, FILE *rp_outfile);
int get_brlens(struct taxon * position, float *total, int *count);
float return_length(char *string);
int untag_taxa(struct taxon *position, int * to_delete, int keep, int count, FILE *rp_outfile);
int print_keep(struct taxon *position, int keep, int count, FILE *rp_outfile);
void resolve_tricotomies_dist (struct taxon * gene_tree, struct taxon *species_tree, int ** scores);
void get_taxa(struct taxon *position, int *presence);
void check_treeisok(struct taxon *position);
void pathmetric_internals(char *string, struct taxon * species_tree, int **scores);
void calculate_withins(struct taxon *position, int **within, int *presence);
long extract_length(char * fullname);
long list_taxa_in_clade(struct taxon * position, int * foundtaxa, struct taxon * longest, long seqlength); /* descend through the tree finding what taxa are there (and putting result into an array) and also identifying the longest sequence (the first number in the <<full>> name of the sequence, after the first "." and before the first "|") */
long list_taxa_above(struct taxon * position, int * foundtaxa, struct taxon * longest, long seqlength);
void identify_species_specific_clades(struct taxon * position);
void prune_monophylies();
void untag_nodes_below(struct  taxon * position);
void untag_nodes_above(struct  taxon * position);
void tips(int num);


int spr_new(struct taxon * master, int maxswaps, int numspectries, int numgenetries);




void controlc1(int signal);
void controlc2(int signal);
void controlc3(int signal);
void controlc4(int signal);
void controlc5(int signal);

/****************** Global variable definitions ****************/
FILE * infile = NULL, *BR_file = NULL, *psfile = NULL, *logfile = NULL, *distributionreconfile = NULL, *onetoonefile = NULL, *strictonetoonefile = NULL;
char **taxa_names = NULL, ***fulltaxanames = NULL, **parsed_command = NULL, **fundamentals = NULL, **stored_funds = NULL, **retained_supers = NULL, **stored_commands = NULL, *tempsuper = NULL, **best_topology = NULL, **tree_names = NULL;
int  *numtaxaintrees = NULL, fullnamesnum = 0, fullnamesassignments = 1, fundamental_assignments = 0, tree_length_assignments = 1, parsed_command_assignments = 1, name_assignments = 0, *taxa_incidence = NULL, number_of_taxa = 0, Total_fund_trees = 0, *same_tree = NULL, **Cooccurrance = NULL, NUMSWAPS = 0;
int ***fund_scores = NULL, ***stored_fund_scores = NULL, **super_scores = NULL, *number_of_comparisons = NULL, *stored_num_comparisons = NULL, **presence_of_taxa = NULL, **stored_presence_of_taxa = NULL, *presenceof_SPRtaxa = NULL;
int num_commands = 0, number_retained_supers = 10, number_of_steps = 3, largest_tree = 0, smallest_tree = 1000000, criterion = 0, parts = 0, **total_coding = NULL, total_nodes = 0, quartet_normalising = 3, splits_weight = 2, dweight =1, *from_tree = NULL, method = 2, tried_regrafts = 0, hsprint = TRUE, max_name_length = NAME_LENGTH, got_weights = FALSE, num_excluded_trees = 0, num_excluded_taxa = 0, calculated_fund_scores = FALSE, select_longest=FALSE;
struct taxon *tree_top = NULL, *temp_top = NULL, *temp_top2 = NULL, *branchpointer = NULL, *longestseq = NULL;
float *scores_retained_supers = NULL, *partition_number = NULL, num_partitions = 0, total_partitions = 0, sprscore = -1, *best_topology_scores = NULL, **weighted_scores = NULL, *sourcetree_scores = NULL, *tree_weights = NULL;
float *score_of_bootstraps = NULL, *yaptp_results = NULL, largest_length = 0, dup_weight = 1, loss_weight = 1, hgt_weight = 1, BESTSCORE = -1;
time_t interval1, interval2;
double sup=1;
char saved_supertree[TREE_LENGTH],  *test_array, inputfilename[100];
int trees_in_memory = 0, *sourcetreetag = NULL, remainingtrees = 0, GC, user_break = FALSE, delimiter = FALSE, print_log = FALSE, num_gene_nodes, testarraypos = 0, taxaorder=0;
int malloc_check =0, count_now = FALSE, another_check =0;


int main(int argc, char *argv[])
    {
	
    int i = 0, j=0, k=0, l=0, m=0, error=FALSE, x, doexecute_command = FALSE, command_line = FALSE, tipnum=0;
    char *command = NULL, HOME[1000], PATH[1000], exefilename[1000];
    time_t time1, time2;
    double diff=0;    
    FILE *tmpclann = NULL;
	
	exefilename[0] = '\0';
	inputfilename[0] = '\0';
    saved_supertree[0] = '\0';
	test_array = malloc(TREE_LENGTH*sizeof(int));
	test_array[0] = '\0';
    
	if(argc > 1)
		{
		for(i=1; i<argc; i++)
			{
			if(strcmp(argv[i], "-n") == 0)
				command_line = TRUE;
			else
				{
				doexecute_command = TRUE;
				strcpy(exefilename, argv[i]);
				}
			}
		}
	
	/********** START OF READLINE STUFF ***********/
	#ifdef HAVE_READLINE    
	if(!command_line)
		{
		PATH[0] ='\0'; HOME[0] = '\0'; 
		
		system("echo $HOME > clann5361temp1023");
		tmpclann = fopen("clann5361temp1023", "r");
		fscanf(tmpclann, "%s", HOME);
		fclose(tmpclann);
		remove("clann5361temp1023");
		strcpy(PATH, HOME);
		strcat(PATH, "/.clannhistory");
		read_history(PATH);
		}
    #endif
    /********** END OF READLINE STUFF *************/
    
    
    time1 = time(NULL);
    /*seed the rand number with the calander time + the PID of the process ---- used for bootstrapping and others */
    #ifdef HAVE_READLINE
    srand((unsigned) (time(NULL)/2)+getpid());
    #else
    srand((unsigned) (time(NULL)/2));
    #endif

    command = malloc(10000*sizeof(char));
    if(!command) memory_error(74);
    command[0] = '\0';

    tempsuper = malloc(TREE_LENGTH*sizeof(char));
    tempsuper[0] = '\0';

    /** allocate the array to store multiple commands that are sepatated by ";" on the commandline */

    stored_commands = malloc(100*sizeof(char *));
    if(!stored_commands) memory_error(62);
    for(i=0; i<100; i++)
        {
        stored_commands[i] = malloc(10000*sizeof(char));
        if(!stored_commands[i]) memory_error(63);
        stored_commands[i][0] = '\0';
        }


    /* allocate the array to hold the retained supertrees */
    retained_supers = malloc(number_retained_supers*sizeof(char *));
    if(!retained_supers) memory_error(45);
    for(i=0; i<number_retained_supers; i++)
        {
        retained_supers[i] = malloc(TREE_LENGTH*sizeof(char));
        if(!retained_supers[i]) memory_error(46);
        retained_supers[i][0] = '\0';
        }
    scores_retained_supers = malloc(number_retained_supers*sizeof(float));
    if(!scores_retained_supers) memory_error(47);
    for(i=0; i<number_retained_supers; i++)
        scores_retained_supers[i] = -1;
    
    /* assign the parsed_command array */
    
    parsed_command = malloc(10000*sizeof(char *));
    for(i=0; i<1000; i++)
        {
        parsed_command[i] = malloc(10000*sizeof(char));
        parsed_command[i][0] = '\0';
        }

    printf("\n\n\n\n\n\t******************************************************************");
    printf("\n\t*                                                                *");
    printf("\n\t*                          Clann  v4.1.3                         *");
    printf("\n\t*                                                                *");
    printf("\n\t*                 web: http://www.creeveylab.org                 *");
    printf("\n\t*                 email: chris.creevey@gmail.com                 *");
    printf("\n\t*                                                                *");
    printf("\n\t*                Copyright Chris Creevey 2003-2016               *");
    printf("\n\t*                                                                *");
    printf("\n\t*         HINT: Type \"help\" to see all available commands        *");
    printf("\n\t******************************************************************\n\n");

     
    if(doexecute_command) /* if the user has specified an input file at the command line */
        {
        execute_command(exefilename, TRUE);
        }
    
    i=0;
    while(strcmp(parsed_command[0], "quit") != 0)  /* while the user has not chosen to quit */
        {
        
        if(num_commands > 0)
            {
            if(strcmp(parsed_command[0], "execute") == 0 || strcmp(parsed_command[0], "exe") == 0)
                {
                if(num_commands == 2 && parsed_command[1][0] == '?')
                    print_commands(1);
                else
                    {
                
                    if(num_commands > 1)
                        {
                        delimiter = FALSE;
                        execute_command(parsed_command[1], TRUE);
							
						/* test fulltaxaname allocations */
						/*for(l=0; l<Total_fund_trees; l++)
							{
							printf("fund tree number: %d\n", l);
							for(m=0; m<numtaxaintrees[l]; m++)
								{
								printf("\ttaxa %d: %s\n", m, fulltaxanames[l][m]);
								}
							}
							printf("here\n"); */
							
						/*spr_dist(); */ 
                        }
                    else
						{
                        printf("Error: no file name specified\n");
						if(print_log) fprintf(logfile, "Error: no file name specified\n");
						}
                    }
                }
            else
                {
                if(strcmp(parsed_command[0], "help") == 0)
                    {
                    print_commands(0);
                    }
                else	
                    {
                    if(strcmp(parsed_command[0], "alltrees") == 0)
                        {
                        if(num_commands == 2 && parsed_command[1][0] == '?')
                            print_commands(2);
                        else
                            {
                            if(number_of_taxa > 0)
                                {                                        
                                if(criterion == 0 || criterion == 2 || criterion == 3)    /** MRP or MSS or QC **/
                                    alltrees_search(TRUE);
                                else
                                    {
                                    BR_file = fopen("coding.nex", "w");
                                    x = coding(0, 0, 0);
                                    fclose(BR_file);
									if(x == FALSE)
										{
										if(system("paup coding.nex") != 0) printf("Error calling paup, please execute the file coding.nex in paup to perform the parsimony step\n");
										}
									if(x == FALSE)
										{
										remove("clanntmp.chr");
										remove("clanntree.chr");
										/*pars("coding.nex", "clanntmp.chr"); */
										remove("clanntmp.chr");
										remove("clanntree.chr");
										}
                                    }
                                }
                            else
                                {
                                printf("Error: You need to load source trees before using this command\n");
                                }
                            }
                        }
                    else
                        {
                        if(strcmp(parsed_command[0], "usertrees") == 0)
                            {
                            if(num_commands == 2 && parsed_command[1][0] == '?')
                                print_commands(3);
                            else
                                {
                                if(number_of_taxa > 0)
                                    {
                                    if(criterion == 1)
                                        printf("Error: You cannont do a user tree search in MRP\n");
                                    else
                                        usertrees_search();
                                    }
                                else
                                    {
                                    printf("Error: You need to load source trees before using this command\n");
                                    }
                                }
                            }
                        else
                            {
                            if(strcmp(parsed_command[0], "hs") == 0 || strcmp(parsed_command[0], "hsearch") == 0)
                                {
                                if(num_commands == 2 && parsed_command[1][0] == '?')
                                    print_commands(4);
                                else
                                    {
                                    if(number_of_taxa > 0)
                                        {
                                        heuristic_search(TRUE, TRUE, 10000, 10);
										
                                        }
                                    else
                                        {
                                        printf("Error: You need to load source trees before using this command\n");
                                        }
                                    }
                                }
                            else
                                {
                                if(strcmp(parsed_command[0], "bootstrap") == 0 || strcmp(parsed_command[0], "boot") == 0)
                                    {
                                    if(num_commands == 2 && parsed_command[1][0] == '?')
                                        print_commands(5);
                                    else
                                        {
                                        if(number_of_taxa > 0)
                                            {
                                            bootstrap_search(); 
                                            }
                                        else
                                            {
                                            printf("Error: You need to load source trees before using this command\n");
                                            }
                                        }
                                    }
                                else
                                    {
                                    if(strcmp(parsed_command[0], "yaptp") == 0)
                                        {
                                        if(num_commands == 2 && parsed_command[1][0] == '?')
                                            print_commands(6);
                                        else	
                                            {
                                            if(number_of_taxa > 0)
                                                {
                                                yaptp_search(); 
                                                }
                                            else
                                                {
                                                printf("Error: You need to load source trees before using this command\n");
                                                }
                                            }
                                        }
                                    else
                                        {
                                        if(strcmp(parsed_command[0], "ideal") == 0)
                                            {
                                            if(num_commands == 2 && parsed_command[1][0] == '?')
                                                print_commands(7);
                                            else
                                                {
                                                if(number_of_taxa > 0)
                                                    {
                                                    /* ideal_search(); */
                                                    }
                                                else
                                                    {
                                                    printf("Error: You need to load source trees before using this command\n");
                                                    }
                                                }
                                            }
                                        else
                                            {
                                            if(strcmp(parsed_command[0], "set") == 0)
                                                {
                                                if(num_commands == 2 && parsed_command[1][0] == '?')
                                                    print_commands(8);
                                                else
                                                    set_parameters(); 
                                                }
                                            else
                                                {
                                                if(strcmp(parsed_command[0], "quit") == 0)
                                                    {
                                                    if(num_commands == 2 && parsed_command[1][0] == '?')
                                                        {
                                                        print_commands(13);
                                                        strcpy(parsed_command[0], "hi");
                                                        }
                                                    }
                                                else
                                                    {
                                                    if(strcmp(parsed_command[0], "deletetrees") == 0)
                                                        {
                                                        if(num_commands == 2 && parsed_command[1][0] == '?')
                                                            print_commands(17);
                                                        else
															{
															if(number_of_taxa > 0)
																exclude(TRUE);
															else
																printf("Error: You need to load source trees before using this command\n");
															}
														}
                                                    else
                                                        {
                                                        if(strcmp(parsed_command[0], "1includetrees") == 0)
                                                            {
                                                            if(num_commands == 2 && parsed_command[1][0] == '?')
                                                                print_commands(18);
                                                            else
																{
																if(number_of_taxa > 0)
																		include(TRUE);
																else
																	printf("Error: You need to load source trees before using this command\n");
																}
                                                            }
                                                        else
                                                            {
                                                            if(strcmp(parsed_command[0], "!") == 0)
                                                                {
                                                                if(num_commands == 2 && parsed_command[1][0] == '?')
                                                                    print_commands(12);
                                                                else
                                                                    {
                                                                    printf("\n\tType exit to return to Clann\n\n");
                                                                    system("csh");
                                                                    }
                                                                }
															else
																{
																if(strcmp(parsed_command[0], "generatetrees") == 0)
																	{
																	if(num_commands == 2 && parsed_command[1][0] == '?')
																			print_commands(14);
																	else
																		{
																		if(number_of_taxa > 0)
																				generatetrees();
																		else
																			printf("Error: You need to load source trees before using this command\n");
																		}
																	}
																else
																	{
																	if(strcmp(parsed_command[0], "showtrees") == 0)
																		{
																		if(num_commands == 2 && parsed_command[1][0] == '?')
																			print_commands(16);
																		else
																			{
																			if(number_of_taxa > 0)
																					showtrees();
																			else
																				printf("Error: You need to load source trees before using this command\n");
																			}
																		}
																	else
																		{
																		if(strcmp(parsed_command[0], "rfdists") == 0)
																			{
																			if(num_commands == 2 && parsed_command[1][0] == '?')
																				print_commands(19);
																			else
																				{
																				if(number_of_taxa > 0)
																					sourcetree_dists();
																				else
																					printf("Error: You need to load source trees before using this command\n");
																				}
																			}
																		else
																			{
																			if(strcmp(parsed_command[0], "deletetaxa") == 0)
																				{
																				if(num_commands == 2 && parsed_command[1][0] == '?')
																					print_commands(20);
																				else
																					{
																					if(number_of_taxa > 0)
																						{
																						exclude_taxa(TRUE);
																						}
																					else
																						printf("Error: You need to load source trees before using this command\n");
																					}
																				}

																			else
																				{
																				if(strcmp(parsed_command[0], "consensus") == 0)
																					{
																					if(num_commands == 2 && parsed_command[1][0] == '?')
																						print_commands(15);
																					else
																						{
																						if(number_of_taxa > 0)
																							do_consensus();
																						else
																							printf("Error: You need to load source trees before using this command\n");
																						}
																					}
																				else
																					{
																					if(strcmp(parsed_command[0], "nj") == 0)
																						{
																						if(num_commands == 2 && parsed_command[1][0] == '?')
																							print_commands(21);
																						else
																							{
																							if(number_of_taxa > 0)
																								nj();
																							else
																								printf("Error: You need to load source trees before using this command\n");
																							}
																						}
																					else
																						{
																						if(strcmp(parsed_command[0], "sprdists") == 0)
																							{
																							if(num_commands == 2 && parsed_command[1][0] == '?')
																								print_commands(22);
																							else
																								{
																								if(number_of_taxa > 0)
																									spr_dist();
																								else
																									printf("Error: You need to load source trees before using this command\n");	
																								}
																							}
																						else
																							{
																							if(strcmp(parsed_command[0], "mapunknowns") == 0)
																								{
																								if(num_commands == 2 && parsed_command[1][0] == '?')
																									print_commands(23);
																								else
																									{
																									if(number_of_taxa > 0)
																										mapunknowns();
																									else
																										printf("Error: You need to load source trees before using this command\n");	
																									}
																								}
																							else
																								{
																								if(strcmp(parsed_command[0], "reconstruct") == 0)
																									{
																									if(num_commands == 2 && parsed_command[1][0] == '?')
																										print_commands(24);
																									else
																										{
																										if(number_of_taxa > 0)
																											reconstruct(TRUE);
																										else
																											printf("Error: You need to load source trees before using this command\n");	
																										}
																									}
																								else
																									{
																									if(strcmp(parsed_command[0], "hgtanalysis") == 0)
																										{
																										if(num_commands == 2 && parsed_command[1][0] == '?')
																											print_commands(23);
																										else
																											{
																											if(number_of_taxa > 0)
																												hgt_reconstruction();
																											else
																												printf("Error: You need to load source trees before using this command\n");	
																											}
																										}
																									else
																										{
																										if(strcmp(parsed_command[0], "exhaustivespr") == 0)
																											{
																											if(num_commands == 2 && parsed_command[1][0] == '?')
																												print_commands(23);
																											else
																												{
																												if(number_of_taxa > 0)
																													exhaustive_SPR(fundamentals[0]);
																												else
																													printf("Error: You need to load source trees before using this command\n");	
																												}
																											}
																										else
																											{
																											if(strcmp(parsed_command[0], "pars") == 0)
																												{
																												if(num_commands == 2 && parsed_command[1][0] == '?')
																													print_commands(23);
																												else
																													{
																													if(number_of_taxa > 0)
																														{
																														remove("clanntmp.chr");
																														remove("clanntree.chr");
																														/*pars();*/
																														remove("clanntmp.chr");
																														remove("clanntree.chr");
																														
																														}
																													else
																														printf("Error: You need to load source trees before using this command\n");	
																													}
																												}
																											else
																												{	
																												if(strcmp(parsed_command[0], "randomisetrees") == 0)
																													{
																													if(num_commands == 2 && parsed_command[1][0] == '?')
																														print_commands(27);
																													else
																														{
																														if(number_of_taxa > 0)
																															{
																															for(j=0; j<Total_fund_trees; j++)
																																{
																																randomise_tree(fundamentals[j]); /*equiprobable */
																																}
																															printf("The source trees have now been randomised\n");
																															}
																														else
																															printf("Error: You need to load source trees before using this command\n");	
																														}
																													}
																												else
																													{
																													if(strcmp(parsed_command[0], "randomprune") == 0)
																														{
																														if(num_commands == 2 && parsed_command[1][0] == '?')
																															print_commands(28);
																														else
																															{
																															if(number_of_taxa > 0)
																																{
																																for(j=0; j<Total_fund_trees; j++)
																																	{
																																	random_prune(fundamentals[j]); /*equiprobable */
																																	}
																																
																																}
																															else
																																printf("Error: You need to load source trees before using this command\n");	
																															}
																														}
																													else
                                                                                                                        {
                                                                                                                        if(strcmp(parsed_command[0], "prunemonophylies") == 0)
                                                                                                                            {
                                                                                                                            if(num_commands == 2 && parsed_command[1][0] == '?')
                                                                                                                                print_commands(25);
                                                                                                                            else
                                                                                                                                {
                                                                                                                                if(number_of_taxa > 0)
                                                                                                                                    {
                                                                                                                                    prune_monophylies();
                                                                                                                                    }
                                                                                                                                else
                                                                                                                                    printf("Error: You need to load source trees before using this command\n"); 
                                                                                                                                }
                                                                                                                            }
                                                                                                                        else
                                                                                                                        	{
                                                                                                                        	if(strcmp(parsed_command[0], "tips") == 0)
                                                                                                                            	{                                     
                                                                                                                            	if(num_commands == 2 && parsed_command[1][0] == '?')
	                                                                                                                                print_commands(26);
	                                                                                                                            else
	                                                                                                                            	{
	                                                                                                                            	if(num_commands == 2 && (tipnum = atoi(parsed_command[1])) > 0 && tipnum < 11)	
	                                                                                                                                	tips(tipnum-1);
	                                                                                                                            	else
	                                                                                                                            		tips((int)fmod(rand(), 10));
	                                                                                                                            	}
	                                                                                                                        	}
	                                                                                                                        else
	                                                                                                                        	printf("Error: command not known.\n\tType help at the prompt to get a list of available commands.\n");
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
                            }
                        }
                    }
                }
            } /* end if*/

        time2 = time(NULL);
        diff = difftime(time2, time1);
        if(diff > 0)
            {
            printf("\nTime taken:"); 
            if(diff > 60)
                {
                if(diff > 3600)
                    {
                    if(diff > 86400)
                        {
                        printf(" %0.0lf Day", diff/86400);
                        if(diff/86400 > 1) printf("s");
                        diff = diff- (86400 * ((int)diff/86400));
                        }
                    printf(" %0.0lf Hour", diff/3600);
                    if(diff/3600 > 1) printf("s");
                    diff = diff - (3600 * ((int)diff/3600));
                    }
                printf(" %0.0lf Minute", diff/60);
                if(diff/60 > 1) printf("s");
                diff = diff - (60 * ((int)diff/60));
                }
            printf(" %0.0lf second", diff);
            if(diff > 1) printf("s");
            printf("\n\n");
            }

        if(i == parts)
            {
            parts =0;
			#ifdef HAVE_READLINE
/****/      free(command);
			#endif
			user_break = FALSE;
			if(signal(SIGINT, controlc4) == SIG_ERR)
				{
				printf("An error occurred while setting a signal handler\n");
				}
		   #ifdef HAVE_READLINE
/*****/    command = readline("clann> "); 
/*****/    command = realloc(command, 10000*sizeof(char));
		   #else
 	        printf("clann> ");
 	        fflush(stdout);
 	        command = xgets(command);
		   #endif
          
           
            /* if the commmand doesn't end in a ";" then put one on */
            k=0;
            while(command[k] != '\0') k++;
            if(command[k-1] != ';')
                strcat(command, ";");  
            
			#ifdef HAVE_READLINE
            if(strcmp(command, ";") != 0 && command_line == TRUE)   /* if the command is not just a blank line */
                {
/******/                add_history(command);
                }
			#endif
            parts = seperate_commands(command); /* if there is more than one command, break each of them up */
            i=0;
			#ifdef HAVE_READLINE
/****/         free(command); 
            command = malloc(10000*sizeof(char));
            command[0] = '\0';
/*****/		#endif


                       
            }
        for(j=0; j<1000; j++)
            {
            parsed_command[j][0] = '\0';
            }
        command[0] = '\0';
        strcpy(command, stored_commands[i]);
        num_commands = parse_command(command);  /* parse each individual part of the command */
        i++;
        time1 = time(NULL);
    
        } /* end while */
    free(command);
    free(tempsuper);
	#ifdef HAVE_READLINE
/****/ if(command_line == TRUE) write_history(PATH);  /* write the history of commands to file */
	#endif
/*	printf("malloc_check = %d x (%d)\n", malloc_check, sizeof(taxon_type));
 */   
	clean_exit(0);
    return(0);
    }


int seperate_commands(char *command)
    {
    int i=0, j=0, numcommands = 0;
    /* each seperate command is to be seperated by a ";" */
    
    /*initialise the stored command array */
    for(i=0; i<100; i++)
        strcpy(stored_commands[i], "");
    
    /* run down the command looking for ";" */
    i=0;
    while(command[i] != '\0')
        {
        j=0;
        while(command[i] != '\0' && command[i] != ';')
            {
            stored_commands[numcommands][j] = command[i];
            i++;j++;
            }
        if(command[i] == ';')   /* this both makes sure 'i' doesn't fall off the array and we disregard any command that doesn't end in a ';' */
            {
            stored_commands[numcommands][j] = ';';
            stored_commands[numcommands][j+1] = '\0';
            numcommands++;
            i++;
            }
        }
    return(numcommands);
    }


 
    
void print_commands(int num)
    {
    
    if(num == 0)
        {
        printf("\nAvailable Commands:\n\n");
        printf("\nThe following commands are always available:\n\n");
        printf("\texecute\t\t- Read in a file of source trees\n");
        printf("\thelp\t\t- Display this message\n");
        printf("\tquit\t\t- Quit Clann\n");
        printf("\tset\t\t- Set the optimality criterion for carrying reconstructing a supertree\n");
        printf("\t!\t\t- Run a shell session, while preserving the current Clann session (type \'exit\' to return)\n");
        printf("\ttips\t\t- Show tips and hints for better use of Clann\n");

        printf("\nThe following commands are only available when there are source trees in memory:\n");

        printf("\nSupertree reconstruction:\n");
        printf("\ths\t\t- Carry out a heuristic search for the best supertree usign the criterion selected\n");
        printf("\tbootstrap\t- Carry out a bootstrap supertree analysis using the criterion selected\n");
        printf("\tnj\t\t- Construct a neighbour-joining supertree\n");
        printf("\talltrees\t- Exhaustively search all possible supertrees\n");
        printf("\tusertrees\t- Assess user-defined supertrees (from seperate file), to find the best scoring\n");
        printf("\tconsensus\t- Calculate a consensus tree of all trees containing all taxa\n");

        printf("\nSource tree selection and modification:\n");
        printf("\tshowtrees\t- Visualise selected source trees in ASCII format (also can save selected trees to file)\n");
       printf("\tdeletetrees\t- Specify source trees to delete from memory (based on a variety of criteria)\n"); 
      /*  printf("\tincludetrees\t- Specify trees for inclusion in the analysis (based on a variety of criteria)\n"); */ /* we are excluding the include command because it is problematic*/
        printf("\tdeletetaxa\t- Specify taxa to delete from all source trees in memory (i.e. prune from the trees while preserving branch lengths)\n");
        printf("\trandomisetrees\t- Randomises the source trees in memory, while preserving taxa composition in each tree\n");

        printf("\nMiscellaneous calculations:\n");
        printf("\trfdists\t\t- Calculate Robinson-Foulds distances between all source trees\n");
        printf("\tgeneratetrees\t- Generate random supertrees & assess  against source trees in memory\n");
        printf("\tyaptp\t\t- \"Yet another permutation-tail-probability\" test - performs a randomisation test\n");



        printf("\nExperimental Options:\n");
        printf("\treconstruct\t- Carry out a gene-tree reconciliation (source trees against a species tree)\n");
        printf("\tprunemonophylies - Prunes clades which consist of multiple sequences from the same species, to a single representative\n");
        printf("\tsprdists\t- Carry out estimation of SPR distances of real data versus ideal and randomised versions of the data\n");

        printf("\n\n\nType a command followed by '?' to get information on the options available i.e.: \"exe ?\"\n");
        printf("Full descriptions of the commands are available in the manual\n\n\n");


        }
   
        

    if(num == 1)
        {
        printf("\nexecute\t<filename>\n");
        printf("\nexe\t<filename>\n\n");
		 printf("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf("\t===========================================================\n");
		printf("\n\tmaxnamelen\t<integer number> | delimited\t*none");
		printf("\n\tsummary\t\tshort | long\t\t\t*long");
        }
    if(num == 2)
        {
        printf("\nalltrees [options]\n\n");
        printf("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf("\t===========================================================\n");
        if(criterion == 0 || criterion == 2 || criterion == 3)
            {
            printf("\n\trange\t\t<treenumber> - <treenumber>\t*all\n\tsavetrees\t<filename>\t\t\ttop_alltrees.txt\n\tcreate\t\tyes | no\t\t\t*no\n");
            if(criterion == 0)
                {
                printf("\tweight\t\tequal | comparisons\t\t");
                if(dweight == 1) printf("comparisons\n");else printf("euqal\n");
                }
            if(criterion == 2)
                {
                printf("\tweight\t\tequal | splits\t\t\t");
                if(splits_weight == 1) printf("equal\n");if(splits_weight == 2) printf("splits\n");
                }
            if(criterion == 3)
                {
                printf("\tweight\t\tequal | taxa | quartets\t\t");
                if(quartet_normalising == 1) printf("equal\n");if(quartet_normalising == 2) printf("taxa\n");if(quartet_normalising == 3) printf("quartets\n");
                }
            }
        
        }
    if(num == 3)
        {
        printf("\nusertrees <filename> [options]\n\n");
        printf("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf("\t===========================================================\n");
        printf("\n\toutfile\t\t<filename>\t\t\tUsertrees_result.txt");
        if(criterion == 0)
            {
            printf("\n\tweight\t\tequal | comparisons\t\t");
            if(dweight == 1) printf("comparisons\n");else printf("splits\n");
            }
        if(criterion == 2)
            {
            printf("\n\tweight\t\tequal | splits\t\t\t");
            if(splits_weight == 1) printf("equal\n");if(splits_weight == 2) printf("splits\n");
            }
        if(criterion == 3)
            {
            printf("\n\tweight\t\tequal | taxa | quartets\t\t");
            if(quartet_normalising == 1) printf("equal\n");if(quartet_normalising == 2) printf("taxa\n");if(quartet_normalising == 3) printf("quartets\n");
            }
		printf("\n\tprintscourcescores\tyes | no\t\t*no\n");
        }
    if(num == 4)
        {
        printf("\nhs (or hsearch) [options]  \n");
        printf("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf("\t===========================================================\n");
        if(criterion == 0 || criterion == 2 || criterion == 3 || criterion==5)
            {
            printf("\n\tsample\t\t<integer number>\t\t*10,000\n\tnreps\t\t<integer number>\t\t*10");
            printf("\n\tswap\t\tnni | spr | tbr\t\t\t");
            if(method == 1) printf("nni");
            if(method == 2) printf("spr"); 
            printf("\n\tnsteps\t\t<integer number>\t\t%d", number_of_steps);
			printf("\n\tstart\t\tnj | random | <filename>\tnj");

			printf("\n\tmaxswaps\t<integer number>\t\t*1,000,000\n\tsavetrees\t<filename>\t\t\tHeuristic_result.txt");
            if(criterion == 0)
                {
                printf("\n\tweight\t\tequal | comparisons\t\t");
                if(dweight == 1) printf("comparisons");else printf("splits");
                }
            if(criterion == 2)
                {
                printf("\n\tweight\t\tequal | splits\t\t\t");
                if(splits_weight == 1) printf("equal");if(splits_weight == 2) printf("splits");
                }
            if(criterion == 3)
                {
                printf("\n\tweight\t\tequal | taxa | quartets\t\t");
                if(quartet_normalising == 1) printf("equal");if(quartet_normalising == 2) printf("taxa");if(quartet_normalising == 3) printf("quartets");
                }
			if(criterion == 0)
				{
				printf("\n\tdrawhistogram\tyes | no\t\t\t*no\n\tnbins\t\t<integer number>\t\t*20\n\thistogramfile\t<filename>\t\t\t*Heuristic_histogram.txt");
				}
            
            }
        if(criterion == 1)
            {
            printf("\n\tanalysis\tparsimony | nj\t\t\t*parsimony\n");
            printf("\n\tParsimony options:\n\tweighted\tyes | no\t\t\t*no\n\tswap\t\tnni | spr | tbr\t\t\t*tbr\n\taddseq\t\tsimple | closest | asis |\n\t\t\trandom | furthest\t\t*random\n\tnreps\t\t<integer number>\t\t*10\n\n\tGeneral Options:\n\tsavetrees\t<filename>\t\t\tMRP.tree\n");
            }
		if(criterion == 4)
			{
			printf("\n\tmissing\t4point | ultrametric\t\t\t*4point\n");
			}
		if(criterion == 5)
			{
			printf("\n\tduplications\t<value>\t\t\t\t*1.0\n\tlosses\t\t<value>\t\t\t\t*1.0");
			printf("\n\tnumspeciesrootings\t<value> | all\t\t*2\n\tnumgenerootings\t\t<value> | all\t\t*2\n");
			}
        }
    if(num == 5)
        {
        printf("\nbootstrap (or boot) [options]  \n\n\n");
        printf("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf("\t===========================================================\n");
        if(criterion != 1)
            {
            printf("\n\tnreps\t\t<integer number>\t\t*100\n\thsreps\t\t<integer number>\t\t*10\n\tsample\t\t<integer number>\t\t*10,000\n\tswap\t\tnni | spr | tbr | all\t\t");
            if(method == 1) printf("nni");
            if(method == 2) printf("spr");
            printf("\n\tstart\t\trandom | <filename>\t\trandom");
            
            printf("\n\tnsteps\t\t<integer number>\t\t%d\n\ttreefile\t<output treefile name>\t\tbootstrap.txt\n\tmaxswaps\t<integer number>\t\t*1,000,000\n", number_of_steps);
            if(criterion == 0)
                {
                printf("\tweight\t\tequal | comparisons\t\t");
                if(dweight == 1) printf("comparisons");else printf("splits");
                }
            if(criterion == 2)
                {
                printf("\tweight\t\tequal | splits\t\t\t");
                if(splits_weight == 1) printf("equal");if(splits_weight == 2) printf("splits");
                }
            if(criterion == 3)
                {
                printf("\tweight\t\tequal | taxa | quartets\t\t");
                if(quartet_normalising == 1) printf("equal");if(quartet_normalising == 2) printf("taxa");if(quartet_normalising == 3) printf("quartets");
                }
			if(criterion == 5)
				{
				printf("\n\tduplications\t<value>\t\t\t\t*1.0\n\tlosses\t\t<value>\t\t\t\t*1.0");
				printf("\n\tnumspeciesrootings\t<value> | all\t\t*2\n\tnumgenerootings\t\t<value> | all\t\t*2\n");
				}

			printf("\n\tconsensus\tstrict | majrule | minor | <proportion>\t*majrule");
			printf("\n\tconsensusfile\t<filename>\t\t\tconsensus.ph\n");
            printf("\n");
            }
        if(criterion == 1)
            printf("\n\tnreps\t\t<integer number>\t\t*100\n\tswap\t\tnni | spr | tbr | all\t\t*tbr\n\taddseq\t\tsimple | closest | asis |\n\t\t\trandom | furthest\t\t*random\n\ttreefile\t<output treefile name>\t\tMRP.tree\n");
        }
    if(num == 6)
        {
        printf("\nyaptp [options]  \n\n\n");
        printf("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf("\t===========================================================\n");
        if(criterion != 1)
            {
            printf("\n\tmethod\t\tequiprobable | markovian\t*equiprobable");
            printf("\n\tnreps\t\t<integer number>\t\t*100\n\thsreps\t\t<integer number>\t\t*10\n\tsample\t\t<integer number>\t\t*10,000\n\tsearch\t\tnni | spr | all\t\t\t");
            if(method == 1) printf("nni");
            if(method == 2) printf("spr");
            printf("\n\tnsteps\t\t<integer number>\t\t%d\n\ttreefile\t<output treefile name>\t\tyaptp.ph\n\tmaxswaps\t<integer number>\t\t*1,000,000\n", number_of_steps);
            if(criterion == 0)
                {
                printf("\tweight\t\tequal | comparisons\t\t");
                if(dweight == 1) printf("comparisons");else printf("splits");
                }
            if(criterion == 2)
                {
                printf("\tweight\t\tequal | splits\t\t\t");
                if(splits_weight == 1) printf("equal");if(splits_weight == 2) printf("splits");
                }
            if(criterion == 3)
                {
                printf("\tweight\t\tequal | taxa | quartets\t\t");
                if(quartet_normalising == 1) printf("equal");if(quartet_normalising == 2) printf("taxa");if(quartet_normalising == 3) printf("quartets");
                }
	if(criterion == 5)
		{
		printf("\n\tduplications\t<value>\t\t\t\t*1.0\n\tlosses\t\t<value>\t\t\t\t*1.0");
		printf("\n\tnumspeciesrootings\t<value> | all\t\t*2\n\tnumgenerootings\t\t<value> | all\t\t*2\n");
		}

            printf("\n");
            }
        else
            {
            printf("\n\tnreps\t\t<integer number>\t\t*100\n");
            }
            
        }
    if(num == 7)
        {
        
        }
    if(num == 8)
        {
        printf("\nset [options]  \n\n\n");
        printf("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf("\t===========================================================\n");
        printf("\n\tcriterion\tdfit | sfit | qfit | mrp | avcon\t");
        if(criterion == 0) printf("dfit");
        if(criterion == 1) printf("mrp");
        if(criterion == 2) printf("sfit");
        if(criterion == 3) printf("qfit");
        if(criterion == 4) printf("avcon");
        printf("\n");
/*        printf("\n\n\t\t\tdfit = best Distance Fit\n\t\t\tsfit = maximum Splits Fit\n\t\t\tqfit = maximum Quartet Fit\n\t\t\tmrp = Matrix representation using parsimony\n");
  */      }
    if(num == 9)
        printf("\nquit \n\tquit from the program and return to the operating system\n");
    if(num == 10)
        {
        printf("\nexclude [options]  \nNOT IMPLEMENTED YET\n\n");
        printf("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf("\t===========================================================\n");

        printf("\n\tall\t\tyes | no\t\t\tno\n\tlessthan\t<integer number>  \n\tgreaterthan\t<integer number>\n\tequalto\t\t<integer number>\n\tcontaining\t<taxanumbers>\n");
        }
    if(num == 11)
        {
        printf("\ninclude [options]  \nNOT IMPLEMENTED YET\n\n");
        printf("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf("\t===========================================================\n");

        printf("\n\tall\t\tyes | no\t\t\tno\n\tlessthan\t<integer number>  \n\tgreaterthan\t<integer number>\n\tequalto\t\t<integer number>\n\tcontaining\t<taxanumbers>\n");
        }
     if(num == 12)
        {
        printf("\n!  \n\n\n");
        printf("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf("\t===========================================================\n");

        printf("\n\tTemporarily returns the user to the operating system by running a C shell\n\tType 'exit' to return to clann");
        }
    if(num == 13)
        printf("\nquit\n\tThis command quits Clann\n");
	
	if(num == 14)
		{
        printf("\ngeneratetrees [options]  \n\n\n");
        printf("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf("\t===========================================================\n");

        printf("\n\tmethod\t\tequiprobable | markovian\t*equiprobable\n\tntrees\t\tall | <integer number>\t\t*100\n\tnbins\t\t<integer number>\t\t*20\n\toutfile\t\t<output file name>\t\t*histogram.txt\n\tsourcedata\treal | randomised | ideal\t*real\n\tsavescores\tyes | no\t\t\t*no\n\tsupertree\tmemory | <supertree file name>\t*memory\n\tsavesourcetrees\tyes | no\t\t\t*no\n\n\tThe option 'supertree' is only used when idealised data are being created\n\n");
		
		
		}
     
	 if(num == 15)
		{
		printf("\nconsensus [options]\n\n\n");
		printf("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf("\t===========================================================\n");
		
		printf("\n\tmethod\tstrict | majrule | minor | <value>\t*minor\n\tguidetree\t<guide tree file name>\t\t*<none>\n\tfilename\t<output file name>\t\t*consensus.ph\n");
		}
	 
	 if(num == 16)
		{
		printf("\nshowtrees\trange | size | namecontains | containstaxa | score \n\n");
		printf("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf("\t===========================================================\n");
		
		printf("\n\trange\t\t<integer value> - <integer value> \t*all\n\tsize\t\tequalto <integer value>\n\t\t\tlessthan <integer value>\n\t\t\tgreaterthan <integer value>\t\t*none\n\tnamecontians\t<character string>\t\t\t*none\n\tcontainstaxa\t<character string>\t\t\t*none\n\tscore\t\t<min score> - <max score>\t\t*none\n\tsavetrees\tyes | no\t\t\t\t*no\n\tfilename\t<output file name>\t\t\t*showtrees.txt\n\tdisplay\t\tyes | no\t\t\t\t*yes\n");
		}

	 if(num == 17)
		{
		printf("\tdeletetrees\trange | size | namecontains | containstaxa | score \n\n");
		printf("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf("\t===========================================================\n");
		
		printf("\n\tsinglecopy\tN/A\n\tmulticopy\tN/A\n\trange\t\t<integer value> - <integer value> \t*all\n\tsize\t\tequalto <integer value>\n\t\t\tlessthan <integer value>\n\t\t\tgreaterthan <integer value>\t\t*none\n\tnamecontains\t<character string>\t\t\t*none\n\tcontainstaxa\t<character string>\t\t\t*none\n\tscore\t\t<min score> - <max score>\t\t*none\n");
		}

	 if(num == 18)
		{
		printf("\nincludetrees\trange | size | namecontains | containstaxa | score \n\n");
		printf("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf("\t===========================================================\n");
		
		printf("\n\trange\t\t<integer value> - <integer value> \t*all\n\tsize\t\tequalto <integer value>\n\t\t\tlessthan <integer value>\n\t\t\tgreaterthan <integer value>\t\t*none\n\tnamecontians\t<character string>\t\t\t*none\n\tcontainstaxa\t<character string>\t\t\t*none\n\tscore\t\t<min score> - <max score>\t\t*none\n");
		}
	
	 if(num == 19)
		{
		printf("\nrfdists [options]\t\n\n");
		printf("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf("\t===========================================================\n");
		printf("\n\tfilename\t<output file name>\t\t*robinson_foulds.txt\n\toutput\t\tmatrix | vector\t\t*matrix\n\tmissing\t\tnone | 4point | ultrametric\t*none\n");
		}
	if(num == 20)
		{
		printf("\ndeletetaxa\t <taxa name> <taxa name> etc...\n\n");
		printf("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf("\t===========================================================\n");
		}
	if(num == 21)
		{
		printf("\nnj [options]\t\n\n");
		printf("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf("\t===========================================================\n");
		printf("\n\tmissing\t\t4point | ultrametric\t\t*4point\n\tsavetrees\t<file name>\t\t\t*NJtree.ph\n");
		}
	 if(num == 22)
		{
		printf("\nsprdists [options]\t\n\n");
		printf("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf("\t===========================================================\n");
		printf("\n\tsupertree\tcreate | memory | <file name>\t*create\n\tdorandomisation\tyes | no\t\t\t*yes\n\toutfile\t\t<file name>\t\t\t*SPRdistances.txt\n");
		}
	 if(num == 23)
		{
		printf("\nmapunknowns \t\n\n");
		printf("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf("\t===========================================================\n");
		printf("\n\t\n");
		}
	 if(num == 24)
		{
		printf("\nreconstruct \t\n\n");
		printf("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf("\t===========================================================\n");
		printf("\n\tduplications\t<value>\t\t\t\t*1.0\n\tlosses\t\t<value>\t\t\t\t*1.0");
		printf("\n\tshowrecon\tyes | no\t\t\t*no\n\tbasescore\t<value>\t\t\t\t*1.0\n\tprintfiles\tyes | no\t\t\t*yes");
		}
     if(num == 25)
        {
        printf("\nprunemonophylies \t\n\n");
        printf("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf("\t===========================================================\n");
        printf("\n\tfilename\t<output file name>\t\t*prunedtrees.txt");
        printf("\n\tselection\trandom | length\t\t\t*random\n\n\tIf \"length\" is chosen, then name MUST have a number directly following the name of the species\n\t representing the sequence length. i.e.: \"Speces.length.XXXXXX\n\n" );
        }
    if(num == 26)
        {
        printf("\ntips \t\n\n");
        printf("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf("\t===========================================================\n");
        printf("\n\tnumber\t<value between 1 and 10>\t\trandom");
        }
    if(num == 27)
        {
        printf("\trandomisetrees \t\n\n");
        printf("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf("\t===========================================================\n");
        printf("\n\tThis command randomises all source trees in memory\n");
        }
    if(num == 28)
        {
        printf("\trandomprune \t\n\n");
        printf("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf("\t===========================================================\n");
        printf("\n\tThis command randomly prunes the source trees in memory\n");
        }
	
	
		
     if(num != 0)
        {
        printf("\n\n\t\t\t\t\t\t\t*Option is nonpersistent\n");
        printf("\t===========================================================\n\n");    
        }
    }
    
    
    
    
    
int parse_command(char *command)
    {
    int i=0, j=0, k=0;
    
    for(i=0; i<1000*parsed_command_assignments; i++)  /* reset the array for the new commands */ 
        {
        strcpy(parsed_command[i], "");
        }
    i=0;
    
    while(command[i] != '\0' && command[i] != ';')
        {
        while(command[i] == ' ' || command[i] == '=' || command[i] == '-')
            i++;
        k = 0;
        while(command[i] != ' ' && command[i] != '=' && command[i] != '-' && command[i] != '\0' && command[i] != ';')
            {
           /* parsed_command[j][k] = tolower(command[i]); */
			if(command[i] == '\'' || command[i] == '"')
				{
				i++;
				while(command[i] != '\'' && command[i] != ';' && command[i] != '\0' &&  command[i] != '"')
					{
					parsed_command[j][k] = command[i];
					i++; k++;
					}
				i++;
				}
			else
				{
				parsed_command[j][k] = command[i];
				i++;k++;
				}
            }
        parsed_command[j][k] = '\0';
        j++;i++;
        while(command[i] == ' ' || command[i] == '=' || command[i] == '-')
            i++;
		if(j+1 == 1000*parsed_command_assignments)
			{
			parsed_command_assignments++;
			parsed_command = realloc(parsed_command, parsed_command_assignments*1000*sizeof(char *));
			k=0;
			for(k=((parsed_command_assignments-1)*1000); k<(parsed_command_assignments*1000); k++)
				{
				parsed_command[k] = malloc(1000*sizeof(char));
				parsed_command[k][0] = '\0';
				}
			}
        }

    return(j);
    }


        
void execute_command(char *commandline, int do_all)
    {
    int i = 0, j=0, k=0, printfundscores = FALSE, error = FALSE;
    char c = '\0', temp[NAME_LENGTH], filename[100], *newbietree, string_num[1000];
	float num = 0;

    for(i=0; i<num_commands; i++)
        {
        if(strcmp(parsed_command[i], "print") == 0)
            printfundscores = TRUE;
        
		if(strcmp(parsed_command[i], "summary") == 0)
			{
			if(strcmp(parsed_command[i+1], "short") == 0)
				do_all = FALSE;
			}
		
		
        if(strcmp(parsed_command[i], "maxnamelen") == 0)
            {
            max_name_length = toint(parsed_command[i+1]);
            if(max_name_length == 0)
                {
				if(strcmp(parsed_command[i+1], "delimited") == 0)
					{
					max_name_length = NAME_LENGTH;
					delimiter = TRUE;
					}
				else
					{
					error = TRUE;
					printf("Error: value %s not valid for maxnamelen\n\n", parsed_command[i+1]);
					}
				}
            } 
        }
    /********** open the input file ************/
    filename[0] = '\0';
 
    if(commandline == '\0')
        {
        strcpy(filename, parsed_command[1]);
        }
    else
        strcpy(filename, commandline);
    
    if((infile = fopen(filename, "r")) == '\0')		/* check to see if the file is there */
        {								/* Open the source tree file */
        printf("Cannot open file %s\n", filename);
        error = TRUE;
        }
	strcpy(inputfilename, filename);
    if(!error)
        {
		largest_tree = 0;
		if(yaptp_results != '\0')
			{
			free(yaptp_results);
			yaptp_results = NULL;
			}
		num_excluded_trees = 0;
		num_excluded_taxa = 0;
		trees_in_memory = 0;
            /************************ Assign the dynamic arrays *************************/
        newbietree = malloc(TREE_LENGTH*sizeof(char));
		if(newbietree == NULL) memory_error(110);
		newbietree[0] = '\0';
		if(weighted_scores != NULL)
			{
			for(i=0; i<number_of_taxa; i++)
				free(weighted_scores[i]);
			free(weighted_scores);
			weighted_scores = NULL;
			}
		
        if(fundamentals != NULL)  /* if we have assigned the fundamental array earlier */ 
            {
            for(i=0; i<fundamental_assignments*FUNDAMENTAL_NUM; i++)
                {
                free(fundamentals[i]);
                }
            free(fundamentals);
			
            fundamentals = NULL;
            }

        
		
		
        if(presence_of_taxa != NULL)  /* if we have already assigned presence_of_taxa */
            {
            for(i=0; i<fundamental_assignments*FUNDAMENTAL_NUM; i++)
                {
                free(presence_of_taxa[i]);
                }
            free(presence_of_taxa);
            presence_of_taxa = NULL;
            }
    
    
        
        fundamental_assignments = 1;
        tree_length_assignments = 1;   
		

		if(fulltaxanames != NULL)
			{
			for(i=0; i<fundamental_assignments*FUNDAMENTAL_NUM; i++)
				{
				if(fulltaxanames[i] != NULL)
					{
					for(j=0; j<numtaxaintrees[i]; j++)
						{
						if(fulltaxanames[i][j] != NULL)
							free(fulltaxanames[i][j]);
						}
					free(fulltaxanames[i]);
					}
				}
			free(fulltaxanames);
			fulltaxanames = NULL;
			}
		
		if(numtaxaintrees != NULL)
			free(numtaxaintrees);		
		fulltaxanames = malloc((fundamental_assignments*FUNDAMENTAL_NUM)*sizeof(char **));
		if(!fulltaxanames) memory_error(106);
		numtaxaintrees = malloc((fundamental_assignments*FUNDAMENTAL_NUM)*sizeof(int));
		if(!numtaxaintrees) memory_error(106);
		tree_names = malloc((fundamental_assignments*FUNDAMENTAL_NUM)*sizeof(char *));
		if(!tree_names) memory_error(106);
		tree_weights = malloc((fundamental_assignments*FUNDAMENTAL_NUM)*sizeof(float));
		if(!tree_weights) memory_error(107);
        fundamentals = malloc((fundamental_assignments*FUNDAMENTAL_NUM)*sizeof(char *));
        if(!fundamentals) memory_error(1);
        for(i=0; i<fundamental_assignments*FUNDAMENTAL_NUM; i++)
            {
			
			fulltaxanames[i] = NULL;

			tree_names[i] = malloc(1000*sizeof(char));
			if(!tree_names) memory_error(108);
			tree_names[i][0] = '\0';
			tree_weights[i] = 1;
			fundamentals[i] = NULL;
         /*   fundamentals[i] = malloc((tree_length_assignments*TREE_LENGTH)*sizeof(char));
            if(!fundamentals[i]) memory_error(2);
                
            fundamentals[i][0] = '\0';
           */ }
    
    
        if(fund_scores != NULL)  /* if we have already assigned fund_scores */
            {
            for(i=0; i<Total_fund_trees; i++)
                {
                for(j=0; j<number_of_taxa; j++)
                    {
                    free(fund_scores[i][j]);
                    }
                free(fund_scores[i]);
                }
            free(fund_scores);
            fund_scores = NULL;
            }

    
        if(taxa_names != NULL) /* if we have assigned taxa names earlier */
            {
            for(i=0; i<TAXA_NUM*name_assignments; i++)
                {
                free(taxa_names[i]);
                }
            free(taxa_names);
            free(taxa_incidence);
            taxa_names = NULL;
            }
    
        if(Cooccurrance != NULL) /* if we have assigned taxa names earlier */
            {
            for(i=0; i<TAXA_NUM*name_assignments; i++)
                {
                free(Cooccurrance[i]);
                }
            free(Cooccurrance);
            free(same_tree);
            Cooccurrance = NULL;
            same_tree = NULL;
            }
            
        name_assignments = 1;
        
        presence_of_taxa = malloc((fundamental_assignments*FUNDAMENTAL_NUM)*sizeof(int *));
        if(!presence_of_taxa) memory_error(3);
           
        for(i=0; i<fundamental_assignments*FUNDAMENTAL_NUM; i++)
            {
            presence_of_taxa[i] = malloc((TAXA_NUM*name_assignments)*sizeof(int));
            if(!presence_of_taxa[i]) memory_error(4);
                
            for(j=0; j<TAXA_NUM*name_assignments; j++)
                {
                presence_of_taxa[i][j] = FALSE;
                }
            }
            


        taxa_names = malloc((TAXA_NUM*name_assignments)*sizeof(char *));
            if(!taxa_names) memory_error(5);
                
            for(j=0; j<TAXA_NUM*name_assignments; j++)
                {
                taxa_names[j] = malloc(NAME_LENGTH*sizeof(char));
                if(!taxa_names[j]) memory_error(6);
                    
                taxa_names[j][0] = '\0';
                }
                
        Cooccurrance = malloc((TAXA_NUM*name_assignments)*sizeof(int *));
        if(!Cooccurrance) memory_error(7);
                
        for(j=0; j<TAXA_NUM*name_assignments; j++)
            {
            Cooccurrance[j] = malloc((TAXA_NUM*name_assignments)*sizeof(int));
            if(!Cooccurrance[j]) memory_error(8);
                
            for(k=0; k<TAXA_NUM*name_assignments; k++)
                Cooccurrance[j][k] = 0;
            }
    
        same_tree = malloc((TAXA_NUM*name_assignments)*sizeof(int));
        if(!same_tree) memory_error(9);
            
        for(j=0; j<TAXA_NUM*name_assignments; j++)
            same_tree[j] = 0;
    
    
        taxa_incidence = malloc((TAXA_NUM*name_assignments)*sizeof(int));
        if(!taxa_incidence) memory_error(10);
            
        for(j=0; j<TAXA_NUM*name_assignments; j++)
            taxa_incidence[j] = 0;
    
        
        if(number_of_comparisons != NULL)
            {
            free(number_of_comparisons);
            number_of_comparisons = NULL;
            }
        
        number_of_comparisons = malloc((fundamental_assignments*FUNDAMENTAL_NUM)*sizeof(int));
        if(!number_of_comparisons) memory_error(11);
            
        for(j=0; j<fundamental_assignments*FUNDAMENTAL_NUM; j++)
            number_of_comparisons[j] =0;

    
		hsprint = TRUE;
        number_of_taxa = 0;
        num_excluded_trees = 0;
        Total_fund_trees = 0;
         /************** next read in all the source trees into memory **************/
        /* we need to determine whether the input file is in NEXUS or in newhapshire (nested parenthsis) format */
        /* If it is in NEXUS format then the first (non redundant) character will be a "#" */ 

        while(((c = getc(infile)) == ' ' || c == '\n' || c == '\r' || c == '[') && !feof(infile)) 
            {
            if(c == '[')
                while(!feof(infile) && (c = getc(infile)) != ']');
            }		/* Skip over nontree characters */
        
        if(c == '#')
            {
            printf("\nReading Nexus format source tree file \n");


			if(nexusparser(infile) == TRUE)
				{
				printf("error reading nexus file\n");
				}


            
            
            
            }  /* End reading NEXUS format tree */
        else
            {
			printf("\nReading Newhampshire (Phylip) format source trees\n");
            
			if(c != '(')
                {
                printf("Error: Treefile not in correct format\n");
                }
            else
                {
                got_weights = FALSE;
                while(!feof(infile)) /* While we haven't reached the end of the file */
                    {
                    for(i=0; i<10; i++) string_num[i] = '\0';
					i=0;
                    while(!feof(infile) && c != ';' )
                        {
						if(c != ' ' && c != '\n' && c != '\r')
							{
							newbietree[i] = c;
							i++;
							}
						c = getc(infile);
						if(c == '[')
							{
							j=0;
							while(c != ']' && !feof(infile))
								{
								c = getc(infile);
								string_num[j] = c; j++;
								}
							string_num[j] = '\0';
							tree_weights[Total_fund_trees] = atof(string_num);
							c = getc(infile);
							}
						}
					
					newbietree[i] = ';';
					newbietree[i+1] = '\0';

					input_fund_tree(newbietree, Total_fund_trees);

					c = getc(infile);
					while((c == ' ' || c == '\t' || c == '\n' || c == '\r') && !feof(infile)) c = getc(infile);
					if(c == '[')
						{
						j=0;
						c = getc(infile);
						while(c != ']' && !feof(infile))
							{

							tree_names[Total_fund_trees-1][j] = c;
							if(j<99) j++;
							c = getc(infile);
							}

						tree_names[Total_fund_trees-1][j] = '\0';
						while((c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == ']') && !feof(infile)) c = getc(infile);

						}
                    }  /* End of reading tree */

                }
            } /* end reading Phylip format trees */
	fclose(infile); 
	infile = NULL;
	remainingtrees = Total_fund_trees;
		if(sourcetreetag != NULL) free(sourcetreetag);
		sourcetreetag = malloc(Total_fund_trees*sizeof(int));
		for(i=0; i<Total_fund_trees; i++) sourcetreetag[i] = TRUE;		
		if(presenceof_SPRtaxa) free(presenceof_SPRtaxa);
		presenceof_SPRtaxa = malloc(number_of_taxa*sizeof(int));
		for(i=0; i<number_of_taxa; i++) presenceof_SPRtaxa[i] = -1;
	weighted_scores = malloc(number_of_taxa * sizeof(float*));
	if(!weighted_scores) memory_error(94);
	for(k=0; k<number_of_taxa; k++)
		{
		weighted_scores[k] = malloc(number_of_taxa*sizeof(float));
		if(!weighted_scores[k]) memory_error(95);
		for(j=0; j<number_of_taxa; j++)
			weighted_scores[k][j] = 0;
		}	
	

        input_file_summary(do_all);
	

		
        if(sourcetree_scores) free(sourcetree_scores);
		sourcetree_scores = malloc(Total_fund_trees*sizeof(float));
		for(i=0; i<Total_fund_trees; i++) sourcetree_scores[i] = -1;

        fund_scores = malloc(Total_fund_trees*sizeof(int**));
        if(fund_scores == NULL) memory_error(35);
            
        for(i=0; i<Total_fund_trees; i++)
            {
            fund_scores[i] = malloc((number_of_taxa)*sizeof(int*));
            if(fund_scores[i] == NULL) memory_error(36);
            else
                {
                for(j=0; j<(number_of_taxa); j++)
                    {
                    fund_scores[i][j] = malloc((number_of_taxa)*sizeof(int));
                    if(fund_scores[i][j] == NULL) memory_error(37);
                    else
                        {
                        for(k=0; k<(number_of_taxa); k++)
                            fund_scores[i][j][k] = 0;
                        }
                    }
                }
            }

		calculated_fund_scores = FALSE;

		
		free(newbietree);

        }
    }
    


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
											while((c = getc(nexusfile)) != ']' &&!feof(nexusfile))
												{
												string[i] = c;
												if(string[i] != ' ') i++;
												}
											string[i] = '\0';
											tree_weights[Total_fund_trees] = tofloat(string);
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
													printf("Error: taxa %s is not in the translation table\n", single);
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
    printf("\n\tSource tree summmary:\n\n");
    printf("\t----------------------------------------------------\n");
    printf("\t\tNumber of input trees: %d\n", Total_fund_trees-num_excluded_trees);
    printf("\t\tNumber of unique taxa: %d\n", number_of_taxa-num_excluded_taxa);
    sup =1;
    for(i=4; i<=number_of_taxa-num_excluded_taxa; i++) sup*=((2*i)-5);
    printf("\t\tTotal unrooted trees in Supertree space?\n\t\t\t%g\n\n", sup);
    

    printf("\tOccurrence summary:\n\n");
    
    printf("\t\tnumber\tTaxa name           \t\tOccurrence\n\n");
    for(i=0; i<number_of_taxa; i++)
        {
        printf("\t\t%-4d\t%-30s\t%d\n", i, taxa_names[i], taxa_incidence[i]);
        }
   
	if(do_all)
		{
		printf("\n\n\tCo-occurrence summary:\n\n");
		printf("\t\tTaxa Number\n      ");


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
				printf("\n\t\t      ");
			else
				printf("\n\tCo-occurence summary continued:\n\t\t      ");
			for(i=previous_limit; i<limit; i++)
				printf("%-5d ", i);
			printf("\n");
			for(i=previous_limit; i<number_of_taxa; i++)
				{
				printf("\t\t%-5d ", i);
				for(j=previous_limit; j<=i; j++)
					{
					if(j<limit)
						{
						if(j== i)
							printf("-     ");
						else
							printf("%-5d ", Cooccurrance[i][j]);
						}
					}
				printf("\n");
				}
			previous_limit = limit;
			limit = limit+12;
			if(limit > number_of_taxa) limit = number_of_taxa;
			}
		}

	printf("\n\n\tSource tree size summary:\n\n  num leaves\n");
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

	printf("\tNumber of single copy trees:\t%d\n", singlecopy);
	printf("\tnumber of multicopy trees:\t%d\n", Total_fund_trees-singlecopy);

    printf("\t----------------------------------------------------\n");
	
	free(size_of_trees);
    }



int assign_taxa_name(char *name,int fund)
	{
	int i =0, j=0, k=0, answer = -1, taxa_on_tree = 0;
	char *delim = NULL;
	
	delim = malloc(10*sizeof(char));
	delim[0] = '.';
	delim[1] = '\0';

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
		while(name[i] != '.' && name[i] != '\0') i++;
		name[i] = '\0';
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

/* A simple version of the standard gets() library function which is unsafe on some systems......... code taken from C: the complete reference
	by Herbert Schildt
*/
char *xgets(char *s)
	{
	char ch, *p = NULL;
	int t = 0;
        
        for(t=0; t<80; t++)
            s[t] = '\0';
	p = s;  /* gets returns a pointer to s */

	for(t=0; t<80; ++t) {
		ch = getchar();

		switch(ch) {
			case '\n':
				s[t] = '\0'; /* terminate the string */
				return(p);
				break;
			case '\r':
				s[t] = '\0'; /* terminate the string */
				return(p);
				break;
			case '\b':
				if(t>0) t--;
				break;
			default:
				s[t] = ch;
				break;
			}
	}
	s[79] = '\0';
	return p;
}



void clean_exit(int error)
    {
    int i=0, j=0;
	if(weighted_scores != NULL)
		{
		for(i=0; i<number_of_taxa; i++)
			free(weighted_scores[i]);
		free(weighted_scores);
		}
	weighted_scores = NULL;
	if(yaptp_results != NULL) free(yaptp_results);
    if(stored_commands[i] != NULL)
		{
		for(i=0; i<100; i++) free(stored_commands[i]);
		free(stored_commands);
		}
	if(parsed_command != NULL)
        {
        for(i=0; i<1000; i++)
			{
			free(parsed_command[i]);
			}
        free(parsed_command);
        }
    if(retained_supers != NULL)
        {
        for(i=0; i<number_retained_supers; i++)
            free(retained_supers[i]);
        free(retained_supers);
        }
	if(fulltaxanames != NULL)
		{
		for(i=0; i<fundamental_assignments*FUNDAMENTAL_NUM; i++)
			{
			if(fulltaxanames[i] != NULL)
				{
				for(j=0; j<numtaxaintrees[i]; j++)
					{
					if(fulltaxanames[i][j] != NULL)
						free(fulltaxanames[i][j]);
					}
				free(fulltaxanames[i]);
				}
			}
		free(fulltaxanames);
		fulltaxanames = NULL;
		}
	if(numtaxaintrees != NULL)
		free(numtaxaintrees);		

    if(fundamentals != NULL)  /* if we have assigned the fundamental array earlier */ 
        {
        for(i=0; i<fundamental_assignments*FUNDAMENTAL_NUM; i++)
            {
            free(fundamentals[i]);
            }
        free(fundamentals);
        fundamentals = NULL;
        }
    if(presence_of_taxa != NULL)  /* if we have already assigned presence_of_taxa */
        {
        for(i=0; i<fundamental_assignments*FUNDAMENTAL_NUM; i++)
            {
            free(presence_of_taxa[i]);
            }
        free(presence_of_taxa);
        presence_of_taxa = NULL;
        }
    if(fund_scores != NULL)  /* if we have already assigned fund_scores */
        {
        for(i=0; i<Total_fund_trees; i++)
            {
            for(j=0; j<number_of_taxa; j++)
                {
                free(fund_scores[i][j]);
                }
            free(fund_scores[i]);
            }
        free(fund_scores);
        fund_scores = NULL;
        }
    if(taxa_names != NULL) /* if we have assigned taxa names earlier */
        {
        for(i=0; i<TAXA_NUM*name_assignments; i++)
            {
            free(taxa_names[i]);
            }
        free(taxa_names);
        free(taxa_incidence);
        taxa_names = NULL;
        }
    if(Cooccurrance != NULL) /* if we have assigned taxa names earlier */
        {
        for(i=0; i<TAXA_NUM*name_assignments; i++)
            {
            free(Cooccurrance[i]);
            }
        free(Cooccurrance);
        free(same_tree);
        Cooccurrance = NULL;
        same_tree = NULL;
        }
    if(number_of_comparisons != NULL)
        {
        free(number_of_comparisons);
        number_of_comparisons = NULL;
        }
    if(tree_top != NULL)
        {
        dismantle_tree(tree_top);
        tree_top = NULL;
        }
    temp_top = NULL;
    if(total_coding != NULL)
        {
        for(i=0; i<total_nodes; i++)
            free(total_coding[i]);
        free(total_coding);
        total_coding = NULL;
        }
    if(partition_number != NULL) free(partition_number);
    if(from_tree != NULL) free(from_tree);

    if(sourcetree_scores != NULL) free(sourcetree_scores);
	if(presenceof_SPRtaxa) free(presenceof_SPRtaxa);
	
	if(sourcetreetag != NULL) free(sourcetreetag);
    if(error)
		exit(1);
    else
		exit(0);
    }
    

/* calculate the path metric for each of the fundamental trees */
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










void unroottree(char * tree)
    {
    int i=0, j=0, k=0, l=0, m=0, basecount = 0, parentheses=0;
    int foundopen = FALSE, foundclose = FALSE;
	float del_nodelen = 0;
	char length[100], restof[TREE_LENGTH];
	
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
    }

    
int texttoint(char c)
    {
    switch(c)
        {
        case '1' :
            return(1);
            break;
        case '2' :
            return(2);
            break;
        case '3':
            return(3);
            break;
        case '4':
            return(4);
            break;
        case '5':
            return(5);
            break;
        case '6':
            return(6);
            break;
        case '7':
            return(7);
            break;
        case '8':
            return(8);
            break;
        case '9':
            return(9);
            break;
        default:
            return(0);
            break;
        }
    }

    
    /* This function returns the number that represents the given tree */
int treeToInt(char *array)
    {
    int i;
    /* The first part is to calculate the path of the tree, given the tree in nested parenthesis format */
    /*  this part of the code assumes that the tree being tested has been created using the intToTree function. The reason for this
        is that that function always puts the position of a new taxon after the parenthesis or taxon it is to be paired with.
        If a tree to be tested is made in any other way, it will lead to these assumptions to be violated and the function may crash or
        the wrong number may be calculated. cc 16/06/2003    */
        
    for(i=number_of_taxa; i>-3; i--)
        {
        
        
        
        }
    
    
    
    
    
    
    
    
    return(1);
    }
    
    
    
    /* this function builds a tree given its number and the number of taxa */
void intTotree(int tree_num, char *array, int num_taxa)
    {
    int  i = 0, j=0, k=0, l=0, exit = 0, bracket_count = 0;
    char tmparray[TREE_LENGTH], *string= NULL;
	double supers = 1, max = 1, min = 0, oldmin = 0, *path = NULL;
	
	for(i=4; i<=num_taxa; i++) supers*=((2*i)-5);

    path = malloc((num_taxa-3)*sizeof(double));  /* this will store the path of the tree */
    string = malloc(30*sizeof(char));
    
    /************* calculate the path of the tree from a triplet to the full size of the tree ***************/
    
    min = 1;
    max = supers; /* max is now the number of possible trees in this tree space */
   /* printf("path of tree number %d with %d taxa is:\n", tree_num, number_of_taxa);
   */
     for(i=4; i<=num_taxa; i++)  /* from a triplet to the full complement of taxa, we figure out the placement of each new taxa at each stage */
        {
        
        path[i-4] = div((tree_num - min), ((max - (min - 1))/((2 * i)- 5))).quot + 1;

   /*     printf("min = %d, max = %d path = %d\n",min, max, path[i-4]);
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
        /*    printf("%s\n", array);  */
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
    
    
    
    free(path);
    free(string);
    }
    
int toint(char *number)
    {
    int charactercount = 0, charactercount1 =0, result = 0, j=0;
    
    while(number[charactercount] != '\0') charactercount ++;
    for(j=charactercount; j>0; j--)
        {
        result += ( pow(10, (charactercount-j)) * (texttoint(number[j-1])));
        }
    return(result);
    }




float tofloat(char *number)
    {
    int charactercount = 0,charactercount1 = 0, j=0;
    float real = 0, fraction = 0, total = 0;



    while(number[charactercount] !=  '.' && number[charactercount] != '\0') charactercount ++;
    for(j=charactercount; j>0; j--)
        {
        real += ( pow(10, (charactercount-j)) * (texttoint(number[j-1])));
        }
    
    
    if(number[charactercount] == '.')
        {
        charactercount1 = charactercount;
        while(number[charactercount1] != '\0') charactercount1++;
        for(j=charactercount+1; j<charactercount1; j++)
            {
            fraction += ( pow(10,(charactercount-j))*(texttoint(number[j])));
            }
        }
	if(real == 0 && fraction == 0) total = 0;
	else total = real+fraction;
    return(total);
    }
     



void alltrees_search(int user)
    {
    int i = 0, j=0, all = TRUE, start = 0, end = 0, error = FALSE, keep = 0;
    char *tree = NULL, *best_tree = NULL, outfilename[100];
    float score = 0, best_score = 0, worst = 0;
    FILE *treesfile = NULL, *userfile = NULL;
    
    
    outfilename[0] = '\0';
    strcpy(outfilename, "alltrees.ph");
    
    if(user)  /*if this was called by the user and not by bootstrap or yaptp */
        {
        for(i=0; i<num_commands; i++)
            {
            if(strcmp(parsed_command[i], "all") == 0)
                all = TRUE;
            else
                {
                if(strcmp(parsed_command[i], "range") == 0)
                    {
                    start = toint(parsed_command[i+1]);
                    end = toint(parsed_command[i+2]);
                    if(end == start) error = TRUE;
                    } 
                }
            if(strcmp(parsed_command[i], "keep") == 0)
                {
                worst = tofloat(parsed_command[i+1]);
                if(worst == 0)
                    {
                    printf("Error: '%s' is an invalid value for keep\n", parsed_command[i+1]);
                    error = TRUE;
                    }
                }
            if(strcmp(parsed_command[i], "nbest") == 0)
                {
                keep = toint(parsed_command[i+1]);
                if(keep == 0)
                    {
                    printf("Error: '%s' is an invalid values for nbest\n", parsed_command[i+1]);
                    error = TRUE;
                    }
                }
            if(strcmp(parsed_command[i], "savetrees") == 0 && criterion != 1)
             {
                if((userfile = fopen(parsed_command[i+1], "w")) == NULL)
                 {
                    printf("Error opening file named %s\n", parsed_command[i+1]);
                    error = TRUE;
                 }
                else
                 {
                    printf("opened output file %s\n", parsed_command[i+1]);
                    strcpy(outfilename, parsed_command[i+1]);
                 }
             }
            
            if(strcmp(parsed_command[i], "create") == 0)
                {
                if((treesfile = fopen("alltrees.ph", "w")) == NULL) 
                    {
                    printf("Error opening file named 'alltrees.ph'\n");
                    error = TRUE;
                    }
                }
            if(criterion == 3)
                {
                if(strcmp(parsed_command[i], "weight") == 0)
                    {
                    if(strcmp(parsed_command[i+1], "equal") == 0)
                        quartet_normalising = 1;
                    else
                        {
                        if(strcmp(parsed_command[i+1], "taxa") == 0)
                            quartet_normalising = 2;
                        else
                            {
                            if(strcmp(parsed_command[i+1], "quartets") == 0)
                                quartet_normalising = 3;
                            else
                                {
                                printf("Error: weight option '%s' is unknown\n", parsed_command[i+1]);
                                error = TRUE;
                                }
                            }
                        }
                    }
                }
            if(criterion == 0)
                {
                if(strcmp(parsed_command[i], "weight") == 0)
                    {
                    if(strcmp(parsed_command[i+1], "equal") == 0)
                        dweight = 0;
                    else
                        {
                        if(strcmp(parsed_command[i+1], "comparisons") == 0)
                            dweight = 1;
                        else
                            {
                            printf("Error: weight option '%s' is unknown\n", parsed_command[i+1]);
                            error = TRUE;
                            }
                        }
                    }
                }            
            if(criterion == 2)
                {
                if(strcmp(parsed_command[i], "weight") == 0)
                    {
                    if(strcmp(parsed_command[i+1], "equal") == 0)
                        splits_weight = 1;
                    else
                        {
                        if(strcmp(parsed_command[i+1], "splits") == 0)
                            splits_weight = 2;
                        else
                            {
                            printf("Error: weight option '%s' is unknown\n", parsed_command[i+1]);
                            error = TRUE;
                            }
                        }
                    }
                }            
            }
        }
    if(!error)
        {
        if(userfile == NULL)
            {
            if((userfile = fopen("top_alltrees.txt", "w")) == NULL) 
                {
                printf("Error opening file named 'alltrees.ph'\n");
                error = TRUE;
                }
            }
        }
        
    if(!error)
        {
        
        if(start == 0 && end == 0)
            {
            start = 1;
            end = sup;
            }
            
        printf("\n\nAlltrees (exhaustive search) settings:\n\trange: tree numbers %d to %d inclusive\n\tOutput file: %s\n\tCreate all trees? ", start, end, outfilename );
        if(treesfile != NULL) printf("Yes\n");
        else printf("No\n");
        printf("\tWeighting Scheme = ");
        if(criterion==0)
            {
            if(dweight == 0) printf("equal\n");
            if(dweight == 1) printf("comparisons\n");
            }
        if(criterion == 2)
            {
            if(splits_weight == 1) printf("equal\n");
            if(splits_weight == 2) printf("splits\n");
            }
        if(criterion == 3)
            {
            if(quartet_normalising == 1) printf("equal\n");
            if(quartet_normalising == 2) printf("taxa\n");
            if(quartet_normalising == 3) printf("quartets\n");
            }
        printf("\n\n");
            
		if(!calculated_fund_scores && criterion == 0)
			{
			cal_fund_scores(FALSE);    /* calculate the path metrics for all the fundamental trees */
			calculated_fund_scores = TRUE;
			}
                        
        psfile = fopen("supertree.ps", "w");

        for(i=0; i<number_of_taxa; i++) presenceof_SPRtaxa[i] = -1;
        /***** define the dynamic arrays  **********/
        
        tree = malloc(TREE_LENGTH*sizeof(char));
        if(!tree) memory_error(25);
            
        tree[0] = '\0';
    
    
        best_tree = malloc(TREE_LENGTH*sizeof(char));
        if(!best_tree)  memory_error(26);
            
        best_tree[0] = '\0';
    
        /* this hold the distances calculated with the pathmetric on the supertree */
        
        if(super_scores == NULL)
            {
            super_scores = malloc(number_of_taxa*sizeof(int *));
            if(!super_scores)  memory_error(27);
                
            for(i=0; i<number_of_taxa; i++)
                {
                super_scores[i] = malloc(number_of_taxa*sizeof(int));
                if(!super_scores[i])  memory_error(28);
                    
                for(j=0; j<number_of_taxa; j++)
                    {
                    super_scores[i][j] = 0;
                    }
                }
            }
        else
            {
            for( i=0; i<number_of_taxa; i++)
                {
                for(j=i; j<number_of_taxa; j++)
                    {
                    super_scores[i][j] = 0;
                    super_scores[j][i] = 0;
                    }
                }
            }
    
        if(user) printf("Progress indicator:");
        
        /************ End assign dynamic arrays **************/
    
        if(criterion == 2 || criterion == 3)
            {
            coding(0,0,0);
            if(criterion == 2) condense_coding();
            }
    
    
        if(start == 0 && end == 0)
            {
            start = 1;
            end = sup;
            }
            
		 if(signal(SIGINT, controlc5) == SIG_ERR)
			{
			printf("An error occurred while setting a signal handler\n");
			}
			
        for(i=start; i<=end; i++)
            {            
            tree[0] = '\0';
            if(user_break)
				{
				printf("%d trees sampled\n", i);
				i = end+1;
				}
            interval2 = time(NULL);
            if(difftime(interval2, interval1) > 5) /* every 10 seconds print a dot to the screen */
                {
                printf("=");
                fflush(stdout);
                interval1 = time(NULL);
                }
    
    
            intTotree(i, tree, number_of_taxa);  /* create the supertree number i */
            if(criterion == 0)  /* if we are using mssa */
                {
                /****** We now need to build the Supertree in memory *******/
                if(tree_top != NULL)
                    {
                    dismantle_tree(tree_top);
                    tree_top = NULL;
                    }
                temp_top = NULL;
                tree_build(1, tree, tree_top, FALSE, -1);
                tree_top = temp_top;
                temp_top = NULL;
                score = compare_trees(FALSE);
                }
            if(criterion == 2)  /* if we are using MRC */
                {
                unroottree(tree);
                score = (float)MRC(tree);
                }
            if(criterion == 3)  /* if we are using QC */
                {
                unroottree(tree);
                score = (float)quartet_compatibility(tree);
                }
                
            if(treesfile != NULL)  /* if the create option was selected */
                {
                strcpy(best_tree, "");
                print_named_tree(tree_top, best_tree);
                fprintf(treesfile,"%s;\t[%f]\n", best_tree, score);
                }
    
    
            
            if(i==1)
                {
				retained_supers[0] = realloc(retained_supers[0], (strlen(tree)+10)*sizeof(char));
                strcpy(retained_supers[0], tree);
                scores_retained_supers[0] = score;
                best_score = score;
                }
            else
                {
                if(score < best_score)
                    {
					retained_supers[0] = realloc(retained_supers[0], (strlen(tree)+10)*sizeof(char));
                    strcpy(retained_supers[0], tree);
                    scores_retained_supers[0] = score;
                    best_score = score;
                    j=1;
                    while(scores_retained_supers[j] != -1 && j < number_retained_supers)
                        {
                        strcpy(retained_supers[j], "");
                        scores_retained_supers[j] = -1;
                        j++;
                        }
                    }
                else
                    {
                    if(score == best_score)
                        {
                        j = 1;
                        while(scores_retained_supers[j] != -1)
                            {
                            j++;
                            if(j+1 == number_retained_supers) reallocate_retained_supers();
                            }
						retained_supers[j] = realloc(retained_supers[j], (strlen(tree)+10)*sizeof(char));
                        strcpy(retained_supers[j], tree);
                        scores_retained_supers[j] = score;
                        }
                    }
                }
            }
        
        /***** Print out the best tree found ******/
        if(user)
            {
            printf("\n");
            i=0;
            
            
           /**** Print out the best trees found *******/
                
            tree[0] = '\0';
            i=0; j=0;
            while(scores_retained_supers[j] != -1)
                {
                j++;
                }
            
            while(scores_retained_supers[i] != -1)
                {
                if(tree_top != NULL)
                    {
                    dismantle_tree(tree_top);
                    tree_top = NULL;
                    }
                temp_top = NULL;
                tree_build(1, retained_supers[i], tree_top, FALSE, -1);
                tree_top = temp_top;
                temp_top = NULL;
            
                strcpy(best_tree, "");
            
                print_named_tree(tree_top, best_tree);
                
                if(userfile != NULL) fprintf(userfile, "%s;\t[%f]\n", best_tree, scores_retained_supers[i] );

                tree_coordinates(best_tree, FALSE, TRUE, FALSE, -1);
                printf("\nSupertree %d of %d score = %f\n", i+1, j, scores_retained_supers[i] );
                i++;
                }

            while(i>=0)
                {
                strcpy(retained_supers[i], "");
                scores_retained_supers[i] = -1;
                i--;
                }
            
          
            }
        free(tree);
        free(best_tree);
        if(treesfile != NULL)
            fclose(treesfile);
        treesfile = NULL;
        fclose(psfile);
        
        fclose(userfile);
        }
    }
    


    
/* This function does the checking of every fundamental tree to the supertree at hand. It returns a float 
	This function also bootstraps the values for the fundamental trees, if the do_bootstrap value is greater than 0 these values are printed to bootstrap.txt*/
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
		if(sourcetreetag[i]) /* if this sourcetree is to be used in the analysis ---- defined by the exclude command */
			{
			
			found = FALSE; here1 = TRUE; here2 = TRUE;
			
			/*** Next check to see if the pruning used affected this fundamental tree ***/
			if(presenceof_SPRtaxa[i] != -1 && sourcetree_scores[i] != -1 && spr)
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
				if(print_pruned_tree(tree_top, 0, pruned_tree, FALSE) >1)
					{
					tmp[0] = '\0';
					strcpy(tmp, "(");
					strcat(tmp, pruned_tree);
					strcat(tmp, ")");
					strcpy(pruned_tree, tmp);
					}
				strcat(pruned_tree, ";");
			
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






/* Tree_build:
	This function reads in a file and from it builds the tree in memory using the taxon_type definition */
	
int tree_build (int c, char *treestring, struct taxon *parent, int fromfile, int fundnum)
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
				c = tree_build(c, treestring, position, fromfile, fundnum);
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
					


				if(fromfile)  /* if the tree has come from a file, then the tree will contain the actual taxa names, unlike if it was created by all_trees, where the tree will contain taxa numbers */
					{
					/* assign the name to the taxon */
					j = number_of_taxa;
					position->name = assign_taxa_name(temp, FALSE);
					if(j != number_of_taxa)
						{
						printf("Error, Taxa name %s is not contained in the source trees\n", temp);
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


/* This makes the taxon structure when we need it so I don't have to keep typing the assignments all the time */
struct taxon * make_taxon(void)
	{
	struct taxon *position = NULL;
	
	if(count_now)malloc_check++;
	
	position = malloc(sizeof(taxon_type));
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

void prune_tree_from_array(struct taxon * super_pos, int * array)
	{
	int i=0, found = FALSE;
	struct taxon *start = super_pos;
	
	while(super_pos != NULL)
		{
		super_pos->tag2 = super_pos->tag;
		super_pos = super_pos->next_sibling;
		}
	super_pos = start;
	while(super_pos != NULL)
		{
		found = FALSE;
		if(super_pos->daughter != NULL)
			{
			prune_tree_from_array(super_pos->daughter, array);  /* If this is pointer sibling, move down the tree */
			}
		if(super_pos->name != -1) /* If there is an actual taxa on this sibling */
			{	
			/* Check to see if that taxa is on the array */
			if(array[super_pos->name] == 0)  	
				{  /* if its not there */
				super_pos->tag = FALSE;
				}
			}
		super_pos = super_pos->next_sibling;
		}	
	
	}

void add_internals_from_array(struct taxon * super_pos, int *array)
	{
	int i=0, found = FALSE;
	while(super_pos != NULL)
		{
		found = FALSE;
		if(super_pos->daughter != NULL)
			{
			if(array[super_pos->tag2] > 0)  	
				{  /* if it's there */
				super_pos->tag = super_pos->tag2;
				}	
			add_internals_from_array(super_pos->daughter, array);  /* If this is pointer sibling, move down the tree */
			}
		super_pos = super_pos->next_sibling;
		}	
	
	}


/* Prune tree: This is a recursive function that is called for every node position of the supertree
	it then checks to see if any of the siblings on this node are not contained in the fundamental tree, these siblings are then turned off.
	This only turns off taxa, pointer siblings will have to be turned off using a separate program */

void prune_tree(struct taxon * super_pos, int fund_num)
	{

	
	/* Traverse the supertree, visiting every taxa and checking if that taxa is on the fundamental tree */
	
	while(super_pos != NULL)
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
	
	}
	






/* This function travels through the tree recursively untagging any pointer siblings that are not being used. 
	this effectively shrinks the tree to the size of the fundamental tree that it is being compared to 
	this recursive function is called by shrink_tree to count how many active taxa there are below any given pointer sibling */
int shrink_tree (struct taxon * position)
	{
	int count = 0, tot = 0, i, j, k, l, havelabel = FALSE;
	float one, two;
	char tempstr[100] , tmplabel[100];
	struct taxon *tmppos = NULL;
	
	tmplabel[0] = '\0';
	
	while(position != NULL)
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
	return(count);
	
	}
	
	

/* This function is only used to print the pruned supertree */
int print_pruned_tree(struct taxon * position, int count, char *pruned_tree, int fullname)
    {
    char *name =NULL, temper[100];
    int i=0;
    
    name = malloc(1000000*sizeof(char));
    if(!name) memory_error(33);
    
    while(position != NULL)
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
                print_pruned_tree(position->daughter, 0, pruned_tree, fullname);
                strcat(pruned_tree, ")");
				}
            else
                {
                count = print_pruned_tree(position->daughter, count, pruned_tree, fullname);
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
    free(name);
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

char inttotext(int c)
    {
   
    	
    switch(c)
        {
        case 1 :
            return('1');
            break;
        case 2 :
            return('2');
            break;
        case 3:
            return('3');
            break;
        case 4:
            return('4');
            break;
        case 5:
            return('5');
            break;
        case 6:
            return('6');
            break;
        case 7:
            return('7');
            break;
        case 8:
            return('8');
            break;
        case 9:
            return('9');
            break;
        default:
            return('0');
            break;
        }
    }



				
void reset_tree(struct	taxon * position)
	{
	
	while(position != NULL)
		{
		position->tag = TRUE;
		if(position->daughter != NULL) reset_tree(position->daughter);
		position = position->next_sibling;
		}
	
	}		
				
int count_taxa(struct taxon * position, int count)
	{
	struct taxon * start = position;
	int i =0;
	
	while(position != NULL)
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
	
	while(position != NULL)
		{
		if(position->daughter != NULL)
			num = number_tree(position->daughter, num);
		position = position->next_sibling;
		}
	position = start;
	while(position != NULL)
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
	
/* this identifies the taxa in a subtree passed to it */
void identify_taxa(struct taxon * position, int *name_array)
	{
	while(position != NULL)
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
	
	while(position != NULL)
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
	
	/* first scan through this level and go down any pointer there are here */
	while(position != NULL)
		{
		if(position->daughter != NULL)
			dismantle_tree(position->daughter);
			
		position = position->next_sibling;
		}
	position = start;
	while(position->next_sibling != NULL) position= position->next_sibling;
	
	/* now remove every thing at this level */
	while(position != NULL)
		{
		
		if(count_now)malloc_check--;
		start = position;
		position = position->prev_sibling;
		if(start->donor != NULL) free(start->donor);
		if(start->fullname != NULL) free(start->fullname);
		free(start);
		}	
	}		


void bootstrap_search(void)
    {
    int i=0, j=0, k=0, l=0, random_num = 0, error = FALSE, *taxa_present = NULL, missing_method = 1;
    int Nreps = 100, search = 1, allpresent = TRUE, num_results = 0;
    char filename[1000], best_tree[1000], **bootstrap_results = NULL, consensusfilename[1000];
    FILE *bootfile = NULL, *temp = NULL, *consensusfile = NULL;
	float percentage = .5;

	consensusfilename[0] = '\0';
    filename[0] = '\0';
    strcpy(filename, "bootstrap.txt");
	strcpy(consensusfilename, "consensus.ph");
    
    for(i=0; i<num_commands; i++)
        {
        if(strcmp(parsed_command[i], "nreps") == 0)
            {
            Nreps = toint(parsed_command[i+1]);
            if(Nreps == 0)
                {
                printf("Error: '%s' is an invalid number of repetitions\n", parsed_command[i+1]);
                error = TRUE;
                }
            }
        if(strcmp(parsed_command[i], "swap") == 0)
            {
            if(strcmp(parsed_command[i+1], "all") == 0) search = 0;
            else
                {
                if(strcmp(parsed_command[i+1], "nni") == 0) search = 1;
                else
                    {
					if(strcmp(parsed_command[i+1], "spr") == 0) search = 1;   /* should be set to 1 when SPR is implemented  NOTE: THE HEURISTIC SEARCH ALGORITHM WILL PICK UP THE SETTINGS AND IMPLEMENT THE CORRECT SEARCH METHOD*/
                    else
                        {
						if(strcmp(parsed_command[i+1], "tbr") == 0) search = 1;   /* should be set to 1 when tbr is implemented */
						else
							{
							printf("Error: swap option '%s' unknown\n", parsed_command[i+1]);
							error = TRUE;
							}
						}
                    }
                }
            }
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
					printf("Error: '%s' is an invalid method to be assigned to missing\n", parsed_command[i+1]);
					missing_method = 1;
					error = TRUE;
					}
				}
			}

			
		if(strcmp(parsed_command[i], "consensus") == 0)
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
							printf("Error: The cut off for a consensus tree must be between .5 and 1.0\n");
							error = TRUE;
							}
						if(percentage == 0)
							{
							printf("Error: consensus option %s unknown\n", parsed_command[i+1]);
							error = TRUE;
							}
						}
					}
				}
			}
		 if(strcmp(parsed_command[i], "consensusfile") == 0)
            {
            strcpy(consensusfilename, parsed_command[i+1]);
            if(filename[0] == '\0')
                {
                printf("Error: Unable to open file named '%s'\n", parsed_command[i+1]);
                error = TRUE;
                strcpy(filename, "consensus.ph");
                }
            }
			
        if(strcmp(parsed_command[i], "treefile") == 0)
            {
            strcpy(filename, parsed_command[i+1]);
            if(filename[0] == '\0')
                {
                printf("Error: Unable to open file named '%s'\n", parsed_command[i+1]);
                error = TRUE;
                strcpy(filename, "bootstrap.txt");
                }
            }
        }


    if(!error)
        {

        printf("\n\nBootstrap Settings:\n\tNumber of Bootstrap replicates = %d\n\tSearching Supertree-space:", Nreps);
        if(search==0)printf(" Exhaustive search\n");
        if(search==1)printf(" Heuristic search ");
	
        printf("\tBootstrap tree file: %s\n", filename);
		if(percentage == 1)
			printf("\tStrict ");
		if(percentage == .5)
			printf("\tMajority rule ");
		if(percentage < 0.5)
			printf("\tMajority rule with minor components ");
		if(percentage > 0.5 && percentage != 1)
			printf("%f cut-off ", percentage);
		printf("consensus tree is to be constructed\n");
		printf("\tConsensus file name = %s\n\n\n", consensusfilename);
		
        taxa_present = malloc(number_of_taxa*sizeof(int));  /** This is so that we can check at each replicate whether or not every taxa is included **/
        for(i=0; i<number_of_taxa; i++)
            taxa_present[i] = FALSE;

        if(!calculated_fund_scores && criterion == 0)
			{
			cal_fund_scores(FALSE);    /* calculate the path metrics for all the fundamental trees */
			calculated_fund_scores = TRUE;
			}
        
        bootfile = fopen(filename, "w");
        /********* initalise the array to store the arrays *********/
        stored_num_comparisons = malloc(Total_fund_trees*sizeof(int));
        if(!stored_num_comparisons) memory_error(60);
        for(i=0; i<Total_fund_trees; i++) stored_num_comparisons[i] = number_of_comparisons[i];
        
        stored_fund_scores = malloc(Total_fund_trees*sizeof(int**));
        if(stored_fund_scores == NULL) memory_error(38);
            
        for(i=0; i<Total_fund_trees; i++)
            {
            stored_fund_scores[i] = malloc((number_of_taxa)*sizeof(int*));
            if(stored_fund_scores[i] == NULL) memory_error(39);
                
            else
                {
                for(j=0; j<(number_of_taxa); j++)
                    {
                    stored_fund_scores[i][j] = malloc((number_of_taxa)*sizeof(int));
                    if(stored_fund_scores[i][j] == NULL) memory_error(40);
                        
                    else
                        {
                        for(k=0; k<(number_of_taxa); k++)
                            stored_fund_scores[i][j][k] = fund_scores[i][j][k];  /* copy the origninal fundamental scores for safe keeping */
                        }
                    }
                }
            }
    
        stored_presence_of_taxa = malloc(Total_fund_trees*sizeof(int *));
        if(!stored_presence_of_taxa) memory_error(41);
            
        for(i=0; i<Total_fund_trees; i++)
            {
            stored_presence_of_taxa[i] = malloc((number_of_taxa)*sizeof(int));
            if(!stored_presence_of_taxa[i]) memory_error(42);
                
            for(j=0; j<number_of_taxa; j++)
                {
                stored_presence_of_taxa[i][j] = presence_of_taxa[i][j];
                }
            }
    
    
        
        stored_funds = malloc(Total_fund_trees*sizeof(char *));
        if(!stored_funds) memory_error(31);
    
        for(i=0; i<Total_fund_trees; i++)
            {
            stored_funds[i] = malloc((strlen(fundamentals[i])+100)*sizeof(char));
            if(!stored_funds) memory_error(32);
            stored_funds[i][0] = '\0';
            strcpy(stored_funds[i], fundamentals[i]);   /* copy the original fundamentals for safe keeping */
    
            }
        
		bootstrap_results = malloc(1*sizeof(char*));
		score_of_bootstraps = malloc(1*sizeof(float));
		bootstrap_results[0] = malloc(TREE_LENGTH*sizeof(char));
		bootstrap_results[0][0] = '\0';

		
        
        if(criterion == 0 || criterion == 2 || criterion == 3)     /* if the criterion used id MSS or MRC **/
            {
            /*****************************************/
            hsprint = TRUE;
            printf("Bootstrap progress indicator: "); 
            for(i=0; i<Nreps; i++)  /* for all reps of the bootstrap algorithm */
                {
				if(!user_break)
					{
					if(i > 0) hsprint = FALSE;
					/*** initialise everything **/
					for(j=0; j<Total_fund_trees; j++)
						{
						strcpy(fundamentals[j], "");
						for(k=0; k<number_of_taxa; k++)
							{
							presence_of_taxa[j][k] = FALSE;
							for(l=0; l<number_of_taxa; l++)
								{
								fund_scores[j][k][l] = 0;
								}
							}
						}
					
					
					/***** Create the bootstrap replicate of fundamental trees  **********/
					for(j=0; j<Total_fund_trees; j++)
						{
						random_num = (int)fmod(rand(), Total_fund_trees);
						fundamentals[j] = realloc(fundamentals[j], (strlen(stored_funds[random_num])+100)*sizeof(char));
						fundamentals[j][0] = '\0';
						strcpy(fundamentals[j],  stored_funds[random_num]);
						number_of_comparisons[j] = stored_num_comparisons[random_num];
						
						for(k=0; k<number_of_taxa; k++)
							{
							presence_of_taxa[j][k] = stored_presence_of_taxa[random_num][k];
							for(l=0; l<number_of_taxa; l++)
								{
								fund_scores[j][k][l] = stored_fund_scores[random_num][k][l];
								}
							}
								
						}
						
					/****** Test whether or not every taxa is represented in this replicate ******/
					
					allpresent = TRUE;
					for(j=0; j<number_of_taxa; j++)
						taxa_present[j] = FALSE;
					for(j=0; j<Total_fund_trees; j++)
						{
						for(k=0; k<number_of_taxa; k++)
							{
							if(presence_of_taxa[j][k] > 0)
								taxa_present[k] = TRUE;
							}
						}
					
					for(j=0; j<number_of_taxa; j++)
						{
						if(taxa_present[j] == FALSE)
							{
							printf("\n\trepetition %d:\tTaxon %s is not present\n", i+1, taxa_names[j]);
							allpresent=FALSE;
							}
						}
					

					
					
					/***** End Test *****/
					if(allpresent)
						{
						/* Now find the best tree for this bootstrapped set of fundamental trees */
						printf("\n\trepetition %d:", i+1);
						fflush(stdout);
						if(search == 0)
							{
							alltrees_search(FALSE);
							}
						else
							{
							heuristic_search(FALSE, FALSE, 10000, 10);
							}
						
						best_tree[0] = '\0';
						k=0; l=0;
						while(scores_retained_supers[l] != -1)
							{
							l++;
							}
						bootstrap_results = realloc(bootstrap_results, (l+num_results)*sizeof(char*));
						score_of_bootstraps = realloc(score_of_bootstraps, (l+num_results)*sizeof(float));
						for(j=num_results; j<l+num_results; j++) 
								{
								bootstrap_results[j] = malloc(TREE_LENGTH*sizeof(char));
								bootstrap_results[j][0] = '\0';
								}
						while(scores_retained_supers[k] != -1)
							{
							if(search ==0)
								{
								if(tree_top != NULL)
									{
									dismantle_tree(tree_top);
									tree_top = NULL;
									}
								temp_top = NULL;
								tree_build(1, retained_supers[k], tree_top, FALSE, -1);
								tree_top = temp_top;
								temp_top = NULL;
							
								strcpy(best_tree, "");
								print_named_tree(tree_top, best_tree);
							   /* printf("\n%s[%f];\t[%f]\n", best_tree, best_tree,(float)((float)1/(float)l), scores_retained_supers[k]);  */
								fprintf(bootfile, "%s [%f];\t[score = %f]\n", best_tree,(float)((float)1/(float)l), scores_retained_supers[k]);
								}
							else
								{
							   /* printf("\n%s\t[%f]\n", retained_supers[k], scores_retained_supers[k]); */
								j=0;
								strcpy(best_tree, "");
								while(retained_supers[k][j] != ';')
									{
									best_tree[j] = retained_supers[k][j];
									j++;
									}
								best_tree[j] = '\0';
								/*printf("\n%s[%f];\t[%f]\n", best_tree, (float)((float)1/(float)l), scores_retained_supers[k]); */
								fprintf(bootfile, "%s [%f];\t[score = %f]\n", best_tree, (float)((float)1/(float)l), scores_retained_supers[k]);
								}
							
							strcpy(bootstrap_results[num_results+k], retained_supers[k]);
							score_of_bootstraps[num_results+k] = (float)((float)1/(float)l);
							
							fflush(bootfile);
				
							k++;
							}
						num_results +=l;
						for(k=0; k<number_retained_supers; k++)
							{
							strcpy(retained_supers[k], "");
							scores_retained_supers[k] = -1;
							}
					   }
					else
						{
						i=Nreps;
						printf("\n\nError: This data may not be suitable for bootstrapping due to the low\noccurrences of some taxa in the source trees. Please check the\noccurance summary to identify problematic taxa\n");
						}
					}
				}
            printf("\n");
            fclose(bootfile);
			
            
            } /** end criterion 0 & 2 & 3 **/
    if(criterion == 1 || criterion == 4)
            {
            /** MRP bootstrap **/
            BR_file = fopen("coding.nex", "w");
            for(i=0; i<Nreps; i++)
                {
                 /*** initialise everything **/
                for(j=0; j<Total_fund_trees; j++)
                    {
                    strcpy(fundamentals[i], "");
                    for(k=0; k<number_of_taxa; k++)
                        {
                        presence_of_taxa[j][k] = FALSE;
                        for(l=0; l<number_of_taxa; l++)
                            {
                            fund_scores[j][k][l] = 0;
                            }
                        }
                    }
                
                /***** Create the bootstrap replicate of fundamental trees  **********/
                for(j=0; j<Total_fund_trees; j++)
                    {
                    random_num = (int)fmod(rand(), Total_fund_trees);
					fundamentals[j] = realloc(fundamentals[j], (strlen(stored_funds[random_num])+100)*sizeof(char));
					fundamentals[j][0] = '\0';
                    strcpy(fundamentals[j],  stored_funds[random_num]);
                    number_of_comparisons[j] = stored_num_comparisons[random_num];
                    
                    for(k=0; k<number_of_taxa; k++)
                        {
                        presence_of_taxa[j][k] = stored_presence_of_taxa[random_num][k];
                        for(l=0; l<number_of_taxa; l++)
                            {
                            fund_scores[j][k][l] = stored_fund_scores[random_num][k][l];
                            }
                        }
                            
                    }
                    
                /****** Test whether or not every taxa is represented in this replicate ******/
                
                allpresent = TRUE;
                for(j=0; j<number_of_taxa; j++)
                    taxa_present[j] = FALSE;
                for(j=0; j<Total_fund_trees; j++)
                    {
                    for(k=0; k<number_of_taxa; k++)
                        {
                        if(presence_of_taxa[j][k] > 0)
                            taxa_present[k] = TRUE;
                        }
                    }
                
                for(j=0; j<number_of_taxa; j++)
                    {
                    if(taxa_present[j] == FALSE)
                        {
                        printf("\n\trepetition %d:\tTaxon %s is not present\n", i+1, taxa_names[j]);
                        allpresent=FALSE;
                        }
                    }
                

                
                /***** End Test *****/
                
                if(allpresent)
                    {
                
                    if(criterion == 1)
						{
						if(Nreps == 1)
							{
							j = coding(0, search, 0);
							if(j == 1) error = TRUE;
							}
						else
							{
							if(i==0)
								{
								j = coding(1, search, 0);
								if(j == 1) error = TRUE;
								}
							if(i>0 && i<Nreps-1)
								{
								j = coding(2, search, 0);
								if(j == 1) error = TRUE;
								}
							if(i== Nreps-1)
								{
								j = coding(3, search, 0);
								if(j == 1) error = TRUE;
								}
							}
						if(error == 1)
							{
							printf("error\n");
							i=Nreps;
							}	
						}
					if(criterion == 4)
						{
						if(Nreps == 1)
							error = average_consensus(0, missing_method, filename, BR_file);
						else
							{
							if(i==0)
								error = average_consensus(1, missing_method, filename, BR_file);
							if(i>0 && i<Nreps-1)
								error = average_consensus(2, missing_method, filename, BR_file);
							if(i== Nreps-1)
								error = average_consensus(3, missing_method, filename, BR_file);
							}
						if(error) i=Nreps;
						}
                    }
                else
                    {
                    i=Nreps;
                    printf("\n\nError: This data may not be suitable for bootstrapping due to the low\noccurrences of some taxa in the source trees. Please check the\nCo-occurance summary to identify problematic taxa\n");
                    error = TRUE;
		    }
                }
            fclose(BR_file);
            if(!error && allpresent)
                {
                if(system("paup coding.nex") != 0) printf("Error calling paup, please execute the file coding.nex in paup to perform the parsimony step\n");
				}
            
            }
		if(!error && criterion != 1 && criterion != 4)
			{
			/***** Do a consensus of the results ******/
			consensusfile = fopen(consensusfilename, "w");
			consensus(num_results, bootstrap_results, Nreps, percentage, consensusfile, NULL); 
			fclose(consensusfile);
			}
			
		/* copy back the arrays to their original places */
		
		for(i=0; i<Total_fund_trees; i++)
			{
			
			fundamentals[i] = realloc(fundamentals[i], (strlen(stored_funds[i])+100)*sizeof(char));
			fundamentals[i][0] = '\0';
			strcpy(fundamentals[i], stored_funds[i]);
			number_of_comparisons[i] = stored_num_comparisons[i];
			for(j=0; j<number_of_taxa; j++)
				{
				presence_of_taxa[i][j] = stored_presence_of_taxa[i][j];
				for(k=0; k<number_of_taxa; k++)
					{
					fund_scores[i][j][k] = stored_fund_scores[i][j][k];
					}
				}
			}
			
        /* free all the memory allocated */
        hsprint = TRUE;
        
        for(i=0; i<Total_fund_trees; i++)
            {
            free(stored_funds[i]);
            free(stored_presence_of_taxa[i]);
            for(j=0; j<number_of_taxa; j++)
                {
                    free(stored_fund_scores[i][j]);
                }
            free(stored_fund_scores[i]);
            }
		if(bootstrap_results)
			{
			for(i=0; i<num_results; i++)
				free(bootstrap_results[i]);
			free(bootstrap_results);
			free(score_of_bootstraps);
			}
        free(stored_funds);
        free(stored_presence_of_taxa);
        free(stored_fund_scores);
        free(stored_num_comparisons);
        free(taxa_present);
        }
    }








void memory_error(int error_num)  /*123 so far*/
    {
    printf("Error: Out of memory @ %d\n", error_num);
    clean_exit(1);
    }




void print_named_tree(struct taxon * position, char *tree)
	{
	struct taxon *place = position;
	int count = 0, j=0;
        char *name = NULL;
	strcat(tree, "(");
        
	name = malloc(30*sizeof(char));
	name[0] = '\0';
	while(position != NULL)
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


void print_tree(struct taxon * position, char *tree)
	{
	struct taxon *place = position;
	int count = 0, j=0;
        char *name = NULL;
	strcat(tree, "(");
        
	name = malloc(30*sizeof(char));
	name[0] = '\0';
	while(position != NULL)
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
		/*	printf("%d\n", name);*/
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



void usertrees_search(void)
    {
    FILE *userfile = NULL, *outfile = NULL, *sourcescoresfile = NULL;
    int keep = 0, nbest = 0, error = FALSE, i=0, tree_number = 0, j=0, k=0, prev = 0, print_source_scores = FALSE;
    char *user_super = NULL, c = '\0', best_tree[TREE_LENGTH], *temp = NULL;
    float score = 0, best_score = 0;
    
    if((userfile = fopen(parsed_command[1], "r")) == NULL)
        {
        printf("Error opening file named %s\n", parsed_command[1]);
        error = TRUE;
        }
    else
        {
        for(i=0; i<num_commands; i++)
            {
			if(strcmp(parsed_command[i], "printsourcescores") == 0)
                {
				if(strcmp(parsed_command[i+1], "yes") == 0)
					{
					print_source_scores = TRUE;
					sourcescoresfile = fopen("sourcetree_scores.txt", "w");
					}
				}
            if(criterion == 3)
                {
                if(strcmp(parsed_command[i], "weight") == 0)
                    {
                    if(strcmp(parsed_command[i+1], "equal") == 0)
                        quartet_normalising = 1;
                    else
                        {
                        if(strcmp(parsed_command[i+1], "taxa") == 0)
                            quartet_normalising = 2;
                        else
                            {
                            if(strcmp(parsed_command[i+1], "quartets") == 0)
                                quartet_normalising = 3;
                            else
                                {
                                printf("Error: weight option '%s' is unknown\n", parsed_command[i+1]);
                                error = TRUE;
                                }
                            }
                        }
                    }
                }  
            if(criterion == 0)
                {
                if(strcmp(parsed_command[i], "weight") == 0)
                    {
                    if(strcmp(parsed_command[i+1], "equal") == 0)
                        dweight = 0;
                    else
                        {
                        if(strcmp(parsed_command[i+1], "comparisons") == 0)
                            dweight = 1;
                        else
                            {
                            printf("Error: weight option '%s' is unknown\n", parsed_command[i+1]);
                            error = TRUE;
                            }
                        }
                    }
                }            
            if(criterion == 2)
                {
                if(strcmp(parsed_command[i], "weight") == 0)
                    {
                    if(strcmp(parsed_command[i+1], "equal") == 0)
                        splits_weight = 1;
                    else
                        {
                        if(strcmp(parsed_command[i+1], "splits") == 0)
                            splits_weight = 2;
                        else
                            {
                            printf("Error: weight option '%s' is unknown\n", parsed_command[i+1]);
                            error = TRUE;
                            }
                        }
                    }
                }            

            if(strcmp(parsed_command[i], "outfile") == 0)
                {
                if((outfile = fopen(parsed_command[i+1], "w")) == NULL)
                    {
                    printf("Error opening output file named: %s\n", parsed_command[i+1]);
                    error = TRUE;
                    }
                }

            }
        if(outfile == NULL)
            {
            if((outfile = fopen("Usertrees_result.txt", "w")) == NULL)
                {
                printf("Error opening output file  Usertrees_result.txt\n");
                error = TRUE;
                }
            }
        }
    if(criterion == 1 || criterion == 4 || criterion == 5)
		{
		printf("ERROR: Usertree search is not available under the criteria ");
		switch(criterion)
			{
			case 1:
				printf("Matrix Representation Using Parsimony (MRP)\n");
				break;
			case 4:
				printf("Average consensus (AVCON)\n");
				break;
			default:
				printf("Reconstruction of duplications and losses (RECON)\n");
				break;
			}
		error = TRUE;
		}
    if(!error)
        {   
		
		printf("\nUsertree Search settings:\n");
		printf("\tCriterion = ");
		switch(criterion)
			{
			case 1:
				printf("Matrix Representation Using Parsimony (MRP)\n");
				break;
			case 0:
				printf("Most Similar Supertree (dfit)\n");
				break;
			case 2:
				printf("Maximum split fit (SFIT)\n");
				break;
			case 4:
				printf("Average consensus (AVCON)\n");
				break;
			case 5:
				printf("Reconstruction of duplications and losses (RECON)\n");
				break;
			default:
				printf("Maximum quartet fit (QFIT)\n");
				break;
				
			}
		if(criterion != 1 && criterion != 4)
			{
			if(criterion != 5)
				printf("\tWeighting Scheme = ");
			if(criterion==0)
				{
				if(dweight == 0) printf("equal\n");
				if(dweight == 1) printf("comparisons\n");
				}
			if(criterion == 2)
				{
				if(splits_weight == 1) printf("equal\n");
				if(splits_weight == 2) printf("splits\n");
				}
			if(criterion == 3)
				{
				if(quartet_normalising == 1) printf("equal\n");
				if(quartet_normalising == 2) printf("taxa\n");
				if(quartet_normalising == 3) printf("quartets\n");
				}
			printf("\tUser defined supertrees from file: %s\n", parsed_command[1]);
			
			}
		
		if(criterion==5) printf("\n\tDuplication weight = %f\n\tLosses weight = %f\n\n", dup_weight, loss_weight);
		
		if(print_source_scores)printf("\nScores of each source tree compared to the best user supertree to be printed to sourcetree_scores.txt\n");
		
		
		
		
		
		
        /* this hold the distances calculated with the pathmetric on the supertree */
        if(!calculated_fund_scores && criterion == 0)
			{
			cal_fund_scores(FALSE);    /* calculate the path metrics for all the fundamental trees */
			calculated_fund_scores = TRUE;
			}
		for(i=0; i<number_of_taxa; i++)presenceof_SPRtaxa[i] = -1;
		
        if(super_scores == NULL)
            {
            super_scores = malloc(number_of_taxa*sizeof(int *));
            if(!super_scores)  memory_error(27);
                
            for(i=0; i<number_of_taxa; i++)
                {
                super_scores[i] = malloc(number_of_taxa*sizeof(int));
                if(!super_scores[i])  memory_error(28);
                    
                for(j=0; j<number_of_taxa; j++)
                    {
                    super_scores[i][j] = 0;
                    }
                }
            }
        else
            {
            for( i=0; i<number_of_taxa; i++)
                {
                for(j=i; j<number_of_taxa; j++)
                    {
                    super_scores[i][j] = 0;
                    super_scores[j][i] = 0;
                    }
                }
            }
    
        psfile = fopen("supertree.ps", "w");

    
        temp = malloc(TREE_LENGTH*sizeof(char));
        if(!temp) memory_error(50);
        temp[0] = '\0';
        
        user_super = malloc(TREE_LENGTH*sizeof(char));
        if(!user_super) memory_error(51);
        user_super[0] = '\0';
    
         for(i=0; i<number_retained_supers; i++)
			{
			strcpy(retained_supers[i], "");
			scores_retained_supers[i] = -1;
			}
    
        if(criterion == 2 || criterion == 3)
            {
            coding(0,0,0);
            if(criterion == 2) condense_coding();
            }
            
        c = getc(userfile);
        while((c == ' ' || c == '\r'  || c == '\n' || c == '[') && !feof(userfile))  /* skip past all the non tree characters */
            {
            if(c == '[')
                {
                while((c = getc(userfile)) != ']' && !feof(userfile));
                }
            c = getc(userfile);
            }
            
		while(!feof(userfile)) /* While we haven't reached the end of the file */
			{
			strcpy(user_super, "");
			i=0;
			while(!feof(userfile) && c != ';' )
				{
				user_super[i] = c;
				i++;
				c = getc(userfile);

				}  /* End of tree */
			user_super[i] = ';';
			
			/* next score this tree */
			/****** We now need to build the Supertree in memory *******/
			
			if(tree_top != NULL)
				{
				dismantle_tree(tree_top);
				tree_top = NULL;
				}
			temp_top = NULL;
			tree_build(1, user_super, tree_top, TRUE, -1);
			tree_top = temp_top;
			temp_top = NULL;
			/*check_tree(tree_top); */
			if(criterion == 0) score = compare_trees(FALSE);
			if(criterion == 2)  /* calculate the distance using the MRC criterion */
				{
				strcpy(temp, "");
				print_tree(tree_top, temp);
				strcat(temp, ";");
				unroottree(temp);
				score = (float)MRC(temp);
				}
			if(criterion == 3)
				{
				strcpy(temp, "");
				print_tree(tree_top, temp);
				strcat(temp, ";");
				unroottree(temp);
				score = (float)quartet_compatibility(temp);
				}
							
			/*printf("%s\t[%f]\n", user_super, score); */
			
			c = getc(userfile);
			while((c == ' ' || c == '\r'  || c == '\n'|| c == '\t' || c == '[') && !feof(userfile))  /* skip past all the non tree characters */
				{
				if(c == '[')
					{
					while((c = getc(userfile)) != ']' && !feof(userfile));
					}
				c = getc(userfile);
				}

			if(tree_number==0)
				{
				tree_number++;
				retained_supers[0] = realloc(retained_supers[0], (strlen(user_super)+10)*sizeof(char));
				strcpy(retained_supers[0], user_super);
				scores_retained_supers[0] = score;
				best_score = score;
				}
			else
				{
				tree_number++;
				if(score < best_score)
					{
					retained_supers[0] = realloc(retained_supers[0], (strlen(user_super)+10)*sizeof(char));
					strcpy(retained_supers[0], user_super);
					scores_retained_supers[0] = score;
					best_score = score;
					j=1;
					while(scores_retained_supers[j] != -1 && j < number_retained_supers)
						{
						strcpy(retained_supers[j], "");
						scores_retained_supers[j] = -1;
						j++;
						}
					}
				else
					{
					if(score == best_score)
						{
						j = 1;
						while(scores_retained_supers[j] != -1)
							{
							j++;
							if(j+1 == number_retained_supers) reallocate_retained_supers();
							}
						retained_supers[j] = realloc(retained_supers[j], (strlen(user_super)+10)*sizeof(char));
						strcpy(retained_supers[j], user_super);
						scores_retained_supers[j] = score;
						}
					}
				}
		   
			} /* end of file */

        /**** Print out the best trees found *******/
            
            i=0; j=0;
			while(scores_retained_supers[j] != -1)
				{
				j++;
				}
			
			while(scores_retained_supers[i] != -1)
				{
				if(print_source_scores && criterion==0)
					{
					fprintf(sourcescoresfile, "Scores of sources trees compared to best User Supertree %d\n\n", i+1);
					if(tree_top != NULL)
						{
						dismantle_tree(tree_top);
						tree_top = NULL;
						}
					temp_top = NULL;

					temp_top = NULL;
					tree_build(1, retained_supers[i], tree_top, TRUE, -1);
					tree_top = temp_top;
					temp_top = NULL;
					/*check_tree(tree_top); */
					if(criterion == 0) score = compare_trees(FALSE);
					for(k=0; k<Total_fund_trees; k++)
						{
						if(strcmp(tree_names[k], "") != 0)
							fprintf(sourcescoresfile, "sourcetree %s:\t", tree_names[k]);
						else
							fprintf(sourcescoresfile, "sourcetree %d:\t", k);
						fprintf(sourcescoresfile, "%f\n", sourcetree_scores[k]);
						}
					}
				
				if(outfile != NULL) fprintf(outfile, "%s\t[%f]\n", retained_supers[i], scores_retained_supers[i] );
				tree_coordinates(retained_supers[i], FALSE, TRUE, FALSE, -1);
				printf("\nSupertree %d of %d score = %f\n", i+1, j, scores_retained_supers[i] );
				i++;
				}

			trees_in_memory = i;
		
        
            fclose(userfile);
            if(outfile != NULL) fclose(outfile);
            free(temp);
            free(user_super);
            fclose(psfile);
			}
	if(sourcescoresfile != NULL) fclose(sourcescoresfile);
    }
        
void controlc1(int signal)
	{
	char *c = NULL;
	
	printf("\n\nCompleted %d random samples before CTRL-C was caught\nDo you want to stop the random sampling and start the heuristuc searches now? (Y/N): \n", GC);
	c = malloc(10000*sizeof(char));
	xgets(c);
	while(strcmp(c, "Y") != 0 && strcmp(c, "y") != 0 && strcmp(c, "N") != 0 && strcmp(c, "n") != 0)
		{
		printf("Error '%s' is not a valid response, please respond Y or N\n", c);
		xgets(c);
		}
	
	if(strcmp(c, "Y") == 0 || strcmp(c, "y") == 0)
		{
		GC = 100000000;
		printf("Starting Heuristic searches\n");
		}

	if(strcmp(c, "N") == 0 || strcmp(c, "n") == 0)
			printf("Continuing random samples.....\n");
	
	free(c);

	}

void controlc2(int signal)
	{
	char *c = NULL;
	
	c = malloc(10000*sizeof(char));
	
	printf("\n\n\n\nCompleted the assessment of %d trees...best tree found so far has a score of %f\nDo you want to stop the heuristic searches now? (Y/N): \n", NUMSWAPS, BESTSCORE);
	xgets(c);
	while(strcmp(c, "Y") != 0 && strcmp(c, "y") != 0 && strcmp(c, "N") != 0 && strcmp(c, "n") != 0)
		{
		printf("Error '%s' is not a valid response, please respond Y or N\n", c);
		xgets(c);
		}
	
	if(strcmp(c, "Y") == 0 || strcmp(c, "y") == 0)
		{
		user_break = TRUE;
		printf("Stopping heuristic search. The resulting tree may be sub-optimal\n\n");
		}


	if(strcmp(c, "N") == 0 || strcmp(c, "n") == 0)
			printf("Continuing heuristic searches.....\n");
	
	free(c);
	}

void controlc3(int signal)
	{
	char *c = NULL;
	
	c = malloc(10000*sizeof(char));
	
	printf("\n\nDo you want to stop the generation of trees now? (Y/N): ");
	xgets(c);
	while(strcmp(c, "Y") != 0 && strcmp(c, "y") != 0 && strcmp(c, "N") != 0 && strcmp(c, "n") != 0)
		{
		printf("Error '%s' is not a valid response, please respond Y or N\n", c);
		xgets(c);
		}
	
	if(strcmp(c, "Y") == 0 || strcmp(c, "y") == 0)
		{
		user_break = TRUE;
		printf("Stopping the generation of treesn\n");
		}


	if(strcmp(c, "N") == 0 || strcmp(c, "n") == 0)
			printf("Continuing the generation of trees.....\n");
	
	free(c);
	}

void controlc4(int signal)
	{
	char *c = NULL;
	
	c = malloc(10000*sizeof(char));
	
	printf("\n\nDo you really wish to quit? (Y/N): ");
	xgets(c);
	while(strcmp(c, "Y") != 0 && strcmp(c, "y") != 0 && strcmp(c, "N") != 0 && strcmp(c, "n") != 0)
		{
		printf("Error '%s' is not a valid response, please respond Y or N\n", c);
		xgets(c);
		}
	
	if(strcmp(c, "Y") == 0 || strcmp(c, "y") == 0)
		{
		user_break = TRUE;
		printf("Bye\n");
		clean_exit(0);
		}
	
	free(c);
	}

void controlc5(int signal)
	{
	char *c = NULL;
	
	c = malloc(10000*sizeof(char));
	
	printf("\n\nDo you want to stop the exhaustive searches now? (Y/N): ");
	xgets(c);
	while(strcmp(c, "Y") != 0 && strcmp(c, "y") != 0 && strcmp(c, "N") != 0 && strcmp(c, "n") != 0)
		{
		printf("Error '%s' is not a valid response, please respond Y or N\n", c);
		xgets(c);
		}
	if(strcmp(c, "Y") == 0 || strcmp(c, "y") == 0)
		{
		user_break = TRUE;
		printf("Stopping exhaustive search. The resulting tree may be sub-optimal\n\n");
		}


	if(strcmp(c, "N") == 0 || strcmp(c, "n") == 0)
			printf("Continuing heuristic searches.....\n");
	
	free(c);
	}


void heuristic_search(int user, int print, int sample, int nreps)
    {
    int i=0, j=0, k=0, l=0, swaps = 0, keep = 0, nbest = 0, start = 2, error = FALSE, numswaps = 1000000, different=TRUE, do_histogram = FALSE, here = FALSE, **taxa_comp = NULL, bins = 20, found = FALSE, missing_method = 1, random_num, numspectries = 2, numgenetries = 2;
    char *tree = NULL, c = '\0', *best_tree = NULL, *temptree = NULL, **starths = NULL, userfilename[10000], useroutfile[10000], histogramfile_name[10000];
    FILE *userfile = NULL, *outfile = NULL, *paupfile = NULL, *histogram_file = NULL;
    float distance=0, number=0, *startscores = NULL, used_weights = 0;
    
	
	for(i=0; i<number_of_taxa; i++) presenceof_SPRtaxa[i] = -1;
	
    best_tree = malloc(TREE_LENGTH*sizeof(char));
    if(!best_tree) memory_error(75);
    best_tree[0] = '\0';
    
    temptree = malloc(TREE_LENGTH*sizeof(char));
    if(!temptree) memory_error(76);
    temptree[0] = '\0';
    
    userfilename[0] = '\0';
    useroutfile[0] = '\0';
	histogramfile_name[0] = '\0';
	strcpy(histogramfile_name, "Heuristic_histogram.txt");
    strcpy(useroutfile, "Heuristic_result.txt");
	
	/*if(criterion == 5) start = 0; */ /* use random starting trees by default for the recon criterion */
	
    /* this hold the distances calculated with the pathmetric on the supertree */
	if(criterion != 5)
		{
		if(!calculated_fund_scores && criterion == 0)
			{
			cal_fund_scores(FALSE);    /* calculate the path metrics for all the fundamental trees */
			calculated_fund_scores = TRUE;
			}
		}


    for(i=0; i<num_commands; i++)
        {
        if(strcmp(parsed_command[i], "numspeciesrootings") == 0)
			{
			if(strcmp(parsed_command[i+1], "all") == 0)
				numspectries=-1;
			else
				{
				numspectries = toint(parsed_command[i+1]);
				if(numspectries == 0)
					{
					printf("Error: '%s' is an invalid value to be assigned to numspeciesrootings\n", parsed_command[i+1]);
					numspectries = 2;
					error = TRUE;
					}
				}
			}
        if(strcmp(parsed_command[i], "numgenerootings") == 0)
			{
			if(strcmp(parsed_command[i+1], "all") == 0)
				numgenetries=-1;
			else
				{
				numgenetries = toint(parsed_command[i+1]);
				if(numgenetries == 0)
					{
					printf("Error: '%s' is an invalid value to be assigned to numgenerootings\n", parsed_command[i+1]);
					numgenetries = 2;
					error = TRUE;
					}
				}
			}
		}



    for(i=0; i<num_commands; i++)
        {
        if(strcmp(parsed_command[i], "missing") == 0 && criterion != 5)
			{
			if(strcmp(parsed_command[i+1], "4point") == 0)
				missing_method = 1;
			else
				{
				if(strcmp(parsed_command[i+1], "ultrametric") == 0)
					missing_method = 0;
				else
					{
					printf("Error: '%s' is an invalid method to be assigned to missing\n", parsed_command[i+1]);
					missing_method = 1;
					error = TRUE;
					}
				}
			}
			
        if(strcmp(parsed_command[i], "drawhistogram") == 0)
            {
            if(strcmp(parsed_command[i+1], "yes") == 0)
				do_histogram = TRUE;
			else
				{
				if(strcmp(parsed_command[i+1], "no") == 0);
				else
					{
					printf("Error: '%s' is an invalid value for drawhistogram\n", parsed_command[i+1]);
					error = TRUE;
					}
				}
            }
		if(strcmp(parsed_command[i], "nbins") == 0)
			{
			bins = toint(parsed_command[i+1]);
			if(bins == 0)
				{
				printf("Error: '%s' is an invalid integer number to be assigned to nbins\n", parsed_command[i+1]);
				bins = 20;
				error = TRUE;
				}
			}
		
		
		if(strcmp(parsed_command[i], "histogramfile") == 0)
            {
			strcpy(histogramfile_name, parsed_command[i]);
			}
			
        if(strcmp(parsed_command[i], "nsteps") == 0)
            {
            number_of_steps = toint(parsed_command[i+1]);
            if(number_of_steps == 0)
                {
                printf("Error: '%s' is an invalid integer number to be assigned to nsteps\n", parsed_command[i+1]);
                number_of_steps = 5;
                error = TRUE;
                }
            }
        if(strcmp(parsed_command[i], "nbest") == 0)
            {
            nbest = toint(parsed_command[i+1]);
            if(nbest == 0)
                {
                printf("Error: '%s' is an invalid integer number to be assigned to nbest\n", parsed_command[i+1]);
                error = TRUE;
                }
            }
        if(strcmp(parsed_command[i], "maxswaps") == 0)
            {
            numswaps = toint(parsed_command[i+1]);
            if(numswaps == 0)
                {
                printf("Error: '%s' is an invalid integer number to be assigned to maxswaps\n", parsed_command[i+1]);                
                error = TRUE;
                }
            }
        

        if((strcmp(parsed_command[i], "hsreps") == 0 )&& user == FALSE)  /* the user can assign the number of hs reps from the boot command */
            {
            nreps = toint(parsed_command[i+1]);
            if(nreps == 0)
                {
                printf("Error: '%s' is an invalid integer number to be assigned to nreps\n", parsed_command[i+1]);
                error = TRUE;
                }
            }
        
        if((strcmp(parsed_command[i], "nreps") == 0) && user == TRUE)
            {
            nreps = toint(parsed_command[i+1]);
            
            if(nreps == 0)
                {
                printf("Error: '%s' is an invalid integer number to be assigned to nreps\n", parsed_command[i+1]);
                error = TRUE;
                }
            }
            
        if(strcmp(parsed_command[i], "swap") == 0 && criterion != 1)
            {
            if(strcmp(parsed_command[i+1], "nni") == 0)
                {
                method = 1;
                }
            else
                {
                if(strcmp(parsed_command[i+1], "spr") == 0)
                    method = 2;
                else
                    {
					if(strcmp(parsed_command[i+1], "tbr") == 0)
						method = 3;
					else
						{
						printf("Error: swap option '%s' not known\n", parsed_command[i+1]);
						error = TRUE;
						}
                    }
                }
            }
        if(strcmp(parsed_command[i], "savetrees") == 0 && criterion != 1)
            {
            if(criterion != 4)
                {
                if((outfile = fopen(parsed_command[i+1], "w")) == NULL)
                    {
                    printf("Error opening file named %s\n", parsed_command[i+1]);
                    error = TRUE;
                    strcpy(useroutfile, parsed_command[i+1]);
                    }
                else
                    {
					strcpy(useroutfile, parsed_command[i+1]);
                    printf("opened output file %s\n", parsed_command[i+1]);
                    }
                }
            else
                strcpy(useroutfile, parsed_command[i+1]);
            }
        if(criterion == 0)
            {
            if(strcmp(parsed_command[i], "weight") == 0)
                {
                if(strcmp(parsed_command[i+1], "equal") == 0)
                    dweight = 0;
                else
                    {
                    if(strcmp(parsed_command[i+1], "comparisons") == 0)
                        dweight = 1;
                    else
                        {
                        printf("Error: weight option '%s' is unknown\n", parsed_command[i+1] );
                        error = TRUE;
                        }
                    }
                }
            }            
        if(criterion == 3)
            {
            if(strcmp(parsed_command[i], "weight") == 0)
                {
                if(strcmp(parsed_command[i+1], "equal") == 0)
                    quartet_normalising = 1;
                else
                    {
                    if(strcmp(parsed_command[i+1], "taxa") == 0)
                        quartet_normalising = 2;
                    else
                        {
                        if(strcmp(parsed_command[i+1], "quartets") == 0)
                            quartet_normalising = 3;
                        else
                            {
                            printf("Error: weight option '%s' is unknown\n", parsed_command[i+1] );
                            error = TRUE;
                            }
                        }
                    }
                }
                        
            }
        if(criterion == 2)
            {
            if(strcmp(parsed_command[i], "weight") == 0)
                {
                if(strcmp(parsed_command[i+1], "equal") == 0)
                    splits_weight = 1;
                else
                    {
                    if(strcmp(parsed_command[i+1], "splits") == 0)
                        splits_weight = 2;
                    else
                        {
                        printf("Error: weight option '%s' is unknown\n", parsed_command[i+1]);
                        error = TRUE;
                        }
                    }
                }
            }            
        if(strcmp(parsed_command[i], "start") == 0 && criterion != 1)
            {
            if(strcmp(parsed_command[i+1], "random") == 0)
                start = 0;
			else
				{
				if(strcmp(parsed_command[i+1], "nj") == 0)
					{
					if(criterion != 5)
						start = 2;
					else
						{
						printf("nj is not a valid option for starting trees in the recon criterion\n");
						error = TRUE;
						}
					}
				else
					{
					start = 1;
					if((userfile = fopen(parsed_command[i+1], "r")) == NULL)
						{
						
						printf("Error opening file named %s\n", parsed_command[i+1]);
						error = TRUE;
						
						}
					else
						{
						printf("opened user file %s\n", parsed_command[i+1]);
						strcpy(userfilename, parsed_command[i+1]);
						}
					} 
				}
			}
		if(criterion==5)
			{
			if(strcmp(parsed_command[i], "duplications") == 0)
				{
				dup_weight = tofloat(parsed_command[i+1]);
				if(dup_weight < 0)
					{
					printf("Error: '%s' is an invalid value for dups\n", parsed_command[i+1]);
					error = TRUE;
					}
				}
			if(strcmp(parsed_command[i], "losses") == 0)
				{
				loss_weight = tofloat(parsed_command[i+1]);
				if(loss_weight < 0)
					{
					printf("Error: '%s' is an invalid value for losses\n", parsed_command[i+1]);
					error = TRUE;
					}
				}
			}
		
    
        }
    for(i=0; i<num_commands; i++)   /* This has to be outside the normal for loop because if the assignment of hsreps (or nreps) is made after the assignment of smaple, then an error can be reported */
        {
 
        if(strcmp(parsed_command[i], "sample") == 0)  /* the user can assign the number of hs reps from the boot command */
            {
            sample = toint(parsed_command[i+1]);
            if(sample < nreps)
                {
                printf("Error: '%s' is an invalid integer number to be assigned to sample\n It must be equal to or larger than", parsed_command[i+1]);
                if(user)printf(" nreps (%d)\n", nreps);
                if(!user)printf(" hsreps (%d)\n", nreps);
                error = TRUE;
                }
            }      
        }

	

    if(!error)
        {
		
		if(sample > sup) sample=sup;
		if(nreps > sup) nreps=sup;
		
        if(hsprint==TRUE)
            {
            /***** Print out summary of settings for this run */
            printf("\n\nHeuristic Search settings:\n");
            printf("\tCriterion = ");
            switch(criterion)
                {
                case 1:
                    printf("Matrix Representation Using Parsimony (MRP)\n");
                    break;
                case 0:
                    printf("Most Similar Supertree (dfit)\n");
                    break;
                case 2:
                    printf("Maximum split fit (SFIT)\n");
                    break;
                case 4:
                    printf("Average consensus (AVCON)\n");
                    break;
				case 5:
					printf("Reconstruction of duplications and losses (RECON)\n");
					break;
                default:
                    printf("Maximum quartet fit (QFIT)\n");
                    break;
                    
                }
            if(criterion != 1 && criterion != 4)
                {
                printf("\tHeuristic search algorithm = ");
                if(method == 1) printf("Nearest Neighbor Interchange (NNI)\n");
                if(method == 2) printf("Sub-tree Pruning and Regrafting (SPR)\n");
				if(method == 3) printf("Tree Bisection and Reconnection (TBR)\n");
                printf("\tMaximum Number of Steps (nsteps) = %d\n", number_of_steps);
                printf("\tMaximum Number of Swaps (maxswaps) = %d\n", numswaps);
                printf("\tNumber of repetitions of Heuristic search = %d\n", nreps);
                if(criterion != 5)
					printf("\tWeighting Scheme = ");
                if(criterion==0)
                    {
                    if(dweight == 0) printf("equal\n");
                    if(dweight == 1) printf("comparisons\n");
                    }
                if(criterion == 2)
                    {
                    if(splits_weight == 1) printf("equal\n");
                    if(splits_weight == 2) printf("splits\n");
                    }
                if(criterion == 3)
                    {
                    if(quartet_normalising == 1) printf("equal\n");
                    if(quartet_normalising == 2) printf("taxa\n");
                    if(quartet_normalising == 3) printf("quartets\n");
                    }
                printf("\tStarting trees = ");
                if(start == 0)
                    printf("Top %d random trees chosen from %d random samples\n",nreps, sample);
                if(start == 1)
                    printf("User defined from file named %s\n", userfilename);
				if(start == 2)
					{
					printf("neighbor-joining tree from Average consensus distances\n\tMissing data estimated using ");
					if(missing_method == 1) printf(" 4 point condition distances\n");
					if(missing_method == 0) printf(" ultrametric distances\n");
					}
				if(user) printf("\tOutput file = %s\n", useroutfile);
                }
            }
		if(criterion == 0 && do_histogram == TRUE)
			{
			printf("\tSource tree scores to be plotted against best supertree(s)\n");
			printf("\tNumber of bins to summarise the source tree scores = %d\n", bins);
			printf("\tHistogram to be written to file '%s'\n", histogramfile_name);
			}
		if(criterion==5 && hsprint == TRUE)
			{
			printf("\n\tDuplication weight = %f\n\tLosses weight = %f\n", dup_weight, loss_weight);
			printf("\tNumber of species tree rootings = ");
			if(numspectries == -1) printf("all possible\n");
			else printf("%d\n", numspectries);
			printf("\tNumber of gene tree rootings = ");
			if(numgenetries == -1) printf("all possible\n");
			else printf("%d\n", numgenetries);
				
			}
		
        for(i=0; i<Total_fund_trees; i++) sourcetree_scores[i] = -1;
        if(criterion == 1 || criterion == 4)
            {
            if(criterion==1)  /* DO MRP */
                {
				BR_file = fopen("coding.nex", "w");
				error = coding(0, 1, 0);
				fclose(BR_file);
				if(error == FALSE)
					{
					if(system("paup coding.nex") != 0) printf("Error calling paup, please execute the file coding.nex in paup to perform the parsimony step\n");
					}
				if(error == FALSE)
					{
					remove("clanntmp.chr");
					remove("clanntree.chr");
					/*pars("coding.nex", "clanntmp.chr"); */
					printf("\n");
					}

                }
            if(criterion == 4)   /* DO AVERAGE CONSENSUS */
                {
					
				if(!got_weights)
					{
					printf("Warning: There were no weights included in the input trees.\nThe analysis will assume that all branch lengths are set to unity (1)\n");
					}
				paupfile = fopen("average_consensus.nex", "w");
				error = average_consensus(0, missing_method,  useroutfile, paupfile);
				fclose(paupfile);
				if(system("paup average_consensus.nex")!= 0) 
					printf("Error calling PAUP*\n\tPlease execute the file average_consensus.nex in PAUP to complete the analysis\n");

                    
                }
            }  /***** FINISH AVERAGE CONSENSUS */
        else
            {
            
            psfile = fopen("supertree.ps", "w");
            
            tree = malloc(TREE_LENGTH*sizeof(char));
            if(!tree) memory_error(53);
            tree[0] = '\0';

            
            if(criterion == 2 || criterion == 3)
                {
                coding(0,0,0);   /** if we are using MRC or QC we need to build the source trees matrix first */
                if(criterion==2) condense_coding();
                }
            
            if(outfile == NULL && print)
                {
                if(outfile == NULL)
                    {
                    if((outfile = fopen("Heuristic_result.txt", "w")) == NULL)
                        {
                        printf("Error opening file named Heuristic_result.txt\n");
                        error = TRUE;
                        }
                    }
                }
                
        
            for(i=0; i<number_retained_supers; i++)
                {
                strcpy(retained_supers[i], "");
                scores_retained_supers[i] = -1;
                }
			

            if(super_scores == NULL)
                {
                super_scores = malloc(number_of_taxa*sizeof(int *));
                if(!super_scores)  memory_error(54);
                    
                for(i=0; i<number_of_taxa; i++)
                    {
                    super_scores[i] = malloc(number_of_taxa*sizeof(int));
                    if(!super_scores[i])  memory_error(55);
                        
                    for(j=0; j<number_of_taxa; j++)
                        {
                        super_scores[i][j] = 0;
                        }
                    }
                }
            else
                {
                for( i=0; i<number_of_taxa; i++)
                    {
                    for(j=i; j<number_of_taxa; j++)
                        {
                        super_scores[i][j] = 0;
                        super_scores[j][i] = 0;
                        }
                    }
                }
			
			
        
            if(start == 0)  /* if we are to use a random starting tree */
                {
                swaps = 0;
                i=0;
                tried_regrafts = 0;
				NUMSWAPS=0;
                interval1 = time(NULL);
                if(print)
                    {
                    printf("\nRandom sampling progress indicator:");
                    fflush(stdout);
                    }
                
                /******* 	HERE WE ARE GOING TO ADD THE ABILITY OF THE HS TO EVALULATE 10,000 RANDOM TREES 
                            TO IDENTIFY THE NECESSARY NUMBER OF STARTING POINTS FOR EACH OF THE REPS OF THE HS.
                            THE NECESSARY NUMBER OF TREES WILL BE THE BEST ONES FOUND FROM THE 10,000 EVALULATED.
                ********/
                starths = malloc(nreps*sizeof(char*));
				startscores = malloc(nreps*sizeof(float));
				if(!startscores) memory_error(82);
                if(!starths) memory_error(82);
                for(j=0; j<nreps; j++)
                    {
                    starths[j] = malloc((number_of_taxa*100)*sizeof(char));
                    startscores[j] = -1;
                    if(!starths[j]) memory_error(83);
                    starths[j][0] = '\0';
                    }
                if(signal(SIGINT, controlc1) == SIG_ERR)
						{
						printf("An error occurred while setting a signal handler\n");
						}
                for(GC=0; GC<sample; GC++)
                    {
					i++;
                    interval2 = time(NULL);
                    if(difftime(interval2, interval1) > 5) /* every 5 seconds print a dot to the screen */
                        {
                        printf("*");
                        fflush(stdout);
                        interval1 = time(NULL);
                        }
                    
                    strcpy(tree, "");
                    random_star_decom(tree);  /* create a random starting tree */
					/*average_consensus(0, missing_method, NULL, NULL);
					neighbor_joining(FALSE, tree, TRUE); */
                    if(tree_top != NULL)
                        {
                        dismantle_tree(tree_top);
                        tree_top = NULL;
                        }
                    temp_top = NULL;
                    tree_build(1, tree, tree_top, TRUE, -1);
                    tree_top = temp_top;
                    temp_top = NULL;
            
                    strcpy(best_tree, "");
                    print_named_tree(tree_top, best_tree);
                    strcat(best_tree, ";");
                    /* evaluate the random tree */
                    if(criterion == 0) distance = compare_trees(FALSE);  /* calculate the distance from the super tree to all the fundamental trees */
                    if(criterion == 2)  /* calculate the distance using the MRC criterion */
                        {
						strcpy(temptree, "");
                        print_tree(tree_top, temptree);
                        strcat(temptree, ";");
                        unroottree(temptree);
                        distance = (float)MRC(temptree);
                        }
                    if(criterion == 3)
                        {
						strcpy(temptree, "");
                        print_tree(tree_top, temptree);
                        strcat(temptree, ";");
                        unroottree(temptree);
                        distance = (float)quartet_compatibility(temptree);
                        }
					if(criterion==5)
						{
						strcpy(temptree, "");
                        print_tree(tree_top, temptree);
                        strcat(temptree, ";");
                        unroottree(temptree);
						distance = get_recon_score(temptree, numspectries, numgenetries);
						}
                    
                    
                    /**** check to see if this is one of the top n trees ***/
					if(distance < startscores[0] || BESTSCORE == -1) BESTSCORE = distance;
                    for(k=0; k<nreps; k++)
                        {
                        if(startscores[k] == -1)
                            {
                            strcpy(starths[k], best_tree);
                            startscores[k] = distance;
                            }
                        else
                            {
                            if(distance <= startscores[k])
                                {
                                for(l=nreps-1; l>k+1; l--)
                                    {
                                    startscores[l] = startscores[l-1];
                                    strcpy(starths[l], starths[l-1]);
                                    }
                                startscores[k] = distance;
                                strcpy(starths[k], best_tree);
                                k = nreps;
                                }
                                
                            }
                        }
                
                    }
                /******* FInished pre-evaulating trees ********/
                swaps+=i;
				NUMSWAPS+= i;
				printf("\nsofar: %d\n", swaps);
				
                if(signal(SIGINT, controlc2) == SIG_ERR)
					{
					printf("An error occurred while setting a signal handler\n");
					}
				for(i=0; i<Total_fund_trees; i++)sourcetree_scores[i] = -1;
				i=0;
				 if(print) printf("\nHeuristic search progress indicator:");
                    fflush(stdout);
                while(i != nreps && !user_break) /* Start searching tree space */
                    {
					if(print)printf("\nRepetition %d of %d:", i+1, nreps);
                    fflush(stdout);
                    swaps += do_search(starths[i], TRUE, print, numswaps, outfile, numspectries, numgenetries);
                    i++;
                    }
                for(i=0; i<nreps; i++)
                    free(starths[i]);
                free(starths);
                free(startscores);

                }
            else   /* if we are to use a nj tree as a starting point */
                {
				if(start == 2)
					{
					average_consensus(0, missing_method, NULL, NULL);
					neighbor_joining(FALSE, temptree, FALSE);
					if(signal(SIGINT, controlc2) == SIG_ERR)
						{
						printf("An error occurred while setting a signal handler\n");
						}
					for(i=0; i<Total_fund_trees; i++)sourcetree_scores[i] = -1;
					 if(print) printf("\nHeuristic search progress indicator:");
					for(i=0; i<nreps; i++)
						{
						if(print) printf("\nRepetition %d of %d:", i+1, nreps);
						if(i != 0)
							{
							random_num =(int)fmod(rand(), 10);
							for(j=0; j<random_num; j++)
								string_SPR(temptree);  /* do a random number of spr permutations to the tree */
							}
						strcpy(tree, temptree);
						returntree(tree);
						swaps += do_search(tree, TRUE, print, numswaps, outfile, numspectries,numgenetries);
						}
					}
				else/* if we are to use a tree (or trees) from a file as the starting tree */
					{
					c = getc(userfile);
					while((c == ' ' || c == '\r'  || c == '\n' || c == '[') && !feof(userfile))  /* skip past all the non tree characters */
						{
						if(c == '[')
							{
							while((c = getc(userfile)) != ']' && !feof(userfile));
							}
						c = getc(userfile);
						}
					k=0; swaps=0;
					while(!feof(userfile))
						{
						if(print)printf("\t\nusertree %d progress indicator: =", k+1);
						i=0;
						while(c != ';')
							{
							if(c != ' ' && c != '\n' && c != '\r')
								{
								tree[i] = c;
								i++;
								}
							tree[i] = ';';
							tree[i+1] = '\0';
							c = getc(userfile);
							}
						i++;

						/* Score the input tree BEFORE doing any swaps !*/
						/****** We now need to build the Supertree in memory *******/
						if(tree_top != NULL)
							{
							dismantle_tree(tree_top);
							tree_top = NULL;
							}
						temp_top = NULL;
						tree_build(1, tree, tree_top, user, -1);
						tree_top = temp_top;
						temp_top = NULL;
						
									
						strcpy(best_tree, "");
						print_named_tree(tree_top, best_tree);
						strcat(best_tree, ";");
						/*  printf("!:assigned tree: %s\n", best_tree); */
								
								
						if(criterion == 0) distance = compare_trees(FALSE);  /* calculate the distance from the super tree to all the fundamental trees */
						if(criterion == 2)  /* calculate the distance using the MRC criterion */
							{
							strcpy(temptree, "");
							print_tree(tree_top, temptree);
							strcat(temptree, ";");
							unroottree(temptree);
							distance = (float)MRC(temptree);
							}
						if(criterion == 3)
							{
							strcpy(temptree, "");
							print_tree(tree_top, temptree);
							strcat(temptree, ";");
							unroottree(temptree);
							distance = (float)quartet_compatibility(temptree);
							}
						if(criterion==5)
							{
							strcpy(temptree, "");
							print_tree(tree_top, temptree);
							strcat(temptree, ";");
							unroottree(temptree);
							distance = get_recon_score(temptree, numspectries, numgenetries);
							}

						if(scores_retained_supers[0] == -1)
							{
							if(BESTSCORE == -1)BESTSCORE = distance;
							retained_supers[0] = realloc(retained_supers[0], (strlen(best_tree)+10)*sizeof(char));
							strcpy(retained_supers[0], best_tree);
							scores_retained_supers[0] = distance;
							number = distance;
							}
						else
							{
							if(distance < number)
								{
								if(distance < BESTSCORE) BESTSCORE = distance;
								j=0;
								while(scores_retained_supers[j] != -1 && j < number_retained_supers)
									{
									strcpy(retained_supers[j], "");
									scores_retained_supers[j] = -1;
									j++;
									}
								number = distance;    /* assign this as the new distance */
								/* we need a string copy of this tree */
								strcpy(best_tree, "");
								print_named_tree(tree_top, best_tree);
								strcat(best_tree, ";");
								retained_supers[0] = realloc(retained_supers[0], (strlen(best_tree)+10)*sizeof(char));
								strcpy(retained_supers[0], best_tree);
								scores_retained_supers[0] = distance;
								
								}
							else
								{
								if(distance == number)
									{

									j = 0;
									strcpy(best_tree, "");
									print_named_tree(tree_top, best_tree);
									strcat(best_tree, ";");
									if(!check_if_diff_tree(best_tree))
										different = FALSE; /* This is here to check if this tree with the same score actually has a different topology from those already stroed */

									if(different)
										{
										while(scores_retained_supers[j] != -1)
											{
											j++;
											if(j+1 == number_retained_supers) reallocate_retained_supers();
											}
										retained_supers[j] = realloc(retained_supers[j], (strlen(best_tree)+10)*sizeof(char));
										strcpy(retained_supers[j], best_tree);
										scores_retained_supers[j] = distance;
										}
									}
								}
							}
					
						for(i=0; i<Total_fund_trees; i++)sourcetree_scores[i] = -1;
						swaps+= do_search(tree, TRUE, print, numswaps, outfile, numspectries, numgenetries);
						strcpy(tree, "");
						c = getc(userfile);
						while((c == ' ' || c == '\r'  || c == '\n' || c == '[') && !feof(userfile))  /* skip past all the non tree characters */
							{
							if(c == '[')
								{
								while((c = getc(userfile)) != ']' && !feof(userfile));
								}
							c = getc(userfile);
							}
						k++;
						}
					}
                }
                if(print)printf("\n");
                if(print)printf("Number of topologies tried: %d\n",swaps);
                
                
                /**** Print out the best trees found *******/
                
                tree[0] = '\0';
                i=0; j=0;
				if(do_histogram && criterion == 0) histogram_file = fopen(histogramfile_name, "w");
                if(print)
                    {
                    while(scores_retained_supers[j] != -1)
                        {
                        j++;
                        }
                    trees_in_memory = j;
					
                    while(scores_retained_supers[i] != -1)
                        {
                        if(outfile != NULL) fprintf(outfile, "%s\t[%f]\n", retained_supers[i], scores_retained_supers[i] );
                        tree_coordinates(retained_supers[i], FALSE, TRUE, FALSE, -1);
                        printf("\nSupertree %d of %d score = %f\n", i+1, j, scores_retained_supers[i] );
						


						if(do_histogram && criterion == 0)  /* If we are going to draw a histogram of the scores of the source trees when compared to the best supertree */
							{		/* right now we can only do this for mssa (dfit) */
							if(tree_top != NULL)
								{
								dismantle_tree(tree_top);
								tree_top = NULL;
								}
							temp_top = NULL;
							
							tree_build(1, retained_supers[i], tree_top, TRUE, -1);
							tree_top = temp_top;
							temp_top = NULL;
							
							/**** evaluate its fit to the source trees in memory *****/
							temptree[0] = '\0';
							for(j=0; j<Total_fund_trees; j++) sourcetree_scores[j] = -1;
							compare_trees(FALSE);  /* calculate the distance from the super tree to all the fundamental trees */
							printf("\nPlot of source tree scores against supertree number %d", i+1);
							draw_histogram(histogram_file, bins, sourcetree_scores, Total_fund_trees);
							
							
							}
							
                        i++;
                        }

                    }
				
                if(do_histogram && criterion == 0) fclose(histogram_file);
                if(outfile != NULL) fclose(outfile);
                if(userfile != NULL) fclose(userfile);
                fclose(psfile);
                free(tree);
             }
        }
	
	
    free(temptree);
    free(best_tree);
    }

int average_consensus(int nrep, int missing_method, char * useroutfile, FILE *paupfile)
	{
	int **taxa_comp = NULL, i, j, k, l, found = FALSE, here = FALSE, error = FALSE;
	float used_weights = 0;
	char *temptree = NULL;
	
	
	temptree = malloc(TREE_LENGTH*sizeof(char));
	if(!temptree) printf("out of memory'n");
	temptree[0] = '\0';
	taxa_comp = malloc(number_of_taxa*sizeof(int*));
	if(!taxa_comp) printf("out of memory'n");
	for(i=0; i<number_of_taxa; i++)
		{
		taxa_comp[i] = malloc(number_of_taxa*sizeof(int));
		if(!taxa_comp[i]) printf("out of memory'n");
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
		printf("\n\nERROR: the overlap in the data is too sparse to calculate missing cells using ");
		if(missing_method == 0) printf("an ultrametric estimate\n");
		if(missing_method == 1) printf("a 4 point condition estimate\n");
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



int do_search(char *tree, int user, int print, int maxswaps, FILE *outfile, int numspectries, int numgenetries)
    {
    int swaps = 0, i=0, better_score = TRUE;
	char temporary_tree[TREE_LENGTH];
	
	
	temporary_tree[0] = '\0';
    unroottree(tree);
        /****** We now need to build the Supertree in memory *******/
        if(tree_top != NULL)
            {
            dismantle_tree(tree_top);
            tree_top = NULL;
            }
        temp_top = NULL;
        tree_build(1, tree, tree_top, user, -1);
        tree_top = temp_top;
        temp_top = NULL;
/*		print_tree(tree_top, temporary_tree);
		printf("built tree = %s\n", temporary_tree);
  */      if(method == 1)
            {
            /* if NNI is to be carried out */       
            swaps = branchswap(maxswaps,-1, numspectries, numgenetries);
            }
        if(method == 2 || method == 3)
            {
            sprscore = -1;
            tried_regrafts = 0;
            /* if SPR is to be carried out */
/*            while(better_score == TRUE && tried_regrafts < maxswaps && !user_break)
                {
                do
                    {
                    branchpointer = NULL;
                    better_score = spr(tree_top, maxswaps, numspectries, numgenetries);
                    }while(better_score == FALSE && remaining_spr(tree_top) > 0 && tried_regrafts < maxswaps && !user_break);
					printf("better_score = %d\t, remaining_spr = %d\ttried_regrafts = %d\tuserbreak = %d\n", better_score, remaining_spr(tree_top), tried_regrafts, user_break);
                reset_spr(tree_top);
                }
            swaps = tried_regrafts;
 */
 

			do/* for a single rep */
				{
				branchpointer = NULL;
				better_score = spr_new(tree_top, maxswaps, numspectries, numgenetries);
				}while(better_score == TRUE && tried_regrafts < maxswaps && !user_break);
				
            swaps = tried_regrafts;
		   

            
            }
    return(swaps);
    }

int remaining_spr(struct taxon *position)
    {
    int count = 0;
    
    while(position != NULL)
        {
        if(position->daughter != NULL)
            {
            if(!position->spr)
                count++;
            else
                count += remaining_spr(position->daughter);
            }
        position = position->next_sibling;
        }
    return(count);
    }





    

/* The NNI branch swapping part of this program will be in several parts:
	part 1: This will call the recursive function to search for the first better score found. It will count how many swaps it makes, and stops the 
		process if it goes through the rescursion without finding a better score, or after a certain number of branch swaps.
		
	part 2: This is the recursive function that goes to the bottom of the tree and picks the structures in memory to swap. If a better score was found this 
		recursive function stops and returns that a better score was found. Otherwise it passes that a worse score was found.
	
	part 3: This actually does the swapping and then calls the scoring algorithm to see what this new topology gets. If the score is less than the last
		score, then keep the new topology, otherwise swap back. This should return whether or not a better score was found.
	
*/				
int branchswap(int number_of_swaps, float score, int numspectries, int numgenetries)
	{  /* This is the function that carries out part 1 of the algorithm */
					
        float number = score, last_number = 0;
        int  i=0, j=0, last = 0;
        
        
        while(i < number_of_swaps && !user_break)
            {
            last_number = number;
            i+= find_swaps(&number, tree_top, number_of_swaps, numspectries, numgenetries);
            last = i;
            if(last_number == number && branchpointer != NULL)
                    {
                    branchpointer = NULL; 
                    }
            else
                {
                
                if(last_number == number && branchpointer == NULL)
                    {
                    last = i;
                    i = number_of_swaps;
                    }
                }
        
            }
			
	return(last);
	}				

				
/* part 2 */				
int find_swaps(float * number, struct taxon * position, int number_of_swaps, int numspectries, int numgenetries)
	{
	struct taxon *start1 = position, *start2 = NULL, *first_swap = NULL, *second_swap = NULL, *tmp = NULL;
	int swaps = 0, better_score = FALSE, *siblings = NULL, i=0, j=0, k=0, count =0;
	float distance = *number;


        /** count how many pointer siblings there are on this level */
        i=0;
        while(position != NULL && !user_break)
            {
            if(position->daughter != NULL) i++;
            position = position->next_sibling;
            }

        siblings = malloc(i*sizeof(int));
        for(j=0; j<i; j++) siblings[j] = FALSE;
        
        position = start1;
	/* Search for pointer siblings and if found, call the function recursively for that pointer sibling */
	while(count < i && swaps < number_of_swaps && !better_score && !user_break )
		{
                /* choose the pointer sibling to travel down */
                j =(int)fmod(rand(), i);  /* the jth pointer sibling is the chosen one */
                while(siblings[j])
                    j =(int)fmod(rand(), i);  /* the jth pointer sibling is the chosen one */

                /* go to the chosen sibling pointer */
                position = start1;
                k=-1;
                while(k != j)
                    {
                    if(position->daughter != NULL) k++;
                    if(k != j)
                        position = position->next_sibling;
                    }
                swaps += find_swaps(number, position->daughter, number_of_swaps, numspectries, numgenetries);
                if(distance != *number) better_score = TRUE;
                count++;
                siblings[j] = TRUE;
		 }		
				
	if(swaps < number_of_swaps && !better_score && !user_break)  /* if we haven't exceeded the total number of swaps allowed */
		{
				
		/* Use the all siblings at this level to swap against the siblings at the chosen level (except the pointer sibling pointing to this level) */	
		/* the pointer first_swap will point to the sibling at this level we want to swap */
		first_swap = start1;
	
		/* The pointer second_swap will point to the sibling at the previous level we want to swap */
		/* This means we have to go up one level and rewind to the start of the list */
	
		if(start1->parent != NULL && !user_break)  /* As long as we are not at the top of the tree */
			{
			
			/** call swapper  **/
                        
			better_score = swapper(start1->parent, start1, 1, first_swap, second_swap, number, &swaps, number_of_swaps, numspectries, numgenetries);
                        
			
			 }
		 
		 }
        free(siblings);
	return(swaps);
	
	}
		 		
		 		
/* This function actually does the swapping */	
/* Part 3 */	
void do_swap(struct taxon * first, struct taxon * second)
	{
	struct taxon *next = NULL, *prev = NULL, *parent = NULL, *tmp = NULL;
	
        interval2 = time(NULL);
        if(difftime(interval2, interval1) > 5) /* every 5 seconds print a dot to the screen */
            {
            printf("=");
            fflush(stdout);
            interval1 = time(NULL);
            }
            
	/* Change everything pointing to the two being swapped */
	if(first->parent != NULL)
		(first->parent)->daughter = second;	
		
	if(second->parent != NULL)
		(second->parent)->daughter = first;
		
	if(first->prev_sibling != NULL)
		(first->prev_sibling)->next_sibling = second;
		
	if(second->prev_sibling != NULL)
		(second->prev_sibling)->next_sibling = first;
		
	if(first->next_sibling != NULL)
		(first->next_sibling)->prev_sibling = second;
		
	if(second->next_sibling != NULL)
		(second->next_sibling)->prev_sibling = first;
		
	
	
	
	/* Change what the two are pointing to */
	next = first->next_sibling;
	prev = first->prev_sibling;
	parent = first->parent;
	
	first->next_sibling = second->next_sibling;
	first->prev_sibling = second->prev_sibling;
	first->parent = second->parent;
	
	second->next_sibling = next;
	second->prev_sibling = prev;
	second->parent = parent;
	
	if(tree_top == first)
		{
		 tree_top = second;
		 }
	else
		{
		if(tree_top == second)
			{
			tree_top = first;
			}
		}
	
	}
	


/* this is the function which will allow the swapping to take place x number of steps away, this will make the branch swapping that much more robust. The function from its start point will travel UP from its present position and then will travel up or down the tree (but never back to its start point) looking for viable swaps x number of steps away */

int swapper(struct taxon * position,struct taxon * prev_pos, int stepstaken, struct taxon * first_swap, struct taxon * second_swap, float * number, int * swaps, int number_of_swaps, int numspectries, int numgenetries)
	{
	
  	float  distance = *number;
	int better_score = FALSE, i=0, j=0, different = TRUE;
	struct taxon *start2 = NULL, *start1 = NULL;
        char *best_tree = NULL, *temptree = NULL;
        
        temptree = malloc(TREE_LENGTH*sizeof(char));
        if(!temptree) memory_error(66);
        temptree[0] = '\0';
        best_tree = malloc(TREE_LENGTH*sizeof(char));
        if(!best_tree) memory_error(52);
        best_tree[0] = '\0';
	
	while(position->prev_sibling != NULL) position = position->prev_sibling;  /* rewind to start of this level */
	second_swap = position;

	/* travel up one if it exists */
	if(stepstaken < number_of_steps && position->parent != NULL && position->parent != prev_pos && !user_break) better_score = swapper(position->parent, position, stepstaken+1, first_swap, second_swap, number, swaps, number_of_swaps, numspectries, numgenetries);

	
	/*   if we aren't the required number steps away, then call another instance of the program    */
	while(stepstaken < number_of_steps && position != NULL && !better_score && !user_break)
		{
		if(position->daughter != NULL && position->daughter != prev_pos)
			{
			better_score = swapper(position->daughter, position, stepstaken+1, first_swap, second_swap, number, swaps, number_of_swaps, numspectries, numgenetries);
			}
		position = position->next_sibling;
		}
	if(!better_score && stepstaken <= number_of_steps && !user_break)
		{
		/* do the swapping thing and check if its better, if it's not then swap it back, otherwise keep it and assign the 
		new score to score, this will act as a signal to the recursive function to stop */

		start1 = first_swap;
			
		if(second_swap == NULL) printf("second_swap is null!\n");
		while(second_swap->prev_sibling != NULL) second_swap = second_swap->prev_sibling;  /* rewinding */
		
		start2 = second_swap; /* start2 now points to the start levle to be swapped against */
	
                /* now first_swap is at the first sibling at this level, and second_swap is at the first sibling at the other level */
                /* The next step is to swap every sibling at this level with every sibling at the level above */
		while(first_swap != NULL && *swaps < number_of_swaps && !better_score && !user_break)
			{
			while(second_swap != NULL && *swaps < number_of_swaps && !better_score && !user_break)
				{
				if(second_swap->daughter != prev_pos )   /* As long as the sibling at the previous level is not the one that defines the clade we are in */
					{
					do_swap(first_swap, second_swap);

                                        if(check_taxa(tree_top) == number_of_taxa)
                                            {
											(*swaps)++;
                                           /* *swaps = *swaps + 1; */ /* count the number of swaps done */
											NUMSWAPS++;
                                            strcpy(best_tree, "");
                                            print_named_tree(tree_top, best_tree);
                                            strcat(best_tree, ";");
                                            /*  printf("!:assigned tree: %s\n", best_tree);  */
                                            
                                            
                                            if(criterion == 0) distance = compare_trees(FALSE);  /* calculate the distance from the super tree to all the fundamental trees */
                                            if(criterion == 2)  /* calculate the distance using the MRC criterion */
                                                {
                                                print_tree(tree_top, temptree);
                                                strcat(temptree, ";");
                                                unroottree(temptree);
                                                distance = (float)MRC(temptree);
                                                }
                                            if(criterion == 3)
                                                {
                                                print_tree(tree_top, temptree);
                                                strcat(temptree, ";");
                                                unroottree(temptree);
                                                distance = (float)quartet_compatibility(temptree);
                                                }

											if(criterion==5)
												{
												strcpy(temptree, "");
												print_tree(tree_top, temptree);
												strcat(temptree, ";");
												unroottree(temptree);
												distance = get_recon_score(temptree, numspectries, numgenetries);
												}


                                            if(scores_retained_supers[0] == -1 || *number == -1)
                                                {
												if(BESTSCORE == -1)BESTSCORE = distance;
                                                if(scores_retained_supers[0] == -1)
                                                    {
													retained_supers[0] = realloc(retained_supers[0], (strlen(best_tree)+10)*sizeof(char));
                                                    strcpy(retained_supers[0], best_tree);
                                                    scores_retained_supers[0] = distance;
                                                    }
                                                *number = distance;
                                                }
                                            else
                                                {
                                                if(distance < scores_retained_supers[0] || distance < *number )
                                                    {
													if(distance < BESTSCORE)BESTSCORE = distance;
                                                    if(distance < *number)
                                                        {
                                                        better_score = TRUE;
                                                        *number = distance;
                                                        }
                                                    if(distance < scores_retained_supers[0] )
                                                        {
                                                        /* we need a string copy of this tree */
                                                        strcpy(best_tree, "");
                                                        print_named_tree(tree_top, best_tree);
                                                        strcat(best_tree, ";");
                                                        /*printf("assigned tree: %s\n", best_tree); */
                                                        retained_supers[0] = realloc(retained_supers[0], (strlen(best_tree)+10)*sizeof(char));
														strcpy(retained_supers[0], best_tree);
                                                        scores_retained_supers[0] = distance;
                                                        j=1;
                                                        while(scores_retained_supers[j] != -1 && j < number_retained_supers)
                                                            {
                                                            strcpy(retained_supers[j], "");
                                                            scores_retained_supers[j] = -1;
                                                            j++;
                                                            }
                                                        
                                                        branchpointer = second_swap;
                                                        }
                                                    }
                                                else
                                                    {
                                                    if(distance == scores_retained_supers[0])
                                                        {
                                                        j = 0;
                                                        strcpy(best_tree, "");
                                                        print_named_tree(tree_top, best_tree);
                                                        strcat(best_tree, ";");
                                                        different = TRUE;
                                                        if(!check_if_diff_tree(best_tree))
                                                            different = FALSE; /* This is here to check if this tree with the same score actually has a different topology from those already stroed */

                                                        if(different)
                                                            {
                                                            while(scores_retained_supers[j] != -1)
                                                                {
                                                                j++;
                                                                if(j+1 == number_retained_supers) reallocate_retained_supers();
                                                                }
															retained_supers[j] = realloc(retained_supers[j], (strlen(best_tree)+10)*sizeof(char));
                                                            strcpy(retained_supers[j], best_tree);
                                                            scores_retained_supers[j] = distance;
                                                            }
                                                        }
                                                    
                                                    /*** swap them back ****/
                                                    if(!better_score) do_swap(first_swap, second_swap);
                                                    }
                                                }
                                            }
                                        else
                                            {
                                            do_swap(first_swap, second_swap);
                                            }
    
					}	 		
			 	second_swap = second_swap->next_sibling;
		 		
			 	}
		 	
			 first_swap = first_swap->next_sibling;
			 second_swap = start2;
		 	
			 }

		}
        free(best_tree);
        free(temptree);
	return(better_score);

	}



void yaptp_search(void)
    {
    int i=0, j=0, k=0, l=0, random_num = 0, error = FALSE, yaptp_method = 1;
    int Nreps = 100, search = 1;
    char filename[1000], best_tree[TREE_LENGTH];
    FILE *yaptpfile = NULL;

    filename[0] = '\0';
    strcpy(filename, "yaptp.txt");
    if(criterion == 5) yaptp_method = 2;    

    for(i=0; i<num_commands; i++)
        {
        if(strcmp(parsed_command[i], "nreps") == 0)
            {
            Nreps = toint(parsed_command[i+1]);
            if(Nreps == 0)
                {
                printf("Error: '%s' is an invalid number of repetitions\n", parsed_command[i+1]);
                error = TRUE;
                }
            }
        if(strcmp(parsed_command[i], "method") == 0)
            {
            if(strcmp(parsed_command[i+1], "equiprobable") == 0)
		{
		if (criterion != 5)
        		yaptp_method =1;
		else
			{
			printf("Using the citerion recon it is only possible to use the markovian method for randomisation\n");
			yaptp_method = 2;
			}
		}
            else
                {
                if(strcmp(parsed_command[i+1], "markovian") == 0)
                    yaptp_method =2;
                else
                    {
                    printf("Error: method option '%s' unknown\n", parsed_command[i+1]);
                        error = TRUE;
                    }
                }
            }
            
        if(strcmp(parsed_command[i], "search") == 0)
            {
            if(strcmp(parsed_command[i+1], "all") == 0) search = 0;
            else
                {
                if(strcmp(parsed_command[i+1], "nni") == 0) search = 1;
                else
                    {
                    if(strcmp(parsed_command[i+1], "spr") == 0) search = 1;   /* should be set to 2 when SPR is implemented */
                    else
                        {
                        printf("Error: search option '%s' unknown\n", parsed_command[i+1]);
                        error = TRUE;
                        }
                    }
                }
            }
        if(strcmp(parsed_command[i], "treefile") == 0)
            {
            strcpy(filename, parsed_command[i+1]);
            if(filename[0] == '\0')
                {
                printf("Error: Unable to open file named '%s'\n", parsed_command[i+1]);
                error = TRUE;
                strcpy(filename, "yaptp.txt");
                }
            }
        }
        
    if(!error)
        {
		
        printf("\n\nYAPTP Settings:\n\tNumber of YAPTP repetitions: %d\n\tSearching Supertree-space: ", Nreps);
        if(search == 0)printf("Exhaustive Search of Supertree-spce\n");
        if(search == 1)printf("Heuristic Search\n");
        printf("\tRandomisation method (method) = ");
        if(yaptp_method==1) printf("equiprobable\n");
        if(yaptp_method==2) printf("markovian\n");
        printf("\tYAPTP output file: yaptp.txt\n\n");
        
		if(!calculated_fund_scores && criterion == 0)
			{
			cal_fund_scores(FALSE);    /* calculate the path metrics for all the fundamental trees */
			calculated_fund_scores = TRUE;
			}
		
		
		if(yaptp_results != NULL) free(yaptp_results);
        yaptp_results = malloc((Nreps+1)*sizeof(float));
		for(i=0; i<Nreps+1; i++)
			yaptp_results[i] = 0;
		
		if(trees_in_memory > 0)
			{
			yaptp_results[0] = scores_retained_supers[0];
			strcpy(saved_supertree, retained_supers[0]);
			}
        else
			{
			heuristic_search(FALSE, TRUE, 1000, 1);
			yaptp_results[0] = scores_retained_supers[0];
			strcpy(saved_supertree, retained_supers[0]);
			}
		
        if(criterion == 1)  /* If we are using MRP and PAUP **/
            {
            BR_file = fopen("coding.nex", "w");
            coding(0, 2, Nreps);
            fclose(BR_file);
            if(system("paup coding.nex") != 0) printf("Error calling paup, please execute the file coding.nex in paup to perform the parsimony step\n");
            }
        else
            {
    
        
            yaptpfile = fopen(filename, "w");
            /********* initalise the array to store the arrays *********/
            
            stored_fund_scores = malloc(Total_fund_trees*sizeof(int**));
            if(stored_fund_scores == NULL) memory_error(38);
                
            for(i=0; i<Total_fund_trees; i++)
                {
                stored_fund_scores[i] = malloc((number_of_taxa)*sizeof(int*));
                if(stored_fund_scores[i] == NULL) memory_error(39);
                    
                else
                    {
                    for(j=0; j<(number_of_taxa); j++)
                        {
                        stored_fund_scores[i][j] = malloc((number_of_taxa)*sizeof(int));
                        if(stored_fund_scores[i][j] == NULL) memory_error(40);
                            
                        else
                            {
                            for(k=0; k<(number_of_taxa); k++)
                                stored_fund_scores[i][j][k] = fund_scores[i][j][k];  /* copy the origninal fundamental scores for safe keeping */
                            }
                        }
                    }
                }
            stored_presence_of_taxa = malloc(Total_fund_trees*sizeof(int *));
            if(!stored_presence_of_taxa) memory_error(41);
                
            for(i=0; i<Total_fund_trees; i++)
                {
                stored_presence_of_taxa[i] = malloc((number_of_taxa)*sizeof(int));
                if(!stored_presence_of_taxa[i]) memory_error(42);
                    
                for(j=0; j<number_of_taxa; j++)
                    {
                    stored_presence_of_taxa[i][j] = presence_of_taxa[i][j];
                    }
                }
        
        
            stored_funds = malloc(Total_fund_trees*sizeof(char *));
            if(!stored_funds) memory_error(31);
        
            for(i=0; i<Total_fund_trees; i++)
                {
                stored_funds[i] = malloc((tree_length_assignments*TREE_LENGTH)*sizeof(char));
                if(!stored_funds) memory_error(32);
        
                stored_funds[i][0] = '\0';
                strcpy(stored_funds[i], fundamentals[i]);   /* copy the original fundamentals for safe keeping */
                }
            
            
        
            /*****************************************/
            hsprint=TRUE;
            printf("YAPTP progress indicator: "); 
			GC = 0;
            for(i=0; i<Nreps; i++)  /* for all reps of the yaptp algorithm */
                {
				if(!user_break)
					{
					if(i>0)hsprint=FALSE;
					printf("\n\trepetition %d:=", i+1);
					fflush(stdout);
					/***** Create the replicate of randomised fundamental trees  **********/
					
					for(j=0; j<Total_fund_trees; j++)
						{
						if(yaptp_method == 1 && criterion != 5) randomise_tree(fundamentals[j]); /*equiprobable */
						if(yaptp_method == 2 || criterion == 5) randomise_taxa(fundamentals[j]);  /*Markovian*/
						}

					fflush(stdout);
					cal_fund_scores(FALSE);
					/* Now find the best tree for this randomised set of fundamental trees */
					
				
					if(search == 0)
						{
						alltrees_search(FALSE);
						}
					else
						{

						heuristic_search(FALSE, FALSE, 100, 1);

						}

					best_tree[0] = '\0';
					k=0;
					yaptp_results[i+1] = scores_retained_supers[0];
					while(scores_retained_supers[k] != -1)
						{
						if(search ==0)
							{
							if(tree_top != NULL)
								{
								dismantle_tree(tree_top);
								tree_top = NULL;
								}
							temp_top = NULL;
							tree_build(1, retained_supers[k], tree_top, FALSE, -1);
							tree_top = temp_top;
							temp_top = NULL;
						
							strcpy(best_tree, "");
							print_named_tree(tree_top, best_tree);
							/*printf("%s;\t[%f]\n", best_tree, scores_retained_supers[k]);  */
							fprintf(yaptpfile, "%s\t[%f %d]\n", best_tree, scores_retained_supers[k], i);
							}
						else
							{
						/* printf("%s;\t[%f]\n", retained_supers[k], scores_retained_supers[k]);  */
							fprintf(yaptpfile, "%s\t[%f %d]\n", retained_supers[k], scores_retained_supers[k], i);
							}
						fflush(yaptpfile);
			
						k++;
						}
					while(k>=0)
						{
						strcpy(retained_supers[k], "");
						scores_retained_supers[k] = -1;
						k--;
						}
					GC++;
					}
                }
            
            
            
            /* copy back the arrays to their original places */
            
            for(i=0; i<Total_fund_trees; i++)
                {
                strcpy(fundamentals[i], stored_funds[i]);
                for(j=0; j<number_of_taxa; j++)
                    {
                    presence_of_taxa[i][j] = stored_presence_of_taxa[i][j];
                    for(k=0; k<number_of_taxa; k++)
                        {
                        fund_scores[i][j][k] = stored_fund_scores[i][j][k];
                        }
                    }
                }
            
			
			/* Draw the histogram as a result */
			
			 draw_histogram(yaptpfile, 10, yaptp_results, GC+1);
			
			
			fclose(yaptpfile);
            /* free all the memory allocated */
            
            for(i=0; i<Total_fund_trees; i++)
                {
                free(stored_funds[i]);
                free(stored_presence_of_taxa[i]);
                for(j=0; j<number_of_taxa; j++)
                    {
                        free(stored_fund_scores[i][j]);
                    }
                free(stored_fund_scores[i]);
                }
            free(stored_funds);
            free(stored_presence_of_taxa);
            free(stored_fund_scores);
            }
        }
    hsprint=TRUE;
    }


void randomise_tree(char *tree)
    {
    int i=0, j=0, k=0, l=0, x=0, y=0, treecount = 0, random=0, supers = 0, actual_num = 0;
    char **array = NULL, temptree[TREE_LENGTH], *newtree = NULL, *tmp;
    /** allocate the array **/
    array = malloc(number_of_taxa*sizeof(char *));
    if(!array) memory_error(56);
    
    for(i=0; i<number_of_taxa; i++)
        {
        array[i] = malloc(100*sizeof(char));
        if(!array[i]) memory_error(57);
        array[i][0] = '\0';
        }
       
    temptree[0] = '\0';
    
    tmp = malloc(10*sizeof(char));
    newtree = malloc(TREE_LENGTH*sizeof(char));
    
    /** run through the tree recording the names of the taxa present into a dynamically allocated array */
    i=0;
    while(tree[i] != ';')
        {
		if(tree[i] == ':')
			{
			while(tree[i] != '(' && tree[i] != ')' && tree[i] != ',' && tree[i] != ';')
				i++;
			}
		if(tree[i] == ')')
			{
			i++;
			while(tree[i] != '(' && tree[i] != ')' && tree[i] != ',' && tree[i] != ';' && tree[i] != ';')
				i++;
			}

		else
			{
			if(tree[i] != '(' && tree[i] != ')' && tree[i] != ',' && tree[i] != ';' && tree[i] != ':')
				{
				k=0;
				while(tree[i] != '(' && tree[i] != ')' && tree[i] != ',' && tree[i] != ';' && tree[i] != ':')
					{
					array[j][k] = tree[i];
					i++;
					k++;
					}
				array[j][k] = '\0';
				j++;
				}
			else
				i++;
			}
        }
    /* now all the names are stored and j is the number of taxa on the tree */
    
    /* now randomly pick a tree with j taxa and put the taxa back on this tree */
    
    /* figure out how many trees there are for j number of taxa */
    actual_num = number_of_taxa;
    number_of_taxa = j;
    supers = 1;
    for(i=4; i<=j; i++) supers*=((2*i)-5);
    /* randomly pick a number between 1 and the number of taxa */
     random = (int)fmod(rand(), supers)+1;
    for(i=0; i<TREE_LENGTH; i++) newtree[i] = '\0';
    /* now build the tree */
    intTotree(random, newtree, number_of_taxa);
	
    /* now put the names back on the tree in place of the numbers returned by intTotree */
    i=0; j=0; l=0;
    temptree[0] = '\0';
    while(newtree[i] != ';')
        {
        if(newtree[i] != '(' && newtree[i] != ')' && newtree[i] != ',' && newtree[i] != ';')
            {
            x=0;
            while(newtree[i] != '(' && newtree[i] != ')' && newtree[i] != ',' && newtree[i] != ';')
                {
                tmp[x] = newtree[i];
                i++; x++;
                }
            tmp[x] = '\0';
            strcat(temptree, array[toint(tmp)]);    
            }
        else
            {
            if(newtree[i] == '(')
                strcat(temptree, "(");
            if(newtree[i] == ')')
                strcat(temptree, ")");
            if(newtree[i] == ',')
                strcat(temptree, ",");
            if(newtree[i] == ';')
                strcat(temptree, ";");

            i++;
            }
        }
    strcat(temptree, ";");
    number_of_taxa = actual_num;
    strcpy(tree, temptree);

    for(i=0; i<number_of_taxa; i++)
        {
        free(array[i]);
        array[i] = NULL;
        }
    free(array);
    array = NULL;
    free(newtree);
    newtree = NULL;
    free(tmp);
    }

void randomise_taxa(char *tree)
    {
    int i=0, j=0, k=0, l=0, x=0, y=0, treecount = 0, random=0, tottax;
    char **array = NULL, temptree[TREE_LENGTH];
    

    /* Start by counting the number of taxa in the tree (there may be more then the variable "number_of_taxa" because of the recon criterion */

    i=0; tottax=1;
    while(tree[i] != ';')
	{
	if(tree[i] == ',')tottax++; /* tottax is the num,ber of taxa */
	i++;
	}

    /** allocate the array **/
    array = malloc(tottax*sizeof(char *));
    if(!array) memory_error(56);
    
    for(i=0; i<tottax; i++)
        {
        array[i] = malloc(100*sizeof(char));
        if(!array[i]) memory_error(57);
        array[i][0] = '\0';
        }
       
    temptree[0] = '\0';
    
    /** run through the tree recording the names of the taxa present into a dynamically allocated array */
    
    i=0;
    while(tree[i] != ';')
        {
		if(tree[i] == ':')
			{
			while(tree[i] != '(' && tree[i] != ')' && tree[i] != ',' && tree[i] != ';')
				i++;
			}
        if(tree[i] != '(' && tree[i] != ')' && tree[i] != ',' && tree[i] != ';')
            {
            k=0;
            while(tree[i] != '(' && tree[i] != ')' && tree[i] != ',' && tree[i] != ';' && tree[i] != ':')
                {
                array[j][k] = tree[i];
                i++;
                k++;
                }
            array[j][k] = '\0';
            j++;
            }
        else
            i++;
        }
        
    /* now all the names are stored and j is the number of taxa on the tree */
    
    /* randomly replace the names of the taxa onto the tree */
    i=0; k=0;
    while(tree[i] != ';')
        {
		if(tree[i] == ':')
			{
			while(tree[i] != '(' && tree[i] != ')' && tree[i] != ',' && tree[i] != ';')
				i++;
			}
        if(tree[i] != '(' && tree[i] != ')' && tree[i] != ',' && tree[i] != ';' && tree[i] != ':')  /* if there is a name here */
            {
            x=0; y=-1;
            treecount = -1;
            /* randomly pick a taxa from those remaining */
            l =(int)fmod(rand(), j);  /* the l^th tree in the array is the chosen one */
            
            /** find the l^th tree in the array  **/
            do
                {
                y++;
                if(array[y][0] != '\0')
                    treecount++;
                }while(treecount != l);
                
            /** put this taxa into the tree **/
            while(array[y][x] != '\0')
                {
                temptree[k] = array[y][x];
                k++;
                x++;
                }
                
            strcpy(array[y], "");  /* we've used this taxa, so delete it */
            
            /** skip forward in the real tree to the end of the taxa name **/
            while(tree[i] != '(' && tree[i] != ')' && tree[i] != ',' && tree[i] != ';' && tree[i] != ':')
                i++;
            j--;  /* one less tree from those remaining */
                    
            }
        else
            {
            temptree[k] = tree[i];
            i++;
            k++;

            }
        }
        temptree[k] = ';';
        temptree[k+1] = '\0';
        strcpy(tree, temptree);

    for(i=0; i<tottax; i++)
        {
        free(array[i]);
        array[i] = NULL;
        }
    free(array);
    array = NULL;
    }


void random_star_decom(char *tree)
    {
    char **treearray = NULL, *tempchar = NULL;
    int i=0, j=0, k=0, l=0, nodes = number_of_taxa, random1 = 0, random2 = 0;
    
    tempchar = malloc((100*number_of_taxa)*sizeof(char));
    tempchar[0] = '\0';
    treearray = malloc(number_of_taxa*sizeof(char *));
    if(!treearray) memory_error(58);
    for(i=0; i<number_of_taxa; i++)
        {
		treearray[i] = NULL;
        treearray[i] = malloc(TREE_LENGTH*sizeof(char));
        if(!treearray[i]) memory_error(59);
        treearray[i][0] = '\0';
        }
    /** copy the names of the taxa into the array **/
    
    for(i=0; i<number_of_taxa; i++)
        {
        strcpy(treearray[i], taxa_names[i]);
        }
    /** at each round, combine two of the nodes into one node **/
    while(nodes > 1)
        {
        /** randomly pick the two nodes to be combined **/
        random1 = (int)fmod(rand(), nodes);
        random2 = (int)fmod(rand(), nodes);
        while(random2 == random1)
            {
            random2 = (int)fmod(rand(), nodes);
            }
        /** these numbers represent the ith numbers in turn left in the array **/
        i=-1;k=0;
        do
            {
            if(strcmp(treearray[k], "") != 0 )
                {
                i++;
                }
            k++;
            }while(i != random1 && k < number_of_taxa);
        k--;
        j=-1;l=0;
        do
            {
            if(strcmp(treearray[l], "") != 0 )
                {
                j++;
                }
            l++;
            }while(j != random2 && k < number_of_taxa);
        l--;
        /** combine the contents of the two nodes **/
        strcpy(tempchar, "(");
        strcat(tempchar, treearray[k]);
        strcat(tempchar, ",");
        strcat(tempchar, treearray[l]);
        strcat(tempchar, ")");
        strcpy(treearray[k], "");
		treearray[k] = realloc(treearray[k], (10)*sizeof(char));
		treearray[l] = realloc(treearray[l], (strlen(tempchar)+10)*sizeof(char));
		strcpy(treearray[l], tempchar);
        
        nodes--;
        }
    /** find the position of the completed tree and copy it into the string passed to the function **/
    i=0;
    while(strcmp(treearray[i], "") == 0 && i < number_of_taxa)
        i++;
    strcpy(tree, treearray[i]);
    strcat(tree, ";");
    
    /** free up allocated memory **/
    for(i=0; i<number_of_taxa; i++)
        {
        free(treearray[i]);
        treearray[i] = NULL;
        }
    free(treearray);
    treearray = NULL;
    free(tempchar);
    tempchar = NULL;
    
    }
    
 
/* This function check a given tree against those stored in retained_tree to make sure that they're not the same **/
/* It does this by checking the differences in path metrics of th etwo trees. It they're the same the path metrics will be the same **/
int check_if_diff_tree(char *tree)
    {
    int i=0, j=0, k=0, l=0, intname = -1, different = TRUE;
    int **scores1 = NULL, **scores2 = NULL, temp = 0;
    char tree1[TREE_LENGTH], tree2[TREE_LENGTH], *name = NULL;
    
    scores1 = malloc(number_of_taxa*sizeof(int *));
    if(!scores1) memory_error(62);
    for(i=0; i<number_of_taxa; i++)
        {
        scores1[i] = malloc(number_of_taxa*sizeof(int));
        if(!scores1[i]) memory_error(63);
        }

    scores2 = malloc(number_of_taxa*sizeof(int *));
    if(!scores2) memory_error(62);
    for(i=0; i<number_of_taxa; i++)
        {
        scores2[i] = malloc(number_of_taxa*sizeof(int));
        if(!scores2[i]) memory_error(63);
        }
    
    name = malloc(1000*sizeof(int));
    if(!name) memory_error(61);
    /** since both the trees being checked will be created using a heuristic search, they will contain the actual taxa names, it is necessary to get their taxan number Equivalent */
    i=0;
    /* scan along the tree and change the names of the taxa to their numbers */
    while(tree[i] != ';')
        {
        if(tree[i] == '(' || tree[i] == ')' || tree[i] == ',' || tree[i] == ';')
            {
            tree1[j] = tree[i];
            i++; j++;
            }
        else
            {
            name[0] = '\0';
            k=0;
            while(tree[i] != '(' && tree[i] != ')' && tree[i] != ',' && tree[i] != ';' )
                {
                name[k] = tree[i];
                i++; k++;
                }
            name[k] = '\0';
            intname = assign_taxa_name(name, FALSE);  /* we now have the taxa number */
            /* next we need this integer written as a string */
            name[0] = '\0';
            totext(intname, name);
            /* now put this string version of the taxa number onto the tree */
            k=0;
            while(name[k] != '\0')
                {
                tree1[j] = name[k];
                j++; k++;
                }
            }
        } /* end while */
    tree1[j] = ';';
    tree1[j+1] = '\0';
    i=0;
    /** now scan along the retain_supers list and convert them to ints and then do the path metric **/
    while(scores_retained_supers[i] != -1)
        {
        j=0;k=0;l=0;
        tree2[0] = '\0';
        /* scan along the tree on retained supers and change the names of the taxa to their numbers */
        while(retained_supers[i][j] != ';')
            {
            if(retained_supers[i][j] == '(' || retained_supers[i][j] == ')' || retained_supers[i][j] == ',' || retained_supers[i][j] == ';')
                {
                tree2[k] = retained_supers[i][j];
                j++; k++;
                }
            else
                {
                name[0] = '\0';
                l=0;
                while(retained_supers[i][j] != '(' && retained_supers[i][j] != ')' && retained_supers[i][j] != ',' && retained_supers[i][j] != ';' )
                    {
                    name[l] = retained_supers[i][j];
                    l++; j++;
                    }
                name[l] = '\0';
                intname = assign_taxa_name(name, FALSE);  /* we now have the taxa number */
                /* next we need this integer written as a string */
                name[0] = '\0';
                totext(intname, name);
                /* now put this string version of the taxa number onto the tree */
                l=0;
                while(name[l] != '\0')
                    {
                    tree2[k] = name[l];
                    l++; k++;
                    }
                }
            } /* end while */
        tree2[k] = ';';
        tree2[k+1] = '\0';
		
        
        /** now do the pathmetric for tree1 and for tree2 and see if they are different */
        /* initialise the super_scores array */
        for(j=0; j<number_of_taxa; j++)
            {
            for(k=j; k<number_of_taxa; k++)
                {
                scores1[j][k] = 0;
                scores1[k][j] = 0;
                scores2[j][k] = 0;
                scores2[k][j] = 0;
                }
            }
        pathmetric(tree1, scores1);  /*calculate the pathmetric for tree1 */        
        pathmetric(tree2, scores2);  /*calculate the pathmetric for tree2 */         
        temp = 0;
        /* score the differences between the two trees */
        for(j=0; j<number_of_taxa; j++)
            {
            for(k=(j+1); k<number_of_taxa; k++)
                {
                if(scores1[j][k] != scores2[j][k])
                    {
                    if(scores1[j][k] > scores2[j][k]) 
                        {
                        temp += scores1[j][k] - scores2[j][k];
                        }
                    else
                        {
                        temp += scores2[j][k] - scores1[j][k];
                        }
                    }
                }
            }
        
        i++;
        if(temp == 0)
            different=FALSE;
        } /* end while */

    for(i=0; i<number_of_taxa; i++)
        {
        free(scores1[i]);
        free(scores2[i]);
        }
    free(scores2);
    free(name);
    return(different);
    }
     
    
/* Next is the code needed for the Baum/ragan coding scheme */
int coding(int nrep, int search, int ptpreps)
    {
    int i=0, j=0, k=0, nreps=10,  parenthesis = 0, count =0, *tracking = NULL, total = 0, **BR_coding = NULL, nodecount = 0, position = 0, split_count = 0, calculate_inhouse = FALSE;
    char number[100], string[TREE_LENGTH], filename[1000], **temptrees = NULL;
    int x=0, one_in_this = FALSE, zero_in_this = FALSE, swap=3, addseq=4, error=FALSE, *num_fund_taxa = NULL, njbuild = FALSE, weighted = FALSE;

    filename[0] = '\0';
    for(i=0; i<num_commands; i++)
        {       
        if(strcmp(parsed_command[i], "nreps") == 0 && (strcmp(parsed_command[0], "boot") != 0 && strcmp(parsed_command[0], "bootstrap") != 0))
            {
            nreps = toint(parsed_command[i+1]);
            
            if(nreps == 0)
                {
                printf("Error: '%s' is an invalid integer number to be assigned to nreps\n", parsed_command[i+1]);
                error = TRUE;
                }
            }

        if(strcmp(parsed_command[i], "hsreps") == 0 )
            {
            nreps = toint(parsed_command[i+1]);
            
            if(nreps == 0)
                {
                printf("Error: '%s' is an invalid integer number to be assigned to nreps\n", parsed_command[i+1]);
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
                    printf("Error: option %s not valid analysis\n", parsed_command[i+1]);
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
                printf("Error addseq option %s not known\n", parsed_command[i+1]);
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
                printf("Error swap option %s not known\n", parsed_command[i+1]);
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
						if(strcmp(filename, "") == 0) fprintf(BR_file,"MRP.tree Format=Phylip");    
						else fprintf(BR_file,"%s Format=Phylip", filename);
						}
					if(search==0) 
						{
						fprintf(BR_file,"set increase=auto notifybeep=no errorbeep=no;\n\talltrees;\n\tshowtrees;\n\tsavetrees FILE=");
						if(strcmp(filename, "") == 0) fprintf(BR_file,"MRP.tree Format=Phylip");    
						else fprintf(BR_file,"%s Format=Phylip", filename);
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
    return(error);
    }

int MRP_matrix(char **trees, int num_trees, int consensus)
	{
	int i=0, j=0, k=0, count =0, x=0, *tracking = NULL, total = 0, **BR_coding = NULL, nodecount = 0, position = 0, split_count = 0;
    char number[100], string[TREE_LENGTH];
	int *num_fund_taxa = NULL, one_in_this = FALSE, zero_in_this = FALSE;
	
	num_fund_taxa = malloc(num_trees*sizeof(int));
	if(!num_fund_taxa) memory_error(68);
	for(i=0; i<num_trees; i++) num_fund_taxa[i] = 0;
	
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
	return(position);	
	}
		

                    
        
 void set_parameters(void)
    {
    int i=0;
    for(i=0; i<num_commands; i++)
        {
        if(strcmp(parsed_command[i], "criterion") == 0)
            {
            if(strcmp(parsed_command[i+1], "mrp") == 0 || strcmp(parsed_command[i+1], "MRP") == 0)
                {
                if(criterion != 1) printf("Scoring criterion set to matrix representation using parsimony (MRP)\n");
                criterion = 1;
                }
            else
                {
                if(strcmp(parsed_command[i+1], "dfit") == 0 || strcmp(parsed_command[i+1], "DFIT") == 0)
                    {
                    if(criterion != 0) printf("Scoring criterion set to best distance fit (DFIT)\n");
                    criterion = 0;
                    }
                else
                    {
                    if(strcmp(parsed_command[i+1], "sfit") == 0 || strcmp(parsed_command[i+1], "SFIT") == 0)
                        {
                        if(criterion != 2) printf("Scoring criterion set to maximum split fit (SFIT)\n");
                        criterion = 2;
                        }
                    else
                        {
                        if(strcmp(parsed_command[i+1], "qfit") == 0 || strcmp(parsed_command[i+1], "QFIT") == 0)
                            {
                            if(criterion != 3) printf("Scoring criterion set to maximum quartet fit (QFIT)\n");
                            criterion = 3;
                            }
                        else
                            {
                            if(strcmp(parsed_command[i+1], "avcon") == 0 || strcmp(parsed_command[i+1], "AVCON") == 0)
                                {
								if(criterion != 4) printf("Scoring criterion set to average consensus (AVCON)\n");
								criterion = 4;
                                }
                            else
								{
								if(strcmp(parsed_command[i+1], "recon") == 0 || strcmp(parsed_command[i+1], "RECON") == 0)
									{
									if(criterion != 5) printf("Scoring criterion set to duplication and loss reconstruction (RECON)\n");
									criterion = 5;
									}
								else
									printf("Error: %s not known as criterion tpye\n", parsed_command[i+1]);
								}
							}
                        }
                    }
                }
            }
        }
    
    }
 

/* this function calculates the score of a supertree to a set of source trees using the criterion of Matrix representation using Compatibility */
/* This is analogous to a clique analysis */
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
                         


/* this function calculates the score of a supertree to a set of source trees using the criterion of Matrix representation using Compatibility */
/* This is analogous to a clique analysis */
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





/* This function is a second method of tree traversal to find the optimum tree. Because the NNI and SPS (subtree pruning and swapping) algorithm
	is such a gentle method, we need something that will change the actual structure of the tree. each node will be pruned out of its present 
	position on the tree and will be grafted into a new position on a new node somewhere else in the tree, this will continue until a better score
	is found. This is called SPR (subtree pruning and regrafting).
   The algorithm is as follows:
		1) if newbie is not already created, create a new node called newbie
   		2) pick a node x, and remove from the present position and shrink the tree
		3) make node x the daughter of the newly created node
		4) traversing through the tree (limited to a certain number of nodes away), for each node (taxa or pointer)
			i)   replace the node with newbie
			ii)  make the node a sibling of node x (which will mean its a daughter of newbie)
			iii) check the score, if better then STOP. If not better then put the node back and move onto the next node to replace with newbie
		5) If the whole tree is traversed  (or we have reached the limit of nodes) and no better score is found then put node x back where it was
*/

int spr_new(struct taxon * master, int maxswaps, int numspectries, int numgenetries)
	{
	struct taxon * start = NULL, * newbie = NULL, *newbie_backup = NULL, * latest = NULL, * tmp = NULL, *copy = NULL, *temper = NULL, *temper1 = NULL, *master_copy = NULL, *position = NULL;
	int better_score = FALSE, lastinline = FALSE, numofsiblings = 0, donenextlevel = FALSE, i=0, j=0, x=0, y=0, q=0, r=0;
	char *debugtree = NULL;
	
	debugtree = malloc(TREE_LENGTH*sizeof(char));
	debugtree[0] = '\0';
    

	start = master;
	print_tree(master, debugtree);
	debugtree[0] = '\0';



	/* As opposed to travelling down the tree in a recursive manner we shall implement a methid which replaces the tree used after each regrafting */
	
	/* 1) First make a copy of the master tree passed to the application, called master_copy */
	
	temp_top = NULL;
	duplicate_tree(master, NULL); /* make a copy of the tree */
	master_copy = temp_top;
	temp_top = NULL;
	/* end 1) */
	/* 2) count the number of nodes that it is possible to break this tree (internal and external branches) */
	
	y = number_tree(master, 0); /* count the parts of newbie */
							debugtree[0] = '\0';
							print_tree_withinternals(master, debugtree);
							debugtree[0] = '\0';

	
	/* end 2) */
	for(x=0; x<y; x++) /* for each branch of the master tree, create an instance where it is bisected at this point */
		{
		/* create a new instance of "master" from master_copy */
		if(master != NULL)
			dismantle_tree(master);
		master = NULL;
		temp_top = NULL;
		duplicate_tree(master_copy, NULL); /* make a copy of the tree */
		master = temp_top;
		tree_top = master;
		temp_top = NULL;
		number_tree(master, 0); 
		
		
		/* find the branch named "j" */
		position = NULL;
		position = get_branch(master, x);
		
		/* count number ofsiblings */
		/* rewind to start */
		tmp = position;
		while(tmp->prev_sibling != NULL) tmp = tmp->prev_sibling;
		/* count number of sibligs */
		numofsiblings = 0;
		while(tmp != NULL) 
			{
			numofsiblings++;
			tmp = tmp->next_sibling;
			}
		
		if(position != master) /* if we are at the top of the tree, don't do anything */
			{
			/* 3) break the tree at this point and create "newbie" */

			/* move position to be a daughter of newbie */
			if(position->next_sibling != NULL) latest = position->next_sibling;
			else
				{ 
				latest = position->prev_sibling;
				lastinline = TRUE;
				}

			if(position->prev_sibling != NULL) (position->prev_sibling)->next_sibling = position->next_sibling;
			if(position->next_sibling != NULL) (position->next_sibling)->prev_sibling = position->prev_sibling;
			if(position->parent != NULL)
				{
				(position->parent)->daughter = position->next_sibling;
				(position->next_sibling)->parent = position->parent;
				}
			position->next_sibling = NULL;
			position->prev_sibling = NULL;
			
			if(numofsiblings == 2 && latest->parent != NULL)  /* put the other sibling on the previous level instead of the pointer sibling that points to this level */
				{
				if((latest->parent)->prev_sibling != NULL) ((latest->parent)->prev_sibling)->next_sibling = latest;
				if((latest->parent)->next_sibling != NULL) ((latest->parent)->next_sibling)->prev_sibling = latest;
				if((latest->parent)->parent != NULL) ((latest->parent)->parent)->daughter = latest;
				latest->prev_sibling = (latest->parent)->prev_sibling;
				latest->next_sibling = (latest->parent)->next_sibling;
				tmp = latest->parent;
				latest->parent = (latest->parent)->parent;
				if(master == tmp)
				 {
				 master = latest;
				 tree_top = master;
				 }
				tmp->next_sibling = NULL;
				tmp->prev_sibling = NULL;
				tmp->daughter = NULL;
				tmp->parent = NULL;
				free(tmp);
				}
			newbie = position;
				print_tree(newbie, debugtree);
				debugtree[0] = '\0';
			
			tmp = latest;
			while(tmp->prev_sibling != NULL) tmp = tmp->prev_sibling;  /* rewinding */
			/**** Now that we have the subtree that we want to regraft, identify the taxa in that subtree for the purposes of speeding up the calculations */
			for(i=0; i<number_of_taxa; i++)presenceof_SPRtaxa[i] = FALSE;
			identify_taxa(newbie, presenceof_SPRtaxa);

			/* end 3) */
			

			/* At this point there are two paths that can be taken. */
			/* If we are doing a spr search, then we use the "regraft" function to find all possible positions to regraft the subtree onto the original tree */
			/* If we are doing a TBR search we need to reroot the subtree "newbie" at all rootings and for each, use "regraft" to try all positions of the rerooted subtree on the original tree */
			/* This makes the search longer, but allows the better penetration of treespace */
			if(method == 2) /* spr method */
				{
				temper = make_taxon();
				temper->daughter = newbie; 
				newbie->parent = temper;
				newbie = temper;
				temper= NULL;
				newbie->spr = TRUE;
				
				print_tree(newbie, debugtree);
				debugtree[0] = '\0';
				
				better_score = regraft(latest, newbie, NULL, 1, maxswaps, numspectries, numgenetries);
				
				}
			if(method == 3) /*tbr method */
				{

				i = number_tree(newbie, 0); /* count the parts of newbie */
				
				if(i > 4) /* no point rerooting a clade with only 2 taxa or a single taxa */
					{
					temper = newbie->daughter;
					temper->parent = NULL;
					free(newbie);
					newbie = temper;
					temper = NULL;

					print_tree(newbie, debugtree);
					strcat(debugtree, ";");
					unroottree(debugtree);
					
					dismantle_tree(newbie);
					newbie = NULL;
					temp_top = NULL;
					tree_build(1, debugtree, newbie, FALSE, -1);
					newbie = temp_top;
					temp_top = NULL;
					
					debugtree[0] = '\0';
					print_tree(newbie, debugtree);
					debugtree[0] = '\0';	
					
					r = number_tree(newbie, 0); /* count the parts of newbie */
					temp_top = NULL;
					duplicate_tree(newbie, NULL); /* make a copy of the subtree */
					copy = temp_top;
					temp_top = NULL;
					
					print_tree_withinternals(copy, debugtree);
					debugtree[0] = '\0';
					
				
					for(q=0; q<r; q++) /* reroot newbie at each of the "r" positions on the subtree */
						{
						temper1 = get_branch(newbie, q);
						
						
						temp_top = newbie;
						reroot_tree(temper1);
						newbie = temp_top;
							print_tree_withinternals(newbie, debugtree);
							debugtree[0] = '\0';

						
						temper = make_taxon();
						temper->daughter = newbie; 
						newbie->parent = temper;
						newbie = temper;
						temper= NULL;
						newbie->spr = TRUE;


						better_score = regraft(tmp, newbie, NULL, 1, maxswaps, numspectries, numgenetries);
						

						if(better_score || user_break) q = r;
						else
							{
							dismantle_tree(newbie);
							newbie = NULL;
							temp_top = NULL;

							duplicate_tree(copy, NULL); /* make a copy of the subtree */
							newbie = temp_top;
							temp_top = NULL;
							number_tree(newbie, 0);
							newbie->spr = TRUE;
							}
						}
					if(copy != NULL)dismantle_tree(copy);
					
					copy = NULL;
					}
				else
					{
					temper = make_taxon();
					temper->daughter = newbie; 
					newbie->parent = temper;
					newbie = temper;
					temper= NULL;
					newbie->spr = TRUE;
					
					print_tree(newbie, debugtree);
					debugtree[0] = '\0';
					
					better_score = regraft(latest, newbie, NULL, 1, maxswaps, numspectries, numgenetries); 
					}
				}
			if(!better_score && newbie != NULL) dismantle_tree(newbie);
			newbie = NULL;
            }
		if(better_score == TRUE || (position == branchpointer) || donenextlevel || tried_regrafts >= maxswaps || user_break) 
			{
			x = y;
			}
		}
		free(debugtree);
		if(master_copy != NULL)
			{
			dismantle_tree(master_copy);
			master_copy = NULL;
			}

        return(better_score);
	}







/* This function is a second method of tree traversal to find the optimum tree. Because the NNI and SPS (subtree pruning and swapping) algorithm
	is such a gentle method, we need something that will change the actual structure of the tree. each node will be pruned out of its present 
	position on the tree and will be grafted into a new position on a new node somewhere else in the tree, this will continue until a better score
	is found. This is called SPR (subtree pruning and regrafting).
   The algorithm is as follows:
		1) if newbie is not already created, create a new node called newbie
   		2) pick a node x, and remove from the present position and shrink the tree
		3) make node x the daughter of the newly created node
		4) traversing through the tree (limited to a certain number of nodes away), for each node (taxa or pointer)
			i)   replace the node with newbie
			ii)  make the node a sibling of node x (which will mean its a daughter of newbie)
			iii) check the score, if better then STOP. If not better then put the node back and move onto the next node to replace with newbie
		5) If the whole tree is traversed  (or we have reached the limit of nodes) and no better score is found then put node x back where it was
*/

int spr(struct taxon * position, int maxswaps, int numspectries, int numgenetries)
	{
	struct taxon * start = NULL, * newbie = NULL, * latest = NULL, * tmp = NULL, *copy = NULL, *temper = NULL, *temper1 = NULL;
	int better_score = FALSE, lastinline = FALSE, numofsiblings = 0, donenextlevel = FALSE, i=0, j=0;
	char *debugtree = NULL;
	
	debugtree = malloc(TREE_LENGTH*sizeof(char));
	debugtree[0] = '\0';
    
		
								start = position;
								print_tree(position, debugtree);
								printf("%s\n", debugtree);
								debugtree[0] = '\0';

	
        /* travel down to the bottom of the tree */
	while(position != NULL && better_score == FALSE && (position != branchpointer) && !donenextlevel && tried_regrafts < maxswaps && !user_break)
		{
		if(position->daughter != NULL && !position->spr && !user_break)
			{
			if(donenextlevel == FALSE)
				{
				better_score = spr(position->daughter, maxswaps, numspectries, numgenetries);
				donenextlevel = TRUE;
				}
			}
		numofsiblings++;
                
		if(!user_break)position = position->next_sibling;
		}
        if(!donenextlevel && !user_break)
            {
          /*  if(position == branchpointer && !better_score && !user_break)
                {
                while(position != NULL && !user_break)
                    {
                    numofsiblings++;
                    position = position->next_sibling;
                    }
                }
	*/		printf("numsiblingd=%d\n", numofsiblings);
            position = start;
            if(!better_score && start->parent != NULL && tried_regrafts < maxswaps && !user_break)
                    {
                    while(position != NULL&& start->parent != NULL && !better_score && tried_regrafts < maxswaps && !user_break)
                            {
                            /*if(newbie == NULL)
                                    newbie = make_taxon();
                            newbie->spr = TRUE; */
                            /* move position to be a daughter of newbie */
                            if(position->next_sibling != NULL) latest = position->next_sibling;
                            else
                                { 
                                latest = position->prev_sibling;
                                lastinline = TRUE;
                                }
                           /* newbie->daughter = position; */
                            if(position->prev_sibling != NULL) (position->prev_sibling)->next_sibling = position->next_sibling;
                            if(position->next_sibling != NULL) (position->next_sibling)->prev_sibling = position->prev_sibling;
                            if(position->parent != NULL)
                                {
                                (position->parent)->daughter = position->next_sibling;
                                (position->next_sibling)->parent = position->parent;
                                }
                            position->next_sibling = NULL;
                            position->prev_sibling = NULL;
                           /* position->parent = newbie;*/
                            
                            if(numofsiblings == 2 && latest->parent != NULL)  /* put the other sibling on the previous level instead of the pointer sibling that points to this level */
                                {
                                if((latest->parent)->prev_sibling != NULL) ((latest->parent)->prev_sibling)->next_sibling = latest;
                                if((latest->parent)->next_sibling != NULL) ((latest->parent)->next_sibling)->prev_sibling = latest;
                                if((latest->parent)->parent != NULL) ((latest->parent)->parent)->daughter = latest;
                                latest->prev_sibling = (latest->parent)->prev_sibling;
                                latest->next_sibling = (latest->parent)->next_sibling;
                                tmp = latest->parent;
                                latest->parent = (latest->parent)->parent;
                                if(tree_top == tmp) tree_top = latest;
                                tmp->next_sibling = NULL;
                                tmp->prev_sibling = NULL;
                                tmp->daughter = NULL;
                                tmp->parent = NULL;
                                free(tmp);
                                }
							newbie = position;
								print_tree(newbie, debugtree);
								printf("NEWBIE=%s\n", debugtree);
								debugtree[0] = '\0';
							
                            tmp = latest;
                            while(tmp->prev_sibling != NULL) tmp = tmp->prev_sibling;  /* rewinding */
							/**** Now that we have the subtree that we want to regraft, identify the taxa in that subtree for the purposes of speeding up the calculations */
							for(i=0; i<number_of_taxa; i++)presenceof_SPRtaxa[i] = FALSE;
							identify_taxa(newbie, presenceof_SPRtaxa);
							
							/* At this point there are two paths that can be taken. */
							/* If we are doing a spr search, then we use the "regraft" function to find all possible positions to regraft the subtree onto the original tree */
							/* If we are doing a TBR search we need to reroot the subtree "newbie" at all rootings and for each, use "regraft" to try all positions of the rerooted subtree on the original tree */
							/* This makes the search longer, but allows the better penetration of treespace */
							if(method == 2) /* spr method */
								{
								temper = make_taxon();
								temper->daughter = newbie; 
								newbie->parent = temper;
								newbie = temper;
								temper= NULL;
								newbie->spr = TRUE;
								
								print_tree(newbie, debugtree);
								printf("%s\n", debugtree);
								debugtree[0] = '\0';
								
								better_score = regraft(tmp, newbie, NULL, 1, maxswaps, numspectries, numgenetries);
								
								}
							if(method == 3) /*tbr method */
								{

								i = number_tree(newbie, 0); /* count the parts of newbie */
								
								if(i > 4) /* no point rerooting a clade with only 2 taxa or a single taxa */
									{
									temper = newbie->daughter;
									temper->parent = NULL;
									free(newbie);
									newbie = temper;
									temper = NULL;

									print_tree(newbie, debugtree);
									strcat(debugtree, ";");
									unroottree(debugtree);
									printf("unrooted = %s\n", debugtree);
									
									dismantle_tree(newbie);
									newbie = NULL;
									temp_top = NULL;
									tree_build(1, debugtree, newbie, FALSE, -1);
									newbie = temp_top;
									temp_top = NULL;
									
									debugtree[0] = '\0';
									print_tree(newbie, debugtree);
									printf("now = %s\n", debugtree);
									debugtree[0] = '\0';	
									
									i = number_tree(newbie, 0); /* count the parts of newbie */
									printf("i=%d\n", i);
									temp_top = NULL;
									duplicate_tree(newbie, NULL); /* make a copy of the subtree */
									copy = temp_top;
									temp_top = NULL;
									
									print_tree_withinternals(copy, debugtree);
									printf("copy = %s\n", debugtree);
									debugtree[0] = '\0';
									
								
									for(j=0; j<i; j++) /* reroot newbie at each of the "i" positions on the subtree */
										{
										printf("j=%d\n", j);
										temper1 = get_branch(newbie, j);
										
										
										printf("I\n");
										temp_top = newbie;
										reroot_tree(temper1);
										newbie = temp_top;
											print_tree_withinternals(newbie, debugtree);
											printf("rerooted: %d\t%s\n", temper1->tag, debugtree);
											debugtree[0] = '\0';

										printf("II\n");
										
										temper = make_taxon();
										temper->daughter = newbie; 
										newbie->parent = temper;
										newbie = temper;
										temper= NULL;
										newbie->spr = TRUE;


										printf("a=\n");
										better_score = regraft(tmp, newbie, NULL, 1, maxswaps, numspectries, numgenetries);
										
										printf("b\n");

										if(better_score || user_break) j = i;
										else
											{
											printf("change\n");
											dismantle_tree(newbie);
											newbie = NULL;
											temp_top = NULL;
											printf("dismantled\n");

											duplicate_tree(copy, NULL); /* make a copy of the subtree */
											newbie = temp_top;
											temp_top = NULL;
											number_tree(newbie, 0);
											newbie->spr = TRUE;
											printf("end\n");
											}
										}
									}
								else
									{
									temper = make_taxon();
									temper->daughter = newbie; 
									newbie->parent = temper;
									newbie = temper;
									temper= NULL;
									newbie->spr = TRUE;
									
									print_tree(newbie, debugtree);
									printf("%s\n", debugtree);
									debugtree[0] = '\0';
									
									better_score = regraft(tmp, newbie, NULL, 1, maxswaps, numspectries, numgenetries);
									}
								}
                            if(!better_score && !user_break)
                                    {
                                    /* we have to put the taxa back where it was before */
                                    if(numofsiblings == 2)
                                            {
                                            newbie->next_sibling = latest->next_sibling;
                                            newbie->prev_sibling = latest->prev_sibling;
                                            newbie->parent = latest->parent;
    
                                            if(newbie->next_sibling != NULL) (newbie->next_sibling)->prev_sibling = newbie;
                                            if(newbie->prev_sibling != NULL) (newbie->prev_sibling)->next_sibling = newbie;
                                            if(newbie->parent != NULL) (newbie->parent)->daughter = newbie;
                                            if(tree_top == latest) tree_top = newbie;
    
                                            if(lastinline)
                                                    {
                                                    latest->next_sibling = position;
                                                    latest->prev_sibling = position->prev_sibling;
                                                    position->prev_sibling = latest;
                                                    position->next_sibling = NULL;
                                                    latest->parent = position->parent;
                                                    if(latest->parent != NULL) (latest->parent)->daughter = latest;
                                                    position->parent = NULL;
    
    
                                                    }
                                            else
                                                    {
                                                    latest->prev_sibling = position;
                                                    latest->next_sibling = position->next_sibling;
                                                    position->next_sibling = latest;
                                                    if(latest->next_sibling != NULL) (latest->next_sibling)->prev_sibling = latest;
                                                    latest->parent = NULL;
                                                    }
                                            newbie = NULL;
                                            }
                                    else
                                            {
                                            if(lastinline)
                                                    {
                                                    latest->next_sibling = position;
                                                    position->prev_sibling = latest;
                                                    position->next_sibling = NULL;
                                                    position->parent = NULL;
                                                    newbie->daughter = NULL;
                                                    }
                                            else
                                                    {
                                                    position->prev_sibling = latest->prev_sibling;
                                                    if(position->prev_sibling != NULL) (position->prev_sibling)->next_sibling = position;
                                                    position->next_sibling = latest;
                                                    latest->prev_sibling = position;
                                                    position->parent = latest->parent;
                                                    if(position->parent != NULL) (position->parent)->daughter = position;
                                                    latest->parent = NULL;
                                                    newbie->daughter = NULL;
                                                    }
                                            }
    
    
                                    }
                            else
                                    {
                                    newbie = NULL;
                                    branchpointer = latest;
                                    }
                            position = position->next_sibling;
                            }
                    }
            }
		free(debugtree);
        return(better_score);
	}




void reset_spr (struct taxon *position)
    {
    
    while(position != NULL)
        {
        position->spr = FALSE;
        if(position->daughter != NULL) reset_spr(position->daughter);
        position = position->next_sibling;
        }
    }
    


int regraft(struct taxon * position, struct taxon * newbie, struct taxon * last, int steps, int maxswaps, int numspectries, int numgenetries)
	{
	struct taxon * start = NULL, * latest = NULL;
	int better_score = FALSE, lastinline = FALSE, j=0, different = TRUE, i=0;
	float tmpscore = 0, *tmp_fund_scores = NULL;
	char *best_tree = NULL, *temptree = NULL;

	tmp_fund_scores = malloc(Total_fund_trees*sizeof(float));
	best_tree = malloc(TREE_LENGTH*sizeof(char));
	if(!best_tree) memory_error(72);
	temptree = malloc(TREE_LENGTH*sizeof(char));
	if(!temptree) memory_error(73);
	
	
	best_tree[0] = '\0';
	temptree[0] = '\0';
            
        
	while(position->prev_sibling != NULL && !user_break) position = position->prev_sibling; /* rewinding */
	start = position;
	
	while(position != NULL && !better_score && tried_regrafts < maxswaps && !user_break)
		{
		if(steps < number_of_steps && position->parent != NULL && position->parent != last && !user_break) better_score = regraft(position->parent, newbie, position, steps+1, maxswaps, numspectries, numgenetries);  /* go up the tree */
		position = position->next_sibling;
		}
	
	position = start;
	if(!better_score && !user_break)
		{
		while(position != NULL && !better_score && tried_regrafts < maxswaps && !user_break)
			{
			if(steps < number_of_steps && !user_break)
				if(position->daughter != NULL && position->daughter != last && !user_break) better_score = regraft(position->daughter, newbie, position, steps+1, maxswaps, numspectries, numgenetries);

			position = position->next_sibling;
			}
		position = start;

		while(position != NULL && !better_score && tried_regrafts < maxswaps && !user_break)
			{
			if(!better_score && steps <= number_of_steps && !user_break)
				{
				/* regraft newbie here */
                                
                                
			/*	if((newbie->daughter)->name != -1)printf("newbie daughter = %d\n", newbie->daughter->name);
				else printf("newbie daughter = pointer\n");
				if(position->name != -1) printf("new sibling daughter = %d\n", position->name);
				else printf("pointer to %d\n", position->daughter->name);
	*/
				newbie->next_sibling = position->next_sibling;
				newbie->prev_sibling = position->prev_sibling;
				newbie->parent = position->parent;
				
				if(newbie->prev_sibling != NULL) (newbie->prev_sibling)->next_sibling = newbie;
				if(newbie->next_sibling != NULL) (newbie->next_sibling)->prev_sibling = newbie;
				if(newbie->parent != NULL) (newbie->parent)->daughter = newbie;
				if(tree_top == position) tree_top = newbie;

				/* put position as a daughter of newbie */
				
				(newbie->daughter)->next_sibling = position;
				position->prev_sibling = newbie->daughter;
				position->next_sibling = NULL;
				position->parent = NULL;

				/* the new tree is complete, so now score this new super tree */

                                for(i=0; i<Total_fund_trees; i++) tmp_fund_scores[i] = sourcetree_scores[i];
                                interval2 = time(NULL);
                                if(difftime(interval2, interval1) > 5) /* every 5 seconds print a dot to the screen */
                                    {
                                    printf("=");
                                    fflush(stdout);
                                    interval1 = time(NULL);
                                    }
                        
                                
                                if(check_taxa(tree_top) == number_of_taxa && !user_break)
                                    {
                                    tried_regrafts++;
									NUMSWAPS++;
                                    
                                    strcpy(best_tree, "");
                                    print_named_tree(tree_top, best_tree);
                                    strcat(best_tree, ";");
                                    if(criterion == 0) tmpscore = compare_trees(TRUE);  /* calculate the distance from the super tree to all the fundamental trees */
                                    if(criterion == 2)  /* calculate the distance using the MRC criterion */
                                        {
										strcpy(temptree, "");
                                        print_tree(tree_top, temptree);
                                        strcat(temptree, ";");
                                        unroottree(temptree);
                                        tmpscore = (float)MRC(temptree);
                                        }
                                    if(criterion == 3)
                                        {
										strcpy(temptree, "");
                                        print_tree(tree_top, temptree);
                                        strcat(temptree, ";");
                                        unroottree(temptree);
                                        tmpscore = (float)quartet_compatibility(temptree);
                                        }
														
									if(criterion==5)
										{
										strcpy(temptree, "");
										print_tree(tree_top, temptree);
										strcat(temptree, ";");
										tmpscore = get_recon_score(temptree, numspectries, numgenetries);
										}

                                    /*printf("scores_retained_supers[0] = %f\n", scores_retained_supers[0]); */
                                    if(scores_retained_supers[0] == -1 || sprscore == -1)  /* Then this is the first tree to be checked */
                                        {
										if(BESTSCORE == -1)BESTSCORE = tmpscore;
                                        if(scores_retained_supers[0] == -1)
                                            {
											retained_supers[0] = realloc(retained_supers[0], (strlen(best_tree)+10)*sizeof(char));
                                            strcpy(retained_supers[0], best_tree);
                                            scores_retained_supers[0] = tmpscore;
                                            }
                                        sprscore = tmpscore;
                                        better_score = TRUE;
                                        }
                                    else
                                        {
                                        if(tmpscore < scores_retained_supers[0] || tmpscore < sprscore)
                                            {
											if(tmpscore < BESTSCORE) BESTSCORE= tmpscore;
                                            if(tmpscore < sprscore)
                                                {
                                                better_score = TRUE;
                                                sprscore = tmpscore;
												
                                                }
                                            if(tmpscore < scores_retained_supers[0])
                                                {
                                                better_score = TRUE;
                                                /* we need a string copy of this tree */
                                                strcpy(best_tree, "");
                                                print_named_tree(tree_top, best_tree);
                                                strcat(best_tree, ";");
												retained_supers[0] = realloc(retained_supers[0], (strlen(best_tree)+10)*sizeof(char));
                                                strcpy(retained_supers[0], best_tree);
                                                scores_retained_supers[0] = tmpscore;
                                                j=1;
                                                while(scores_retained_supers[j] != -1 && j < number_retained_supers)
                                                    {
                                                    strcpy(retained_supers[j], "");
                                                    scores_retained_supers[j] = -1;
                                                    j++;
                                                    }
                                                }
                                            }
                                        else
                                            {
                                            if(tmpscore == scores_retained_supers[0])
                                                {
                                                j = 0;
                                                strcpy(best_tree, "");
                                                print_named_tree(tree_top, best_tree);
                                                strcat(best_tree, ";");
                                                different = TRUE;
                                                if(!check_if_diff_tree(best_tree))
													{
                                                    different = FALSE; /* This is here to check if this tree with the same score actually has a different topology from those already stroed */
													}
                                                if(different)
                                                    {
                                                    while(scores_retained_supers[j] != -1)
                                                        {
                                                        j++;
                                                        if(j+1 == number_retained_supers) reallocate_retained_supers();
                                                        }
													retained_supers[j] = realloc(retained_supers[j], (strlen(best_tree)+10)*sizeof(char));
                                                    strcpy(retained_supers[j], best_tree);
                                                    scores_retained_supers[j] = tmpscore;
                                                    }
                                                }
            
											
                                            /* put the node back */
                                            position->next_sibling = newbie->next_sibling;
                                            position->prev_sibling = newbie->prev_sibling;
                                            position->parent = newbie->parent;
                                            if(tree_top == newbie) tree_top = position;
                                            
                                            if(position->prev_sibling != NULL) (position->prev_sibling)->next_sibling = position;
                                            if(position->next_sibling != NULL) (position->next_sibling)->prev_sibling = position;
                                            if(position->parent != NULL) (position->parent)->daughter = position;
                                            
                                            newbie->next_sibling = NULL;
                                            newbie->prev_sibling = NULL;
                                            newbie->parent = NULL;
                                            (newbie->daughter)->next_sibling = NULL;
                                            
                                            /* the node is now back in position */
                                            for(i=0; i<Total_fund_trees; i++) sourcetree_scores[i] = tmp_fund_scores[i];
                                            
                                            }  /*else */
                                        } /* else */
                                    } /* if */
                                else /* the new supertree is not valid */
                                    {
									if(!user_break)
										{
										/* put the node back */
										position->next_sibling = newbie->next_sibling;
										position->prev_sibling = newbie->prev_sibling;
										position->parent = newbie->parent;
										if(tree_top == newbie) tree_top = position;

										if(position->prev_sibling != NULL) (position->prev_sibling)->next_sibling = position;
										if(position->next_sibling != NULL) (position->next_sibling)->prev_sibling = position;
										if(position->parent != NULL) (position->parent)->daughter = position;

										newbie->next_sibling = NULL;
										newbie->prev_sibling = NULL;
										newbie->parent = NULL;
										(newbie->daughter)->next_sibling = NULL;
										for(i=0; i<Total_fund_trees; i++) sourcetree_scores[i] = tmp_fund_scores[i];
										}
                                    /* the node is now back in position */
                                    }
                                } /* if */
                            
			position = position->next_sibling;
			} /* while */
		}
		free(tmp_fund_scores);
        free(temptree);
        free(best_tree);
	return(better_score);
	}

 void get_lengths(struct taxon *position)
    {
    while(position != NULL)
        {
        if(position->daughter != NULL)
            {
            get_lengths(position->daughter);
			printf("%f\n", position->length);
            if(position->length > largest_length) largest_length = position->length;
            }
        position = position->next_sibling;
        }
    }   
    
       
       
          

/*** The next three functions is designed to produce x  co-ordinates for each node in the tree for graphical representation ***/

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


void tree_coordinates(char *tree, int bootstrap, int build, int mapping, int fundnum)
    {
    char **treearray = NULL;
    int i=0, j=0, taxa_count = 0, deepest = 0, acrosssize = 20, upsize = 5;
    
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
		taxaorder = 0;
		tree_build(1, tree, tree_top, 1, fundnum);
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
            printf("%c", treearray[i][j]);
            }
        printf("\n");
        }
	
	if(bootstrap)
		fclose(psfile);
		
    for(i=0; i<taxa_count*2+1; i++)
        free(treearray[i]);
    free(treearray);
    }
    
void generatetrees(void)
	{
	int i, j, k, ntrees = 100, error = FALSE, gen_method = 1, random = TRUE, n = 20, data = 1, tree_rand_method = 1, super = 1, print_all_scores = FALSE, saveideal = FALSE;
	float *results = NULL;
	char *temptree = NULL, *rand_tree = NULL, filename[100], *pruned_tree = NULL, tmp[TREE_LENGTH], superfilename[1000], c;
	double min = -1, max = -1, a = 0.0, b = 1.0;
	FILE *outfile = NULL, *superfile = NULL, *allscores = NULL, *idealfile = NULL;
	
	filename[0] = '\0';
	strcpy(filename, "histogram.txt");
	for(i=0; i<num_commands; i++)
		{
		if(strcmp(parsed_command[i], "savescores") == 0)
			{
			if(strcmp(parsed_command[i+1], "yes") == 0)
				print_all_scores = TRUE;
			}
		if(strcmp(parsed_command[i], "method") == 0)
			{
			if(strcmp(parsed_command[i+1], "equiprobable") == 0)
				{
				gen_method = 1;
				}
			else
				{
				if(strcmp(parsed_command[i+1], "markovian") == 0)
					{
					gen_method = 2;
					}
				else
					{
					printf("Error: %s is not a valid method\n\n", parsed_command[i+1]);
					error = TRUE;
					}
				}
			}
		
		if(strcmp(parsed_command[i], "ntrees") == 0)
			{
			if(strcmp(parsed_command[i+1], "all") == 0)
				{
				random = FALSE;
				ntrees = sup;
				}
			else
				{
				ntrees = toint(parsed_command[i+1]);
				if(ntrees == 0) 
					{
					error = TRUE;
					printf("%d is not a valid value for ntrees\n", ntrees);
					}
				}
			}
						
		if(strcmp(parsed_command[i], "outfile") == 0)
			{
			strcpy(filename, parsed_command[i+1]);
			}
				
		if(strcmp(parsed_command[i], "nbins") == 0)
			{
			n = toint(parsed_command[i+1]);
			if(n == 0) 
				{
				error = TRUE;
				printf("%d is not a valid value for nbins\n", ntrees);
				}
			}
					
		if(strcmp(parsed_command[i], "sourcedata") == 0)
			{	
			if(strcmp(parsed_command[i+1], "real") == 0)
				data = 1;
			else
				{
				if(strcmp(parsed_command[i+1], "randomised") == 0)
					data = 2;
				else
					{
					if(strcmp(parsed_command[i+1], "ideal") == 0)
						data = 3;
					else
						{
						printf("Error: %s not valid sourcedata option\n", parsed_command[i+1]);
						error = TRUE;
						}
					}
				}
			}

		if(strcmp(parsed_command[i], "supertree") == 0)
			{
			if(strcmp(parsed_command[i], "memory") == 0)
				super = 1;
			else
				{
				superfilename[0] = '\0';
				strcpy(superfilename, parsed_command[i+1]);
				super = 2;
				}
			}
						
		if(strcmp(parsed_command[i], "savesourcetrees") == 0)
			{
			if(strcmp(parsed_command[i+1], "yes") == 0 || strcmp(parsed_command[i+1], "YES") == 0 )
				saveideal = TRUE;
			else
				{
				if(strcmp(parsed_command[i+1], "no") == 0 || strcmp(parsed_command[i+1], "NO") == 0 )
					saveideal = FALSE;
				else
					{
					printf("Error: %s not valid savesourcetrees option\n", parsed_command[i+1]);
					error = TRUE;
					}
				}
			}
		}
		

	if(criterion == 4 || criterion == 1)
		{
		printf("Error: This command is not currently implemented for MRP or Average consensus\n\tplease select another criterion to continue\n");
		error = TRUE;
		}
	if(data == 3 && super == 1 && trees_in_memory == 0)
		{
		printf("\n\nError: there are no trees in memory, a supertree is required to create ideal data\n\tYou have two choices:\n\t1) choose to run a Heuristic search to create a supertree\n\t2) specify a file containing the supertree to use\n\n\n");
		error = TRUE;
		}
	if(!error)
		{
		printf("\n\nGeneratetrees Settings:\n");
		printf("\tGenerating");
		if(random) printf(" %d random trees\n", ntrees);
		else printf(" all %.0f possible trees\n", sup);
		if(random)
			{
			printf("\tRandom supertree generator used: ");
			if(gen_method == 1) printf("Equiprobable\n");
			else printf("Markovian\n");
			}
		printf("\tNumber of bins = %d\n", n);
		printf("\tCalculating the tree scores using ");
		if(data == 1) printf("the real source trees\n");
		if(data == 2) printf("randomised versions of the source trees\n");
		if(data == 3) 
			{
			printf("idealised versions of the source trees\n");
			printf("\tIdeal data created from first supertree in ");
			if(super == 1)printf("memory\n");
			if(super == 2)printf("file named %s\n", superfilename);
			if(saveideal) printf("\tIdeal data to be saved to the file 'idealtrees.ph'\n");
			}
		printf("\tCriterion used to calculate supertree scores = ");
		if(criterion == 0) printf("dfit\n");
		if(criterion == 2) printf("sfit\n");
		if(criterion == 3) printf("qfit\n");
		printf("\tOutput file = %s\n", filename);
		if(print_all_scores)printf("\tSaving all %d supertree scores to file 'allscores.txt'\n", ntrees);
		
		outfile = fopen(filename, "w");
		
		if(!random) ntrees = sup;
		results = malloc(ntrees*sizeof(float));
		if(!results) memory_error(96);
		
		temptree = malloc(TREE_LENGTH*sizeof(char));
		if(!temptree) memory_error(97);
		temptree[0] = '\0';
		
		
		
		rand_tree = malloc(TREE_LENGTH*(sizeof(char)));
		if(!rand_tree) memory_error(98);
		rand_tree[0] = '\0';
		
		if(signal(SIGINT, controlc3) == SIG_ERR)
			{
			printf("An error occurred while setting a signal handler\n");
			}
		for(i=0; i<number_of_taxa; i++) presenceof_SPRtaxa[i] = '\0';
		
		if(super_scores == NULL)
			{
			super_scores = malloc(number_of_taxa*sizeof(int *));
			if(!super_scores)  memory_error(54);
				
			for(i=0; i<number_of_taxa; i++)
				{
				super_scores[i] = malloc(number_of_taxa*sizeof(int));
				if(!super_scores[i])  memory_error(55);
					
				for(j=0; j<number_of_taxa; j++)
					{
					super_scores[i][j] = 0;
					}
				}
			}
		else
			{
			for(i=0; i<number_of_taxa; i++)
				{
				for(j=i; j<number_of_taxa; j++)
					{
					super_scores[i][j] = 0;
					super_scores[j][i] = 0;
					}
				}
			}
		
		/***** source tree data has to be changed if the ideal of randomised options have been chosen */
		if(data != 1)
			{
			/** create somewhere to store the source trees during the operation **/
			
			stored_fund_scores = malloc(Total_fund_trees*sizeof(int**));
			if(stored_fund_scores == NULL) memory_error(38);
				
			for(i=0; i<Total_fund_trees; i++)
				{
				stored_fund_scores[i] = malloc((number_of_taxa)*sizeof(int*));
				if(stored_fund_scores[i] == NULL) memory_error(39);
					
				else
					{
					for(j=0; j<(number_of_taxa); j++)
						{
						stored_fund_scores[i][j] = malloc((number_of_taxa)*sizeof(int));
						if(stored_fund_scores[i][j] == NULL) memory_error(40);
							
						else
							{
							for(k=0; k<(number_of_taxa); k++)
								stored_fund_scores[i][j][k] = fund_scores[i][j][k];  /* copy the origninal fundamental scores for safe keeping */
							}
						}
					}
				}
			stored_presence_of_taxa = malloc(Total_fund_trees*sizeof(int *));
			if(!stored_presence_of_taxa) memory_error(41);
				
			for(i=0; i<Total_fund_trees; i++)
				{
				stored_presence_of_taxa[i] = malloc((number_of_taxa)*sizeof(int));
				if(!stored_presence_of_taxa[i]) memory_error(42);
					
				for(j=0; j<number_of_taxa; j++)
					{
					stored_presence_of_taxa[i][j] = presence_of_taxa[i][j];
					}
				}
		
		
			stored_funds = malloc(Total_fund_trees*sizeof(char *));
			if(!stored_funds) memory_error(31);
		
			for(i=0; i<Total_fund_trees; i++)
				{
				stored_funds[i] = malloc((strlen(fundamentals[i])+10)*sizeof(char));
				if(!stored_funds) memory_error(32);
		
				stored_funds[i][0] = '\0';
				strcpy(stored_funds[i], fundamentals[i]);   /* copy the original fundamentals for safe keeping */
				}

			
			if(data == 2)  /* randomise the source trees */
				{
                /***** Create the replicate of randomised fundamental trees  **********/
                for(j=0; j<Total_fund_trees; j++)
                    {
					
                    if(tree_rand_method == 1) randomise_tree(fundamentals[j]); /*equiprobable */
                    if(tree_rand_method == 2) randomise_taxa(fundamentals[j]);  /*Markovian*/
                    }
				if(saveideal)
					{
					idealfile = fopen("idealtrees.ph", "w");
					for(i=0; i<Total_fund_trees; i++)
						{
						strcpy(temptree, fundamentals[i]);
						returntree(temptree);
						fprintf(idealfile, "%s [%s]\n", temptree, tree_names[i]);
						}
					fclose(idealfile);
					printf("finished printing out trees\n");
					}	
				
                cal_fund_scores(FALSE);
				}
			if(data == 3) /* make an ideal dataset from the best supertree in memory */
				{
				if(super == 1 && trees_in_memory == 0)
					{
					
					}
				else
					{
					pruned_tree = malloc(TREE_LENGTH*sizeof(char));
					
					/****** We now need to build the Supertree in memory *******/
					if(tree_top != NULL)
						{
						dismantle_tree(tree_top);
						tree_top = NULL;
						}
					temp_top = NULL;
					if(super == 1)
						tree_build(1, retained_supers[0], tree_top, TRUE, -1);
					else
						{
						if((superfile = fopen(superfilename, "r")) == NULL)
							{
							printf("Error opening supertree file %s\n", superfilename);
							error = TRUE;
							}
						else
							{
							strcpy(tmp, "");
							i=0;
							c =fgetc(superfile);
							while(!feof(superfile) && c != ';')
								{
								tmp[i] = c;
								c = fgetc(superfile);
								i++;
								}
							tmp[i] = ';'; tmp[i+1] = '\0';
							fclose(superfile);
							tree_build(1, tmp, tree_top, TRUE, -1);
							}
						}
					tree_top = temp_top;
					temp_top = NULL;
					if(!error)
						{
						for(i=0; i< Total_fund_trees; i++)
							{
							prune_tree(tree_top, i);  /* Prune the supertree so that it has the same taxa as the fundamental tree i */
							shrink_tree(tree_top);    /* Shrink the pruned tree by switching off any internal nodes that are not needed */
							for(j=0; j<TREE_LENGTH; j++) pruned_tree[j] = '\0';
							pruned_tree[0] = '\0'; /* initialise the string */
							if(print_pruned_tree(tree_top, 0, pruned_tree, FALSE) >1)
								{
								tmp[0] = '\0';
								strcpy(tmp, "(");
								strcat(tmp, pruned_tree);
								strcat(tmp, ")");
								strcpy(pruned_tree, tmp);
								}
							strcat(pruned_tree, ";");

							strcpy(fundamentals[i], pruned_tree);
							reset_tree(tree_top);  /* reset the supertree for the next comparison */
							}
						if(saveideal)
							{
							idealfile = fopen("sourcetrees.ph", "w");
							for(i=0; i<Total_fund_trees; i++)
								{
								strcpy(temptree, fundamentals[i]);
								returntree(temptree);
								fprintf(idealfile, "%s [%s]\n", temptree, tree_names[i]);
								}
							fclose(idealfile);
							}
						 cal_fund_scores(FALSE);				
						}
					
					
					}
				
				
				}			
				
            
				
				
			}
		
		
		if(!error)
			{
		
			if(criterion == 2 || criterion == 3)
				{
				coding(0,0,0);   /** if we are using MRC or QC we need to build the source trees matrix first */
				if(criterion==2) condense_coding();
				}
				
			interval1 = time(NULL);
			printf("\nProgress Indicator:");
			fflush(stdout);
			if(print_all_scores) allscores = fopen("allscores.txt", "w");

			GC = 0;
			for(i=0; i<ntrees; i++)
				{
				if(!user_break)
					{
					interval2 = time(NULL);
					if(difftime(interval2, interval1) > 5) /* every 5 seconds print a dot to the screen */
						{
						printf("=");
						fflush(stdout);
						interval1 = time(NULL);
						}
					/**** generate the tree ****/
					if(random)
						{
						if(gen_method == 1)
							{
							strcpy(rand_tree, "");
							intTotree((int)fmod(rand(), sup)+1, rand_tree, number_of_taxa);
							}
						if(gen_method == 2)
							{
							random_star_decom(rand_tree);
							}
						
						}
					else
						intTotree( i+1, rand_tree, number_of_taxa);
					
					if(tree_top != NULL)
						{
						dismantle_tree(tree_top);
						tree_top = NULL;
						}
					temp_top = NULL;
					
					if(gen_method == 1)tree_build(1, rand_tree, tree_top, FALSE, -1);
					if(gen_method == 2)tree_build(1, rand_tree, tree_top, TRUE, -1);
					tree_top = temp_top;
					temp_top = NULL;
			/**** evaluate its fit to the source trees in memory *****/
					temptree[0] = '\0';
					if(criterion == 0) results[i] = compare_trees(FALSE);  /* calculate the distance from the super tree to all the fundamental trees */
					if(criterion == 2)  /* calculate the distance using the MRC criterion */
						{
						strcpy(temptree, "");
						print_tree(tree_top, temptree);
						strcat(temptree, ";");
						unroottree(temptree);
						results[i] = (float)MRC(temptree);
						}
					if(criterion == 3)
						{
						strcpy(temptree, "");
						print_tree(tree_top, temptree);
						strcat(temptree, ";");
						unroottree(temptree);
						results[i] = (float)quartet_compatibility(temptree);
						}
					/* Display results */
					if(print_all_scores) fprintf(allscores, "%f\n", results[i]);
					GC++;
					}
				}
			draw_histogram(outfile, n, results, GC);
			}
		if(data != 1)
			{
			/* copy back the arrays to their original places */
			
			for(i=0; i<Total_fund_trees; i++)
				{
				strcpy(fundamentals[i], stored_funds[i]);
				for(j=0; j<number_of_taxa; j++)
					{
					presence_of_taxa[i][j] = stored_presence_of_taxa[i][j];
					for(k=0; k<number_of_taxa; k++)
						{
						fund_scores[i][j][k] = stored_fund_scores[i][j][k];
						}
					}
				}
			
			/* free all the memory allocated */
			
			for(i=0; i<Total_fund_trees; i++)
				{
				free(stored_funds[i]);
				free(stored_presence_of_taxa[i]);
				for(j=0; j<number_of_taxa; j++)
					{
						free(stored_fund_scores[i][j]);
					}
				free(stored_fund_scores[i]);
				}
			free(stored_funds);
			free(stored_presence_of_taxa);
			free(stored_fund_scores);



				
				}
			
		free(temptree);
		free(rand_tree);
		free(results);
		fclose(outfile);
		if(print_all_scores) fclose(allscores);
		}
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
		
	if(outfile != NULL)printf("\nResults as follows:\n\n\n");
	for(i=0; i<bins; i++)
		{
		 x = (hist[i]/max1)*40;
		if(outfile != NULL)fprintf(outfile, "%.2f-%.2f\t%d\n", min+(i*k), (min+((i+1)*k))-0.01, hist[i]);
		if(outfile != NULL)printf("%.2f - %.2f\t|", min+(i*k), (min+((i+1)*k))-0.01);
		else printf("\t%-4.0f|", min+(i*k));
		for(j=0; j<x; j++)
			printf("=");
		printf(" (%d)\n",hist[i]);
		y=i;
		z =0;
		while(y < bins && hist[y] == 0)
			{
			y++;
			z++;
			}
		if(z > 15)
			{
			printf("\t    .\n\t    .\n\t    .\n\t    \n");
			i = y-3;
			}
		
		}
		printf("\n\n\n");

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
	
	
	if(outfile != NULL)printf("Moments of the Distribution:\n\n\tMean = %f\n\tVariance = %f\n\tStandard Deviation = %f\n\tSkewness = %f\n\tStandard deviation of skewness = %f\n\n", mean, variance, std_dev, skewness, sqrt(((double)6)/(double)num_results));
	
	
	
	
	
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
				printf("Cannot open file %s\n", guidetreename);
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
							printf("Error: The cut off for a consensus tree must be between .5 and 1.0\n");
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
						printf("Error: %s is an invalid option for data\n", parsed_command[i+1]);
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
				printf("\nConsensus settings:\n");
				printf("\tConsensus of %d universally distributed source trees\n", numtrees);
				if(useguide)
					{
					printf("\tOnly relationships as defined by the guidetree in %s ", guidetreename);
					}
				else
					{
					printf("\tOnly relationships with ");
					if(percentage == 1) printf("100%% support ");
					if(percentage == 0) printf("50%% support or greater (including congruent minor components) ");
					if(percentage != 0 && percentage != 1) printf("%0.0f%% support or greater ", percentage*100);
					}
				printf("are included in the consensus\n");
				printf("\tConsensus file = %s\n", consensusfilename);
				
				consensusfile = fopen(consensusfilename, "w");
				
				consensus(numtrees, temptrees, numtrees, percentage, consensusfile, guidetreefile);
				
				fclose(consensusfile);
				}
			else
				printf("There are no tress that contain all the taxa, unable to construct a consensus tree\n");
					
			}
		}
	
	}

    
/*** consensus is a function to carry out a majority-rule consensus on the results of a bootstrap analysis ****/
void consensus(int num_trees, char **trees, int num_reps, float percentage, FILE *outfile, FILE *guidetreefile)
    {
    int i, j, k, q, r, same1 = FALSE, same2 = FALSE, same3 = FALSE, same4 = FALSE,same5 = FALSE,same6 = FALSE,same7 = FALSE,same8 = FALSE, l, **sets, *in, *tmpcoding = NULL, end = FALSE, found = -1;
	/* The first thing needed is to create a Baum-Ragan coding scheme holding all the information from the bootstrapped trees */
	char **string, *tmp = NULL, name[100], rest[TREE_LENGTH], value[100];
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
					if(k != 0) printf("%d ", shorthand[k-1]);
					}
				shorthand[i][k] += (((int)(pow(2, j%16))) * total_coding[i][j]);
				}
			 printf("%d ", shorthand[i][k]);
			printf("\n");
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

		
		printf("\n\n\nSets included in the consensus tree\n");
		for(i=0; i<num_partitions; i++)
			{
			if(in[i] && (float)(partition_number[i]/num_reps) >= percentage && (float)(partition_number[i]/num_reps) >= .5)
				{
				printf("\t");
				for(j=0; j<number_of_taxa; j++)
					{
					if(total_coding[i][j] == 1)
						printf("*");
					else
						printf(".");
					}
				printf("\t%.2f\n", (float)(partition_number[i]/num_reps));
				}
			}
		
		printf("\n\n\nMinor Components included in the consensus tree\n");
		for(i=0; i<num_partitions; i++)
			{
			if(in[i] && (float)(partition_number[i]/num_reps) < .5)
				{
				printf("\t");
				for(j=0; j<number_of_taxa; j++)
					{
					if(total_coding[i][j] == 1)
						printf("*");
					else
						printf(".");
					}
				printf("\t%.2f\n", (float)(partition_number[i]/num_reps));
				}
			}
		
		printf("\n\nSets not included in the consensus tree\n");
		for(i=0; i<num_partitions; i++)
			{
			if(!in[i])
				{
				printf("\t");
				for(j=0; j<number_of_taxa; j++)
					{
					if(total_coding[i][j] == 1)
						printf("*");
					else
						printf(".");
					}
				printf("\t%.2f\n", (float)(partition_number[i]/num_reps));
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
	free(tmp);
    }



void showtrees(void)
	{
	int worst = -2, best = -2,savetrees = FALSE, found = TRUE, taxachosen = 0, counter = 0, mode[5] = {TRUE, FALSE, FALSE, FALSE, FALSE}, start = 0, end = Total_fund_trees, error = FALSE, i=0, j=0, k=0, l=0, num=0, equalto = -1, greaterthan =0, lessthan = 1000000000, taxa_count = 0;
	char *temptree, string_num[10], namecontains[100], **containstaxa = NULL, savedfile[100], temptree1[TREE_LENGTH], tmp[TREE_LENGTH];
	FILE *showfile = NULL;
	float bestscore =10000000, worstscore = 0, **tempscores = NULL;
	int *tempsourcetreetag = NULL, display = TRUE, best_total = -1, total = 0;
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
	strcpy(savedfile, "showtrees.txt");
	containstaxa = malloc(number_of_taxa*sizeof(char*));
	for(i=0; i<number_of_taxa; i++)
		{
		containstaxa[i] = malloc(1000*sizeof(char));
		containstaxa[i][0] = '\0';
		}
	temptree = malloc(TREE_LENGTH*sizeof(char));
	temptree[0] = '\0';
	temptree1[0] = '\0';
	for(i=0; i<num_commands; i++)
		{
		if(strcmp(parsed_command[i], "display") == 0)
			{
			if(strcmp(parsed_command[i+1], "no") == 0)
				display = FALSE;
			else
				{
				if(strcmp(parsed_command[i+1], "yes") == 0)
					{
					display = TRUE;
					savetrees = TRUE;
					}
				else
					{
					printf("ERROR: %s not valid modifier of \"display\"\n", parsed_command[i+1]);
					error = TRUE;
					}
				}
			}
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
					printf("Error: %s not a valid modifier of \"savetrees\"\n", parsed_command[i+1]);
					error = TRUE;
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
				printf("Error: the start of the range must not be less than 1 and the end of the range must no be greater than the number of source trees\n");
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
				printf("Error: the number of taxa chosen os greater than the number of taxa in memory\n");
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
					printf("Error in size \"equalto\"\n\n");
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
						printf("Error in size \"greaterthan\"\n\n");
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
							printf("Error in size \"lessthan\"\n\n");
							}
						}
					else
						{
						printf("Error: %s not valid option for \"size\"\n", parsed_command[i+1]);
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
			printf("Error: There are no saved supertrees in memory from which to calculate scores\n");
			error = TRUE;
			}
		}
			
	if(savetrees)
		{
		if((showfile = fopen(savedfile, "w")) == NULL)		/* check to see if the file is there */
			{								/* Open the source tree file */
			printf("Cannot open file %s\n", savedfile);
			error = TRUE;
			}
		}
	if(!error)
		{
		printf("showtrees settings:\n\n");
		if(savetrees)
			printf("\tsaving selection of trees in phylip format to file: %s\n", savedfile);
		
		
		
		if(mode[0])
			{
			for(i=0; i<start; i++)
				tempsourcetreetag[i] = FALSE;
			for(i=end+1; i<Total_fund_trees; i++)
				tempsourcetreetag[i] = FALSE;
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
				tree_build(1, temptree, tree_top, 1, -1);

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
				tree_build(1, temptree, tree_top, 1, -1);

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
			        tree_build(1, fundamentals[j], tree_top, 0, j); /* build the tree passed to the function */

			        tree_top = temp_top;
			        temp_top = NULL;
			        reset_tree(tree_top);

			        identify_species_specific_clades(tree_top);  /* Call recursive function to travel down the tree looking for species-specific clades */
			        shrink_tree(tree_top);    /* Shrink the pruned tree by switching off any internal nodes that are not needed */
			        temptree[0] = '\0'; /* initialise the string */
			        if(print_pruned_tree(tree_top, 0, temptree, TRUE) >1)
			            {
			            tmp[0] = '\0';
			            strcpy(tmp, "(");
			            strcat(tmp, temptree);
			            strcat(tmp, ")");
			            strcpy(temptree, tmp);
			            }
			        strcat(temptree, ";");

			        if(strcmp(tree_names[j], "") != 0)
			        	fprintf(showfile, "%s[%s", temptree, tree_names[j]);
			        else
			        	fprintf(showfile, "%s[%d", temptree, j);
			        if(trees_in_memory > 0)
						fprintf(showfile, " %f]\n", sourcetree_scores[i]);
					else
						fprintf(showfile, "]\n");
					}
				if(display)
					{
					temptree[0] = '\0';
					strcpy(temptree, fundamentals[i]);
					returntree(temptree);
					printf("\n\n\nTree number %d\nTree name = %s\n", i+1, tree_names[i]);
					printf("Weight = %f\n", tree_weights[i]);
					if(trees_in_memory > 0)printf("Score = %f\n", sourcetree_scores[i]);
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
					printf("\n\n\nTree number %d\nTree name = %s\n", i+1, tree_names[i]);
					printf("Weight = %f\n", tree_weights[i]);
					if(trees_in_memory > 0)printf("Score = %f\n", sourcetree_scores[i]);
					tree_coordinates(temptree, TRUE, TRUE, FALSE, -1);
					}
				counter++;
				}
			} */ /* Old save trees before names were delimited */
		
		printf("\n%d source trees met with the criteria specified\n", counter);
		}
	if(savetrees) fclose(showfile);
	free(temptree);
	for(i=0; i<Total_fund_trees; i++)
		free(tempscores[i]);
	for(i=0; i<number_of_taxa; i++)
		free(containstaxa[i]);
	free(tempscores);
	free(containstaxa);
	free(tempsourcetreetag);
	}

/* quick sort */
void quick(float ** items, int count)
	{
	qs(items, 0, count-1);
	}

/* the quick sort */
void qs(float **items, int left, int right)
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
	int worst = -2, best = -2,savetrees = FALSE, found = TRUE, taxachosen = 0, counter = 0, mode[10] = {FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE}, start = 0, end = Total_fund_trees, error = FALSE, i=0, j=0, k=0, l=0, num=0, equalto = -1, greaterthan =1000000000, lessthan = , taxa_count = 0;
	char *temptree, string_num[10], namecontains[100], **containstaxa = NULL, savedfile[100], *command = NULL, tmp[TREE_LENGTH];
	FILE *showfile = NULL, *tempfile = NULL;
	float bestscore =10000000, worstscore = 0, **tempscores = NULL;
	int *tempsourcetreetag = NULL, countedout =0, *temp_incidence = NULL;
	
	
	
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
				printf("Error: the start of the range must not be less than 1 and the end of the range must no be greater than the number of source trees\n");
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
				printf("Error: the number of taxa chosen os greater than the number of taxa in memory\n");
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
				printf("Error: the range must end with a larger or equal score to the start of the range\n");
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
					printf("Error in size \"equalto\"\n\n");
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
						printf("Error: %s not valid option for \"size\"\n", parsed_command[i+1]);
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
			printf("Error: There are no saved supertrees in memory from which to calculate scores\n");
			error = TRUE;
			}
		}
	if(mode[4])
		{
		if(trees_in_memory == 0)
			{
			printf("Error: There are no saved supertrees in memory from which to calculate scores\n");
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
				tree_build(1, temptree, tree_top, 1, -1);

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
				tree_build(1, temptree, tree_top, 1, -1);

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

		if(mode[4])
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
		if(mode[5] || mode[6])
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



		if(mode[0] || mode[1] || mode[2] || mode[3] || mode[4] || mode[5] || mode[6])
			{
			for(i=0; i<number_of_taxa; i++) temp_incidence[i] = 0;
			for(i=0; i<Total_fund_trees; i++)
				{
				if(tempsourcetreetag[i])
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
					tree_build(1, temptree, tree_top, 1, -1);

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
			l=1;
			if(l> 0)
				{
				/*printf("\nWarning: %d Taxa are no longer represented in the included source trees\nThese taxa are as follows:\n", l);
				for(i=0; i<number_of_taxa; i++)
					{
					if(temp_incidence[i] == 0)
						printf("\t%s\n", taxa_names[i]);
					}
				command = malloc(10000*sizeof(char));
				printf("This will permenantly remove these trees from memory\nAre you sure you wish to continue: (yes/no) ");
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
						        tree_build(1, fundamentals[j], tree_top, 0, j); /* build the tree passed to the function */

						        tree_top = temp_top;
						        temp_top = NULL;
						        reset_tree(tree_top);

						        identify_species_specific_clades(tree_top);  /* Call recursive function to travel down the tree looking for species-specific clades */
						        shrink_tree(tree_top);    /* Shrink the pruned tree by switching off any internal nodes that are not needed */
						        temptree[0] = '\0'; /* initialise the string */
						        if(print_pruned_tree(tree_top, 0, temptree, TRUE) >1)
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
					printf("\nAction aborted\n");
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
				printf("\n%d source trees were excluded, %d trees remain in memory\n", counter, countedout );	
				}
			}
		}
	
	free(temptree);
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

void returntree(char *temptree)
	{
	char string_num[10], string[TREE_LENGTH];
	int i=0, j=0, k=0, l=0, num;
	
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
				printf("Error: the start of the range must not be less than 1 and the end of the range must no be greater than the number of source trees\n");
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
				printf("Error: the number of taxa chosen os greater than the number of taxa in memory\n");
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
				printf("Error: the range must end with a larger score than the start of the range\n");
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
					printf("Error in size \"equalto\"\n\n");
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
						printf("Error in size \"greaterthan\"\n\n");
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
							printf("Error in size \"lessthan\"\n\n");
							}
						}
					else
						{
						printf("Error: %s not valid option for \"size\"\n", parsed_command[i+1]);
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
			printf("Error: There are no saved supertrees in memory from which to calculate scores\n");
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
				tree_build(1, temptree, tree_top, 1, -1);

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
		if(mode[2])
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
				tree_build(1, temptree, tree_top, 1, -1);

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
		if(mode[4])
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
		printf("\n%d source trees were re-included, %d trees remain in memory\n", counter, countedout );	
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

void sourcetree_dists(void)
	{
	int i=0, j=0, k=0, l=0, m=0, y=0, ***trees_coding = NULL, **included = NULL, *tracking = NULL, count = 0, x, total, *tree1 = NULL, *tree2 = NULL, *t1tag = NULL, *t2tag = NULL, *t1score = NULL, *t2score = NULL, same, same1, same2, same3;
	char number[100], string[TREE_LENGTH], RFfilename[100];
	int p1 = 0, p2 = 0, counter = 0, remaining_taxa = 0, num_falsed = 0, error = FALSE, r = 0, **shared_taxa = NULL, output_format = 0, missing_method = 2, here = TRUE, found = TRUE;
	float **results = NULL;
	FILE *RFfile = NULL;

	RFfilename[0] = '\0';
	strcpy(RFfilename, "robinson-foulds.txt");
	 for(i=0; i<num_commands; i++)
        {
		if(strcmp(parsed_command[i], "filename") == 0)
			{
			strcpy(RFfilename, parsed_command[i+1]);
			}
			
		if(strcmp(parsed_command[i], "output") == 0)
			{
			if(strcmp(parsed_command[i+1], "matrix") == 0)
				output_format = 0;
			else
				{
				if(strcmp(parsed_command[i+1], "vector") == 0)
					output_format = 1;
				else
					{
					printf("Error: '%s' not valid option for output\n", parsed_command[i+1]);
					error = TRUE;
					}
				}
			}
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
					if(strcmp(parsed_command[i+1], "none") == 0)
						missing_method = 2;
					else
						{
						printf("Error: '%s' not valid option for missing\n", parsed_command[i+1]);
						error = TRUE;
						}
					}
				}
			}
		}
	
	if((RFfile = fopen(RFfilename, "w")) == NULL)		/* check to see if the file is there */
			{								/* Open the source tree file */
			printf("Cannot open file %s\n", RFfilename);
			error = TRUE;
			}
	if(!error)
		{
		
		results = malloc((Total_fund_trees-num_excluded_trees)*sizeof(float *));
		if(!results) memory_error(111);
		shared_taxa = malloc((Total_fund_trees-num_excluded_trees)*sizeof(float *));
		if(!shared_taxa) memory_error(112);
		tree1 = malloc(number_of_taxa*sizeof(int));
		if(!tree1) memory_error(113);
		tree2 = malloc(number_of_taxa*sizeof(int));
		if(!tree2) memory_error(114);
		t1tag = malloc(number_of_taxa*sizeof(int));
		if(!t1tag) memory_error(115);
		t2tag = malloc(number_of_taxa*sizeof(int));
		if(!t2tag) memory_error(116);
		t1score = malloc(number_of_taxa*sizeof(int));
		if(!t1score) memory_error(117);
		t2score = malloc(number_of_taxa*sizeof(int));
		if(!t2score) memory_error(118);
		
		trees_coding = malloc((Total_fund_trees - num_excluded_trees)*sizeof(int**));
		if(!trees_coding) memory_error(119);
		for(i=0; i<Total_fund_trees - num_excluded_trees; i++)
			{
			shared_taxa[i] = malloc((Total_fund_trees-num_excluded_trees)*sizeof(float));
			if(!shared_taxa[i]) memory_error(120);
			results[i] = malloc((Total_fund_trees-num_excluded_trees)*sizeof(float));
			if(!results[i]) memory_error(121);
			for(j=0; j<Total_fund_trees - num_excluded_trees; j++)
				{
				results[i][j] = 0;
				shared_taxa[i][j] = 0;
				}
			trees_coding[i] = malloc(number_of_taxa*sizeof(int*));
			if(!trees_coding[i]) memory_error(122);
			for(j=0; j<number_of_taxa; j++)
				{
				trees_coding[i][j] = malloc(number_of_taxa*sizeof(int));
				if(!trees_coding[i][j]) memory_error(123);
				for(k=0; k<number_of_taxa; k++)
					trees_coding[i][j][k] = 0;
				}
			}
		included = malloc((Total_fund_trees - num_excluded_trees)*sizeof(int *));
		if(!included) memory_error(124);
		for(i=0; i<Total_fund_trees - num_excluded_trees; i++)
			{
			included[i] = malloc(number_of_taxa*sizeof(int));
			if(!included[i]) memory_error(125);
			for(j=0; j<number_of_taxa; j++)
				included[i][j] = TRUE;
			}
		tracking = malloc((2*number_of_taxa)*sizeof(int));
			if(!tracking) memory_error(79);
			for(i=0; i<(2*number_of_taxa); i++) tracking[i] = FALSE;

		printf("Calculating the Robinson-Foulds distance between the trees .... \n\n");
		/**** calculate a BR coding for each sourcetree in memory ****/ 
		y=0;
		for(x=0; x<Total_fund_trees; x++)
			{
			unroottree(fundamentals[x]);
			if(sourcetreetag[x])
				{
				count=0;
				i=0;
				/* we scan the nested parenthesis string once from left to right calculating the baum-ragan coding scheme */
				while(fundamentals[x][i] != ';')
					{
					
					switch(fundamentals[x][i])
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
						
						case ':':
							while(fundamentals[x][i] != '(' && fundamentals[x][i] != ')' && fundamentals[x][i] != ',' && fundamentals[x][i] != ';')
								i++;
							break;
			
						default :
							for(j=0; j<100; j++)
								number[j] = '\0';
							j=0;
							while(fundamentals[x][i] != '(' && fundamentals[x][i] != ')' && fundamentals[x][i] != ',' && fundamentals[x][i] != ':')
								{
								number[j] = fundamentals[x][i];
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
									trees_coding[y][j][total] = 1;
									}
								}
							break;
						}
					}

				/*** Put ? at those taxa that are not present **/
				for(i=0; i<number_of_taxa; i++)
					{
					if(!presence_of_taxa[x][i])
						{
						for(j=0; j<number_of_taxa; j++)
							{
							trees_coding[y][j][i] = 3;
							}
						}
					}

				
				
				y++;
				}
			}
		x = y = 0;
		for(i=0; i<Total_fund_trees; i++)
			{
			if(sourcetreetag[i])
				{
				
				y=x+1;
				for(j=i+1; j<Total_fund_trees; j++)
					{
					if(sourcetreetag[j])
						{
						num_falsed = 0;
						for(k=0; k<number_of_taxa; k++)
							{
							if(presence_of_taxa[i][k] > 0) tree1[k] = TRUE;
							else
								{
								tree1[k] = FALSE;
								num_falsed ++;
								}
							}
						for(k=0; k<number_of_taxa; k++)
							{
							t1tag[k] = TRUE;
							t2tag[k] = TRUE;
							t1score[k] = 0;
							t2score[k] = 0;
							if(!presence_of_taxa[j][k])
								{
								tree1[k] = FALSE;
								if(presence_of_taxa[i][k] > 0)num_falsed ++;
								}
							}
						/**** check to make sure that there are at least 4 taxa in common between the two trees *****/
						remaining_taxa = 0;
						for(k=0; k<number_of_taxa; k++)
							{
							
							if(presence_of_taxa[i][k] > 0  && presence_of_taxa[j][k] > 0) remaining_taxa++;
							}
						shared_taxa[x][y] = shared_taxa[y][x] = remaining_taxa;
						if(remaining_taxa >3)
							{
							/** step 1 calculate the binary score of each column **/
							for(k=1; k<number_of_taxa; k++)
								{
								r = 0;
								p1 = p2 = 0;
								for(l=0; l<number_of_taxa; l++)
									{
									if(tree1[l])
										{
										p1 += trees_coding[x][k][l];
										p2 += trees_coding[y][k][l];
										t1score[k] += (trees_coding[x][k][l]);
										
										t2score[k] += (trees_coding[y][k][l]);
										r++;
										}
									}
								if(p1 < 2 || p1+num_falsed > number_of_taxa-2) t1score[k] = -1;
								if(p2 < 2 || p2+num_falsed > number_of_taxa-2) t2score[k] = -1;							
								}
							/** step 2 remove any duplicated columns from each tree **/
							for(k=0; k<number_of_taxa; k++)
								{
								for(l=k+1; l<number_of_taxa; l++)
									{
									same = TRUE; same1 = TRUE; same2 = TRUE, same3 = TRUE;
									for(m=0; m<number_of_taxa; m++)
										{
										if(tree1[m])
											{
											if(trees_coding[x][k][m] != trees_coding[x][l][m]) same = FALSE;
											else same2 = FALSE;
											if(trees_coding[y][k][m] != trees_coding[y][l][m]) same1 = FALSE;
											else same3 = FALSE;
											}
										}
									if(same || same2) t1score[l] = -1;
									if(same1 || same3) t2score[l] = -1;
									}
								}
							
							/****** check to see if the scores of any of the columns are found in the other tree ***/
							p1 = p2 = 0;
							counter = 0;
							for(k=1; k<number_of_taxa; k++)
								{
								if(t1score[k] != -1 && t1score[k] != 0)
									p2++;
								if(t2score[k] != -1 && t2score[k] != 0) counter++;
								if(t1score[k] != -1 && t1score[k] != 0)
									{
									for(l=1; l<number_of_taxa; l++)
										{
										if(t2score[l] != -1 && t2score[l] != 0)
											{
											same = TRUE; same1 = TRUE;
											for(m=0; m<number_of_taxa; m++)
												{
												if(tree1[m])
													{
													if(trees_coding[x][k][m] != trees_coding[y][l][m]) same = FALSE;
													else same1 = FALSE;
													}
												}
											if(same || same1)
												{
											/*	if(t1score[k] == t2score[l])
											*/		p1++;
												}
											}
										}
									}
								}
							
							results[x][y] = (float)((float)((p2+counter)-(2*p1))/(float)((2*(number_of_taxa-num_falsed))-6));
							results[y][x] = results[x][y];
							}
						else
							{
							results[x][y] = -1;
							results[y][x] = -1;
							}
						y++;
						}
					}
				x++;
				}
				
			}
		x=0;
		
		/*** if we want to try to calculate the missing cells ***/
		
		if(missing_method != 2)
			{
			printf("Estimating the missing data using ");
			if(missing_method == 0) printf("ultrametric distances\n");
			else printf("4 point condition distances\n");
			while(here && found)
				{
				for(i=0; i<Total_fund_trees-num_excluded_trees; i++)
					{
					for(j=0; j<Total_fund_trees-num_excluded_trees; j++)
						{
						if(results[i][j] == -1)
							{
							/** Calculate an approximation for the missing cell **/
							if(missing_method == 0) /*** Do the ultrametric calculation ***/
								{
								/** First find a cell shard by i and j that has a distace for both **/
								found = FALSE;
								for(k=0; k<Total_fund_trees-num_excluded_trees; k++)
									{
									if(k != i && k != j && results[i][k] != -1 && results[j][k] != -1 && !found)
										{
										found = TRUE;
										if(results[i][k] > results[j][k])
											results[i][j] = results[j][i] = results[i][k];
										else
											results[i][j] = results[j][i] = results[j][k];
										
										i = j = Total_fund_trees-num_excluded_trees;
										k = l = Total_fund_trees-num_excluded_trees;
										}
									}
								}
							if(missing_method == 1) /*** Do the 4 point condition calculation ***/
								{
								
								/** this method requires that i find two cells that have values for both i and j and for each other **/
								found = FALSE;
								for(k=0; k<Total_fund_trees-num_excluded_trees; k++)
									{
									for(l=0; l<Total_fund_trees-num_excluded_trees; l++)
										{
										if(k != i && k != j && l != i && l != j && l != k && results[i][k] != -1 && results[i][l] != -1 && results[j][k] != -1 && results[j][l] != -1 && results[k][l] != -1 && !found)
											{
											found = TRUE;
											if((results[i][l] + results[j][k]) > (results[i][k] + results[j][l]))
												results[i][j] = results[j][i] = (results[i][l] + results[j][k]) - results[k][l];
											else
												results[i][j] = results[j][i] = (results[i][k] + results[j][l]) - results[k][l];
											
											i = j = Total_fund_trees-num_excluded_trees;
											k = l = Total_fund_trees-num_excluded_trees;
											}
										}
									}
								}
							}
						}
					}
				
		
				here = FALSE;
				
				for(i=0; i< Total_fund_trees-num_excluded_trees; i++)
					{
					for(j=0; j< Total_fund_trees-num_excluded_trees; j++)
						{
						if(results[i][j] == -1) here = TRUE;
						}
					}
				}
			
		
		
			if(!found)
				{
				printf("\n\nERROR: the overlap in the data is too sparse to calculate missing cells using ");
				if(missing_method == 0) printf("an ultrametric estimate\n");
				if(missing_method == 1) printf("a 4 point condition estimate\n");
				error = TRUE;
				}
			}
				
		
		/**** END missing cells ***/
		
		if(output_format == 1)
			{
			fprintf(RFfile, "Tree number\tSize\tTree name\n");
			
			for(i=0; i<Total_fund_trees-num_excluded_trees; i++)
				{
				if(sourcetreetag[i])
					{
					j=0;
					for(k=0; k<number_of_taxa; k++)
						{
						if(presence_of_taxa[x][k] > 0) j++;
						}
					fprintf(RFfile, "%d        \t%d \t%s\n", x, j, tree_names[x]);
					x++;
					}
				}
			fprintf(RFfile, "\n\n\nCompared trees\tDistance\tshared leaves\n");
			}
			
		for(i=0; i<Total_fund_trees-num_excluded_trees; i++)
			{
			for(j=0; j<i; j++)
				{
				if(output_format == 0)
					{
					if(results[i][j] == -1)
						fprintf(RFfile,"NA        \t");
					else
						fprintf(RFfile, "%-10f\t", results[i][j]);
					}
				else
					{
					if(results[i][j] == -1)
						fprintf(RFfile, "%d Vs. %d\t\tNA        \t%d\n", i , j, shared_taxa[i][j]);
					else
						fprintf(RFfile, "%d Vs. %d\t\t%f\t%d\n", i , j, results[i][j], shared_taxa[i][j]);
					
					}
				}
			fprintf(RFfile, "\n");
			}
			
		printf("\nRobinson-Foulds distances of the source trees have been writen to file named %s\n\n", RFfilename);
		
		for(i=0; i<Total_fund_trees-num_excluded_trees; i++)
			{
			free(results[i]);
			for(j=0; j<number_of_taxa; j++)
				free(trees_coding[i][j]);
			free(trees_coding[i]);
			free(included[i]);
			}
		free(results);
		free(trees_coding);
		free(included);
		free(tree1);
		free(tree2);
		free(t1tag);
		free(t2tag);
		free(t1score);
		free(t2score);
		}
	}

void exclude_taxa(int do_all)
	{
	char  *pruned_tree = NULL, tmp[TREE_LENGTH], *command = NULL, tmpfilename[10000], previnputfilename[10000];
	int i=0, j=0, q=0, error = FALSE, taxachosen = 0, found = FALSE, *tobeexcluded = NULL, k=0, l=0, done = FALSE, num_left = 0, min_taxa = 4;
	FILE *tempfile = NULL;
	
	tmp[0] = '\0';
	tmpfilename[0] = '\0';
	previnputfilename[0] = '\0';
	pruned_tree = malloc(TREE_LENGTH*sizeof(char));

	tobeexcluded = malloc(number_of_taxa*sizeof(int));


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
						printf("Error the minimum number of taxa cannot be set below 1\n");
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
						printf("Error: Cannot find any taxa that match the string \"%s\"\n", parsed_command[k]);
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
		printf("Error: you must specify the name of at least one taxa to delete\n");
		error = TRUE;
		}
	else
		printf("\n%d taxa will be permanently deleted from the source trees\n", q);


		/*** go through the names given and identify the taxa ***/
	if(!error)
		{
		strcpy(previnputfilename, inputfilename);
		strcpy(tmpfilename, inputfilename);
		strcat(tmpfilename, ".tempclann.chr");
		tempfile = fopen(tmpfilename, "w");
		command = malloc(10000*sizeof(char));
		command[0] = '\0';
	/*	printf("\nWarning: This command will permenantly delete any chosen taxa from each tree in memory\nThis will also delete any trees that contain less than 4 taxa after the pruning\n\nAre you sure you wish to continue? (yes/no): ");
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
					/***  build the sourcetree in memory ***/
					if(tree_top != NULL)
                        {
                        dismantle_tree(tree_top);
                        tree_top = NULL;
                        }
                    temp_top = NULL;
                    tree_build(1, fundamentals[i], tree_top, FALSE, -1);
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
					if(print_pruned_tree(tree_top, 0, pruned_tree, FALSE) >1)
						{
						tmp[0] = '\0';
						strcpy(tmp, "(");
						strcat(tmp, pruned_tree);
						strcat(tmp, ")");
						strcpy(pruned_tree, tmp);
						}
					strcat(pruned_tree, ";");
					/*** print the pruned sourcetree to the output file ***/
					returntree(pruned_tree);
					strncpy(tmp, pruned_tree, strlen(pruned_tree)-1);
					strcpy(pruned_tree, tmp);
					fprintf(tempfile, "%s", pruned_tree);
					fprintf(tempfile, " [%f]; [ %s]\n", tree_weights[i], tree_names[i]);
					}
				else
					{
					if(!done)printf("\n\nThe following trees have been deleted from memory because they now contain less than 4 taxa\nTree Number\tTree name\n");
					done = TRUE;
					printf("%-10d\t%s\n", i+1, tree_names[i]);
					}
				}  
			}
		}

		
	free(tobeexcluded);
	free(pruned_tree);
	
	if(!error)
		{
		fclose(tempfile);
		if(num_left == 0)
			printf("\nError: This will delete all trees from memory..... aborting deletetaxa\n");
		else
			execute_command(tmpfilename, do_all);
		strcpy(inputfilename, previnputfilename);
		remove(tmpfilename);
		}
	
	}
	
	/* Prune tree: This is a recursive function that is called for every node position of the supertree
	it then checks to see if any of the siblings on this node are not contained in the fundamental tree, these siblings are then turned off.
	This only turns off taxa, pointer siblings will have to be turned off using a separate program */

void prune_taxa_for_exclude(struct taxon * super_pos, int *tobeexcluded)
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

void spr_dist(void)
	{
	float real_score = 0, sprscore = 0, bestreal = 0,  amountspr, previous, totalnow, bestfake = 0, *results = NULL;
	int i=0, j=0, k=0, l=0, x=0, y=0, error = FALSE, diff, *originaldiff = NULL, numbersprs =0, nreps=100, bestnumSPR = 0, minsprs = 0, best = FALSE, now = 0, bestscore = -1, ***scores_original = NULL, **scores_changed = NULL;
	char *pruned_tree = NULL, tmp[TREE_LENGTH], ideal[TREE_LENGTH], userinfile[1000], c, outputfile[1000], inputtree[TREE_LENGTH];
	int starting_super = 0, randomisation = TRUE;
	FILE *outfile = NULL, *infile = NULL;
	
	userinfile[0] = '\0'; outputfile[0] = '\0'; inputtree[0] = '\0';
	strcpy(outputfile, "SPRdistances.txt");
	for(i=0; i<num_commands; i++)
        {
		if(strcmp(parsed_command[i], "supertree") == 0)
			{
			if(strcmp(parsed_command[i+1], "create") == 0)
				starting_super = 0;
			else
				{
				if(strcmp(parsed_command[i+1], "memory") == 0)
					{
					starting_super = 1;
					if(trees_in_memory == 0)
						{
						printf("Error there are no trees in memory\n");
						error = TRUE;
						}
					}
				else
					{
					if((infile = fopen(parsed_command[i+1], "r")) == NULL)
						{
						printf("Error opening file named %s\n", parsed_command[i+1]);
						error = TRUE;
						strcpy(userinfile, parsed_command[i+1]);
						}
					else
						{
						printf("opened input file %s\n", parsed_command[i+1]);
						starting_super = 2;
						}
					}
				}
			}
		if(strcmp(parsed_command[i], "dorandomisation") == 0)
			{
			if(strcmp(parsed_command[i+1], "no") == 0)
				randomisation = FALSE;
			else
				{
				if(strcmp(parsed_command[i+1], "yes") == 0)
					randomisation = TRUE;
				else
					{
					printf("Error: %s not valid option for dorandomisation\n", parsed_command[i+1]);
					error = TRUE;
					}
				}
			}
		if(strcmp(parsed_command[i], "outfile") == 0)
			{
			strcpy(outputfile, parsed_command[i+1]);
			}
			
		}
		
	if((outfile = fopen(outputfile, "w")) == NULL)
		{
		printf("Error opening file named %s\n", outputfile);
		error = TRUE;
		}
	else
		{
		printf("opened output file %s\n", outputfile);
		}
	
	if(!error)
		{
		originaldiff = malloc(Total_fund_trees*sizeof(int));
		scores_original = malloc(Total_fund_trees*sizeof(int**));
		for(i=0; i<Total_fund_trees; i++)
			{
			scores_original[i] = malloc(number_of_taxa*sizeof(int*));
			originaldiff[i] = 0;
			for(j=0; j<number_of_taxa; j++)
				{
				scores_original[i][j] = malloc(number_of_taxa*sizeof(int));
				for(k=0; k<number_of_taxa; k++)
					scores_original[i][j][k] = 0;
				}
			}
		scores_changed = malloc(number_of_taxa*sizeof(int*));
		for(i=0; i<number_of_taxa; i++)
			{
			scores_changed[i] = malloc(number_of_taxa*sizeof(int));
			for(j=0; j<number_of_taxa; j++)
				scores_changed[i][j] = 0;
			}

		
			/** create somewhere to store the source trees during the operation **/
		tmp[0] = '\0';
		ideal[0] = '\0';
		
		stored_fund_scores = malloc(Total_fund_trees*sizeof(int**));
		if(stored_fund_scores == NULL) memory_error(38);
			
		for(i=0; i<Total_fund_trees; i++)
			{
			stored_fund_scores[i] = malloc((number_of_taxa)*sizeof(int*));
			if(stored_fund_scores[i] == NULL) memory_error(39);
				
			else
				{
				for(j=0; j<(number_of_taxa); j++)
					{
					stored_fund_scores[i][j] = malloc((number_of_taxa)*sizeof(int));
					if(stored_fund_scores[i][j] == NULL) memory_error(40);
						
					else
						{
						for(k=0; k<(number_of_taxa); k++)
							stored_fund_scores[i][j][k] = fund_scores[i][j][k];  /* copy the origninal fundamental scores for safe keeping */
						}
					}
				}
			}
		stored_presence_of_taxa = malloc(Total_fund_trees*sizeof(int *));
		if(!stored_presence_of_taxa) memory_error(41);
			
		for(i=0; i<Total_fund_trees; i++)
			{
			stored_presence_of_taxa[i] = malloc((number_of_taxa)*sizeof(int));
			if(!stored_presence_of_taxa[i]) memory_error(42);
				
			for(j=0; j<number_of_taxa; j++)
				{
				stored_presence_of_taxa[i][j] = presence_of_taxa[i][j];
				}
			}


		stored_funds = malloc(Total_fund_trees*sizeof(char *));
		if(!stored_funds) memory_error(31);

		for(i=0; i<Total_fund_trees; i++)
			{
			stored_funds[i] = malloc((tree_length_assignments*TREE_LENGTH)*sizeof(char));
			if(!stored_funds) memory_error(32);

			stored_funds[i][0] = '\0';
			strcpy(stored_funds[i], fundamentals[i]);   /* copy the original fundamentals for safe keeping */
			}


		
		/**** 1. Perform a heuristic search of tree space to identify the score of the best tree **/
		
		if(starting_super == 0) /* do a heuristic search to find the best supertree */
			{
			heuristic_search(FALSE, FALSE, 1000, 1);
			bestreal = scores_retained_supers[0];
			printf("\nbestreal = %f\n", bestreal);
			hsprint = FALSE;
			if(signal(SIGINT, controlc4) == SIG_ERR)
				printf("An error occurred while setting a signal handler\n");
			}
		if(starting_super == 2) /* read in the tree from file */
			{
			c = getc(infile);
			i=0;
			while(c != ';' && !feof(infile))
				{
				tmp[i] = c;
				c = getc(infile);
				i++;
				}
			fclose(infile);
			}
		strcpy(inputtree, tmp);
		/**** 3. Create ideal source trees ***/
		pruned_tree = malloc(TREE_LENGTH*sizeof(char));
		
		/****** We now need to build the Supertree in memory *******/

		/**** 4. Indroduce increasing levels of SPRs per tree until the score of the best tree is the same (within reason) as the score for the best tree from real data ****/
		/*** start with  min % SPR and a max % SPR **/
		/*** By default t hey will be 0.01 for min and 0.8 for max ****/
		/*** Calculate the score of the best tree with 0.01 SPR and if the score is less than the real score **/
		/*** assign  */
		for(y=0; y<2; y++)
			{		
			if(tree_top != NULL)
				{
				dismantle_tree(tree_top);
				tree_top = NULL;
				}
			temp_top = NULL;
			

			if(starting_super != 2) tree_build(1, retained_supers[0], tree_top, TRUE, -1);
			else
				tree_build(1, inputtree, tree_top, TRUE, -1);

			tree_top = temp_top;
			temp_top = NULL;
			strcpy(tmp, "");

			if(!error)
				{
				
				for(i=0; i< Total_fund_trees; i++)
					{
					for(k=0; k<number_of_taxa; k++)
						{
						for(j=0; j<number_of_taxa; j++)
							scores_changed[k][j]  = 0;
						}
					
					if(y== 1)  /* if this is the randomisation part */
						{
						/*randomise_tree(fundamentals[i]); */ /*equiprobable */
						randomise_taxa(fundamentals[i]);  /*markovian*/
						}
					/** calculate the pathmetric for the real tree */
					pathmetric(fundamentals[i], scores_changed);

					prune_tree(tree_top, i);  /* Prune the supertree so that it has the same taxa as the fundamental tree i */
					shrink_tree(tree_top);    /* Shrink the pruned tree by switching off any internal nodes that are not needed */
					for(j=0; j<TREE_LENGTH; j++) pruned_tree[j] = '\0';
					pruned_tree[0] = '\0'; /* initialise the string */
					if(print_pruned_tree(tree_top, 0, pruned_tree, FALSE) >1)
						{
						tmp[0] = '\0';
						strcpy(tmp, "(");
						strcat(tmp, pruned_tree);
						strcat(tmp, ")");
						strcpy(pruned_tree, tmp);
						}
					strcat(pruned_tree, ";");

					strcpy(fundamentals[i], pruned_tree);
					/** calculate the pathmetric for the ideal tree */
					pathmetric(fundamentals[i], scores_original[i]);

					reset_tree(tree_top);  /* reset the supertree for the next comparison */
					originaldiff[i] = 0;
					for(k=0; k<number_of_taxa; k++)
						{
						for(j=k+1; j<number_of_taxa; j++)
							{
							if(scores_original[i][k][j] > scores_changed[k][j])
								originaldiff[i] += (scores_original[i][k][j] - scores_changed[k][j]);
							else
								originaldiff[i] += (scores_changed[k][j] - scores_original[i][k][j]);
							}
						}

					}
				 cal_fund_scores(FALSE);				
				}
				
				
			/**** 5. Continue increasing the level of SPRs per tree until the score of the best tree is the same as the mean of the distribution of best trees from the Yaptp test ****/
			/*** begin by making a dataset with 50% spr and evaluate if this has a score better of lessthan the real data */

			totalnow = 0;
			if(y == 0) fprintf(outfile, "Number of SPRs needed to make Ideal data as bad as the real data\n\n");
			else fprintf(outfile, "Number of SPRs needed to make Ideal data as bad as randomised data\n\n");
			fprintf(outfile, "num SPRs\tTree number\tNum Leaves\tTree name\n");
			for(l=0; l<Total_fund_trees; l++)
				{
				if(sourcetreetag[l])
					{

					/**** For each source tree, take its ideal equivalent ant perform SPRs until its score is comparable to the real version ***/
					/** remember the original ideal version ***/
					strcpy(ideal, fundamentals[l]);
					/*** make 100 attempts to get a score the same (or as close as possible) to the real score ***/
					bestnumSPR = -1;
					bestscore = -1;
					for(x=0; x<100; x++)
						{
						
						previous = 1;
						diff = 0;

						/** copy the version of the changed ideal tree we have now **/
						strcpy(fundamentals[l], ideal);
						numbersprs = 0;
						strcpy(tmp, fundamentals[l]);
						/** continuously change this tree until we reach a score that is worse than the difference between the previous score and the real score ***/
						while(diff < originaldiff[l])
							{
							previous = diff;
							/** perform a random SPR on the tree as long as it is not the first rep, this allows us to assess if the tree is the same as ideal **/

							string_SPR(fundamentals[l]);
							numbersprs++;

							for(i=0; i<number_of_taxa; i++)
								{
								for(j=0; j<number_of_taxa; j++)
									scores_changed[i][j] = 0;
								}
							
							
							
							pathmetric(fundamentals[l], scores_changed);
							strcpy(tmp, fundamentals[l]);
							diff=0;
							for(i=0; i<number_of_taxa; i++)
								{
								for(j=i+1; j<number_of_taxa; j++)
									{
									if(scores_original[l][i][j] > scores_changed[i][j])
										diff+= scores_original[l][i][j] - scores_changed[i][j];
									else
										diff+= scores_changed[i][j] - scores_original[l][i][j];
									}
								}

							} /* end while */
						
						if(bestscore == -1)
							{
							bestscore = abs(diff-originaldiff[l]);
							now = diff;
							bestnumSPR = numbersprs;
							}
						else
							{
							if(abs(diff-originaldiff[l]) < bestscore)
								{
								bestscore = abs(diff-originaldiff[l]);
								now = diff;
								bestnumSPR = numbersprs;
								}
							if(abs(diff-originaldiff[l]) == bestscore)
								{
								if(numbersprs < bestnumSPR)
									bestnumSPR = numbersprs;
								}
							}	
						} /** end for **/
					j=0;
					for(i=0; i<number_of_taxa; i++)
						{
						if(presence_of_taxa[l][i] > 0) j++;
						}
					fprintf(outfile, "%8d\t%11d\t%10d\t%s\n", bestnumSPR, l, j, tree_names[l]);
					totalnow += bestnumSPR;
					bestfake += now;
					strcpy(fundamentals[l], tmp);
					}
				
				}
			
			
			/**** Report the amount of SPR per tree between ideal and real and the the same for real and random ***/
			if(y==0)printf("amount SPR per tree required to make Ideal as incongruent as real data = %f\n", totalnow/(Total_fund_trees-num_excluded_trees));
			else printf("amount SPR per tree required to make Ideal as incongruent as random data = %f\n", totalnow/(Total_fund_trees-num_excluded_trees));
			if(!randomisation) y = 2;
			}
		/**** Finish! ****/
		
		
		fclose(outfile);
			/* copy back the arrays to their original places */
		
		for(i=0; i<Total_fund_trees; i++)
			{
			strcpy(fundamentals[i], stored_funds[i]);
			for(j=0; j<number_of_taxa; j++)
				{
				presence_of_taxa[i][j] = stored_presence_of_taxa[i][j];
				for(k=0; k<number_of_taxa; k++)
					{
					fund_scores[i][j][k] = stored_fund_scores[i][j][k];
					}
				}
			}
		cal_fund_scores(FALSE);
		
			/* free all the memory allocated */
		
		for(i=0; i<Total_fund_trees; i++)
			{
			free(stored_funds[i]);
			free(stored_presence_of_taxa[i]);
			for(j=0; j<number_of_taxa; j++)
				{
					free(stored_fund_scores[i][j]);
				}
			free(stored_fund_scores[i]);
			}
		free(stored_funds);
		free(stored_presence_of_taxa);
		free(stored_fund_scores);
		for(i=0; i<number_of_taxa; i++)
			free(scores_changed[i]);
		
		for(i=0; i<Total_fund_trees; i++)
			{
			for(j=0; j<number_of_taxa; j++)
				free(scores_original[i][j]);
			free(scores_original[i]);
			}
		
		free(scores_original);
		free(scores_changed);
		free(originaldiff);
		}

	}


int string_SPR(char * string)
	{
	int i=0, j=0, k=0, l=0, components = 0, random_num = 0, done = FALSE, found = FALSE, **scores_original = NULL, **scores_changed = NULL, attempts = 0;
	char extracted[TREE_LENGTH], *temptree = NULL, *tmp = NULL, original[TREE_LENGTH], *string1 = NULL;
	
	string1=malloc(TREE_LENGTH*sizeof(char));
	string1[0] = '\0';
	strcpy(string1, string);

	scores_original = malloc(number_of_taxa*sizeof(int*));
	scores_changed = malloc(number_of_taxa*sizeof(int*));
	for(i=0; i<number_of_taxa; i++)
		{
		scores_original[i] = malloc(number_of_taxa*sizeof(int));
		scores_changed[i] = malloc(number_of_taxa*sizeof(int));
		for(j=0; j<number_of_taxa; j++)
			scores_original[i][j] = scores_changed[i][j] = 0;
		}
	
	temptree = malloc(TREE_LENGTH*sizeof(char));
	temptree[0] = '\0';
	tmp = malloc(TREE_LENGTH*sizeof(char));
	while(k == 0)
		{
		unroottree(string1);
		attempts++;
		tmp[0] = '\0';
		extracted[0] = '\0';
		/**** Step 1, count the number of components (internal beanches + taxa) in the tree */
		original[0] = '\0';
		strcpy(original, string1);
		i=1;
		while(string1[i] != ';')
			{
			switch(string1[i])
				{
				case '(':
					if(i!= 0)components++;
					i++;
					break;
					
				case ')':
				case ':':
					i++;
					while(string1[i] != ')' && string1[i] != ':' && string1[i] != ',' && string1[i] != '(' && string1[i] != ';')
						i++;
					break;
					
				case ',':
					i++;
					break;
				
				default:
					components++;
					while(string1[i] != ')' && string1[i] != ':' && string1[i] != ',' && string1[i] != '(' && string1[i] != ';')
						i++;
					break;
				}
			}
		/** now randomly choose a number between 1 and the number of components */
		 random_num = (int)fmod(rand(), components)+1;
		 /**** find the component with that number **/
		 i=1; components = 0;
		while(string1[i] != ';' && !found)
			{
			switch(string1[i])
				{
				case '(':
					if(i!= 0)
						components++;
					if(random_num == components) found = TRUE;
					else i++;
					break;
					
				case ')':
				case ':':
					i++;
					while(string1[i] != ')' && string1[i] != ':' && string1[i] != ',' && string1[i] != '(' && string1[i] != ';')
						i++;
					break;
					
				case ',':
					i++;
					break;
				
				default:
					components++;
					if(random_num == components) found = TRUE;
					else
						{
						while(string1[i] != ')' && string1[i] != ':' && string1[i] != ',' && string1[i] != '(' && string1[i] != ';')
							i++;
						}
					break;
				}
			}
		 /*** Now that we have the start of the component, extract it entirely */
		j=0; k=0;
		 if(string1[i] == '(')
			{
			if(string1[i-1] == ',') 
				{
				string1[i-1] = '~';
				done = TRUE;
				}
			extracted[k] = string1[i];
			string1[i] = '~';
			i++; k++;
			while(string1[i] != ')' || j != 0 )
				{
				if(string1[i] == '(') j++;
				if(string1[i] == ')') j--;
				extracted[k] = string1[i];
				string1[i] = '~';
				i++; k++;
				}
			extracted[k] = string1[i];
			k++;
			string1[i] = '~';
			if(!done) string1[i+1] = '~';
			}
		else
			{
			if(string1[i-1] == ',')
				{
				string1[i-1] = '~';
				done = TRUE;
				}
			while(string1[i] != ')' && string1[i] != ':' && string1[i] != ',' && string1[i] != '(' && string1[i] != ';')
				{
				extracted[k] = string1[i];
				string1[i] = '~';
				k++; i++;
				}
			if(!done) string1[i] = '~';
			}
		
		extracted[k] = '\0';
		/*** now delete the '~' characters from the original string */
		k=0; l=0; i=0; j=0;
		
		while(string1[i] != ';' && string1[i] != '~')
			i++;
		j=i;
		while(string1[j] == '~')
			j++;
		
		while(string1[j] != ';')
			{
			string1[i] = string1[j];
			i++;j++;
			}
		string1[i] = ';';
		string1[i+1] = '\0';
		
		/*** delete unwanted parts of the original tree ***/
		
		/*** build the post-extraction tree in memory ****/
		if(tree_top != NULL)
			{
			dismantle_tree(tree_top);
			tree_top = NULL;
			}
		temp_top = NULL;
		
		tree_build(1, string1, tree_top, FALSE, -1);

		tree_top = temp_top;
		temp_top = NULL;

		/**** shrink the tree ***/
		shrink_tree(tree_top);
		temptree[0] = '\0';
		/** printf the result **/
		if(print_pruned_tree(tree_top, 0, temptree, FALSE) >1)
			{
			tmp[0] = '\0';
			strcpy(tmp, "(");
			strcat(tmp, temptree);
			strcat(tmp, ")");
			strcpy(temptree, tmp);
			}
		strcat(temptree, ";");
		unroottree(temptree);
		strcpy(string1, temptree);

		/** choose a position to put the extracted part back into */
			/****  count the number of components (internal beanches + taxa) in the tree */
		components = 0; i=1;
		while(string1[i] != ';')
			{
			switch(string1[i])
				{
				case '(':
					if(i!= 0)components++;
					i++;
					break;
					
				case ')':
				case ':':
					i++;
					while(string1[i] != ')' && string1[i] != ':' && string1[i] != ',' && string1[i] != '(' && string1[i] != ';')
						i++;
					break;
					
				case ',':
					i++;
					break;
				
				default:
					components++;
					while(string1[i] != ')' && string1[i] != ':' && string1[i] != ',' && string1[i] != '(' && string1[i] != ';')
						i++;
					break;
				}
			}
		/** now randomly choose a number between 1 and the number of components */
		 random_num = (int)fmod(rand(), components)+1;
		 /**** find the component with that number **/
		 i=1; components = 0; found = FALSE;
		while(string1[i] != ';' && !found)
			{
			switch(string1[i])
				{
				case '(':
					if(i!= 0)
						components++;
					if(random_num == components) found = TRUE;
					else i++;
					break;
					
				case ')':
				case ':':
					i++;
					while(string1[i] != ')' && string1[i] != ':' && string1[i] != ',' && string1[i] != '(' && string1[i] != ';')
						i++;
					break;
					
				case ',':
					i++;
					break;
				
				default:
					components++;
					if(random_num == components) found = TRUE;
					else
						{
						while(string1[i] != ')' && string1[i] != ':' && string1[i] != ',' && string1[i] != '(' && string1[i] != ';')
							i++;
						}
					break;
				}
			}
		/** i is now pointing at the position in which we need to place the extracted subtree back into ***/
		/** We will do this by making the extracted part a sister group of the component */
		
		if(string1[i] == '(')
			{
			/*** find the end of this component **/
			k=i; l=0;
			k++;
			while(string1[k] != ')' || l != 0)
				{
				if(string1[k] == '(') l++;
				if(string1[k] == ')') l--;
				k++;
				}
			}
		else
			{
			k=i+1;
			while(string1[k] != '(' && string1[k] != ')' && string1[k] != ',' && string1[k] != ';')
				k++;
			k--;
			}
		
	
		/* k now has the position of the end of this component */
		/* put the close bracket in after the component */
		l = strlen(string1);
		for(j=strlen(string1)+1; j>k; j--)
			string1[j] = string1[j-1];
		string1[k+1] = ')';
		string1[l+1] = '\0';
		
		/* now put the extracted part before the component */
		strcpy(tmp, "(");
		strcat(tmp, extracted);
		strcat(tmp, ",");
		strcpy(extracted, tmp);
		l = strlen(extracted)+strlen(string1);
		for(j=l; j>=i; j--) {
			string1[j] = string1[j-strlen(extracted)];
			}
		for(j=0; j<strlen(extracted); j++)
			string1[i+j] = extracted[j];
			
		string1[l+1] = '\0';
		
		
		strcpy(string, string1);
		/*** lastly check to make sure that the new tree is actually different to the original tree */
		pathmetric(original, scores_original);
		pathmetric(string, scores_changed);
		
		k=0;
		for(i=0; i<number_of_taxa; i++)
			{
			for(j=i+1; j<number_of_taxa; j++)
				{
				if(scores_original[i][j] > scores_changed[i][j])
					k+= scores_original[i][j] - scores_changed[i][j];
				else
					k+= scores_changed[i][j] - scores_original[i][j];
				}
			}
		}
	
	for(i=0; i<number_of_taxa; i++)
		{
		free(scores_original[i]);
		free(scores_changed[i]);
		}
	free(scores_original);
	free(scores_changed);
	free(temptree);
	free(tmp);
	free(string1);
	
	
	return(attempts);
	}





void exhaustive_SPR(char * string)
	{
	int i, j, k, l, x, y, q, r=-1, labelonly = FALSE, exnum, components =0, pruned_components =0, num, *component_index = NULL, *pruned_component_index = NULL, **scores_original = NULL, **scores_changed = NULL;
	char labeledtree[TREE_LENGTH], taxaname[100], tmp_labeledtree[TREE_LENGTH], filename[100], filename1[100], pasted_name[10000], cut_name[100], extractedpart[TREE_LENGTH], tmp_tree[TREE_LENGTH], pruned_tree[TREE_LENGTH], tmp[TREE_LENGTH];
	FILE *sproutfile = NULL, *sprdescriptor = NULL, *labeledtreefile = NULL;
	
	 for(i=0; i<num_commands; i++)
        {
        if(strcmp(parsed_command[i], "labelonly") == 0)
            labelonly = TRUE;
		}
	component_index = malloc(2*number_of_taxa*sizeof(int));
	for(i=0; i<2*number_of_taxa; i++) component_index[i] = 0;   /* This records the index position of the beginning of each of the components... so we don't have to look for them later */
	pruned_component_index = malloc(2*number_of_taxa*sizeof(int));
	for(i=0; i<2*number_of_taxa; i++) pruned_component_index[i] = 0;   /* This records the index position of the beginning of each of the components... so we don't have to look for them later */

	scores_original = malloc(number_of_taxa*sizeof(int*));
	scores_changed = malloc(number_of_taxa*sizeof(int*));
	for(i=0; i<number_of_taxa; i++)
		{
		scores_original[i] = malloc(number_of_taxa*sizeof(int));
		scores_changed[i] = malloc(number_of_taxa*sizeof(int));
		for(j=0; j<number_of_taxa; j++)
			scores_original[i][j] = scores_changed[i][j] = 0;
		}
	
	/*** first count the number of components in the tree and label them */
	taxaname[0] ='\0';
	pasted_name[0] = '\0';
	cut_name[0] = '\0';
	pruned_tree[0] = '\0';
	tmp[0] = '\0';
	tmp_labeledtree[0] = '\0';
	tmp_tree[0] = '\0';
	labeledtree[0] = '\0';
	filename[0] = '\0';
	filename1[0] = '\0';
	extractedpart[0] = '\0';
	unroottree(string);
	strcpy(labeledtree, string);
	i=1; components=0;
	/* NB this assumes that there are no bootstrap supports on the tree, otherwise it will crash */
	/** Add labels to the tree if there are none */
	while(labeledtree[i] != ';')
		{
		switch(labeledtree[i])
			{
			case '(':
				if(i!= 0)
					{
					components++;
					component_index[components] = i;
					/** find the end of this part and add a label **/
					y=i; l=0;
					y++;
					strcpy(tmp_labeledtree, "");
					while(labeledtree[y] != ')' || l != 0)
						{
						if(labeledtree[y] == '(') l++;
						if(labeledtree[y] == ')') l--;
						y++;
						}
					y++;
					if(labeledtree[y] == ':' || labeledtree[y] == ')' || labeledtree[y] == '(' || labeledtree[y] == ',' || labeledtree[y] == ';')
						{
						x = y;
						q=0;
						while(labeledtree[x] != '\0')
							{
							tmp_labeledtree[q] = labeledtree[x];
							x++;
							q++;
							}
						tmp_labeledtree[q] = '\0';
						labeledtree[y] = '\0';
						strcpy(taxaname, "");
						sprintf(taxaname, "I%d", components);
						strcat(labeledtree, taxaname);
						strcat(labeledtree, tmp_labeledtree);
						}
					}
				i++;
				break;
				
			case ')':
			case ':':
				i++;
				while(labeledtree[i] != ')' && labeledtree[i] != ':' && labeledtree[i] != ',' && labeledtree[i] != '(' && labeledtree[i] != ';')
					i++;
				break;
				
			case ',':
				i++;
				break;
			
			default:
				components++;
				component_index[components] = i;
				while(labeledtree[i] != ')' && labeledtree[i] != ':' && labeledtree[i] != ',' && labeledtree[i] != '(' && labeledtree[i] != ';')
					{
					i++;
					}
				break;
			}
		}
	strcpy(tmp_labeledtree, labeledtree);
	returntree(tmp_labeledtree);
	labeledtreefile = fopen("labelledtree.ph", "w");
	fprintf(labeledtreefile, "%s\n", tmp_labeledtree);
	fclose(labeledtreefile);
	printf("written labelledtree to file\n");
	
	if(!labelonly)
		{
		
		/* now for each component, extract it (along with the label) */
		for(num=1; num<=components; num++)
			{
			strcpy(tmp_tree, labeledtree);
			/* sort out files for this component */
			
			strcpy(filename, "");
			sprintf(filename, "sprpart%d.ph", num);
			if(sproutfile != NULL) fclose(sproutfile);
			sproutfile = fopen(filename, "w");
			
			strcpy(filename1, "");
			sprintf(filename1, "sprdescr%d.txt", num);
			if(sprdescriptor != NULL) fclose(sprdescriptor);
			sprdescriptor = fopen(filename1, "w");
			
			
			/* find the num^th component */
			i = component_index[num];

			j=0; x=FALSE;
			/* extract it (along with the label) */
			if(tmp_tree[i] == '(')  /* it's an internal branch */
				{
				k=0; l=0;
				do
					{
					if(tmp_tree[i] == '(') k++;
					if(tmp_tree[i] == ')') k--;
					extractedpart[l] = tmp_tree[i];
					tmp_tree[i] = '~';
					l++; i++;
					} while(k != 0);
				while(tmp_tree[i] != '(' && tmp_tree[i] != ')' && tmp_tree[i] != ','  && tmp_tree[i] != ';' )  /* extract the label too */
					{
					if(tmp_tree[i] == ':')
						{
						x = TRUE;
						cut_name[j] = '\0';
						}
					if(!x) cut_name[j] = tmp_tree[i];
					extractedpart[l] = tmp_tree[i];
					tmp_tree[i] = '~';
					l++; i++; j++;
					}
				fprintf(sprdescriptor,"Name of cut branch: %s\nNames of pasted positions:\n", cut_name);
				}
			else /* it's a taxa */
				{
				l=0;
				while(tmp_tree[i] != '(' && tmp_tree[i] != ')' && tmp_tree[i] != ',' &&  tmp_tree[i] != ';' )
					{
					if(tmp_tree[i] == ':')
						{
						x = TRUE;
						cut_name[j] = '\0';
						}
					if(!x) cut_name[j] = tmp_tree[i];
					extractedpart[l] = tmp_tree[i];
					tmp_tree[i] = '~';
					l++; i++; j++;
					}
				fprintf(sprdescriptor,"Name of cut branch: %s\nNames of pasted positions:\n", taxa_names[atoi(cut_name)]);
				}
			
			extractedpart[l] = '\0';
			/* Delete the '~'s  (i is just after the end of the component)*/
			j = component_index[num];
			
			
			while(tmp_tree[i] != ';')
				{
				tmp_tree[j] = tmp_tree[i];
				i++; j++;
				}
			tmp_tree[j] = ';';
			tmp_tree[j+1] = '\0';
			
			/*** build the post-extraction tree in memory ****/
			if(tree_top != NULL)
				{
				dismantle_tree(tree_top);
				tree_top = NULL;
				}
			temp_top = NULL;
			tree_build(1, tmp_tree, tree_top, FALSE, -1);

			tree_top = temp_top;
			temp_top = NULL;
			/**** shrink the tree ***/
			shrink_tree(tree_top);
			pruned_tree[0] = '\0';
			
			/** print the result **/
			if(print_pruned_tree(tree_top, 0, pruned_tree, FALSE) >1)
				{
				tmp[0] = '\0';
				strcpy(tmp, "(");
				strcat(tmp, pruned_tree);
				strcat(tmp, ")");
				strcpy(pruned_tree, tmp);
				}
			strcat(pruned_tree, ";");
			unroottree(pruned_tree);
			pruned_components = 0;
			for(i=0; i<2*number_of_taxa; i++) pruned_component_index[i] = 0;
			i=0;
			/* count the number of components in the pruned tree */
			while(pruned_tree[i] != ';')
			{
			switch(pruned_tree[i])
				{
				case '(':
					if(i!= 0)
						{
						pruned_components++;
						pruned_component_index[pruned_components] = i;
						}
					i++;
					break;
					
				case ')':
				case ':':
					i++;
					while(pruned_tree[i] != ')'  && pruned_tree[i] != ',' && pruned_tree[i] != '(' && pruned_tree[i] != ';')
						i++;
					break;
					
				case ',':
					i++;
					break;
				
				default:
					pruned_components++;
					pruned_component_index[pruned_components] = i;
					while(pruned_tree[i] != ')' && pruned_tree[i] != ',' && pruned_tree[i] != '(' && pruned_tree[i] != ';')
						{
						i++;
						}
					break;
				}
			}

			tmp_labeledtree[0] = '\0';
			strcpy(tmp_labeledtree, pruned_tree);
			/** choose a position to put the extracted part back into */
			tmp_tree[0] = '\0';
			strcpy(tmp_tree, "(");
			strcat(tmp_tree, extractedpart);
			strcat(tmp_tree, ",");
			strcpy(extractedpart, tmp_tree);
			r=-1;
			for(exnum=1; exnum<=pruned_components; exnum++)
				{
				strcpy(pruned_tree, tmp_labeledtree);
				i = pruned_component_index[exnum];  /* i now points to the start of the component in question */
				
				
				/** i is now pointing at the position in which we need to place the extracted subtree back into ***/
				/** We will do this by making the extracted part a sister group of the component */
				x= FALSE;
				if(pruned_tree[i] == '(')
					{
					/*** find the end of this component **/
					k=i; l=0;
					k++;
					while(pruned_tree[k] != ')' || l != 0)
						{
						if(pruned_tree[k] == '(') l++;
						if(pruned_tree[k] == ')') l--;
						k++;
						}
					k++; l=0;
					while(pruned_tree[k] != '(' && pruned_tree[k] != ')' && pruned_tree[k] != ',' && pruned_tree[k] != ';') /** go to the end of the label on the branch */
						{
						if(pruned_tree[k] == ':')
							{
							x = TRUE;
							pasted_name[l] = '\0';
							}
						if(!x) pasted_name[l] = pruned_tree[k];
						k++; l++;
						}
					k--;
					fprintf(sprdescriptor, "%s\n", pasted_name);
					}
				else
					{
					k=i; l=0;
					strcpy(pasted_name, "");
					while(pruned_tree[k] != '(' && pruned_tree[k] != ')' && pruned_tree[k] != ',' && pruned_tree[k] != ';')
						{
						if(pruned_tree[k] == ':')
							{
							x = TRUE;
							pasted_name[l] = '\0';
							}
						if(!x) pasted_name[l] = pruned_tree[k];
						k++; l++;
						}
					k--;
					fprintf(sprdescriptor, "%s\n", taxa_names[atoi(pasted_name)]);
					}
				/* k now has the position of the end of this component */
				/* put the close bracket in after the component */
				l = strlen(pruned_tree);
				for(j=strlen(pruned_tree)+1; j>k; j--)
					pruned_tree[j] = pruned_tree[j-1];
				pruned_tree[k+1] = ')';
				pruned_tree[l+1] = '\0';
				/* now put the extracted part before the component */
				l = strlen(extractedpart)+strlen(pruned_tree);
				for(j=l; j>=i; j--)
					pruned_tree[j] = pruned_tree[j-strlen(extractedpart)];
				
				for(j=0; j<strlen(extractedpart); j++)
					pruned_tree[i+j] = extractedpart[j];
				pruned_tree[l+1] = '\0';
				/*** lastly check to make sure that the new tree is actually different to the original tree */
				pathmetric(labeledtree, scores_original);
				pathmetric(pruned_tree, scores_changed);
				k=0;
				for(i=0; i<number_of_taxa; i++)
					{
					for(j=i+1; j<number_of_taxa; j++)
						{
						if(scores_original[i][j] > scores_changed[i][j])
							k+= scores_original[i][j] - scores_changed[i][j];
						else
							k+= scores_changed[i][j] - scores_original[i][j];
						}
					}
				returntree(pruned_tree);
				fprintf(sproutfile, "%s\n", pruned_tree);
				if(k==0) r=exnum;
				}
			
			fprintf(sprdescriptor, "tree num %d is the original\n", r);
			
			}
		
		fclose(sproutfile);	
		fclose(sprdescriptor);
		}
	for(i=0; i<number_of_taxa; i++)
		{
		free(scores_original[i]);
		free(scores_changed[i]);
		}
	free(scores_original);
	free(scores_changed);
	free(pruned_component_index);
	free(component_index);


	}	

void neighbor_joining(int brlens, char *tree, int names)
	{
	int num_nodes = number_of_taxa, *deleted = NULL, i, j, smallest_i = -1, smallest_j = -1;
	float *transformed = NULL, smallest = 0, vi, vj, vtmp;
	char **tree_structure = NULL, string[TREE_LENGTH], tmp[TREE_LENGTH], c;
	

	string[0] = '\0'; tmp[0] = '\0';
	/* Assume that the distances array is number_of_taxa * number_of_taxa in size **/
	deleted = malloc(number_of_taxa*sizeof(int));
	if(!deleted) memory_error(99);
	
	transformed = malloc(number_of_taxa*sizeof(float));
	if(!transformed) memory_error(100);
	for(i=0; i<number_of_taxa; i++)
		{
		deleted[i] = FALSE;
		transformed[i] = 0;
		}
	
	tree_structure = malloc(number_of_taxa*sizeof(char *));
	if(!tree_structure) memory_error(101);
	for(i=0; i<number_of_taxa; i++)
		{
		tree_structure[i] = malloc(TREE_LENGTH*sizeof(char));
		if(!tree_structure[i]) memory_error(102);
		tree_structure[i][0] = '\0';
		if(names)strcpy(tree_structure[i], taxa_names[i]);
		else sprintf(tree_structure[i], "%d", i);
		}
	
	while(num_nodes > 2)
		{
		
		/** calculate Ui for each taxa **/
		for(i=0; i<number_of_taxa; i++)
			{
			for(j=0; j<number_of_taxa; j++)
				{
				if(i != j && !deleted[i] && !deleted[j])
					transformed[i] += weighted_scores[i][j];
				}
			transformed[i] = transformed[i]/(num_nodes - 2);
			}
		
		/* find the smallest Dij - Ui - Uj **/
		smallest_i = -1;
		for(i=0; i<number_of_taxa; i++)
			{
			if(!deleted[i])
				{
				for(j=i+1; j<number_of_taxa; j++)
					{
					if(!deleted[j])
						{
						if(smallest_i == -1)
							{
							smallest_i = i; smallest_j = j; smallest = (weighted_scores[i][j] - transformed[i] - transformed[j]);
							}
						else
							{
							if((weighted_scores[i][j] - transformed[i] - transformed[j]) < smallest)
								{
								smallest_i = i; smallest_j = j; smallest = (weighted_scores[i][j] - transformed[i] - transformed[j]);
								}
							}
						}
					}
				}
			}
		/** join i and j into a new node **/
		strcpy(string, "(");
		strcat(string, tree_structure[smallest_i]);
		if(brlens)
			{
			vi = (weighted_scores[smallest_i][smallest_j]*0.5)+((transformed[smallest_i] - transformed[smallest_j])*0.5);
			sprintf(tmp, ":%f", vi);
			strcat(string, tmp);
			}
			
		strcat(string, ",");
		strcat(string, tree_structure[smallest_j]);
		if(brlens)
			{
			vj = (weighted_scores[smallest_i][smallest_j]*0.5)+((transformed[smallest_j] - transformed[smallest_i])*0.5);
			sprintf(tmp, ":%f", vj);
			strcat(string, tmp);
			}
		strcat(string, ")");
		strcpy(tree_structure[smallest_i], string);
		strcpy(tree_structure[smallest_j], "");
		deleted[smallest_j] = TRUE;
		
		/*** Compute the distance between the new node and each of the remaining tips **/
		
		for(i=0; i<number_of_taxa; i++)
			{
			vtmp = weighted_scores[smallest_i][i];
			if(!deleted[i] && i != smallest_i)
				{
				weighted_scores[smallest_i][i] = weighted_scores[i][smallest_i] = ((vtmp + weighted_scores[smallest_j][i]) - weighted_scores[smallest_i][smallest_j])/2;
				}
			}
		num_nodes--;
		
		for(i=0; i<number_of_taxa; i++)
			{
			transformed[i] = 0;
			}
	
		}
	/** now join the last two nodes **/
	/** find i and j that are left **/
	smallest_i = -1; smallest_j = -1;
	for(i=0; i<number_of_taxa; i++)
		{
		if(!deleted[i] && smallest_i == -1) smallest_i = i;
		else
			{
			if(!deleted[i]) smallest_j = i;
			}
		}
	strcpy(string, "(");
	strcat(string, tree_structure[smallest_i]);
	if(brlens)
		{
		vi = (weighted_scores[smallest_i][smallest_j]);
		sprintf(tmp, ":%f", vi);
		strcat(string, tmp);
		}
		
	strcat(string, ",");
	strcat(string, tree_structure[smallest_j]);
	strcat(string, ");");
	strcpy(tree_structure[smallest_i], string);
	strcpy(tree_structure[smallest_j], "");
	deleted[smallest_j] = TRUE;
	
	/** copy the resulting tree into "tree" **/
	strcpy(tree, tree_structure[smallest_i]);
	/** clean up the memory **/
	
	for(i=0; i<number_of_taxa; i++)
		free(tree_structure[i]);
	free(tree_structure);
	free(transformed);
	free(deleted);
	}
	
void nj(void)
	{
	int i, j, missing_method = 1, error = FALSE;
	char *tree = NULL, useroutfile[100], *fakefilename = NULL;
	FILE *outfile = NULL;
	
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
					printf("Error: '%s' is an invalid method to be assigned to missing\n", parsed_command[i+1]);
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
			printf("Error opening file named %s\n", useroutfile);
			error = TRUE;
			}
		}

	if(!error)
		{
		fakefilename = malloc(100*sizeof(char));
		fakefilename[0] = '\0';
		tree = malloc(TREE_LENGTH*sizeof(char));
		printf("\n\nNeighbor-joining settings:\n\tDistance matrix generation by average consensus method\n\tEstimation of missing data using ");
		if(missing_method == 1)
			printf("4 point condition distances\n");
		if(missing_method == 0)
			printf("ultrametric distances\n");
		printf("\tresulting tree saved to file %s\n\n\n", useroutfile);
		
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
		fclose(outfile);
		free(tree);
		}
	
	}
	
	
void reroot_tree(struct taxon *outgroup)
	{
	struct taxon *newbie = NULL, *temp_treetop = NULL, *position = NULL, *temp = NULL, *start = NULL, *parent_start = NULL, *next = NULL, *parent = NULL, *pointer = NULL;
	char *temptree = NULL;
	temptree = malloc(TREE_LENGTH*sizeof(char));
	temptree[0] = '\0';
	newbie = make_taxon();
	temp_treetop = newbie;
	/* every sibling must become a descendent and every descentdent as sibling of that level */
	/*** Assign newbie to point to any siblings of outgroup as a daughter **/
	position = outgroup;
	if(outgroup->prev_sibling != NULL)
		{
		while(position->prev_sibling != NULL) position = position->prev_sibling;
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
	
	parent = start->parent;
	
	/** create new pointer on first level */
	pointer = make_taxon();
	position = start;
	while(position->next_sibling != NULL) position = position->next_sibling;
	position->next_sibling = pointer;
	pointer->prev_sibling = position;
	
	while(parent != NULL)
		{
		position = parent;
		while(position->prev_sibling != NULL) position = position->prev_sibling;
		parent_start = position;
		pointer->daughter = parent_start;
		next = parent_start->parent;
		parent_start->parent = pointer;
		pointer = parent;
		start = parent_start;
		parent = next;
		
		}
	pointer->daughter = NULL;
	
	(newbie->daughter)->parent = newbie;
	
	/* now make outgroup a siblig to newbie */
	newbie->next_sibling = outgroup;
	outgroup->prev_sibling = newbie;
	
	temp_top = temp_treetop;
	/**** Now we need to remove any pointer taxa that are nolonger pointing to anything **/
	clean_pointer_taxa(temp_top);
	/*** Finally we need to put an extra node at the top */
	newbie = make_taxon();
	newbie->daughter = temp_top;
	temp_top->parent = newbie;
	temp_top = newbie;
	free(temptree);
	}
	
void clean_pointer_taxa(struct taxon *position)
	{
	
	struct taxon *start = position, *tmp = NULL;
	int i=0;
	
	while(position != NULL)
		{
		if(position->daughter != NULL) clean_pointer_taxa(position->daughter);
		position = position->next_sibling;
		}
	position = start;
	
	/*** Check to see if there are any pointer taxa not going anywhere ***/
	
	while(position != NULL)
		{
		if(position->name == -1 && position->daughter == NULL)
			{
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
			if(position == tree_top)
				 tree_top = position->next_sibling;
			if(position == start) start = position->next_sibling;
			free(position);
			position = tmp;
			}
		else
			if(position != NULL) position = position->next_sibling;
		}
	position = start;
	
	/*** Count the number of children of each pointer taxa, if any only have 1 then delete that pointer taxa and replace it with its daughter ***/
	
	while(position != NULL)
		{
		if(position->daughter != NULL)
			{
			i=0;
			tmp = position->daughter;
			while(tmp != NULL)
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
				if(position == tree_top) tree_top = position->daughter;
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
	
	
	while(position != NULL && answer == NULL)
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
		while(position != NULL && answer == NULL)
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

void isittagged(struct taxon * position)
	{
	while(position != NULL)
		{
		if(position->tag2 == TRUE) printf("found\n");
		if(position->daughter != NULL)
			{
			isittagged(position->daughter);
			}
		position = position->next_sibling;
		}
	}

void duplicate_tree(struct taxon * orig_pos, struct taxon * prev_dup_pos)
	{
	struct taxon *sibling = NULL, *dup_pos = NULL;
	int done = FALSE, i;
	
	if(orig_pos->parent == NULL) done = TRUE;
	
	while(orig_pos != NULL)
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
			printf("3\n");

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


void printnamesandtags(struct taxon *position)
	{
	while(position != NULL)
		{
		if(position->daughter != NULL) 
			{
			printf(" ->%d", position->tag);
			printnamesandtags(position->daughter);
			printf(":");
			}
		else printf(" %d", position->name);
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



void print_tree_labels(struct taxon *position, int **results, int treenum, struct taxon *species_tree)
	{
	int onetoone;
	while(position != NULL)
		{
		if(position->loss >= 1) results[1][position->tag]++;
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


			
void label_gene_tree(struct taxon * gene_position, struct taxon * species_top, int *presence, int xnum)
	{
	struct taxon * position = gene_position, *tmp = NULL;
	int i =0, j=0, latest = -1;
/*	printf("in Label_gene_tree\n");*/
	while(position != NULL)
		{
		if(position->daughter != NULL) 
			{
	/*		printf("daughter\n"); */
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
	for(i=0; i<number_of_taxa; i++) presence[i] = FALSE;
	
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
			while(position->tag != spec_pos->tag)
					position = position->next_sibling;
					
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

void find(struct taxon * position)
	{
	while(position != NULL)
		{
		if(position->daughter != NULL) find(position->daughter);
		if(position->tag2 == TRUE)
			printf("!%d\n", position->tag);
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


void mapunknowns()
	{
	struct taxon *position = NULL, *species_tree = NULL, *gene_tree = NULL, *best_mapping = NULL, *unknown_fund = NULL, *pos = NULL,*copy = NULL;
	int i, j, k, l,  *presence = NULL, basescore = 1;
	float *overall_placements = NULL, biggest = -1, total, best_total = -1;
	char *temptree, temptree1[TREE_LENGTH];
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
	tree_build(1, temptree, species_tree, 1, -1);
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
		tree_build(1, temptree, gene_tree, 1, -1);
		gene_tree = temp_top;
		temp_top = NULL;
		strcpy(temptree1, temptree);
		
		k=0;
		while(presence_of_taxa[0][k] >0 || presence_of_taxa[l][k] == FALSE) k++;  /* fundamental [0] is the species tree, so whatever is not in there is an unknown */
		/* k is now the number of the unkown taxa */
		printf("unknown to be mapped: %s\n", taxa_names[k]);
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
			printf("a==\n");
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
			
		printf("best reconstruction with a score of %f\n", best_total);
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
	free(presence);
	free(overall_placements);
	}

float get_recon_score(char *giventree, int numspectries, int numgenetries)
	{
	struct taxon *position = NULL, *species_tree = NULL, *gene_tree = NULL, *best_mapping = NULL, *copy = NULL, *temp_top1 = NULL, *temp_top2 = NULL, *spec_copy = NULL;
	int i, j, k, l, m, q, r, spec_start=0, spec_end, gene_start, gene_end, num_species_internal = 0, error = FALSE, num_species_roots = 0, basescore = 1, rand1=0, rand2=0, dospecrand = 1, dogenerand=1;
	float *overall_placements = NULL, biggest = -1, total, best_total = -1, sum_of_totals = 0, rooting_score = -1;
	char *temptree, temptree1[TREE_LENGTH];

	temptree = malloc(TREE_LENGTH*sizeof(char));
	temptree[0] = '\0';
	temptree1[0] = '\0';


	
		/**** BUILD THE SPECIES TREE *****/
		strcpy(temptree,giventree);
		returntree(temptree);
		/* build the tree in memory */
		/****** We now need to build the Supertree in memory *******/
		temp_top = NULL;
		tree_build(1, temptree, species_tree, 1, -1);
		species_tree = temp_top;
		temp_top = NULL;
				printf("4\n");

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
					tree_build(1, temptree, gene_tree, 1, l);
					gene_tree = temp_top;
					temp_top = NULL;
					strcpy(temptree1, temptree);

					if(presence_of_trichotomies(gene_tree)) gene_tree = do_resolve_tricotomies(gene_tree, species_tree, basescore);	
				/*	else printf("no resolving needed\n"); */


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
							/*printf("trying rooting %d of %d in the genetree\t score:", j, i); */
							position = get_branch(gene_tree, j);
							temp_top = gene_tree;
							/*printf("2\n"); */
							reroot_tree(position);
							gene_tree = temp_top;
							temp_top = NULL;
							printf("b\n");
							total = tree_map(gene_tree, species_tree,0);
						/*	printf("%f\t",j, i, total); */
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
	return(rooting_score);
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
				printf("5\n");

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
			printf("\t%s", taxa_names[i]);
		else
			printf("\t%d", i);
		}
	printf("\n");
	for(i=0; i<2*number_of_taxa; i++)	
		{
		if(i<number_of_taxa)
			printf("%s\t", taxa_names[i]);
		else
			printf("%d\t", i);
		for(j=0; j<2*number_of_taxa; j++)
			printf("%d\t", scores[i][j]);
		printf("\n");
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

void reconstruct(int print_settings)
	{
	struct taxon *position = NULL, *species_tree = NULL, *gene_tree = NULL, *best_mapping = NULL, *unknown_fund = NULL, *pos = NULL, *copy = NULL, *newbie = NULL;
	int i, j, k, l, xnum=0, *presence = NULL, **label_results = NULL, num_species_internal = 0, error = FALSE, printfiles = TRUE, how_many = 0, diff_overall =0, dorecon = FALSE, basescore = 1;
	float *overall_placements = NULL, biggest = -1, total, best_total = -1;
	char *temptree, temptree1[TREE_LENGTH], reconfilename[100], otherfilename[100], *tmp1 = NULL, c = '\0';
	FILE *reconstructionfile = NULL, *descendentsfile = NULL, *genebirthfile = NULL;
	
	temptree = malloc(TREE_LENGTH*sizeof(char));
	temptree[0] = '\0';
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
			if(strcmp(parsed_command[i+1], "no") == 0)
				printfiles = FALSE;
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
	
	label_results = malloc(5*sizeof(int*));
	for(i=0; i<5; i++)
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
				printf("Error: '%s' is an invalid value for dups\n", parsed_command[i+1]);
				error = TRUE;
				}
			}
		if(strcmp(parsed_command[i], "losses") == 0)
			{
			loss_weight = tofloat(parsed_command[i+1]);
			if(loss_weight < 0)
				{
				printf("Error: '%s' is an invalid value for losses\n", parsed_command[i+1]);
				error = TRUE;
				}
			}
		}
	if(!error)
		{
		if(print_settings)printf("\nReconstruction of Most likely Duplications and Losses\n\tDuplication weight = %f\n\tLosses weight = %f\n\nTree name:\tReconstruction score\n", dup_weight, loss_weight);
		/**** BUILD THE SPECIES TREE *****/
		strcpy(temptree, fundamentals[0]);
		returntree(temptree);
		/* build the tree in memory */
		/****** We now need to build the Species tree in memory *******/
		temp_top = NULL;
		tree_build(1, temptree, species_tree, 1, -1);
		species_tree = temp_top;
		temp_top = NULL;
		/** add an extra node to the top of the tree */
		temp_top = make_taxon();
		temp_top->daughter = species_tree;
		species_tree->parent = temp_top;
		species_tree = temp_top;
		temp_top = NULL;
		printf("1\n");
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
		
		
		
		for(l=1; l<Total_fund_trees; l++)
			{
			temp_top = NULL;
			strcpy(temptree, "");
			strcpy(temptree, fundamentals[l]);
			unroottree(temptree);
			
			returntree(temptree);
			
			
			/* build the tree in memory */
			/****** We now need to build the gene tree in memory *******/
			temp_top = NULL;
			how_many++;
			taxaorder=0;
			tree_build(1, temptree, gene_tree, 1, l);
			gene_tree = temp_top;
			temp_top = NULL;
			strcpy(temptree1, temptree);
			i = count_taxa(gene_tree, 0);
			
			if(presence_of_trichotomies(gene_tree))
				{
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
			i = number_tree(gene_tree, 0);
			best_total = -1;
			if(i>2)
				{
				for(j=0; j<i; j++)
					{
					malloc_check=0;
					position = get_branch(gene_tree, j);
					temp_top = gene_tree;
					reroot_tree(position);
					gene_tree = temp_top;
					
					diff_overall -= malloc_check;
					printf("c\n");
					total = tree_map(gene_tree, species_tree,1); 
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
				printf("d\n");
				total = tree_map(gene_tree, species_tree,1);				
				best_total = total;
				best_mapping = gene_tree;
				gene_tree = NULL;
				}
			if(strcmp(tree_names[l], "") == 0)
				printf("Tree number: %d\t%f\n", l, best_total);
			else
				printf("%s\t%f\n", tree_names[l], best_total);
			tree_top = best_mapping;
			if(dorecon == TRUE)
				{
				tree_coordinates(temptree, TRUE, FALSE, FALSE, l);
				}
			/******* PRINT_TREE_LABELS TEST *******/	
			/****	label_results[0] = number of copies of the gene at this internal branch
					label_results[1] = number of duplications at this internal branch
					label_results[2] = number of losses at this internal branch
					label_results[3] = number of 1:1 orthologs after this internal branch (allowing duplications and losses in external taxa)
					label_results[4] = number of strict 1:1 orthologs after this internal branch (not allowing any duplications and losses)
			****/

			if(printfiles)
				{	
				for(i=0; i<5; i++)
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
					
					fprintf(reconstructionfile, "%5d (%5d+ %5d- %5d 1:1 %5d 1:1S)\t", label_results[0][j], label_results[1][j], label_results[2][j], label_results[3][j], label_results[4][j]);
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
	free(tmp1);
	free(presence);
	free(overall_placements);
	if(label_results != NULL)
		{
		for(i=0; i<5; i++)
			{
			free(label_results[i]);
			label_results[i] = NULL;
			}
		free(label_results);
		label_results = NULL;
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

void reset_tag2(struct taxon * position)
	{
	while(position != NULL)
		{
		position->tag2 = FALSE;
		if(position->daughter != NULL) reset_tag2(position->daughter);
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

void hgt_reconstruction()
	{
	struct taxon *position = NULL, *species_tree = NULL, *gene_tree = NULL, *best_mapping = NULL, *best_mapping1 = NULL, *best_mapping2 = NULL, *unknown_fund = NULL, *posit = NULL,*copy = NULL, *copy1 = NULL, **parts = NULL, *test_part = NULL, *pos = NULL, *best_donor = NULL, *best_HGT = NULL, *attached = NULL;
	int i, j, k, l,  *presence = NULL,*presence1 = NULL, *presence2 = NULL, hgt_receipient1, hgt_receipient2, **overall_presence = NULL, *overall_reconstruction = NULL, *overall_receptor = NULL, receptor, **tmp_presence1 = NULL, **tmp_presence2 = NULL, *before1 = NULL, *before2 = NULL, *after1 = NULL, *after2 = NULL, *temporary = NULL;
	float *overall_placements = NULL, biggest = -1,  total, best_total = -1,  best_total1 = -1, best_total2 = -1, HGT1 = 0, HGT2 = 0, original = 0, best_HGT_recon = -1, best_reconstruction = -1, sum, HGT_score = -1, donor_score = -1, tmp_allow = FALSE;
	char *temptree = NULL, temptree1[TREE_LENGTH];
	int **species_allowed = NULL, **dependent_species_allowed = NULL, *previous = NULL, xnum =0, x, y, z, partA, partB, q, r, s, allow_HGT1 = TRUE, allow_HGT2 = TRUE, numparts = 1,  place_marker = 1, found_better = FALSE, error = FALSE; 
	int basescore = 1; /** see reconstrution command **/
	
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
	
	printf("Calculating best reconstruction of Duplications, losses and Horizontal gene transfers (HGT)\n\nDuplication weight set to %f\nLoss weight set to %f\nHGT weight set to %f\n\n", dup_weight, loss_weight, hgt_weight);
	
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
	tree_build(1, temptree, species_tree, 1, -1);
	species_tree = temp_top;
	temp_top = NULL;
	/** add an extra node to the top of the tree */
	temp_top = make_taxon();
	temp_top->daughter = species_tree;
	species_tree->parent = temp_top;
	species_tree = temp_top;
	temp_top = NULL;
			printf("2\n");

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
		tree_build(1, temptree, gene_tree, 1, l);
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
			printf("e\n");
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
			printf("ON part %d out of %d parts\n", k, numparts);
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
							printf("f\n");
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
						printf("g\n");
						best_total1 = tree_map(copy1, species_tree,0);
						find_tagged(copy1, presence1);
						best_mapping1 = copy1; /* this will be used later for checking HGT compatibilities (( this version is only used if k != 0 ))*/
						copy1 = NULL;
						
						tree_top = NULL;
						duplicate_tree(test_part, NULL);
						copy1 = tree_top;
						printf("h\n");
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
										
										printf("I\n");
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
									printf("j\n");
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
				printf("HGT%d\n", numparts);
				overall_receptor[place_marker] = receptor;
				overall_placements[k] = donor_score;
				overall_placements[place_marker] = HGT_score;
				if(parts[k] != NULL) dismantle_tree(parts[k]);
				if(parts[place_marker] != NULL) dismantle_tree(parts[place_marker]);
				parts[k] = best_donor;
				best_donor = NULL;
				strcpy(temptree, "");
				print_tree(parts[k], temptree);
				printf("donor:\n%s\n", temptree);
				parts[place_marker] = best_HGT;
				strcpy(temptree, "");
				print_tree(parts[place_marker], temptree);
				printf("hgt:\n%s\n", temptree);

				best_HGT = NULL;
				assign_hgtdonors(parts[k], num_gene_nodes, place_marker);
				for(q=0; q<2*number_of_taxa; q++) overall_presence[place_marker][q] = presence[q]; /* record the possible donor nodes for this HGT */
				
				numparts++;
				k--; /* we need to re-evaluate the donor again to see if removing this HGT means that another HGT becomes apparant */
				}
			
			}
		printf("\n\nprinting results\n");
			
		/**** print out the results */
		i=0;
		while(parts[i] != NULL && i < num_gene_nodes)
			{
			printf("part %d\nScore = %f\nReceptor Node:%d\nPossible Donors: ",i,overall_placements[i], overall_receptor[i] );
			for(j=0; j<2*number_of_taxa; j++) printf("%d,", overall_presence[i][j]);
			printf("\n");
			strcpy(temptree, "");
			strcpy(temptree, "");
			print_tree(parts[i], temptree);
			printf("parts[%d] %s\n",i, temptree);
			printf("\n");
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



/* function for Konrad to automatically collapse clades that have an average branch length of less than X (user defined), it is to leave one randomly chosen taxa to represent the clade. */

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
	tree_build(1, fund_tree, tree_top, 0, -1);
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
		showtrees();
		printf("\n\tSummary of pruned taxa writtin to file prunedtaxa.txt\n");
		}
	else
		{
		printf("0 clades met the criteria specified\n");
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
	char flt_length[TREE_LENGTH];
	
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
	
	if(position->parent != NULL) printf("^^ position->parent %d\n", position->parent->tag);
	while(position != NULL)
		{
		if(position->prev_sibling != NULL) printf("[prev %d\t", position->prev_sibling->tag);
		else printf("[no prev\t");
		if(position->parent != NULL) printf("^^ parent = %d ^^\t ", position->parent->tag);
		if(position->daughter != NULL)
			printf("pointer vv%d (daughter = %d)\t", position->tag, position->daughter->tag);
		else
			printf("taxa %d\t", position->name);
		if(position->next_sibling != NULL) printf("\tnext %d]\t", position->next_sibling->tag);
		else printf("no next]\n");
		position = position->next_sibling;
		}
	position = start;
	while(position != NULL)
		{
		if(position->daughter != NULL) check_treeisok(position->daughter);
		position = position->next_sibling;
		}
	}



/* function for Karen to automatically collapse clades that have all of the same taxa in them, keeping only the one with the longest sequence (found in the full name) */

void prune_monophylies(void)
    {
    int i=0, j=0;
    char *pruned_tree = NULL, *tmp = NULL, filename2[1000];
    FILE *pm_outfile = NULL; 
    select_longest=FALSE;
    filename2[0]='\0';
    strcpy(filename2, "prunedtrees.txt");
    
    tmp = malloc(TREE_LENGTH*sizeof(char));
    tmp[0] = '\0';
    pruned_tree = malloc(TREE_LENGTH*sizeof(char));
    pruned_tree[0] = '\0';

    for(i=0; i<num_commands; i++)
        {
        if(strcmp(parsed_command[i], "selection") == 0)
            {       
            if(strcmp(parsed_command[i+1], "length")==0) select_longest=TRUE; 
            }
        if(strcmp(parsed_command[i], "filename") == 0)
            {       
            strcpy(filename2, parsed_command[i+1]); 
            }

        }

    pm_outfile = fopen(filename2, "w");


    printf("input trees will be pruned where clans of single species exist\nOne remaining representative will be chosen by ");
    if(select_longest == TRUE) printf ("the length of the sequence (in the name, after the species (delimited with a \".\")\n");
        else printf("random\n");

    for(j=0; j<Total_fund_trees; j++)
        {
        if(tree_top != NULL) dismantle_tree(tree_top);  /* Dismantle any trees already in memory */
        tree_top = NULL;
        
        temp_top = NULL;
        taxaorder=0;
        tree_build(1, fundamentals[j], tree_top, 0, j); /* build the tree passed to the function */
        tree_top = temp_top;
        temp_top = NULL;
        reset_tree(tree_top);

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

        if(strcmp(tree_names[j], "") != 0)
        	fprintf(pm_outfile, "%s[%s]\n", pruned_tree, tree_names[j]);
        else
        	fprintf(pm_outfile, "%s[%d]\n", pruned_tree, j);
        }
    
    printf("\nPruning finished. All pruned trees written to the file \"%s\"\n", filename2);
    free(tmp);
    free(pruned_tree);
    fclose(pm_outfile);
    }


void identify_species_specific_clades(struct taxon * position)
    {
    float total = 0;
    int count = 0, taxa_count = 0, all_same_taxon = -1, *foundtaxa = NULL, i;
    long seqlength = 0;
    struct taxon * starting = position, *longest=NULL;
    foundtaxa = malloc(number_of_taxa*sizeof(int));

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
                identify_species_specific_clades(position->daughter); /* if this node has a daughter, then call new instance of the function on the daughter */
            position = position->next_sibling;
            }
        }
    else /* if this was a species-specific clade, then untag everything, except for the longest */
        {
        untag_nodes_below(position);
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
            untag_nodes_above(position->parent);
            }

    free(foundtaxa);

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


void untag_nodes_below(struct  taxon * position)
    {
    
    while(position != NULL)
        {
        if(position != longestseq) position->tag = FALSE;
        if(position->daughter != NULL) untag_nodes_below(position->daughter);
        position = position->next_sibling;
        }
    
    }
    
void untag_nodes_above(struct  taxon * position)
    {
    struct taxon * origin = position;
    while(position->prev_sibling != NULL) position = position->prev_sibling;

    while(position != NULL)
        {
        if(position != longestseq) position->tag = FALSE;
        if(position->daughter != NULL && position != origin) untag_nodes_below(position->daughter);
        if(position->parent != NULL) untag_nodes_above(position->parent);
        position = position->next_sibling;
        }
    
    }


void tips(int num)
	{
	switch(num)
		{
		case(0):
			printf("\n\t1. Clann can be used to transform nexus formatted tree files into newick formatted files.\n\tThis is done by executing the nexus file as normal and then using the command: \"showtrees savetrees=yes\"\n\tIt is also possible to set the name of the file to which the trees are saved, and to stop clann from displaying a graphical representation of each source tree while this is done\n\n");
			break;
		
		case(1):
			printf("\n\t2. Clann can be told only to read the first few characters of each taxa name whenreading the source trees into memory\n\tThis is useful when it is necessary to have unique identifiers (for instance gene IDs) on the source trees\n\tThe option \'maxnamelen\' in the \"exe\" command sets this value\n\tIf the names are not fixed widths, the option \'maxnamelen=delimited\' tells Clann to look for the fist dot \".\" which will specifying the end of the taxon ID in the trees\n\n\t\tFor instance using \"exe maxnamelen=delimited\" on the following tree:\n\t\t(apple.00121,(orange.1435,lemon.3421),pear.1032);\n\t\tResults in clann ignoring the numbers after the dots in the taxa names\n");
			break;

		case(2):
			printf("\n\t3. The equals sign (=), hyphen (-) and space ( ) are special characters in Clann and by default cannot be used in filenames to be read by clann\n\tIf a filename contain some of these characters Clann can only read the name of the file properly by putting the name in inverted commas.\n\t\tFor example: exe \"my-file.txt\"\n");
			break;

		case(3):
			printf("\n\t4. The first command that you should run if you don’t know what to do is \"help\".\n\tThis will display the list of the commands that are available.\n\tCalling any of the commands followed by a question mark (for instance \"hs ?\"), will display the options and defaults associated with that command\n");
			break;
		
		case(4):
			printf("\n\t5. The command \"!\" runs a shell terminal on Unix and Mac operating systems allowing system commands can be run without having to quit Clann\n");
			break;

		case(5):
			printf("\n\t6. Clann can assess supertrees created using other programs\n\tUsing the \"usertrees\" command, clann will read in the file specified and assess all the trees it contains\n\tThe best supertree found in the file is displayed along with its score\n");
			break;

		case(6):
			printf("\n\t7. All commands in Clann should be written completely in lowercase, typing the command \"boot\" is not the same as \"Boot\" and only the first will be recognised as a valid command\n");
			break;

		case(7):
			printf("\n\t8. Heuristic and exhaustive searches of supertree space can be interrupted using the key combination \"control-c\"\n\tThis dispays the score if the best tree found so far and give the user the option to stop the search now or continue.\n\tIf this is done during the random sampling phase of a heuristic search, it will allow the user to move straight to the heuristic search without completing the random sampling\n");
			break;

		case(8):
			printf("\n\t9. Users can assess different configurations of their data by excluding (or including) certain source trees from subsequent commands using the \"excludetrees\" and \"includetrees\" commands\n\tSource trees can be selected based on their name,the taxa they contain, their size (number of taxa they contain) or their score when compared to a supertree\n");
			break;

		case(9):
			printf("\n\t10. Individual (or multiple) taxa can be pruned from the source trees using the command \"deletetaxa\"\n\tBranch lengths are adjusted to take the deletion of thetaxa into account\n\tIf the deletion of taxa from a source tree means that there are less than 4 taxa remaining, that source tree is removed from the analysis\n\tClann will display the names of the source trees removed if this occurs\n");
			break;

		default:
			printf("\n\t1. Clann can be used to transform nexus formatted tree files into newick formatted files.\n\tThis is done by executing the nexus file as normal and then using the command: \"showtrees savetrees=yes\"\n\tIt is also possible to set the name of the file to which the trees are saved, and to stop clann from displaying a graphical representation of each source tree while this is done\n\n");
			printf("\n\t2. Clann can be told only to read the first few characters of each taxa name when reading the source trees into memory\n\tThis is useful when it is necessary to have unique identifiers (for instance gene IDs) on the source trees\n\tThe option \'maxnamelen\' in the \"exe\" command sets this value\n\tIf the names are not fixed widths, the option \'maxnamelen=delimited\' tells Clann to look for the fist dot \".\" which will specify the end of the taxon ID in the trees\n\n\t\tFor instance using \"exe maxnamelen=delimited\" on the following tree:\n\t\t(apple.00121,(orange.1435,lemon.3421),pear.1032);\n\t\tResults in clann ignoring the numbers after the dots in the taxa names\n");
			printf("\n\t3. The equals sign (=), hyphen (-) and space ( ) are special characters in Clann and by default cannot be used in filenames to be read by clann\n\tIf a filename contain some of these characters Clann can only read the name of the file properly by putting the name in inverted commas.\n\t\tFor example: exe \"my-file.txt\"\n");
			printf("\n\t4. The first command that you should run if you don’t know what to do is \"help\".\n\tThis will display the list of the commands that are available.\n\tCalling any of the commands followed by a question mark (for instance \"hs ?\"), will display the options and defaults associated with that command\n");
			printf("\n\t5. The command \"!\" runs a shell terminal on Unix and Mac operating systems allowing system commands can be run without having to quit Clann\n");
			printf("\n\t6. Clann can assess supertrees created using other programs\n\tUsing the \"usertrees\" command, clann will read in the file specified and assess all the trees it contains\n\tThe best supertree found in the file is displayed along with its score\n");
			printf("\n\t7. All commands in Clann should be written completely in lowercase, typing the command \"boot\" is not the same as \"Boot\" and only the first will be recognised as a valid command\n");
			printf("\n\t8. Heuristic and exhaustive searches of supertree space can be interrupted using the key combination \"control-c\"\n\tThis dispays the score if the best tree found so far and give the user the option to stop the search now or continue.\n\tIf this is done during the random sampling phase of a heuristic search, it will allow the user to move straight to the heuristic search without completing the random sampling\n");
			printf("\n\t9. Users can assess different configurations of their data by deleting certain source trees using the \"deletetrees\" command\n\tSource trees can be selected based on their name,the taxa they contain, their size (number of taxa they contain) or their score when compared to a supertree\n");
			printf("\n\t10. Individual (or multiple) taxa can be pruned from the source trees using the command \"deletetaxa\"\n\tBranch lengths are adjusted to take the deletion of thetaxa into account\n\tIf the deletion of taxa from a source tree means that there are less than 4 taxa remaining, that source tree is removed from the analysis\n\tClann will display the names of the source trees removed if this occurs\n");
			break;

		}

	}



