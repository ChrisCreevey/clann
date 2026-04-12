/*
 *  treecompare2.c — Clann v5.0.0
 *  Investigating phylogenetic information through supertree analyses
 *
 *  Created by Chris Creevey, 2003.
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

#include "clann.h"
#include "utils.h"
#include "topology.h"
#include "viz.h"
#include "tree_io.h"
#include "tree_ops.h"
#include "scoring.h"
#include "consensus.h"
#include "reconcile.h"
#include "prune.h"
#include "main.h"
#include "spr_tree.h"
#include "treecluster.h"

BipartSet *fund_bipart_sets = NULL;  /* [Total_fund_trees], precomputed once per analysis */


/* LandscapeMap functions: declared in topology.h */

/*********** Function definitions ********************/

void find(struct taxon * position);
int run_main(int argc, char *argv[]);
/* moved to tree_io.h */ /* input_file_summary */
void clean_exit(int error);
void totext(int c, char *array);
/* moved to tree_io.h */ /* assign_taxa_name */
void execute_command(char *filename, int do_all);
int seperate_commands(char *command);
int parse_command(char *command);
void print_commands(int num);
void cal_fund_scores(int printfundscores);
void pathmetric(char *string, int **scores);
void weighted_pathmetric(char *string, float **scores, int fund_num);
int unroottree(char * tree);
void alltrees_search(int user);
int * apply_singlecopy_filter(void);
void restore_singlecopy_filter(int *saved);
#ifdef _OPENMP
static void hs_alloc_thread_state(void);
static void hs_free_thread_state(void);
static void hs_merge_results(char ***par_retained, float **par_scores, int *par_n, float *par_best, int *par_NUMSWAPS, VisitedSet *par_visited);
static int  hs_same_topology(char *t1, char *t2);
#endif
/* vs_*, tree_topo_hash: declared in topology.h */
float compare_trees(int spr);
void  rf_precompute_fund_biparts(void);
float compare_trees_rf(int spr);
float compare_trees_ml(int spr);
static float       ml_display_score(float s);                    /* total score → lnL  */
static float       ml_display_source_score(float raw, int tidx); /* per-tree score → lnL */
static const char *ml_score_label(void);                         /* "lnL" or "score"   */
struct taxon * make_taxon(void);
void intTotree(int tree_num, char *array, int num_taxa);
int tree_build (int c, char *treestring, struct taxon *parent, int fromfile, int fund_num, int *taxaorder);
void prune_tree(struct taxon * super_pos, int fund_num);
int treeToInt(char *array);
int shrink_tree (struct taxon * position);
int print_pruned_tree(struct taxon * position, int count, char *pruned_tree, int fullname, int treenum);
void reset_tree(struct	taxon * position);
int count_taxa(struct taxon * position, int count);
void check_tree(struct taxon * position, int tag_id, FILE *reconstructionfile);
int check_taxa(struct taxon * position);
int find_taxa(struct taxon * position, char *query);
int number_tree(struct taxon * position, int num);
void dismantle_tree(struct taxon * position);
void bootstrap_search(void);
void mlscores(void);
void memory_error(int error_num);
void print_named_tree(struct taxon * position, char *tree);
void print_fullnamed_tree(struct taxon * position, char *tree, int fundtreenum);
void print_tree(struct taxon * position, char *tree);
void reallocate_retained_supers(void);
void usertrees_search(void);
void heuristic_search(int user, int print, int sample, int nreps);
int average_consensus(int nrep, int missing_method, char * useroutfile, FILE *paupfile);
int do_search(char *tree, int user, int print, int maxswaps, FILE *outfile, int numspectries, int numgenetries);
int branchswap(int number_of_swaps, float score, int numspectries, int numgenetries);
static void fix_parent_pointers(struct taxon *pos, struct taxon *parent);
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
/* xposition1/2, middle_number, yposition0/1/2, print_coordinates,
 * tree_coordinates: declared in viz.h */
void generatetrees(void);
/* draw_histogram: declared in viz.h */
void consensus(int num_trees, char **trees, int num_reps, float percentage, FILE *outfile, FILE *guidetreefile);
/* moved to tree_io.h */ /* input_fund_tree */
/* moved to tree_io.h */ /* nexusparser */
void do_consensus(void);
/* moved to tree_io.h */ /* comment */
/* moved to tree_io.h */ /* showtrees */
/* moved to tree_io.h */ /* exclude */
/* moved to tree_io.h */ /* returntree */
/* moved to tree_io.h */ /* returntree_fullnames */
/* moved to tree_io.h */ /* quick (now static in tree_io.c) */
/* moved to tree_io.h */ /* qs (now static in tree_io.c) */
/* moved to tree_io.h */ /* include */
/* moved to tree_io.h */ /* exclude_taxa */
void sourcetree_dists();
/* moved to tree_io.h */ /* prune_taxa_for_exclude (now static in tree_io.c) */
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
int number_tree2(struct taxon * position, int num);
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
int identify_species_specific_clades(struct taxon * position, int numt, int *taxa_fate, int clannID);
void prune_monophylies();
void untag_nodes_below(struct  taxon * position, int * taxa_fate, int clannID);
void untag_nodes_above(struct  taxon * position, int * taxa_fate, int clannID);
void tips(int num);
void get_taxa_details(struct taxon *position);
void get_taxa_names(struct taxon *position, char **taxa_fate_names);
int basic_tree_build (int c, char *treestring, struct taxon *parent, int fullnames);
int sort_tree(struct taxon *position);
int spr_new(struct taxon * master, int maxswaps, int numspectries, int numgenetries);
static int spr_new2(struct taxon *master, int maxswaps, int numspectries, int numgenetries);
static int tbr_new(struct taxon *master, int maxswaps, int numspectries, int numgenetries);
static int evaluate_candidate(const char *candidate_nwk, float *tmp_fund_scores, int numspectries, int numgenetries);
static int spr_new3(struct taxon *master, int maxswaps, int numspectries, int numgenetries);
static int tbr_new2(struct taxon *master, int maxswaps, int numspectries, int numgenetries);
void do_log(void);
void print_splash(void);




void controlc1(int signal);
void controlc2(int signal);
void controlc3(int signal);
void controlc4(int signal);
void controlc5(int signal);

/****************** Global variable definitions ****************/
FILE * infile = NULL, *BR_file = NULL, *commands_file=NULL, *psfile = NULL, *logfile = NULL, *distributionreconfile = NULL, *onetoonefile = NULL, *strictonetoonefile = NULL, *tempoutfile = NULL;
char **taxa_names = NULL, *commands_filename = NULL, ***fulltaxanames = NULL, **parsed_command = NULL, **fundamentals = NULL, **stored_funds = NULL, **retained_supers = NULL, **stored_commands = NULL, *tempsuper = NULL, **best_topology = NULL, **tree_names = NULL;
char **original_fundamentals = NULL;   /* originals preserved for reconstruct when autoprunemono is active */
int   autoprunemono_active   = 0;      /* set to TRUE when autoprunemono=yes was used at load time */
int  *numtaxaintrees = NULL, fullnamesnum = 0, fullnamesassignments = 1, fundamental_assignments = 0, tree_length_assignments = 1, parsed_command_assignments = 1, name_assignments = 0, *taxa_incidence = NULL, number_of_taxa = 0, Total_fund_trees = 0, *same_tree = NULL, **Cooccurrance = NULL, NUMSWAPS = 0;
int ***fund_scores = NULL, ***stored_fund_scores = NULL, **super_scores = NULL, *number_of_comparisons = NULL, *stored_num_comparisons = NULL, **presence_of_taxa = NULL, **stored_presence_of_taxa = NULL, *presenceof_SPRtaxa = NULL;
int seed, num_commands = 0, number_retained_supers = 10, number_of_steps = 99999, largest_tree = 0, smallest_tree = 1000000, criterion = 0, parts = 0, **total_coding = NULL, *coding_from_tree = NULL, total_nodes = 0, quartet_normalising = 3, splits_weight = 2, dweight =1, *from_tree = NULL, method = 3, tried_regrafts = 0, hsprint = TRUE, max_name_length = NAME_LENGTH, got_weights = FALSE, num_excluded_trees = 0, num_excluded_taxa = 0, calculated_fund_scores = FALSE, select_longest=FALSE;
struct taxon *tree_top = NULL, *temp_top = NULL, *temp_top2 = NULL, *branchpointer = NULL, *longestseq = NULL;
float *scores_retained_supers = NULL, *partition_number = NULL, num_partitions = 0, total_partitions = 0, sprscore = -1, *best_topology_scores = NULL, **weighted_scores = NULL, *sourcetree_scores = NULL, *tree_weights = NULL;
float *score_of_bootstraps = NULL, *yaptp_results = NULL, largest_length = 0, dup_weight = 1, loss_weight = 1, hgt_weight = 1, BESTSCORE = -1;
float ml_beta = 1.0f;   /* L.U.st exponential slope parameter (default 1.0) */
int   ml_scale = 2;     /* ML score scale: 0=paper (raw sum), 1=lust (log10), 2=lnl (default) */
double ml_eta = 0.0;    /* [experimental] tree-size scaling exponent: 0=Steel 2008, 1=normalised, >1=downweight large trees */
double *ml_norm_logZ  = NULL;  /* log Z_{T_i|X_i} per source tree; populated by compare_trees_ml when ml_do_normcorr=1 */
int     ml_do_normcorr = 0;    /* 0=off, 1=apply Bryant & Steel (2008) normalising constant correction in usertrees */
int   bsweight = 0;     /* use per-split BS support as weights in sfit/qfit (0=off, 1=on) */
time_t interval1, interval2;
double sup=1;
char saved_supertree[TREE_LENGTH],  *test_array, inputfilename[10000], delimiter_char = '.', logfile_name[10000], system_call[100000];
volatile int user_break = FALSE;
int trees_in_memory = 0, *sourcetreetag = NULL, remainingtrees = 0, GC, delimiter = TRUE, print_log = FALSE, num_gene_nodes, testarraypos = 0;
int malloc_check =0, count_now = FALSE, another_check =0;
unsigned int thread_seed = 0;   /* per-thread random seed for rand_r(); see threadprivate block below */
uint64_t *taxon_hash_vals = NULL; /* shared read-only: splitmix64 weight per taxon; set at taxa load time */
VisitedSet *visited_set   = NULL; /* threadprivate: per-replicate visited topology hash set */
VisitedSet *thread_visited_acc = NULL; /* threadprivate: accumulates all unique topos across reps in one thread */
LandscapeMap *landscape_map = NULL; /* threadprivate: per-thread visited-topology landscape accumulator */
/* Shared landscape globals (not threadprivate) */
static char         g_landscape_file[4096] = ""; /* filename; empty = feature disabled */
static LandscapeMap *g_landscape_map = NULL;      /* global accumulator across all threads */
/* Landscape clustering options (set by hs option parsing) */
static int   g_cluster_enabled   = 0;                   /* 0=off, 1=on (requires visitedtrees=) */
static char  g_cluster_output[4096] = "treeclusters.tsv"; /* output TSV filename */
static float g_cluster_threshold = 0.2f;                /* max normalized RF distance [0,1] */
static int   g_cluster_orderby   = 0;                   /* 0=score ascending, 1=visits descending */
time_t  rep_start_time    = 0;    /* threadprivate: wall-clock time when current do_search() began */
int     hs_do_print       = 0;    /* threadprivate: mirrors the 'print' param of do_search() */
float   last_status_score = -1.0f;/* threadprivate: sprscore at last periodic status line (for improvement marker) */
/* Shared (non-threadprivate): parallel-mode global progress state, written under omp critical */
float  par_progress_best   = -1.0f; /* best score seen so far across ALL parallel threads */
float  par_last_print_score= -1.0f; /* par_progress_best at last status line (thread-0 only) */
time_t par_search_start    = 0;     /* wall-clock time when the parallel region began */
int    skip_streak         = 0;     /* threadprivate: consecutive already-visited skips since last new topology */
int    nni_swaps           = 0;     /* threadprivate: NNI refinement swaps performed after SPR/TBR in last do_search() rep */
int    hs_maxskips         = -1;    /* shared: stop replicate when skip_streak reaches this (0=disabled, -1=auto: N*N) */
int    hs_maxskips_is_auto = 1;    /* 1=auto-scale to N*N at search start; 0=user has set an explicit value */
int    hs_strategy         = 0;    /* 0=first-improvement (depth-first); 1=best-improvement (breadth-first) */
int    hs_progress_interval= 5;    /* shared: parallel progress print interval in seconds (0=every improvement) */
time_t par_last_progress_time = 0; /* shared: wall-clock time of last parallel progress line printed */
float  hs_droprep          = 0.0f; /* shared: abandon rep if its score is >droprep fraction above par_progress_best; 0=disabled */
int    rep_abandon         = 0;    /* threadprivate: set when this rep should be abandoned due to droprep */
int    hs_par_rep          = 0;    /* threadprivate: 1-based rep number for the rep currently running on this thread */
int    hs_thread_report_interval = 0; /* shared: per-thread status report interval in seconds (0=disabled) */
time_t thread_report_last = 0;        /* threadprivate: wall-clock time of last per-thread status print */

/****** OpenMP thread-private state: one independent copy per thread in parallel regions ******/
#ifdef _OPENMP
#pragma omp threadprivate( \
    tree_top, temp_top, temp_top2, branchpointer, \
    super_scores, sourcetree_scores, presenceof_SPRtaxa, \
    sprscore, tried_regrafts, \
    retained_supers, scores_retained_supers, \
    best_topology, best_topology_scores, number_retained_supers, \
    BESTSCORE, NUMSWAPS, \
    thread_seed, \
    visited_set, thread_visited_acc, landscape_map, \
    rep_start_time, hs_do_print, last_status_score, \
    skip_streak, nni_swaps, rep_abandon, hs_par_rep, thread_report_last, \
    fundamentals, presence_of_taxa, fund_scores, number_of_comparisons, \
    fund_bipart_sets \
)
#endif

/* Snapshots of the master thread's input-data pointers, set just before each
 * heuristic_search parallel region so that worker threads can share them. */
#ifdef _OPENMP
static char **g_hs_fundamentals_snap    = NULL;
static int  **g_hs_presence_snap        = NULL;
static int ***g_hs_fund_scores_snap     = NULL;
static int   *g_hs_num_comp_snap        = NULL;
static BipartSet  *g_hs_fund_bipart_snap = NULL;  /* snapshot of master's fund_bipart_sets for HS worker threads */
#endif



/* CLI helpers, main, seperate_commands, print_splash, print_commands,
 * parse_command, recount_from_tree, autoprunemono_apply, execute_command:

    





















    

/* calculate the path metric for each of the fundamental trees */




    
    















 
 int sort_tree(struct taxon *position)   
 	{
    	/* This function will take as input a tree built in memory */
 		/* It will return the same tree, but where the compoenents have been ordered so that they are in the numberical order that would have arose from the intTOtree */

 	struct taxon *start = position;
 	int mintaxon = number_of_taxa;

 	while(position != NULL)
 		{	
 		if(position -> daughter != NULL)
 			{
 			position->tag = sort_tree(position->daughter);
 			}
 		else
 			{	
 			position->tag = position->name;
 			}
 		position = position->next_sibling;
 		}
 	position = start;

	while(position != NULL)
 		{
 		if(position->tag < mintaxon) mintaxon = position->tag;
 		}

 	/* now sort the sibliing on this level based on the number in the tags */
 	return(1);
    }
    
    
    /* this function builds a tree given its number and the number of taxa */

/* toint: moved to utils.c */




/* tofloat: moved to utils.c */
     



void alltrees_search(int user)
    {
    int i = 0, j=0, all = TRUE, start = 0, end = 0, error = FALSE, keep = 0;
    char *tree = NULL, *best_tree = NULL, outfilename[100];
    float score = 0, best_score = 0, worst = 0;
    FILE *treesfile = NULL, *userfile = NULL;
    int *saved_tags = NULL;  /* for single-copy auto-filter */
    
    
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
                    printf2("Error: '%s' is an invalid value for keep\n", parsed_command[i+1]);
                    error = TRUE;
                    }
                }
            if(strcmp(parsed_command[i], "nbest") == 0)
                {
                keep = toint(parsed_command[i+1]);
                if(keep == 0)
                    {
                    printf2("Error: '%s' is an invalid values for nbest\n", parsed_command[i+1]);
                    error = TRUE;
                    }
                }
            if(strcmp(parsed_command[i], "savetrees") == 0 && criterion != 1)
             {
                if((userfile = fopen(parsed_command[i+1], "w")) == NULL)
                 {
                    printf2("Error opening file named %s\n", parsed_command[i+1]);
                    error = TRUE;
                 }
                else
                 {
                    printf2("opened output file %s\n", parsed_command[i+1]);
                    strcpy(outfilename, parsed_command[i+1]);
                 }
             }
            
            if(strcmp(parsed_command[i], "create") == 0)
                {
                if((treesfile = fopen("alltrees.ph", "w")) == NULL) 
                    {
                    printf2("Error opening file named 'alltrees.ph'\n");
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
                                printf2("Error: weight option '%s' is unknown\n", parsed_command[i+1]);
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
                            printf2("Error: weight option '%s' is unknown\n", parsed_command[i+1]);
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
                            printf2("Error: weight option '%s' is unknown\n", parsed_command[i+1]);
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
                printf2("Error opening file named 'alltrees.ph'\n");
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
            
        printf2("\n\nAlltrees (exhaustive search) settings:\n\trange: tree numbers %d to %d inclusive\n\tOutput file: %s\n\tCreate all trees? ", start, end, outfilename );
        if(treesfile != NULL) printf2("Yes\n");
        else printf2("No\n");
        printf2("\tCriterion = ");
        
        if(criterion==0)
            {
            printf2("DFIT\n\tWeighting Scheme = ");
            if(dweight == 0) printf2("equal\n");
            if(dweight == 1) printf2("comparisons\n");
            }
        if(criterion == 2)
            {
            printf2("SFIT\n\tWeighting Scheme = ");
            if(splits_weight == 1) printf2("equal\n");
            if(splits_weight == 2) printf2("splits\n");
            }
        if(criterion == 3)
            {
            printf2("QFIT\n\tWeighting Scheme = ");
            if(quartet_normalising == 1) printf2("equal\n");
            if(quartet_normalising == 2) printf2("taxa\n");
            if(quartet_normalising == 3) printf2("quartets\n");
            }
        printf2("\n\n");
            
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
    
        if(user) printf2("Progress indicator:");
        
        /************ End assign dynamic arrays **************/
    
        if(criterion == 2 || criterion == 3 || criterion == 6 || criterion == 7)
            rf_precompute_fund_biparts();


        if(start == 0 && end == 0)
            {
            start = 1;
            end = sup;
            }
            
		 if(signal(SIGINT, controlc5) == SIG_ERR)
			{
			printf2("An error occurred while setting a signal handler\n");
			}

        saved_tags = apply_singlecopy_filter();

        for(i=start; i<=end; i++)
            {            
            tree[0] = '\0';
            if(user_break)
				{
				printf2("%d trees sampled\n", i);
				i = end+1;
				}
            interval2 = time(NULL);
            if(difftime(interval2, interval1) > 5) /* every 10 seconds print a dot to the screen */
                {
              /*  printf2("="); */
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
                { int _to = 0; tree_build(1, tree, tree_top, FALSE, -1, &_to); }
                tree_top = temp_top;
                temp_top = NULL;
                score = compare_trees(FALSE);
                }
            if(criterion == 2 || criterion == 3 || criterion == 6 || criterion == 7)
                {
                if(tree_top != NULL) { dismantle_tree(tree_top); tree_top = NULL; }
                temp_top = NULL;
                { int _to = 0; tree_build(1, tree, tree_top, FALSE, -1, &_to); }
                tree_top = temp_top;
                temp_top = NULL;
                if(criterion == 2)      score = compare_trees_sfit(FALSE);
                else if(criterion == 3) score = compare_trees_qfit(FALSE);
                else if(criterion == 6) score = compare_trees_rf(FALSE);
                else                    score = compare_trees_ml(FALSE);
                }

            /* Always populate best_tree with named tree for retained_supers storage.
               For criterion==0/6/7 tree_top is already built; for MRC/QC we convert
               the integer-indexed tree string directly using returntree(). This ensures
               retained_supers[] always holds actual taxon names so that
               'reconstruct speciestree memory' works correctly after alltrees. */
            if((criterion == 0 || criterion == 2 || criterion == 3 || criterion == 6 || criterion == 7) && tree_top != NULL)
                {
                strcpy(best_tree, "");
                print_named_tree(tree_top, best_tree);
                }
            else
                {
                strcpy(best_tree, tree);
                returntree(best_tree);
                }

            if(treesfile != NULL)  /* if the create option was selected */
                fprintf(treesfile,"%s;\t[%f]\n", best_tree, score);



            if(i==1)
                {
				retained_supers[0] = realloc(retained_supers[0], (strlen(best_tree)+10)*sizeof(char));
                strcpy(retained_supers[0], best_tree);
                scores_retained_supers[0] = score;
                best_score = score;
                }
            else
                {
                if(score < best_score)
                    {
					retained_supers[0] = realloc(retained_supers[0], (strlen(best_tree)+10)*sizeof(char));
                    strcpy(retained_supers[0], best_tree);
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
						retained_supers[j] = realloc(retained_supers[j], (strlen(best_tree)+10)*sizeof(char));
                        strcpy(retained_supers[j], best_tree);
                        scores_retained_supers[j] = score;
                        }
                    }
                }
            }
        
        /***** Print out the best tree found ******/
        if(user)
            {
            printf2("\n");
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
                /* retained_supers[] now stores named trees (actual taxon names), so use TRUE */
                { int _to = 0; tree_build(1, retained_supers[i], tree_top, TRUE, -1, &_to); }
                tree_top = temp_top;
                temp_top = NULL;

                strcpy(best_tree, "");

                print_named_tree(tree_top, best_tree);

                if(userfile != NULL) fprintf(userfile, "%s;\t[%f]\n", best_tree, scores_retained_supers[i] );

                tree_coordinates(best_tree, FALSE, TRUE, FALSE, -1);
                printf2("\nSupertree %d of %d %s = %f\n", i+1, j, ml_score_label(), ml_display_score(scores_retained_supers[i]) );
                i++;
                }

            /* Keep retained_supers[] in memory (like hs/nj) so 'reconstruct speciestree memory'
               can use the alltrees result. Set trees_in_memory to the count of best trees found. */
            trees_in_memory = j;


            }
        restore_singlecopy_filter(saved_tags);
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






/* Tree_build:
	This function reads in a file and from it builds the tree in memory using the taxon_type definition */
	

/* Basic_Tree_build:
	This function builds a tree in memory withouth incrementing the number of taxa etc in, bascially this is for building a tree where do don;t needall the bells and whistles. This gets used in Exclude taxa to help deal with gene names */





/* This makes the taxon structure when we need it so I don't have to keep typing the assignments all the time */


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

	






/* This function travels through the tree recursively untagging any pointer siblings that are not being used. 
	this effectively shrinks the tree to the size of the fundamental tree that it is being compared to 
	this recursive function is called by shrink_tree to count how many active taxa there are below any given pointer sibling */
	
	

/* This function is only used to print the pruned supertree */





				




	
	
/* this identifies the taxa in a subtree passed to it */

	




#ifdef _OPENMP
/*  boot_alloc_thread_state -----------------------------------------------
 *  Allocate per-thread copies of the fundamental input-data arrays and
 *  initialise search state for a bootstrap parallel region.
 *  Must be called from inside the bootstrap #pragma omp parallel region.
 */
static void boot_alloc_thread_state(int ntrees, int ntaxa)
    {
    int i, j, is_master;
    is_master = (omp_get_thread_num() == 0);

    /* Allocate thread-private copies of all four input-data arrays.
     * The master thread's copies replace its pre-parallel pointers;
     * worker threads get fresh independent allocations. */
    fundamentals          = malloc(ntrees * sizeof(char *));
    presence_of_taxa      = malloc(ntrees * sizeof(int *));
    fund_scores           = malloc(ntrees * sizeof(int **));
    number_of_comparisons = malloc(ntrees * sizeof(int));
    if(!fundamentals || !presence_of_taxa || !fund_scores || !number_of_comparisons)
        memory_error(200);
    for(i = 0; i < ntrees; i++)
        {
        fundamentals[i]     = malloc(TREE_LENGTH * sizeof(char));
        fundamentals[i][0]  = '\0';
        presence_of_taxa[i] = malloc(ntaxa * sizeof(int));
        fund_scores[i]      = malloc(ntaxa * sizeof(int *));
        if(!fundamentals[i] || !presence_of_taxa[i] || !fund_scores[i]) memory_error(201);
        for(j = 0; j < ntaxa; j++)
            {
            fund_scores[i][j] = malloc(ntaxa * sizeof(int));
            if(!fund_scores[i][j]) memory_error(202);
            }
        }

    /* Allocate / reset search-result state (mirrors hs_alloc_thread_state). */
    if(!is_master)
        {
        int k, init_n = 10;
        retained_supers        = malloc(init_n * sizeof(char *));
        scores_retained_supers = malloc(init_n * sizeof(float));
        best_topology          = malloc(init_n * sizeof(char *));
        best_topology_scores   = malloc(init_n * sizeof(float));
        for(k = 0; k < init_n; k++)
            {
            retained_supers[k]        = malloc(TREE_LENGTH * sizeof(char));
            retained_supers[k][0]     = '\0';
            scores_retained_supers[k] = -1;
            best_topology[k]          = malloc(TREE_LENGTH * sizeof(char));
            best_topology[k][0]       = '\0';
            best_topology_scores[k]   = -1;
            }
        number_retained_supers = init_n;
        super_scores           = malloc(ntaxa * sizeof(int *));
        for(k = 0; k < ntaxa; k++)
            super_scores[k] = malloc(ntaxa * sizeof(int));
        sourcetree_scores  = malloc(ntrees * sizeof(float));
        presenceof_SPRtaxa = malloc(ntaxa  * sizeof(int));
        }
    else
        {
        /* Master: search-state pointers are valid; just reset content. */
        int k;
        for(k = 0; k < number_retained_supers; k++)
            {
            retained_supers[k][0]     = '\0';
            scores_retained_supers[k] = -1;
            if(best_topology && best_topology[k]) best_topology[k][0] = '\0';
            if(best_topology_scores) best_topology_scores[k] = -1;
            }
        if(super_scores)
            {
            for(k = 0; k < ntaxa; k++) free(super_scores[k]);
            free(super_scores);
            }
        super_scores = malloc(ntaxa * sizeof(int *));
        for(k = 0; k < ntaxa; k++)
            super_scores[k] = malloc(ntaxa * sizeof(int));
        }

    /* Common init for all threads */
    {
    int k;
    for(k = 0; k < ntrees; k++) sourcetree_scores[k]  = -1;
    for(k = 0; k < ntaxa;  k++) presenceof_SPRtaxa[k] = -1;
    /* Each bootstrap thread builds its own bipart sets from its resampled fundamentals. */
    fund_bipart_sets = NULL;
    }
    BESTSCORE    = -1;
    NUMSWAPS     = 0;
    sprscore     = -1;
    tried_regrafts = 0;
    tree_top     = NULL;
    temp_top     = NULL;
    temp_top2    = NULL;
    branchpointer = NULL;
    thread_seed  = (unsigned int)(seed + (unsigned int)omp_get_thread_num() * 1000003u);
    }


/*  boot_free_thread_state ------------------------------------------------
 *  Free per-thread fundamental data and search state allocated by
 *  boot_alloc_thread_state.  Must be called inside the parallel region.
 */
static void boot_free_thread_state(int ntrees, int ntaxa)
    {
    int i, j, is_master;
    is_master = (omp_get_thread_num() == 0);
    /* Free search-result state (mirrors hs_free_thread_state). */
    if(tree_top  != NULL) { dismantle_tree(tree_top);  tree_top  = NULL; }
    if(temp_top  != NULL) { dismantle_tree(temp_top);  temp_top  = NULL; }
    if(super_scores)
        {
        for(i = 0; i < ntaxa; i++) free(super_scores[i]);
        free(super_scores); super_scores = NULL;
        }
    if(!is_master)
        {
        if(sourcetree_scores)  { free(sourcetree_scores);  sourcetree_scores  = NULL; }
        if(presenceof_SPRtaxa) { free(presenceof_SPRtaxa); presenceof_SPRtaxa = NULL; }
        for(i = 0; i < number_retained_supers; i++)
            {
            if(retained_supers && retained_supers[i]) free(retained_supers[i]);
            if(best_topology   && best_topology[i])   free(best_topology[i]);
            }
        if(retained_supers)        { free(retained_supers);        retained_supers        = NULL; }
        if(scores_retained_supers) { free(scores_retained_supers); scores_retained_supers = NULL; }
        if(best_topology)          { free(best_topology);          best_topology          = NULL; }
        if(best_topology_scores)   { free(best_topology_scores);   best_topology_scores   = NULL; }
        }

    /* Free per-thread bipart sets built from resampled fundamentals (all threads). */
    if(fund_bipart_sets != NULL)
        {
        int _fi;
        for(_fi = 0; _fi < ntrees; _fi++)
            {
            if(fund_bipart_sets[_fi].hashes)   free(fund_bipart_sets[_fi].hashes);
            if(fund_bipart_sets[_fi].sizes)    free(fund_bipart_sets[_fi].sizes);
            if(fund_bipart_sets[_fi].supports) free(fund_bipart_sets[_fi].supports);
            }
        free(fund_bipart_sets);
        fund_bipart_sets = NULL;
        }

    /* Free the thread-private fundamental data arrays.
     * NON-MASTER threads: free everything.
     * MASTER thread: leave the arrays intact — the restore block in
     * bootstrap_search (after the parallel region) reuses them via
     * realloc/strcpy to restore the original gene trees. */
    if(!is_master)
        {
        for(i = 0; i < ntrees; i++)
            {
            if(fundamentals[i])     free(fundamentals[i]);
            if(presence_of_taxa[i]) free(presence_of_taxa[i]);
            for(j = 0; j < ntaxa; j++)
                if(fund_scores[i][j]) free(fund_scores[i][j]);
            free(fund_scores[i]);
            }
        free(fundamentals);          fundamentals          = NULL;
        free(presence_of_taxa);      presence_of_taxa      = NULL;
        free(fund_scores);           fund_scores           = NULL;
        free(number_of_comparisons); number_of_comparisons = NULL;
        }
    }
#endif /* _OPENMP */


/* -----------------------------------------------------------------------
 * mlscores: estimate the Steel (2008) ML beta parameter for the current
 * supertree and source trees using the closed-form MLE beta = W / WD,
 * where W = sum of tree weights and WD = weighted sum of RF distances.
 * Updates ml_beta so subsequent hs/boot runs use the estimated value.
 *
 * Options:
 *   outfile=<file>   write beta log-likelihood profile (log-uniform scan)
 *   scan=<n>         number of points in profile (default 100)
 *   scanmin=<f>      lower bound of beta scan (default beta/100)
 *   scanmax=<f>      upper bound of beta scan (default beta*10)
 * ----------------------------------------------------------------------- */
void mlscores(void)
    {
    int i, k;
    char outfile_name[10000];
    int    scan = 100, user_set_scanmin = FALSE, user_set_scanmax = FALSE;
    int    do_eta = FALSE, escan = 50, fix_beta = FALSE;
    char   sourcescores_name[10000];
    sourcescores_name[0] = '\0';
    double eta_max = 3.0;
    float  scanmin = 0.0f, scanmax = 0.0f;
    outfile_name[0] = '\0';

    for(i = 1; i < num_commands; i++)
        {
        if(strcmp(parsed_command[i], "outfile") == 0 && i + 1 < num_commands)
            { strncpy(outfile_name, parsed_command[i+1], 9999); outfile_name[9999] = '\0'; i++; }
        else if(strcmp(parsed_command[i], "scan") == 0 && i + 1 < num_commands)
            { scan = atoi(parsed_command[i+1]); i++; }
        else if(strcmp(parsed_command[i], "scanmin") == 0 && i + 1 < num_commands)
            { scanmin = (float)atof(parsed_command[i+1]); user_set_scanmin = TRUE; i++; }
        else if(strcmp(parsed_command[i], "scanmax") == 0 && i + 1 < num_commands)
            { scanmax = (float)atof(parsed_command[i+1]); user_set_scanmax = TRUE; i++; }
        else if(strcmp(parsed_command[i], "eta") == 0 && i + 1 < num_commands)
            { if(strcmp(parsed_command[i+1], "auto") == 0) do_eta = TRUE; i++; }
        else if(strcmp(parsed_command[i], "escan") == 0 && i + 1 < num_commands)
            { escan = atoi(parsed_command[i+1]); i++; }
        else if(strcmp(parsed_command[i], "etamax") == 0 && i + 1 < num_commands)
            { eta_max = atof(parsed_command[i+1]); i++; }
        else if(strcmp(parsed_command[i], "fixbeta") == 0)
            { fix_beta = TRUE; }
        else if(strcmp(parsed_command[i], "sourcescores") == 0 && i + 1 < num_commands)
            { strncpy(sourcescores_name, parsed_command[i+1], 9999); sourcescores_name[9999] = '\0'; i++; }
        }

    if(criterion != 7)
        { printf2("  mlscores: criterion must be set to ml (use 'set criterion ml')\n"); return; }
    if(number_of_taxa == 0 || Total_fund_trees == 0)
        { printf2("Error: You need to load source trees before using this command\n"); return; }
    if(tree_top == NULL)
        { printf2("  mlscores: no supertree in memory (run hs or nj first)\n"); return; }

    if(fund_bipart_sets == NULL) rf_precompute_fund_biparts();

    /* Compute raw (unscaled) RF distances once; reused for both beta and alpha estimation */
    float *dists_raw = malloc(Total_fund_trees * sizeof(float));
    if(!dists_raw) { printf2("mlscores: out of memory\n"); return; }
    compute_raw_rf_dists(dists_raw);

    /* Apply current global ml_eta scaling for the standard beta MLE */
    float *dists = malloc(Total_fund_trees * sizeof(float));
    if(!dists) { free(dists_raw); printf2("mlscores: out of memory\n"); return; }
    for(i = 0; i < Total_fund_trees; i++)
        {
        if(!sourcetreetag[i]) { dists[i] = 0.0f; continue; }
        int ki = fund_bipart_sets[i].count;
        dists[i] = (ml_eta != 0.0 && ki > 0)
                   ? (float)(dists_raw[i] / pow((double)ki, ml_eta))
                   : dists_raw[i];
        }

    double W = 0.0, WD = 0.0;
    for(i = 0; i < Total_fund_trees; i++)
        {
        if(!sourcetreetag[i]) continue;
        double w = (double)tree_weights[i];
        W  += w;
        WD += w * (double)dists[i];
        }
    free(dists);

    if(WD <= 0.0)
        {
        free(dists_raw);
        printf2("  mlscores: all source trees have zero RF distance to supertree; beta is undefined\n");
        return;
        }

    double beta_hat = fix_beta ? (double)ml_beta : W / WD;
    /* Full joint lnL = W*log(beta) + sum(w_i*log(w_i)) - beta*WD - eta*P
     * where P = sum_i log(k_i) for active trees with k_i > 0.
     * sum(w_i*log(w_i)) accounts for individual beta_i = beta*w_i per Steel (2008).
     * When all weights=1 this term is zero.  When eta=0, reduces to Steel (2008). */
    double lnl_penalty = 0.0;
    if(ml_eta != 0.0 && fund_bipart_sets != NULL)
        for(i = 0; i < Total_fund_trees; i++)
            if(sourcetreetag[i] && fund_bipart_sets[i].count > 0)
                lnl_penalty += log((double)fund_bipart_sets[i].count);
    double w_entropy = 0.0;
    for(i = 0; i < Total_fund_trees; i++)
        if(sourcetreetag[i] && tree_weights[i] > 0.0f)
            w_entropy += (double)tree_weights[i] * log((double)tree_weights[i]);
    double logL_hat = W * log(beta_hat) + w_entropy - beta_hat * WD - ml_eta * lnl_penalty;
    if(!fix_beta) ml_beta = (float)beta_hat;

    printf2("  Source trees     : %d   Weighted RF sum : %.4f\n", Total_fund_trees, (float)WD);
    if(ml_eta != 0.0)
        printf2("  ml_eta (active): %.4f  (RF distances scaled by k_i^%.4f)\n", ml_eta, ml_eta);
    if(fix_beta)
        {
        printf2("  Beta (fixed)     : %.6f\n", beta_hat);
        printf2("  Log-likelihood   : %.4f\n", (float)logL_hat);
        printf2("  ml_beta unchanged (fixbeta active)\n");
        }
    else
        {
        printf2("  Optimal beta MLE : %.6f\n", (float)beta_hat);
        printf2("  Log-likelihood   : %.4f\n", (float)logL_hat);
        printf2("  ml_beta updated to %.6f (used for all subsequent hs/boot runs)\n", ml_beta);
        }

    /* --- Beta profile outfile ------------------------------------------- */
    if(outfile_name[0] != '\0')
        {
        if(!user_set_scanmin) scanmin = (float)(beta_hat / 100.0);
        if(!user_set_scanmax) scanmax = (float)(beta_hat * 10.0);
        if(scanmin <= 0.0f) scanmin = 1e-6f;
        if(scanmax <= scanmin) scanmax = scanmin * 1000.0f;
        if(scan < 2) scan = 2;

        FILE *fp = fopen(outfile_name, "w");
        if(!fp)
            { printf2("mlscores: cannot open '%s' for writing\n", outfile_name); }
        else
            {
            fprintf(fp, "# mlscores beta profile log-likelihood\n");
            fprintf(fp, "# Source trees: %d   Weighted RF sum: %.6f\n", Total_fund_trees, WD);
            if(ml_eta != 0.0)
                fprintf(fp, "# ml_eta: %.6f (RF distances scaled by k_i^eta)\n", ml_eta);
            fprintf(fp, "# Optimal beta (MLE): %.6f   Log-likelihood: %.6f\n", beta_hat, logL_hat);
            if(ml_eta != 0.0)
                fprintf(fp, "# eta penalty (eta*sum_log_k): %.6f  (constant offset included in logL)\n",
                        ml_eta * lnl_penalty);
            fprintf(fp, "# beta\tlogL\n");
            double log_min = log((double)scanmin);
            double log_max = log((double)scanmax);
            double step    = (log_max - log_min) / (scan - 1);
            for(i = 0; i < scan; i++)
                {
                double beta_k = exp(log_min + i * step);
                double logL_k = W * log(beta_k) + w_entropy - beta_k * WD - ml_eta * lnl_penalty;
                fprintf(fp, "%.6f\t%.6f\n", beta_k, logL_k);
                }
            fclose(fp);
            printf2("  Beta profile written to: %s\n", outfile_name);
            }
        }

    /* --- Eta grid search -------------------------------------------------
     * Model: P(d_i|S,beta,eta) = (beta/k_i^eta) * exp(-(beta/k_i^eta)*d_i)
     * log L(beta,eta) = n*log(beta) - eta*penalty - beta*WD(eta)
     * penalty = sum_i log(k_i)   [from normalisation constant of each exponential]
     * Profiling beta: beta_hat(eta) = W / WD(eta)
     * log L(eta)   = W*log(W/WD(eta)) - W - eta*penalty
     * The penalty term prevents eta running to infinity.
     * --------------------------------------------------------------------- */
    if(do_eta)
        {
        if(escan < 2) escan = 2;
        if(eta_max <= 0.0) eta_max = 3.0;

        /* Penalty: sum of log(k_i) over tagged trees with splits */
        double penalty = 0.0;
        for(i = 0; i < Total_fund_trees; i++)
            if(sourcetreetag[i] && fund_bipart_sets[i].count > 0)
                penalty += log((double)fund_bipart_sets[i].count);

        double best_logL_a  = -1e300;
        double best_eta     = 0.0;
        double best_beta_a  = beta_hat;  /* stays fixed throughout if fix_beta */

        /* Store grid for output */
        double *g_eta  = malloc(escan * sizeof(double));
        double *g_beta = malloc(escan * sizeof(double));
        double *g_logL = malloc(escan * sizeof(double));
        int     g_n    = 0;

        if(!g_eta || !g_beta || !g_logL)
            { printf2("mlscores eta: out of memory\n"); }
        else
            {
            for(k = 0; k < escan; k++)
                {
                double eta_k = (double)k * eta_max / (escan - 1);
                double WD_k = 0.0;
                for(i = 0; i < Total_fund_trees; i++)
                    {
                    if(!sourcetreetag[i]) continue;
                    int ki = fund_bipart_sets[i].count;
                    double d_k = (eta_k != 0.0 && ki > 0)
                                 ? dists_raw[i] / pow((double)ki, eta_k)
                                 : dists_raw[i];
                    WD_k += (double)tree_weights[i] * d_k;
                    }
                if(WD_k <= 0.0) continue;
                double beta_k = fix_beta ? beta_hat : W / WD_k;
                double logL_k = W * log(beta_k) + w_entropy - beta_k * WD_k - eta_k * penalty;
                g_eta[g_n]  = eta_k;
                g_beta[g_n] = beta_k;
                g_logL[g_n] = logL_k;
                g_n++;
                if(logL_k > best_logL_a)
                    { best_logL_a = logL_k; best_eta = eta_k; best_beta_a = beta_k; }
                }

            ml_eta = best_eta;
            if(!fix_beta) ml_beta = (float)best_beta_a;

            printf2("\n  --- Tree-size scaling exponent (eta) estimation ---\n");
            printf2("  Penalty term (sum log k_i) : %.4f\n", penalty);
            printf2("  Eta grid  : 0 to %.2f  (%d points)\n", eta_max, escan);
            printf2("  Optimal eta : %.4f\n", best_eta);
            if(fix_beta)
                printf2("  Beta (fixed)  : %.6f\n", best_beta_a);
            else
                printf2("  Optimal beta  : %.6f\n", best_beta_a);
            printf2("  Log-likelihood: %.4f\n", best_logL_a);
            printf2("  ml_eta updated to %.4f\n", ml_eta);
            if(fix_beta)
                printf2("  ml_beta unchanged (fixbeta active)\n");
            else
                printf2("  ml_beta  updated to %.6f\n", ml_beta);

            /* Append eta profile to outfile if requested */
            if(outfile_name[0] != '\0' && g_n > 0)
                {
                FILE *fp2 = fopen(outfile_name, "a");
                if(!fp2)
                    { printf2("mlscores: cannot reopen '%s' for eta append\n", outfile_name); }
                else
                    {
                    fprintf(fp2, "\n# --- eta profile ---\n");
                    fprintf(fp2, "# Model: log L(eta) = W*log(W/WD(eta)) - W - eta*penalty\n");
                    fprintf(fp2, "# Penalty (sum log k_i): %.6f\n", penalty);
                    fprintf(fp2, "# Optimal eta (MLE): %.6f\n", best_eta);
                    fprintf(fp2, "# Optimal beta  (MLE): %.6f\n", best_beta_a);
                    fprintf(fp2, "# Joint log-likelihood: %.6f\n", best_logL_a);
                    fprintf(fp2, "# eta\tbeta\tlogL\n");
                    for(k = 0; k < g_n; k++)
                        fprintf(fp2, "%.6f\t%.6f\t%.6f\n", g_eta[k], g_beta[k], g_logL[k]);
                    fclose(fp2);
                    printf2("  Eta profile appended to: %s\n", outfile_name);
                    }
                }
            }

        free(g_eta);
        free(g_beta);
        free(g_logL);
        }

    /* --- Source tree scores file ------------------------------------------ */
    if(sourcescores_name[0] != '\0')
        {
        FILE *ssf = fopen(sourcescores_name, "w");
        if(!ssf)
            { printf2("mlscores: cannot open '%s' for writing\n", sourcescores_name); }
        else
            {
            fprintf(ssf, "# mlscores source tree lnL contributions\n");
            fprintf(ssf, "# input file: %s\n", inputfilename[0] ? inputfilename : "(unknown)");
            fprintf(ssf, "# beta=%.6f  eta=%.6f\n", beta_hat, ml_eta);
            fprintf(ssf, "# name\tweight\tlnL\n");
            for(i = 0; i < Total_fund_trees; i++)
                {
                if(!sourcetreetag[i]) continue;
                double w_k = (double)tree_weights[i];
                int    ki  = (fund_bipart_sets != NULL) ? fund_bipart_sets[i].count : 0;
                double d_k = (ml_eta != 0.0 && ki > 0)
                             ? dists_raw[i] / pow((double)ki, ml_eta)
                             : (double)dists_raw[i];
                double lnL_k = w_k * log(beta_hat);
                if(w_k > 0.0) lnL_k += w_k * log(w_k);
                if(ml_eta != 0.0 && ki > 0)
                    lnL_k -= ml_eta * w_k * log((double)ki);
                lnL_k -= beta_hat * w_k * d_k;
                if(strcmp(tree_names[i], "") != 0)
                    fprintf(ssf, "%s\t%.6f\t%.6f\n", tree_names[i], w_k, lnL_k);
                else
                    fprintf(ssf, "%d\t%.6f\t%.6f\n", i+1, w_k, lnL_k);
                }
            fclose(ssf);
            printf2("  Source tree scores written to: %s\n", sourcescores_name);
            }
        }

    free(dists_raw);
    }

void bootstrap_search(void)
    {
    int i=0, j=0, k=0, l=0, random_num = 0, error = FALSE, *taxa_present = NULL, missing_method = 1;
    int Nreps = 100, search = 1, allpresent = TRUE, num_results = 0;
    int nthreads = 1;
    char filename[1000], **bootstrap_results = NULL, consensusfilename[1000];
    FILE *bootfile = NULL, *temp = NULL, *consensusfile = NULL;
	float percentage = .5;
#ifdef _OPENMP
    nthreads = omp_get_num_procs();
#endif

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
                printf2("Error: '%s' is an invalid number of repetitions\n", parsed_command[i+1]);
                error = TRUE;
                }
            }
        if(strcmp(parsed_command[i], "swap") == 0)
            {
            if(strcmp(parsed_command[i+1], "all") == 0) search = 0;
            else
                {
                if(strcmp(parsed_command[i+1], "nni") == 0) { search = 1; method = 1; }
                else
                    {
					if(strcmp(parsed_command[i+1], "spr") == 0) { search = 1; method = 2; }
                    else
                        {
						if(strcmp(parsed_command[i+1], "tbr") == 0) { search = 1; method = 3; }
						else
							{
							printf2("Error: swap option '%s' unknown\n", parsed_command[i+1]);
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
					printf2("Error: '%s' is an invalid method to be assigned to missing\n", parsed_command[i+1]);
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
							printf2("Error: The cut off for a consensus tree must be between .5 and 1.0\n");
							error = TRUE;
							}
						if(percentage == 0)
							{
							printf2("Error: consensus option %s unknown\n", parsed_command[i+1]);
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
                printf2("Error: Unable to open file named '%s'\n", parsed_command[i+1]);
                error = TRUE;
                strcpy(filename, "consensus.ph");
                }
            }
			
        if(strcmp(parsed_command[i], "treefile") == 0)
            {
            strcpy(filename, parsed_command[i+1]);
            if(filename[0] == '\0')
                {
                printf2("Error: Unable to open file named '%s'\n", parsed_command[i+1]);
                error = TRUE;
                strcpy(filename, "bootstrap.txt");
                }
            }
#ifdef _OPENMP
        if(strcmp(parsed_command[i], "nthreads") == 0)
            {
            nthreads = toint(parsed_command[i+1]);
            if(nthreads < 1)
                {
                printf2("Error: '%s' is an invalid value for nthreads\n", parsed_command[i+1]);
                nthreads = 1;
                }
            }
#endif
        }


    if(!error)
        {

        printf2("\n\nBootstrap Settings:\n\tNumber of Bootstrap replicates = %d\n\tSearching Supertree-space:", Nreps);
        if(search==0)printf2(" Exhaustive search\n");
        if(search==1)printf2(" Heuristic search ");
	
        printf2("\tBootstrap tree file: %s\n", filename);
		if(percentage == 1)
			printf2("\tStrict ");
		if(percentage == .5)
			printf2("\tMajority rule ");
		if(percentage < 0.5)
			printf2("\tMajority rule with minor components ");
		if(percentage > 0.5 && percentage != 1)
			printf2("%f cut-off ", percentage);
		printf2("consensus tree is to be constructed\n");
		printf2("\tConsensus file name = %s\n\n\n", consensusfilename);
		
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

		
        
        if(criterion == 0 || criterion == 2 || criterion == 3 || criterion == 6 || criterion == 7)
            {
            /******* Parallel bootstrap over replicates *******/
            hsprint = FALSE;  /* suppress per-replicate HS settings banner during bootstrap */

            /* fund_bipart_sets is now threadprivate: each bootstrap thread builds its own
             * copy from its resampled fundamentals via rf_precompute_fund_biparts().
             * Parallel bootstrap is therefore safe for all criteria. */

#ifdef _OPENMP
            if(nthreads > 1)
                printf2("Bootstrap: %d replicates using %d threads\n", Nreps, nthreads);
            else
#endif
                printf2("Bootstrap progress indicator:\n");

            /* Shared accumulators — protected by boot_merge critical section */
            {
            int    boot_reps_done = 0;
            int    boot_reps_completed = 0;  /* actual replicate count, updated once per rep */
            int    boot_reps_running = 0;   /* reps currently executing (across all threads) */
            char **boot_results_acc = malloc(1*sizeof(char*));
            float *boot_scores_acc  = malloc(1*sizeof(float));
            boot_results_acc[0] = malloc(TREE_LENGTH*sizeof(char));
            boot_results_acc[0][0] = '\0';
            int    boot_alloc_n = 1;

            /* Save original threadprivate pointers before the parallel region.
             * boot_alloc_thread_state() will overwrite the master thread's
             * threadprivate copies; boot_free_thread_state() will free the
             * per-bootstrap allocations but leaves dangling pointers.  We
             * restore the originals immediately after the parallel region so
             * that the restore-original-data block and clean_exit() see valid
             * pointers to the data loaded by exe. */
            char **boot_save_fundamentals    = fundamentals;
            int  **boot_save_presence        = presence_of_taxa;
            int ***boot_save_fund_scores     = fund_scores;
            int   *boot_save_num_comparisons = number_of_comparisons;
            BipartSet *boot_save_fund_biparts = fund_bipart_sets;

#ifdef _OPENMP
            #pragma omp parallel num_threads(nthreads) default(shared) \
                    private(i,j,k,l,random_num,allpresent)
#endif
                {
                /* Per-thread local result list */
                int    thr_n     = 0;
                int    thr_alloc = 1;
                char **thr_results = malloc(1*sizeof(char*));
                float *thr_scores  = malloc(1*sizeof(float));
                thr_results[0] = malloc(TREE_LENGTH*sizeof(char));
                thr_results[0][0] = '\0';
                int *thr_taxa_present = malloc(number_of_taxa*sizeof(int));

#ifdef _OPENMP
                boot_alloc_thread_state(Total_fund_trees, number_of_taxa);
                #pragma omp for schedule(dynamic,1)
#endif
                for(i=0; i<Nreps; i++)
                    {
                    if(user_break) continue;

                    /* ---- Resample fundamentals into thread-private copies ---- */
                    for(j=0; j<Total_fund_trees; j++)
                        {
#ifdef _OPENMP
                        random_num = (int)fmod((double)rand_r(&thread_seed), (double)Total_fund_trees);
#else
                        random_num = (int)fmod((double)rand(), (double)Total_fund_trees);
#endif
                        fundamentals[j] = realloc(fundamentals[j],
                                           (strlen(stored_funds[random_num])+100)*sizeof(char));
                        fundamentals[j][0] = '\0';
                        strcpy(fundamentals[j], stored_funds[random_num]);
                        number_of_comparisons[j] = stored_num_comparisons[random_num];
                        for(k=0; k<number_of_taxa; k++)
                            {
                            presence_of_taxa[j][k] = stored_presence_of_taxa[random_num][k];
                            for(l=0; l<number_of_taxa; l++)
                                fund_scores[j][k][l] = stored_fund_scores[random_num][k][l];
                            }
                        }

                    /* ---- Refresh RF/ML bipartitions from resampled source trees ---- */
                    /* fund_bipart_sets is precomputed from the original trees; after
                     * resampling we must rebuild it so compare_trees_rf/ml score
                     * against the bootstrap replicate, not the original data. */
                    if(criterion == 2 || criterion == 3 || criterion == 6 || criterion == 7) rf_precompute_fund_biparts();

                    /* ---- Check all taxa present ---- */
                    allpresent = TRUE;
                    for(j=0; j<number_of_taxa; j++) thr_taxa_present[j] = FALSE;
                    for(j=0; j<Total_fund_trees; j++)
                        for(k=0; k<number_of_taxa; k++)
                            if(presence_of_taxa[j][k] > 0) thr_taxa_present[k] = TRUE;
                    for(j=0; j<number_of_taxa; j++)
                        if(thr_taxa_present[j] == FALSE) { allpresent = FALSE; break; }

                    if(!allpresent)
                        {
#ifdef _OPENMP
                        #pragma omp critical (boot_merge)
#endif
                            {
                            printf2("\n\trepetition %d: taxon missing — skipped\n", i+1);
                            boot_reps_completed++;
                            printf2("\r\t%d running, %d done / %d",
                                    boot_reps_running, boot_reps_completed, Nreps);
                            fflush(stdout);
                            }
                        continue;
                        }

                    /* ---- Run the search ---- */
#ifdef _OPENMP
                    #pragma omp critical (boot_merge)
#endif
                        {
                        boot_reps_running++;
                        printf2("\r\t%d running, %d done / %d", boot_reps_running, boot_reps_completed, Nreps);
                        fflush(stdout);
                        }
                    if(search == 0)
                        alltrees_search(FALSE);
                    else
                        heuristic_search(FALSE, FALSE, 10000, 1);

                    /* ---- Collect results from (threadprivate) retained_supers ---- */
                    {
                    int rep_l = 0;
                    char *local_best = NULL;  /* heap-allocated to avoid stack overflow in parallel region */
                    while(scores_retained_supers[rep_l] != -1) rep_l++;
                    if(rep_l == 0) goto next_replicate;
                    local_best = malloc(TREE_LENGTH * sizeof(char));
                    if(!local_best) memory_error(203);
                    local_best[0] = '\0';

                    /* Grow thread-local accumulators */
                    thr_results = realloc(thr_results, (thr_n+rep_l)*sizeof(char*));
                    thr_scores  = realloc(thr_scores,  (thr_n+rep_l)*sizeof(float));
                    for(j=thr_alloc; j<thr_n+rep_l; j++)
                        {
                        thr_results[j] = malloc(TREE_LENGTH*sizeof(char));
                        thr_results[j][0] = '\0';
                        }
                    if(thr_n+rep_l > thr_alloc) thr_alloc = thr_n+rep_l;

                    for(k=0; k<rep_l; k++)
                        {
                        /* Build printable tree string */
                        if(search == 0)
                            {
                            if(tree_top != NULL) { dismantle_tree(tree_top); tree_top = NULL; }
                            temp_top = NULL;
                            { int _to = 0; tree_build(1, retained_supers[k], tree_top, FALSE, -1, &_to); }
                            tree_top = temp_top; temp_top = NULL;
                            local_best[0] = '\0';
                            print_named_tree(tree_top, local_best);
                            }
                        else
                            {
                            int jj=0;
                            local_best[0] = '\0';
                            while(retained_supers[k][jj] != ';')
                                { local_best[jj] = retained_supers[k][jj]; jj++; }
                            local_best[jj] = '\0';
                            }
                        strcpy(thr_results[thr_n+k], retained_supers[k]);
                        thr_scores[thr_n+k] = (float)((float)1/(float)rep_l);
#ifdef _OPENMP
                        #pragma omp critical (boot_merge)
#endif
                            {
                            fprintf(bootfile, "%s [%f];\t[score = %f]\n",
                                    local_best, (float)((float)1/(float)rep_l),
                                    scores_retained_supers[k]);
                            fflush(bootfile);
                            boot_reps_done++;
                            fflush(stdout);
                            }
                        }
                    free(local_best);
                    thr_n += rep_l;

                    /* Reset search results for next replicate */
                    for(k=0; k<number_retained_supers; k++)
                        {
                        retained_supers[k][0] = '\0';
                        scores_retained_supers[k] = -1;
                        }
                    } /* end collect */
                    next_replicate: ;
#ifdef _OPENMP
                    #pragma omp critical (boot_merge)
#endif
                        {
                        boot_reps_running--;
                        boot_reps_completed++;
                        printf2("\r\t%d running, %d done / %d",
                                boot_reps_running, boot_reps_completed, Nreps);
                        fflush(stdout);
                        }
                    } /* end loop body (i) */

#ifdef _OPENMP
                #pragma omp critical (boot_merge)
#endif
                    {
                    /* Merge thread-local results into shared accumulator.
                     * Guard against realloc(ptr,0) which frees the buffer on macOS. */
                    if(thr_n > 0)
                        {
                        boot_results_acc = realloc(boot_results_acc,
                                                   (num_results+thr_n)*sizeof(char*));
                        boot_scores_acc  = realloc(boot_scores_acc,
                                                   (num_results+thr_n)*sizeof(float));
                        for(j=boot_alloc_n; j<num_results+thr_n; j++)
                            {
                            boot_results_acc[j] = malloc(TREE_LENGTH*sizeof(char));
                            boot_results_acc[j][0] = '\0';
                            }
                        if(num_results+thr_n > boot_alloc_n)
                            boot_alloc_n = num_results+thr_n;
                        for(j=0; j<thr_n; j++)
                            {
                            strcpy(boot_results_acc[num_results+j], thr_results[j]);
                            boot_scores_acc[num_results+j] = thr_scores[j];
                            }
                        num_results += thr_n;
                        }
                    }

#ifdef _OPENMP
                boot_free_thread_state(Total_fund_trees, number_of_taxa);
#endif
                for(j=0; j<thr_alloc; j++) free(thr_results[j]);
                free(thr_results);
                free(thr_scores);
                free(thr_taxa_present);
                } /* end parallel region */

            printf2("\n");   /* move off the \r progress line before consensus output */
            hsprint = TRUE;  /* restore after bootstrap */

            /* Restore original threadprivate pointers (boot copies already freed) */
            fundamentals          = boot_save_fundamentals;
            presence_of_taxa      = boot_save_presence;
            fund_scores           = boot_save_fund_scores;
            number_of_comparisons = boot_save_num_comparisons;
            fund_bipart_sets      = boot_save_fund_biparts;

            /* Install merged results into the globals consensus() expects */
            if(bootstrap_results)
                {
                for(j=0; j<1; j++) free(bootstrap_results[j]);
                free(bootstrap_results);
                free(score_of_bootstraps);
                }
            bootstrap_results   = boot_results_acc;
            score_of_bootstraps = boot_scores_acc;

            printf2("\n");
            fclose(bootfile);
            } /* end accumulators block */
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
                        printf2("\n\trepetition %d:\tTaxon %s is not present\n", i+1, taxa_names[j]);
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
							printf2("error\n");
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
                    printf2("\n\nError: This data may not be suitable for bootstrapping due to the low\noccurrences of some taxa in the source trees. Please check the\nCo-occurance summary to identify problematic taxa\n");
                    error = TRUE;
		    }
                }
            fclose(BR_file);
            if(!error && allpresent)
                {
                if(print_log == FALSE) 
					strcpy(system_call, "paup coding.nex");
				else
					{
					strcpy(system_call, "paup coding.nex");
					strcat(system_call, " | tee -a ");
					strcat(system_call, logfile_name);
					}

				if(system(system_call) != 0) printf2("Error calling paup, please execute the file coding.nex in paup to perform the parsimony step\n");

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
    printf2("Error: Out of memory @ %d\n", error_num);
    clean_exit(1);
    }















/* --- ML topology test helpers ------------------------------------------- */

/* Exact two-sided binomial p-value using log-space arithmetic (avoids overflow). */
static double ml_binomial_p2(int n_plus, int n_minus)
    {
    int n = n_plus + n_minus;
    if(n == 0) return 1.0;
    int k_obs = (n_plus < n_minus) ? n_plus : n_minus;
    double lh = -n * log(2.0);
    double p = 0.0, log_binom = 0.0;
    for(int k = 0; k <= k_obs; k++)
        {
        if(k > 0) log_binom += log((double)(n-k+1)) - log((double)k);
        p += exp(log_binom + lh);
        }
    p *= 2.0;
    return (p > 1.0) ? 1.0 : p;
    }

/* One-sided normal p-value: P(Z > z). */
static double ml_normal_p(double z)
    { return erfc(z / sqrt(2.0)) / 2.0; }

/* Significance asterisks. */
static const char *ml_sig_label(double p, float alpha)
    {
    if(p < 0.001) return "***";
    if(p < 0.01)  return "**";
    if(p < (double)alpha) return "*";
    return "ns";
    }

/* Run Winning-Sites, KH and SH tests comparing each candidate tree to the best.
 * all_gene_scores[t][k] = sourcetree_scores[k] for candidate tree t.
 * all_total_scores[t]   = total ML score (positive, minimised) for tree t.
 */
static void run_usertrees_ml_tests(
        float **all_gene_scores, float *all_total_scores,
        int n_all, int nboot, float alpha, const char *testsfile_name,
        int nthreads, char **tree_names, int tree_names_alloc)
    {
    int i, b, k;
    /* Helper: build a display label for tree i (1-based: i+1) */
    #define TREE_LABEL_BUF_SZ 32
    char _tlabel[TREE_LABEL_BUF_SZ];
    #define TREE_LABEL(idx) ( \
        (tree_names && (idx) < tree_names_alloc && tree_names[(idx)]) \
            ? (snprintf(_tlabel, TREE_LABEL_BUF_SZ, "%d:%.*s", (idx)+1, 20, tree_names[(idx)]), _tlabel) \
            : (snprintf(_tlabel, TREE_LABEL_BUF_SZ, "T%d", (idx)+1), _tlabel) \
        )

    /* Find best (minimum score = maximum lnL). */
    int best_idx = 0;
    for(i = 1; i < n_all; i++)
        if(all_total_scores[i] < all_total_scores[best_idx]) best_idx = i;
    float best_score_ml = all_total_scores[best_idx];

    printf2("\nML topology tests  (best = %s,  lnL = %.4f,  beta = %.2f)\n",
            TREE_LABEL(best_idx), ml_display_score(best_score_ml), ml_beta);
    printf2("=============================================================================\n");
    printf2("  %-9s  %-10s  %-12s  %-14s  %-12s\n",
            "Tree", "lnL", "WinSites p", "KH z / p", "SH p");
    printf2("-----------------------------------------------------------------------------\n");

    FILE *tfile = fopen(testsfile_name, "w");
    if(tfile)
        {
        fprintf(tfile, "# ML supertree topology tests\n");
        fprintf(tfile, "# Best tree: %s  lnL=%.6f  beta=%.4f\n",
                TREE_LABEL(best_idx), (double)ml_display_score(best_score_ml), ml_beta);
        fprintf(tfile, "# tree\tdelta_lnL\tn_inform\tn_plus\tn_minus\tWS_p\tKH_z\tKH_p\tSH_p\n");
        }

    for(i = 0; i < n_all; i++)
        {
        if(i == best_idx) continue;

        int n_plus = 0, n_minus = 0, n_inform = 0;
        double delta_total = 0.0;
        float *delta = malloc(Total_fund_trees * sizeof(float));

        for(k = 0; k < Total_fund_trees; k++)
            {
            if(!sourcetreetag[k]) { delta[k] = 0.0f; continue; }
            float d = all_gene_scores[i][k] - all_gene_scores[best_idx][k];
            delta[k] = d;
            if(all_gene_scores[best_idx][k] == 0.0f && all_gene_scores[i][k] == 0.0f)
                continue;  /* uninformative for both trees */
            n_inform++;
            delta_total += d;
            if(d > 0.0f) n_plus++;
            else if(d < 0.0f) n_minus++;
            }

        if(n_inform < 4)
            printf2("  WARNING: %s has only %d informative gene trees — results unreliable\n",
                    TREE_LABEL(i), n_inform);

        /* Winning Sites */
        double ws_p = ml_binomial_p2(n_plus, n_minus);

        /* KH parametric */
        double kh_z = 0.0, kh_p = 0.5;
        if(n_inform > 1)
            {
            double c_bar = delta_total / n_inform;
            double var = 0.0;
            for(k = 0; k < Total_fund_trees; k++)
                {
                if(!sourcetreetag[k]) continue;
                if(all_gene_scores[best_idx][k]==0.0f && all_gene_scores[i][k]==0.0f) continue;
                double diff = delta[k] - c_bar;
                var += diff * diff;
                }
            var /= (n_inform - 1);
            double se = sqrt((double)n_inform * var);
            if(se > 0.0) { kh_z = delta_total / se; kh_p = ml_normal_p(kh_z); }
            }

        /* SH bootstrap */
        double sh_p = 1.0;
        if(n_inform > 1 && nboot > 0)
            {
            double c_bar = delta_total / n_inform;
            float *c_star = malloc(Total_fund_trees * sizeof(float));
            int *inform_idx = malloc(n_inform * sizeof(int));
            int n_inf2 = 0;
            for(k = 0; k < Total_fund_trees; k++)
                {
                if(!sourcetreetag[k]) continue;
                if(all_gene_scores[best_idx][k]==0.0f && all_gene_scores[i][k]==0.0f) continue;
                c_star[k] = (float)(delta[k] - c_bar);
                inform_idx[n_inf2++] = k;
                }
            int exceed = 0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:exceed) num_threads(nthreads) schedule(static)
#endif
            for(b = 0; b < nboot; b++)
                {
                /* per-replicate seed keeps threads independent and avoids
                 * the non-thread-safe rand(); result is statistically
                 * equivalent to a single stream for bootstrap purposes */
                unsigned int local_seed = (unsigned int)((b + 1) * 1000003u);
                double boot_sum = 0.0;
                int s;
                for(s = 0; s < n_inf2; s++)
                    boot_sum += c_star[inform_idx[rand_r(&local_seed) % n_inf2]];
                if(boot_sum >= delta_total) exceed++;
                }
            sh_p = (double)exceed / nboot;
            free(c_star);
            free(inform_idx);
            }
        else if(nboot == 0)
            sh_p = -1.0;  /* sentinel: skipped */

        double lnL_i     = (double)ml_display_score(all_total_scores[i]);
        double delta_lnL = lnL_i - (double)ml_display_score(all_total_scores[best_idx]);

        printf2("  %-9s  lnL=%-6.3f  p=%-8.4f  z=%-5.2f p=%-5.3f %-3s  p=%-6.4f %-3s\n",
                TREE_LABEL(i), lnL_i,
                ws_p,
                kh_z, kh_p, ml_sig_label(kh_p, alpha),
                sh_p >= 0 ? sh_p : 9.9999, sh_p >= 0 ? ml_sig_label(sh_p, alpha) : "--");

        if(tfile)
            fprintf(tfile, "%s\t%.6f\t%d\t%d\t%d\t%.6f\t%.4f\t%.6f\t%.6f\n",
                    TREE_LABEL(i), delta_lnL,
                    n_inform, n_plus, n_minus,
                    ws_p, kh_z, kh_p, sh_p >= 0 ? sh_p : -1.0);

        free(delta);
        }

    printf2("-----------------------------------------------------------------------------\n");
    printf2("  %s p<alpha=%g; ** p<0.01; *** p<0.001; ns=not significant\n", "*", (double)alpha);
    if(nboot > 0)
        printf2("  SH test: %d bootstrap replicates\n", nboot);
    else
        printf2("  SH test: skipped (nboot=0)\n");
    printf2("  AU test: not yet implemented\n");
    printf2("  Ref: Steel & Rodrigo (2008) Syst Biol 57:243; Kishino & Hasegawa (1989)\n");
    printf2("       J Mol Evol 29:170; Shimodaira & Hasegawa (1999) Mol Biol Evol 16:1114\n");
    if(tfile)
        {
        printf2("  Per-gene-tree details written to: %s\n", testsfile_name);
        fclose(tfile);
        }
    #undef TREE_LABEL
    #undef TREE_LABEL_BUF_SZ
    }

/* --- end ML topology test helpers --------------------------------------- */


void usertrees_search(void)
    {
    FILE *userfile = NULL, *outfile = NULL, *sourcescoresfile = NULL;
    int keep = 0, nbest = 0, error = FALSE, i=0, tree_number = 0, j=0, k=0, prev = 0, print_source_scores = FALSE;
    int    do_normcorr_local = 0;       /* 0=off, 1=approx, 2=exact (normcorrect option) */
    char  **user_tree_names     = NULL;   /* [tree_number]: name from [..] bracket for each input tree */
    int     user_tree_names_alloc = 0;
    int    *retained_orig_idx   = NULL;   /* [j]: 0-based input tree index for retained_supers[j] */
    char    utn_buf[NAME_LENGTH + 2];     /* temp buffer for capturing [name] from input */
    /* ML topology test state */
    int    run_ml_tests = FALSE, nboot_tests = 1000;
    float  alpha_tests  = 0.05f;
#ifdef _OPENMP
    int    nthreads_tests = omp_get_num_procs();
#else
    int    nthreads_tests = 1;
#endif
    char   testsfile_name[10000];
    int    n_all = 0, n_all_alloc = 0;
    float **all_gene_scores = NULL;
    float  *all_total_scores = NULL;
    strcpy(testsfile_name, "mltest_results.txt");
    /* Score matrix state */
    int    do_scorematrix    = FALSE;
    char   scorematrixfile_name[10000];
    float **matrix_scores    = NULL;   /* [tree_idx][Total_fund_trees] per-source scores */
    float  *matrix_total     = NULL;   /* [tree_idx] total score per input tree           */
    int    matrix_alloc      = 0;
    strcpy(scorematrixfile_name, "usertrees_scorematrix.txt");
    char *user_super = NULL, c = '\0', *best_tree = malloc(TREE_LENGTH * sizeof(char)), *temp = NULL;
    if(!best_tree) { printf2("Error: out of memory in usertrees_search\n"); return; }
    float score = 0, best_score = 0;
    
    if((userfile = fopen(parsed_command[1], "r")) == NULL)
        {
        printf2("Error opening file named %s\n", parsed_command[1]);
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
            if(strcmp(parsed_command[i], "scorematrix") == 0)
                {
                if(strcmp(parsed_command[i+1], "yes") == 0)
                    do_scorematrix = TRUE;
                }
            if(strcmp(parsed_command[i], "scorematrixfile") == 0)
                strncpy(scorematrixfile_name, parsed_command[i+1], 9999);
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
                                printf2("Error: weight option '%s' is unknown\n", parsed_command[i+1]);
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
                            printf2("Error: weight option '%s' is unknown\n", parsed_command[i+1]);
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
                            printf2("Error: weight option '%s' is unknown\n", parsed_command[i+1]);
                            error = TRUE;
                            }
                        }
                    }
                }            

            if(strcmp(parsed_command[i], "outfile") == 0)
                {
                if((outfile = fopen(parsed_command[i+1], "w")) == NULL)
                    {
                    printf2("Error opening output file named: %s\n", parsed_command[i+1]);
                    error = TRUE;
                    }
                }

            if(strcmp(parsed_command[i], "tests") == 0)
                {
                if(strcmp(parsed_command[i+1], "yes") == 0)
                    run_ml_tests = TRUE;
                }
            if(strcmp(parsed_command[i], "normcorrect") == 0)
                {
                if(parsed_command[i+1] && strcmp(parsed_command[i+1], "exact") == 0)
                    do_normcorr_local = 2;
                else
                    do_normcorr_local = 1;
                }
            if(strcmp(parsed_command[i], "nboot") == 0)
                nboot_tests = atoi(parsed_command[i+1]);
            if(strcmp(parsed_command[i], "testsfile") == 0)
                strcpy(testsfile_name, parsed_command[i+1]);
            if(strcmp(parsed_command[i], "nthreads") == 0)
                {
                nthreads_tests = atoi(parsed_command[i+1]);
                if(nthreads_tests < 1) nthreads_tests = 1;
                }

            }
        if(outfile == NULL)
            {
            if((outfile = fopen("Usertrees_result.txt", "w")) == NULL)
                {
                printf2("Error opening output file  Usertrees_result.txt\n");
                error = TRUE;
                }
            }
        }
    if(criterion == 1 || criterion == 4 || criterion == 5)
		{
		printf2("ERROR: Usertree search is not available under the criteria ");
		switch(criterion)
			{
			case 1:
				printf2("Matrix Representation Using Parsimony (MRP)\n");
				break;
			case 4:
				printf2("Average consensus (AVCON)\n");
				break;
			default:
				printf2("Reconstruction of duplications and losses (RECON)\n");
				break;
			}
		error = TRUE;
		}
    if(!error)
        {   
		
		printf2("\nUsertree Search settings:\n");
		printf2("\tCriterion = ");
		switch(criterion)
			{
			case 1:
				printf2("Matrix Representation Using Parsimony (MRP)\n");
				break;
			case 0:
				printf2("Most Similar Supertree (dfit)\n");
				break;
			case 2:
				printf2("Maximum split fit (SFIT)\n");
				break;
			case 4:
				printf2("Average consensus (AVCON)\n");
				break;
			case 5:
				printf2("Reconstruction of duplications and losses (RECON)\n");
				break;
			case 6:
				printf2("Robinson-Foulds distance (RF)\n");
				break;
			case 7:
				printf2("Maximum likelihood supertree (beta=%.2f, scale=%s)\n",
						ml_beta, ml_scale==1?"lust":ml_scale==2?"lnl":"Steel & Rodrigo 2008");
				break;
			default:
				printf2("Maximum quartet fit (QFIT)\n");
				break;

			}
		if(criterion != 1 && criterion != 4)
			{
			if(criterion != 5)
				printf2("\tWeighting Scheme = ");
			if(criterion==0)
				{
				if(dweight == 0) printf2("equal\n");
				if(dweight == 1) printf2("comparisons\n");
				}
			if(criterion == 2)
				{
				if(splits_weight == 1) printf2("equal\n");
				if(splits_weight == 2) printf2("splits\n");
				}
			if(criterion == 3)
				{
				if(quartet_normalising == 1) printf2("equal\n");
				if(quartet_normalising == 2) printf2("taxa\n");
				if(quartet_normalising == 3) printf2("quartets\n");
				}
			if(criterion == 2 || criterion == 3)
				printf2("\tBS split weights (bsweight) = %s\n", bsweight ? "on" : "off");
			printf2("\tUser defined supertrees from file: %s\n", parsed_command[1]);
			
			}
		
		if(criterion==5) printf2("\n\tDuplication weight = %f\n\tLosses weight = %f\n\n", dup_weight, loss_weight);
		
		if(print_source_scores)printf2("\nScores of each source tree compared to the best user supertree to be printed to sourcetree_scores.txt\n");
		
		
		
		
		
		
        /* this hold the distances calculated with the pathmetric on the supertree */
        if(!calculated_fund_scores && criterion == 0)
			{
			cal_fund_scores(FALSE);    /* calculate the path metrics for all the fundamental trees */
			calculated_fund_scores = TRUE;
			}
        if(criterion == 6 || criterion == 7)
            rf_precompute_fund_biparts();
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

		/* Allocate name-tracking arrays */
		user_tree_names_alloc = 64;
		user_tree_names = calloc(user_tree_names_alloc, sizeof(char *));
		retained_orig_idx = malloc(number_retained_supers * sizeof(int));
		if(retained_orig_idx) { for(i = 0; i < number_retained_supers; i++) retained_orig_idx[i] = -1; }

        if(criterion == 2 || criterion == 3 || criterion == 6 || criterion == 7)
            rf_precompute_fund_biparts();

        /* Enable Bryant & Steel normalising constant correction if requested */
        if(do_normcorr_local && criterion == 7 && ml_scale == 2)
            {
            ml_norm_logZ = calloc(Total_fund_trees, sizeof(double));
            if(ml_norm_logZ)
                {
                ml_do_normcorr = do_normcorr_local;
                if(do_normcorr_local == 2)
                    printf2("  Normalising constant correction: exact subset enumeration (n<=20; fallback to approx for larger pruned trees)\n");
                else
                    printf2("  Normalising constant correction: enabled (Bryant & Steel 2008, large-beta approx)\n");
                }
            else printf2("  Warning: normcorrect allocation failed; correction disabled\n");
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
			{ int _to = 0; tree_build(1, user_super, tree_top, TRUE, -1, &_to); }
			tree_top = temp_top;
			temp_top = NULL;
			/*check_tree(tree_top); */
			if(criterion == 0)      score = compare_trees(FALSE);
			if(criterion == 2)      score = compare_trees_sfit(FALSE);
			if(criterion == 3)      score = compare_trees_qfit(FALSE);
			if(criterion == 6)      score = compare_trees_rf(FALSE);
			if(criterion == 7)      score = compare_trees_ml(FALSE);

			/* Collect per-gene-tree scores for ML topology tests */
			if(run_ml_tests && criterion == 7)
				{
				if(n_all == n_all_alloc)
					{
					n_all_alloc += 64;
					all_gene_scores  = realloc(all_gene_scores,  n_all_alloc * sizeof(float *));
					all_total_scores = realloc(all_total_scores, n_all_alloc * sizeof(float));
					}
				all_gene_scores[n_all] = malloc(Total_fund_trees * sizeof(float));
				for(k = 0; k < Total_fund_trees; k++)
					all_gene_scores[n_all][k] = sourcetree_scores[k];
				all_total_scores[n_all] = score;
				n_all++;
				}

			/*printf("%s\t[%f]\n", user_super, score); */

			/* Capture per-source scores for score matrix (all criteria) */
			if(do_scorematrix && sourcetree_scores)
				{
				if(tree_number == matrix_alloc)
					{
					matrix_alloc += 64;
					matrix_scores = realloc(matrix_scores, matrix_alloc * sizeof(float *));
					matrix_total  = realloc(matrix_total,  matrix_alloc * sizeof(float));
					}
				if(matrix_scores && matrix_total)
					{
					matrix_scores[tree_number] = malloc(Total_fund_trees * sizeof(float));
					if(matrix_scores[tree_number])
						for(k = 0; k < Total_fund_trees; k++)
							matrix_scores[tree_number][k] = sourcetree_scores[k];
					matrix_total[tree_number] = score;
					}
				}

			c = getc(userfile);
			utn_buf[0] = '\0';
			while((c == ' ' || c == '\r'  || c == '\n'|| c == '\t' || c == '[') && !feof(userfile))  /* skip past all the non tree characters */
				{
				if(c == '[')
					{
					int _ni = 0;
					while((c = getc(userfile)) != ']' && !feof(userfile))
						{ if(_ni < NAME_LENGTH) utn_buf[_ni++] = (char)c; }
					utn_buf[_ni] = '\0';
					}
				c = getc(userfile);
				}
			/* Store name for this tree (tree_number is the 0-based index, not yet incremented) */
			if(user_tree_names != NULL)
				{
				if(tree_number >= user_tree_names_alloc)
					{
					int _na = user_tree_names_alloc * 2;
					user_tree_names = realloc(user_tree_names, _na * sizeof(char *));
					if(user_tree_names)
						{ int _ii; for(_ii = user_tree_names_alloc; _ii < _na; _ii++) user_tree_names[_ii] = NULL; }
					user_tree_names_alloc = _na;
					}
				if(user_tree_names && tree_number < user_tree_names_alloc)
					user_tree_names[tree_number] = (utn_buf[0] != '\0') ? strdup(utn_buf) : NULL;
				}

			if(tree_number==0)
				{
				tree_number++;
				retained_supers[0] = realloc(retained_supers[0], (strlen(user_super)+10)*sizeof(char));
				strcpy(retained_supers[0], user_super);
				scores_retained_supers[0] = score;
				best_score = score;
				if(retained_orig_idx) retained_orig_idx[0] = tree_number - 1;
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
					if(retained_orig_idx) retained_orig_idx[0] = tree_number - 1;
					j=1;
					while(scores_retained_supers[j] != -1 && j < number_retained_supers)
						{
						strcpy(retained_supers[j], "");
						scores_retained_supers[j] = -1;
						if(retained_orig_idx) retained_orig_idx[j] = -1;
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
							if(j+1 == number_retained_supers)
								{
								int _old_nrs = number_retained_supers;
								reallocate_retained_supers();
								if(retained_orig_idx)
									{
									retained_orig_idx = realloc(retained_orig_idx, number_retained_supers * sizeof(int));
									int _rr; for(_rr = _old_nrs; _rr < number_retained_supers; _rr++) retained_orig_idx[_rr] = -1;
									}
								}
							}
						retained_supers[j] = realloc(retained_supers[j], (strlen(user_super)+10)*sizeof(char));
						strcpy(retained_supers[j], user_super);
						scores_retained_supers[j] = score;
						if(retained_orig_idx) retained_orig_idx[j] = tree_number - 1;
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

			printf2("\n%d input tree(s) evaluated:\n", tree_number);

			while(scores_retained_supers[i] != -1)
				{
				/*if(print_source_scores && criterion==0) */
				if(print_source_scores)
					{
					fprintf(sourcescoresfile, "Scores of sources trees compared to best User Supertree %d\n\n", i+1);
					if(tree_top != NULL)
						{
						dismantle_tree(tree_top);
						tree_top = NULL;
						}
					temp_top = NULL;

					temp_top = NULL;
					{ int _to = 0; tree_build(1, retained_supers[i], tree_top, TRUE, -1, &_to); }
					tree_top = temp_top;
					temp_top = NULL;
					/*check_tree(tree_top); */
					if(criterion == 0) score = compare_trees(FALSE);
					if(criterion == 6) score = compare_trees_rf(FALSE);
					if(criterion == 7) score = compare_trees_ml(FALSE);
					for(k=0; k<Total_fund_trees; k++)
						{
						if(strcmp(tree_names[k], "") != 0)
							fprintf(sourcescoresfile, "sourcetree %s:\t", tree_names[k]);
						else
							fprintf(sourcescoresfile, "sourcetree %d:\t", k);
						fprintf(sourcescoresfile, "%f\n",
							ml_display_source_score(sourcetree_scores[k], k));
						}
					}

				if(outfile != NULL) fprintf(outfile, "%s\t[%f]\n", retained_supers[i], ml_display_score(scores_retained_supers[i]) );
				tree_coordinates(retained_supers[i], FALSE, TRUE, FALSE, -1);
				{
				int _orig = (retained_orig_idx && retained_orig_idx[i] >= 0) ? retained_orig_idx[i] : i;
				const char *_nm = (user_tree_names && _orig < user_tree_names_alloc && user_tree_names[_orig]) ? user_tree_names[_orig] : NULL;
				if(_nm)
					printf2("\nInput tree %d (%s)  %s = %f\n", _orig+1, _nm, ml_score_label(), ml_display_score(scores_retained_supers[i]) );
				else
					printf2("\nInput tree %d  %s = %f\n", _orig+1, ml_score_label(), ml_display_score(scores_retained_supers[i]) );
				}
				i++;
				}

			trees_in_memory = i;

			/* Write per-source-tree score matrix if requested */
			if(do_scorematrix && matrix_scores && tree_number > 0)
				{
				int _active = 0;
				FILE *_mf = NULL;
				for(k = 0; k < Total_fund_trees; k++)
					if(sourcetreetag[k]) _active++;
				_mf = fopen(scorematrixfile_name, "w");
				if(!_mf)
					printf2("Warning: could not open score matrix file '%s'\n", scorematrixfile_name);
				else
					{
					/* ---- comment header ---- */
					fprintf(_mf, "# usertrees score matrix\n");
					fprintf(_mf, "# criterion=%s  source_trees=%d  input_trees=%d\n",
						ml_score_label(), _active, tree_number);
					fprintf(_mf, "# Total %s per input tree:\n", ml_score_label());
					for(i = 0; i < tree_number; i++)
						{
						const char *_cn = (user_tree_names && i < user_tree_names_alloc
						                   && user_tree_names[i]) ? user_tree_names[i] : NULL;
						if(_cn)
							fprintf(_mf, "#   Tree_%d(%s)\t%s=%g\n",
								i+1, _cn, ml_score_label(),
								(double)ml_display_score(matrix_total[i]));
						else
							fprintf(_mf, "#   Tree_%d\t%s=%g\n",
								i+1, ml_score_label(),
								(double)ml_display_score(matrix_total[i]));
						}
					fprintf(_mf, "#\n");

					/* ---- column header row ---- */
					fprintf(_mf, "SourceTree\tSize\tWeight");
					for(i = 0; i < tree_number; i++)
						{
						const char *_cn = (user_tree_names && i < user_tree_names_alloc
						                   && user_tree_names[i]) ? user_tree_names[i] : NULL;
						if(_cn)
							fprintf(_mf, "\tTree_%d(%s)", i+1, _cn);
						else
							fprintf(_mf, "\tTree_%d", i+1);
						}
					fprintf(_mf, "\n");

					/* ---- data rows: one per active source tree ---- */
					for(k = 0; k < Total_fund_trees; k++)
						{
						if(!sourcetreetag[k]) continue;
						/* row label: number, with name if available */
						if(tree_names && tree_names[k] && strcmp(tree_names[k], "") != 0)
							fprintf(_mf, "%d(%s)", k+1, tree_names[k]);
						else
							fprintf(_mf, "%d", k+1);
						/* tree size: count of taxa present in this source tree */
						{ int _sz = 0, _t;
						  if(presence_of_taxa && presence_of_taxa[k])
						      for(_t = 0; _t < number_of_taxa; _t++)
						          if(presence_of_taxa[k][_t]) _sz++;
						  fprintf(_mf, "\t%d", _sz); }
						/* tree weight */
						fprintf(_mf, "\t%g", (double)(tree_weights ? tree_weights[k] : 1.0f));
						/* score columns */
						for(i = 0; i < tree_number; i++)
							{
							float _v = (matrix_scores[i]) ? matrix_scores[i][k] : -1.0f;
							fprintf(_mf, "\t%g",
								(double)ml_display_source_score(_v, k));
							}
						fprintf(_mf, "\n");
						}
					fclose(_mf);
					printf2("\nScore matrix (%d source trees × %d input trees) written to %s\n",
						_active, tree_number, scorematrixfile_name);
					}
				}

			/* Cleanup score matrix arrays */
			if(matrix_scores)
				{
				for(i = 0; i < tree_number; i++)
					if(matrix_scores[i]) free(matrix_scores[i]);
				free(matrix_scores);
				matrix_scores = NULL;
				}
			if(matrix_total) { free(matrix_total); matrix_total = NULL; }

			/* Run ML topology tests if requested */
			if(run_ml_tests)
				{
				if(criterion != 7)
					printf2("WARNING: tests=yes requires criterion=ml; topology tests skipped.\n");
				else if(n_all < 2)
					printf2("NOTE: only %d tree(s) scored; topology tests require at least 2 trees.\n", n_all);
				else
					run_usertrees_ml_tests(all_gene_scores, all_total_scores,
										   n_all, nboot_tests, alpha_tests, testsfile_name,
										   nthreads_tests, user_tree_names, user_tree_names_alloc);
				}

			/* Cleanup ML test arrays */
			if(all_gene_scores != NULL)
				{
				for(i = 0; i < n_all; i++)
					free(all_gene_scores[i]);
				free(all_gene_scores);
				}
			if(all_total_scores != NULL) free(all_total_scores);

			/* Free name-tracking arrays */
			if(user_tree_names)
				{
				int _fi;
				for(_fi = 0; _fi < user_tree_names_alloc; _fi++)
					if(user_tree_names[_fi]) free(user_tree_names[_fi]);
				free(user_tree_names);
				user_tree_names = NULL;
				}
			if(retained_orig_idx) { free(retained_orig_idx); retained_orig_idx = NULL; }

            fclose(userfile);
            if(outfile != NULL) fclose(outfile);
            free(temp);
            free(user_super);
            fclose(psfile);
			}
	if(sourcescoresfile != NULL) fclose(sourcescoresfile);
	free(best_tree);
    /* Clean up normalising constant correction state */
    ml_do_normcorr = 0;
    if(ml_norm_logZ) { free(ml_norm_logZ); ml_norm_logZ = NULL; }
    }




/* ---------- Single-copy auto-filter helpers ----------
   When delimiter mode is on, supertree searches (hs/nj/alltrees) should
   only score against single-copy gene trees.  These two helpers save and
   restore sourcetreetag[] around a search so the change is invisible to
   all other commands (reconstruct, showtrees, etc.).                     */

int * apply_singlecopy_filter(void)
	{
	int i, j, multicopy, n_single = 0, n_multi = 0;
	int *saved = NULL;

	if(!delimiter) return NULL;   /* only filter when delimiter mode is on */

	saved = malloc(Total_fund_trees * sizeof(int));
	memcpy(saved, sourcetreetag, Total_fund_trees * sizeof(int));

	for(i=0; i<Total_fund_trees; i++)
		{
		if(sourcetreetag[i])
			{
			multicopy = FALSE;
			for(j=0; j<number_of_taxa; j++)
				{
				if(presence_of_taxa[i][j] > 1) { multicopy = TRUE; break; }
				}
			if(multicopy) { sourcetreetag[i] = FALSE; n_multi++; }
			else n_single++;
			}
		}

	if(n_multi > 0)
		printf2("Single-copy filter: using %d single-copy trees for supertree search "
		        "(%d multicopy trees reserved for 'reconstruct').\n", n_single, n_multi);
	else
		{
		/* All active trees are already single-copy -- nothing to do */
		free(saved);
		return NULL;
		}
	return saved;
	}

void restore_singlecopy_filter(int *saved)
	{
	if(saved == NULL) return;
	memcpy(sourcetreetag, saved, Total_fund_trees * sizeof(int));
	free(saved);
	}
/* ------------------------------------------------------ */

/* vs_create, vs_free, vs_clear, vs_contains, vs_insert,
 * lm_create, lm_free, lm_record, lm_merge, lm_write,
 * sth_aux, tree_topo_hash: moved to topology.c */


/* --- ML display helpers --------------------------------------------------
 * ml_display_score: convert internal minimisation score to lnL at the
 * globally fixed ml_beta and ml_eta.
 *
 *   Internal score s = ml_beta * WD(eta)  (positive, minimised during search).
 *   WD = s / ml_beta
 *
 *   lnL at fixed beta:
 *     lnL = W * log(beta) - beta * WD - eta * P
 *         = W * log(ml_beta) - s - eta * P
 *
 *   This value changes when ml_beta or ml_eta change, so it reflects the
 *   parameters shown in the run header.  When ml_beta = beta_hat (MLE), this
 *   equals the profile MLE lnL.
 *
 *   For paper/lust scales the score is returned as-is (no transformation).
 *
 * ml_score_label:  returns column header "lnL" or "score".
 */
static float ml_display_score(float s)
    {
    if(!(criterion == 7 && ml_scale == 2)) return s;
    double W = 0.0, penalty = 0.0, w_entropy = 0.0;
    int _i;
    for(_i = 0; _i < Total_fund_trees; _i++)
        if(sourcetreetag[_i])
            {
            double w = (double)tree_weights[_i];
            W += w;
            if(w > 0.0) w_entropy += w * log(w);
            }
    if(ml_eta != 0.0 && fund_bipart_sets != NULL)
        for(_i = 0; _i < Total_fund_trees; _i++)
            if(sourcetreetag[_i] && fund_bipart_sets[_i].count > 0)
                penalty += log((double)fund_bipart_sets[_i].count);
    /* lnL at fixed ml_beta: W*log(ml_beta) + sum(w_i*log(w_i)) - s - eta*P
     * The sum(w_i*log(w_i)) term accounts for individual beta_i = ml_beta*w_i.
     * When all weights are 1 this term is zero. */
    double lnL;
    double sd = (double)s;
    if(ml_beta > 0.0 && sd >= 0.0)
        lnL = W * log((double)ml_beta) + w_entropy - sd - ml_eta * penalty;
    else
        lnL = w_entropy - sd - ml_eta * penalty;   /* fallback */
    /* Bryant & Steel (2008) normalising constant correction: subtract sum(w_i * log Z_i) */
    if(ml_do_normcorr && ml_norm_logZ != NULL)
        for(_i = 0; _i < Total_fund_trees; _i++)
            if(sourcetreetag[_i] && tree_weights[_i] > 0.0f)
                lnL -= (double)tree_weights[_i] * ml_norm_logZ[_i];
    return (float)lnL;
    }

static const char *ml_score_label(void)
    { return (criterion == 7 && ml_scale == 2) ? "lnL" : "score"; }

/* ml_display_source_score: convert a single source-tree's raw ML score to
 * its individual lnL contribution.
 *
 *   sourcetree_scores[k] = ml_beta * w_k * d_k * k_k^{-eta}
 *                        = the beta*d term for tree k (positive, minimised).
 *
 *   The per-tree lnL contribution from the exponential model is:
 *     lnL_k = w_k * log(beta / k_k^eta) - beta * w_k * d_k * k_k^{-eta}
 *           = w_k * (log(beta) - eta*log(k_k)) - sourcetree_scores[k]
 *
 *   tree_idx is the source tree index (used to look up w_k and k_k).
 */
static float ml_display_source_score(float raw_score, int tree_idx)
    {
    if(!(criterion == 7 && ml_scale == 2)) return raw_score;
    if(ml_beta <= 0.0f) return raw_score;
    double w_k  = (double)tree_weights[tree_idx];
    /* full per-tree lnL: w_k*log(beta_i) where beta_i = ml_beta*w_k
     * = w_k*(log(ml_beta) + log(w_k)) - eta*w_k*log(k_k) - raw_score */
    double lnL_k = w_k * log((double)ml_beta);
    if(w_k > 0.0) lnL_k += w_k * log(w_k);
    if(ml_eta != 0.0 && fund_bipart_sets != NULL
       && fund_bipart_sets[tree_idx].count > 0)
        lnL_k -= ml_eta * w_k * log((double)fund_bipart_sets[tree_idx].count);
    lnL_k -= (double)raw_score;
    /* Bryant & Steel (2008) normalising constant correction for this source tree */
    if(ml_do_normcorr && ml_norm_logZ != NULL && tree_idx >= 0 && tree_idx < Total_fund_trees)
        lnL_k -= w_k * ml_norm_logZ[tree_idx];
    return (float)lnL_k;
    }

/* --- end RF/ML functions ------------------------------------------------- */


#ifdef _OPENMP
/*  hs_alloc_thread_state --------------------------------------------------
 *  Allocate (or reset) the per-thread search globals for the current thread.
 *  Must be called from inside a #pragma omp parallel region.
 */
static void hs_alloc_thread_state(void)
    {
    int k, init_n, is_master;
    init_n = 10;
    is_master = (omp_get_thread_num() == 0);

    if(!is_master)
        {
        /* Non-master threads: all threadprivate vars are undefined -- allocate fresh. */
        retained_supers = malloc(init_n * sizeof(char *));
        scores_retained_supers = malloc(init_n * sizeof(float));
        best_topology = malloc(init_n * sizeof(char *));
        best_topology_scores = malloc(init_n * sizeof(float));
        for(k=0; k<init_n; k++)
            {
            retained_supers[k] = malloc(TREE_LENGTH * sizeof(char));
            retained_supers[k][0] = '\0';
            scores_retained_supers[k] = -1;
            best_topology[k] = malloc(TREE_LENGTH * sizeof(char));
            best_topology[k][0] = '\0';
            best_topology_scores[k] = -1;
            }
        number_retained_supers = init_n;
        super_scores = malloc(number_of_taxa * sizeof(int *));
        for(k=0; k<number_of_taxa; k++)
            super_scores[k] = malloc(number_of_taxa * sizeof(int));
        sourcetree_scores = malloc(Total_fund_trees * sizeof(float));
        presenceof_SPRtaxa = malloc(number_of_taxa * sizeof(int));
        /* Point to master's read-only input data for hs parallel regions.
         * (For bootstrap parallel regions these pointers are set by
         *  boot_alloc_thread_state instead, overriding these assignments.) */
        fundamentals          = g_hs_fundamentals_snap;
        presence_of_taxa      = g_hs_presence_snap;
        fund_scores           = g_hs_fund_scores_snap;
        number_of_comparisons = g_hs_num_comp_snap;
        /* Share master's precomputed bipart sets (read-only during HS parallel reps). */
        fund_bipart_sets = g_hs_fund_bipart_snap;
        }
    else
        {
        /* Master thread: valid pointers from before the parallel region -- reset content. */
        for(k=0; k<number_retained_supers; k++)
            {
            retained_supers[k][0] = '\0';
            scores_retained_supers[k] = -1;
            if(best_topology && best_topology[k]) best_topology[k][0] = '\0';
            if(best_topology_scores) best_topology_scores[k] = -1;
            }
        /* Refresh super_scores (may contain stale path-metric data). */
        if(super_scores)
            {
            for(k=0; k<number_of_taxa; k++) free(super_scores[k]);
            free(super_scores);
            }
        super_scores = malloc(number_of_taxa * sizeof(int *));
        for(k=0; k<number_of_taxa; k++)
            super_scores[k] = malloc(number_of_taxa * sizeof(int));
        }

    /* Common to all threads: */
    thread_visited_acc = vs_create(1024);
    if(g_landscape_file[0])
        landscape_map = lm_create(8192);
    else
        landscape_map = NULL;
    for(k=0; k<Total_fund_trees; k++) sourcetree_scores[k] = -1;
    for(k=0; k<number_of_taxa; k++) presenceof_SPRtaxa[k] = -1;
    BESTSCORE = -1;
    NUMSWAPS = 0;
    sprscore = -1;
    tried_regrafts = 0;
    tree_top = NULL;
    temp_top = NULL;
    temp_top2 = NULL;
    branchpointer = NULL;
    thread_seed = (unsigned int)(seed + (unsigned int)omp_get_thread_num() * 1000003u);
    }


/*  hs_free_thread_state ---------------------------------------------------
 *  Free per-thread search globals.  Must be called inside a parallel region.
 *  The master thread does NOT free retained_supers / scores_retained_supers /
 *  best_topology / sourcetree_scores / presenceof_SPRtaxa -- these are freed
 *  by the caller after the parallel region ends.
 */
static void hs_free_thread_state(void)
    {
    int k, is_master;
    is_master = (omp_get_thread_num() == 0);

    /* Dismantle any tree left in working memory by the last do_search call. */
    if(tree_top != NULL) { dismantle_tree(tree_top); tree_top = NULL; }
    if(temp_top != NULL) { dismantle_tree(temp_top); temp_top = NULL; }

    /* Free super_scores (freshly allocated by hs_alloc_thread_state for all threads). */
    if(super_scores)
        {
        for(k=0; k<number_of_taxa; k++) free(super_scores[k]);
        free(super_scores); super_scores = NULL;
        }

    if(thread_visited_acc) { vs_free(thread_visited_acc); thread_visited_acc = NULL; }
    /* landscape_map should have been merged and freed in hs_merge_results; safety-free if not */
    if(landscape_map) { lm_free(landscape_map); landscape_map = NULL; }

    if(!is_master)
        {
        /* Non-master: free everything that was allocated per-thread. */
        /* fund_bipart_sets points to master's data (shared read-only) -- NULL the pointer, do not free. */
        fund_bipart_sets = NULL;
        if(sourcetree_scores) { free(sourcetree_scores); sourcetree_scores = NULL; }
        if(presenceof_SPRtaxa) { free(presenceof_SPRtaxa); presenceof_SPRtaxa = NULL; }
        for(k=0; k<number_retained_supers; k++)
            {
            if(retained_supers && retained_supers[k]) free(retained_supers[k]);
            if(best_topology && best_topology[k]) free(best_topology[k]);
            }
        if(retained_supers) { free(retained_supers); retained_supers = NULL; }
        if(scores_retained_supers) { free(scores_retained_supers); scores_retained_supers = NULL; }
        if(best_topology) { free(best_topology); best_topology = NULL; }
        if(best_topology_scores) { free(best_topology_scores); best_topology_scores = NULL; }
        }
    /* Master thread: retained_supers etc. are freed + replaced by caller after parallel ends. */
    }




/*  hs_same_topology -------------------------------------------------------
 *  Returns TRUE if the two named-taxon Newick strings t1 and t2 represent
 *  the same unrooted topology, using the path-metric comparison (identical
 *  logic to the inner loop of check_if_diff_tree).
 */
static int hs_same_topology(char *t1, char *t2)
    {
    int i=0, j=0, k=0, l=0, intname=-1, temp=0;
    int **scores1=NULL, **scores2=NULL;
    char *tree1 = malloc(TREE_LENGTH * sizeof(char));
    char *tree2 = malloc(TREE_LENGTH * sizeof(char));
    char name[1000];

    if(!tree1 || !tree2) { free(tree1); free(tree2); memory_error(62); }

    scores1 = malloc(number_of_taxa * sizeof(int *));
    if(!scores1) memory_error(62);
    for(i=0; i<number_of_taxa; i++)
        {
        scores1[i] = malloc(number_of_taxa * sizeof(int));
        if(!scores1[i]) memory_error(63);
        }

    scores2 = malloc(number_of_taxa * sizeof(int *));
    if(!scores2) memory_error(62);
    for(i=0; i<number_of_taxa; i++)
        {
        scores2[i] = malloc(number_of_taxa * sizeof(int));
        if(!scores2[i]) memory_error(63);
        }

    /* Convert t1 taxon names to integers -> tree1 */
    i=0; j=0;
    while(t1[i] != ';')
        {
        if(t1[i]=='(' || t1[i]==')' || t1[i]==',' || t1[i]==';')
            { tree1[j]=t1[i]; i++; j++; }
        else
            {
            k=0;
            while(t1[i]!='(' && t1[i]!=')' && t1[i]!=',' && t1[i]!=';')
                { name[k]=t1[i]; i++; k++; }
            name[k]=' ';
            intname=assign_taxa_name(name, FALSE);
            name[0]=' '; totext(intname, name);
            k=0;
            while(name[k]!=' ') { tree1[j]=name[k]; j++; k++; }
            }
        }
    tree1[j]=';'; tree1[j+1]=' ';

    /* Convert t2 taxon names to integers -> tree2 */
    i=0; j=0;
    while(t2[i] != ';')
        {
        if(t2[i]=='(' || t2[i]==')' || t2[i]==',' || t2[i]==';')
            { tree2[j]=t2[i]; i++; j++; }
        else
            {
            k=0;
            while(t2[i]!='(' && t2[i]!=')' && t2[i]!=',' && t2[i]!=';')
                { name[k]=t2[i]; i++; k++; }
            name[k]=' ';
            intname=assign_taxa_name(name, FALSE);
            name[0]=' '; totext(intname, name);
            k=0;
            while(name[k]!=' ') { tree2[j]=name[k]; j++; k++; }
            }
        }
    tree2[j]=';'; tree2[j+1]=' ';

    /* Initialise score matrices and compute path metrics */
    for(j=0; j<number_of_taxa; j++)
        for(k=j; k<number_of_taxa; k++)
            { scores1[j][k]=0; scores1[k][j]=0; scores2[j][k]=0; scores2[k][j]=0; }
    pathmetric(tree1, scores1);
    pathmetric(tree2, scores2);

    /* Sum absolute differences */
    temp=0;
    for(j=0; j<number_of_taxa; j++)
        for(k=j+1; k<number_of_taxa; k++)
            if(scores1[j][k] != scores2[j][k])
                temp += (scores1[j][k] > scores2[j][k])
                        ? scores1[j][k] - scores2[j][k]
                        : scores2[j][k] - scores1[j][k];

    for(i=0; i<number_of_taxa; i++) { free(scores1[i]); free(scores2[i]); }
    free(scores1); free(scores2);
    free(tree1); free(tree2);

    return(temp == 0);   /* TRUE if same topology */
    }


/*  hs_merge_results -------------------------------------------------------
 *  Merge the current thread's search results (threadprivate retained_supers /
 *  scores_retained_supers / NUMSWAPS) into the shared merge buffers.
 *  Must be called inside #pragma omp critical.
 */
static void hs_merge_results(char ***par_retained, float **par_scores, int *par_n,
                              float *par_best, int *par_NUMSWAPS, VisitedSet *par_visited)
    {
    int mi, mj, ml, is_dup, old_n;
    float msc;
    char *mtr;

    for(mi=0; mi<number_retained_supers; mi++)
        {
        if(scores_retained_supers[mi] < 0) break;
        msc = scores_retained_supers[mi];
        mtr = retained_supers[mi];

        if(*par_best < 0 || msc < *par_best)
            {
            /* New global best: clear merge buffers. */
            *par_best = msc;
            for(mj=0; mj<*par_n; mj++)
                {
                (*par_retained)[mj][0] = '\0';
                (*par_scores)[mj] = -1;
                }
            }

        /* Add tree if it matches the current best score (relative tolerance scales with score magnitude). */
        if(fabsf(msc - *par_best) <= 1e-5f * (fabsf(*par_best) + 1e-6f))
            {
            is_dup = FALSE;
            for(mj=0; mj<*par_n; mj++)
                {
                if((*par_scores)[mj] < 0) break;
                if(hs_same_topology((*par_retained)[mj], mtr)) { is_dup = TRUE; break; }
                }
            if(!is_dup)
                {
                /* Find first empty slot (score == -1). */
                for(mj=0; mj<*par_n && (*par_scores)[mj] >= 0; mj++);
                if(mj >= *par_n)
                    {
                    /* Grow merge buffers. */
                    old_n = *par_n;
                    *par_n += 10;
                    *par_retained = realloc(*par_retained, *par_n * sizeof(char *));
                    *par_scores   = realloc(*par_scores,   *par_n * sizeof(float));
                    for(ml=old_n; ml<*par_n; ml++)
                        {
                        (*par_retained)[ml] = malloc(TREE_LENGTH * sizeof(char));
                        (*par_retained)[ml][0] = '\0';
                        (*par_scores)[ml] = -1;
                        }
                    }
                (*par_retained)[mj] = realloc((*par_retained)[mj],
                                               (strlen(mtr)+10)*sizeof(char));
                strcpy((*par_retained)[mj], mtr);
                (*par_scores)[mj] = *par_best;
                }
            }
        }
    *par_NUMSWAPS += NUMSWAPS;
    if(par_visited && thread_visited_acc) vs_merge(par_visited, thread_visited_acc);
    /* Merge this thread's landscape map into the global accumulator */
    if(g_landscape_file[0] && landscape_map)
        {
        lm_merge(g_landscape_map, landscape_map);
        lm_free(landscape_map);
        landscape_map = NULL;
        }
    }
#endif /* _OPENMP */


/* =======================================================================
 * ensure_stored_source_trees: populate stored_* snapshot arrays from the
 * current active source-tree data if they have not already been filled by
 * a prior bootstrap_search() call.
 *
 * This snapshot is the reference set used by build_bootstrap_nj() to
 * resample starting trees for heuristic search reps.
 * ======================================================================= */
static void ensure_stored_source_trees(void)
	{
	int i, j, k;
	if(stored_funds != NULL) return;   /* already populated (e.g. by boot) */

	stored_num_comparisons = malloc(Total_fund_trees * sizeof(int));
	if(!stored_num_comparisons) memory_error(60);
	for(i = 0; i < Total_fund_trees; i++)
		stored_num_comparisons[i] = number_of_comparisons[i];

	stored_fund_scores = malloc(Total_fund_trees * sizeof(int **));
	if(!stored_fund_scores) memory_error(38);
	for(i = 0; i < Total_fund_trees; i++)
		{
		stored_fund_scores[i] = malloc(number_of_taxa * sizeof(int *));
		if(!stored_fund_scores[i]) memory_error(39);
		for(j = 0; j < number_of_taxa; j++)
			{
			stored_fund_scores[i][j] = malloc(number_of_taxa * sizeof(int));
			if(!stored_fund_scores[i][j]) memory_error(40);
			for(k = 0; k < number_of_taxa; k++)
				stored_fund_scores[i][j][k] = fund_scores[i][j][k];
			}
		}

	stored_presence_of_taxa = malloc(Total_fund_trees * sizeof(int *));
	if(!stored_presence_of_taxa) memory_error(41);
	for(i = 0; i < Total_fund_trees; i++)
		{
		stored_presence_of_taxa[i] = malloc(number_of_taxa * sizeof(int));
		if(!stored_presence_of_taxa[i]) memory_error(42);
		for(j = 0; j < number_of_taxa; j++)
			stored_presence_of_taxa[i][j] = presence_of_taxa[i][j];
		}

	stored_funds = malloc(Total_fund_trees * sizeof(char *));
	if(!stored_funds) memory_error(31);
	for(i = 0; i < Total_fund_trees; i++)
		{
		stored_funds[i] = malloc((strlen(fundamentals[i]) + 100) * sizeof(char));
		if(!stored_funds[i]) memory_error(32);
		strcpy(stored_funds[i], fundamentals[i]);
		}
	}


/* =======================================================================
 * build_bootstrap_nj: generate one bootstrap-resampled NJ starting tree.
 *
 * Resamples Total_fund_trees source trees with replacement from stored_*,
 * verifies every taxon is present, calls average_consensus + neighbor_joining
 * to produce a fresh NJ tree, then restores the active arrays to originals.
 *
 * Up to BOOT_NJ_MAX_ATTEMPTS resamples are tried.  Returns TRUE and writes
 * the Newick to out_tree on success; returns FALSE if every attempt produced
 * a sample missing at least one taxon (caller should fall back to SPR
 * perturbation of the original NJ tree).
 * ======================================================================= */
#define BOOT_NJ_MAX_ATTEMPTS 20
static int build_bootstrap_nj(char *out_tree, int missing_method)
	{
	int attempt, j, k, l, r, allpresent;

	for(attempt = 0; attempt < BOOT_NJ_MAX_ATTEMPTS; attempt++)
		{
		/* Resample Total_fund_trees source trees with replacement */
		for(j = 0; j < Total_fund_trees; j++)
			{
			r = (int)fmod(rand(), Total_fund_trees);
			fundamentals[j] = realloc(fundamentals[j],
			                          (strlen(stored_funds[r]) + 100) * sizeof(char));
			strcpy(fundamentals[j], stored_funds[r]);
			number_of_comparisons[j] = stored_num_comparisons[r];
			for(k = 0; k < number_of_taxa; k++)
				{
				presence_of_taxa[j][k] = stored_presence_of_taxa[r][k];
				for(l = 0; l < number_of_taxa; l++)
					fund_scores[j][k][l] = stored_fund_scores[r][k][l];
				}
			}

		/* Check every taxon appears in at least one resampled source tree */
		allpresent = TRUE;
		for(k = 0; k < number_of_taxa && allpresent; k++)
			{
			int found = FALSE;
			for(j = 0; j < Total_fund_trees && !found; j++)
				if(presence_of_taxa[j][k] > 0) found = TRUE;
			if(!found) allpresent = FALSE;
			}

		if(allpresent)
			{
			average_consensus(0, missing_method, NULL, NULL);
			neighbor_joining(FALSE, out_tree, FALSE);

			/* Restore active arrays to originals */
			for(j = 0; j < Total_fund_trees; j++)
				{
				fundamentals[j] = realloc(fundamentals[j],
				                          (strlen(stored_funds[j]) + 100) * sizeof(char));
				strcpy(fundamentals[j], stored_funds[j]);
				number_of_comparisons[j] = stored_num_comparisons[j];
				for(k = 0; k < number_of_taxa; k++)
					{
					presence_of_taxa[j][k] = stored_presence_of_taxa[j][k];
					for(l = 0; l < number_of_taxa; l++)
						fund_scores[j][k][l] = stored_fund_scores[j][k][l];
					}
				}
			return TRUE;
			}
		}

	/* All attempts failed — restore originals before returning */
	for(j = 0; j < Total_fund_trees; j++)
		{
		fundamentals[j] = realloc(fundamentals[j],
		                          (strlen(stored_funds[j]) + 100) * sizeof(char));
		strcpy(fundamentals[j], stored_funds[j]);
		number_of_comparisons[j] = stored_num_comparisons[j];
		for(k = 0; k < number_of_taxa; k++)
			{
			presence_of_taxa[j][k] = stored_presence_of_taxa[j][k];
			for(l = 0; l < number_of_taxa; l++)
				fund_scores[j][k][l] = stored_fund_scores[j][k][l];
			}
		}
	return FALSE;
	}


void heuristic_search(int user, int print, int sample, int nreps)
    {
    int i=0, j=0, k=0, l=0, swaps = 0, keep = 0, nbest = 0, start = 2, error = FALSE, numswaps = 1000000, different=TRUE, do_histogram = FALSE, here = FALSE, **taxa_comp = NULL, bins = 20, found = FALSE, missing_method = 1, random_num, numspectries = 2, numgenetries = 2, nthreads = 1;
    hs_maxskips_is_auto = 1;   /* reset: auto-scale N² unless user sets maxskips= this call */
    char *tree = NULL, c = '\0', *best_tree = NULL, *temptree = NULL, **starths = NULL, userfilename[10000], useroutfile[10000], histogramfile_name[10000];
    char *memory_start_tree = NULL;  /* saved copy of retained_supers[0] for start=memory (captured at parse time, before the retained_supers reset) */
    FILE *userfile = NULL, *outfile = NULL, *paupfile = NULL, *histogram_file = NULL;
    float distance=0, number=0, *startscores = NULL, used_weights = 0;
    int *saved_tags = NULL;  /* for single-copy auto-filter */
#ifdef _OPENMP
    /* Default to all available logical CPUs; user can override with nthreads= */
    nthreads = omp_get_num_procs();
    /* If called from inside a bootstrap parallel region, run single-threaded
     * to avoid nested parallelism (bootstrap threads already provide the
     * outer-level parallelism). */
    if(omp_in_parallel()) nthreads = 1;
#endif


	for(i=0; i<number_of_taxa; i++) presenceof_SPRtaxa[i] = -1;

    /* Reset landscape recording state for this hs call */
    g_landscape_file[0] = '\0';
    if(g_landscape_map) { lm_free(g_landscape_map); g_landscape_map = NULL; }
    g_cluster_enabled   = 0;
    strcpy(g_cluster_output, "treeclusters.tsv");
    g_cluster_threshold = 0.2f;
    g_cluster_orderby   = 0;

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
					printf2("Error: '%s' is an invalid value to be assigned to numspeciesrootings\n", parsed_command[i+1]);
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
					printf2("Error: '%s' is an invalid value to be assigned to numgenerootings\n", parsed_command[i+1]);
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
					printf2("Error: '%s' is an invalid method to be assigned to missing\n", parsed_command[i+1]);
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
					printf2("Error: '%s' is an invalid value for drawhistogram\n", parsed_command[i+1]);
					error = TRUE;
					}
				}
            }
		if(strcmp(parsed_command[i], "nbins") == 0)
			{
			bins = toint(parsed_command[i+1]);
			if(bins == 0)
				{
				printf2("Error: '%s' is an invalid integer number to be assigned to nbins\n", parsed_command[i+1]);
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
                printf2("Error: '%s' is an invalid integer number to be assigned to nsteps\n", parsed_command[i+1]);
                number_of_steps = 5;
                error = TRUE;
                }
            }
        if(strcmp(parsed_command[i], "nbest") == 0)
            {
            nbest = toint(parsed_command[i+1]);
            if(nbest == 0)
                {
                printf2("Error: '%s' is an invalid integer number to be assigned to nbest\n", parsed_command[i+1]);
                error = TRUE;
                }
            }
        if(strcmp(parsed_command[i], "maxswaps") == 0)
            {
            numswaps = toint(parsed_command[i+1]);
            if(numswaps == 0)
                {
                printf2("Error: '%s' is an invalid integer number to be assigned to maxswaps\n", parsed_command[i+1]);                
                error = TRUE;
                }
            }
        

        if((strcmp(parsed_command[i], "hsreps") == 0 )&& user == FALSE)  /* the user can assign the number of hs reps from the boot command */
            {
            nreps = toint(parsed_command[i+1]);
            if(nreps == 0)
                {
                printf2("Error: '%s' is an invalid integer number to be assigned to nreps\n", parsed_command[i+1]);
                error = TRUE;
                }
            }
        
        if((strcmp(parsed_command[i], "nreps") == 0) && user == TRUE)
            {
            nreps = toint(parsed_command[i+1]);
            
            if(nreps == 0)
                {
                printf2("Error: '%s' is an invalid integer number to be assigned to nreps\n", parsed_command[i+1]);
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
						printf2("Error: swap option '%s' not known\n", parsed_command[i+1]);
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
                    printf2("Error opening file named %s\n", parsed_command[i+1]);
                    error = TRUE;
                    strcpy(useroutfile, parsed_command[i+1]);
                    }
                else
                    {
					strcpy(useroutfile, parsed_command[i+1]);
                    printf2("opened output file %s\n", parsed_command[i+1]);
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
                        printf2("Error: weight option '%s' is unknown\n", parsed_command[i+1] );
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
                            printf2("Error: weight option '%s' is unknown\n", parsed_command[i+1] );
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
                        printf2("Error: weight option '%s' is unknown\n", parsed_command[i+1]);
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
						printf2("nj is not a valid option for starting trees in the recon criterion\n");
						error = TRUE;
						}
					}
				else if(strcmp(parsed_command[i+1], "memory") == 0)
					{
					if(trees_in_memory > 0 && retained_supers[0][0] != '\0')
						{
						start = 3;
						/* Capture NOW — retained_supers[0] is cleared before the search loop */
						memory_start_tree = malloc(TREE_LENGTH * sizeof(char));
						if(!memory_start_tree) memory_error(75);
						strcpy(memory_start_tree, retained_supers[0]);
						printf2("Starting tree: best tree currently in memory\n");
						}
					else
						{
						printf2("Warning: no tree currently in memory — falling back to NJ starting tree\n");
						start = 2;
						}
					}
				else
					{
					start = 1;
					if((userfile = fopen(parsed_command[i+1], "r")) == NULL)
						{

						printf2("Error opening file named %s\n", parsed_command[i+1]);
						error = TRUE;

						}
					else
						{
						printf2("opened user file %s\n", parsed_command[i+1]);
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
		if(strcmp(parsed_command[i], "nthreads") == 0)
			{
			nthreads = toint(parsed_command[i+1]);
			if(nthreads < 1)
				{
				printf2("Error: '%s' is an invalid value for nthreads (must be >= 1)\n", parsed_command[i+1]);
				nthreads = 1;
				error = TRUE;
				}
			}
		if(strcmp(parsed_command[i], "maxskips") == 0)
			{
			hs_maxskips = toint(parsed_command[i+1]);
			if(hs_maxskips < 0)
				{
				printf2("Error: maxskips must be >= 0\n");
				hs_maxskips = -1;
				error = TRUE;
				}
			hs_maxskips_is_auto = 0;   /* user has set an explicit value — don't auto-scale */
			}
		if(strcmp(parsed_command[i], "strategy") == 0)
			{
			if(strcmp(parsed_command[i+1], "first") == 0)
				hs_strategy = 0;
			else if(strcmp(parsed_command[i+1], "best") == 0)
				hs_strategy = 1;
			else
				{
				printf2("Error: strategy must be 'first' or 'best'\n");
				error = TRUE;
				}
			}
		if(strcmp(parsed_command[i], "progress") == 0)
			{
			hs_progress_interval = toint(parsed_command[i+1]);
			if(hs_progress_interval < 0)
				{
				printf2("Error: progress must be >= 0 (0=every improvement, N=every N seconds)\n");
				hs_progress_interval = 5;
				error = TRUE;
				}
			}
		if(strcmp(parsed_command[i], "droprep") == 0)
			{
			hs_droprep = (float)atof(parsed_command[i+1]);
			if(hs_droprep < 0.0f)
				{
				printf2("Error: droprep must be >= 0 (0=disabled, e.g. 0.1 = abandon if >10%% above global best)\n");
				hs_droprep = 0.0f;
				error = TRUE;
				}
			}
		if(strcmp(parsed_command[i], "threadreport") == 0)
			{
			hs_thread_report_interval = atoi(parsed_command[i+1]);
			if(hs_thread_report_interval < 0)
				{
				printf2("Error: threadreport must be >= 0 (0=disabled, e.g. 10 = each thread reports every 10s)\n");
				hs_thread_report_interval = 0;
				error = TRUE;
				}
			}
		if(strcmp(parsed_command[i], "mlbeta") == 0)
			{
			ml_beta = atof(parsed_command[i+1]);
			if(ml_beta <= 0)
				{
				printf2("Error: mlbeta must be > 0\n");
				ml_beta = 1.0f;
				error = TRUE;
				}
			else
				printf2("ML beta parameter set to %.4f\n", ml_beta);
			}
		if(strcmp(parsed_command[i], "mlscale") == 0)
			{
			if(strcmp(parsed_command[i+1], "lnl") == 0 || strcmp(parsed_command[i+1], "lnL") == 0)
				{
				ml_scale = 2;
				printf2("ML scoring: lnL display (reported as true Steel lnL = W*log(beta) - beta*WD)\n");
				}
			else if(strcmp(parsed_command[i+1], "lust") == 0 || strcmp(parsed_command[i+1], "LUSt") == 0)
				{
				ml_scale = 1;
				printf2("ML scoring: L.U.st-compatible scaling / Akanni et al. 2014 (beta * d * log10(e))\n");
				}
			else if(strcmp(parsed_command[i+1], "paper") == 0)
				{
				ml_scale = 0;
				printf2("ML scoring: Steel & Rodrigo (2008) formula (beta * d, minimised directly)\n");
				}
			else
				printf2("Error: mlscale must be 'paper', 'lust', or 'lnl'\n");
			}
		if(strcmp(parsed_command[i], "mleta") == 0)
			{
			double a = atof(parsed_command[i+1]);
			if(a < 0.0)
				printf2("Error: mleta must be >= 0 (0=Steel 2008, 1=normalised, >1=downweight large trees)\n");
			else
				{
				ml_eta = a;
				printf2("[Experimental] ml_eta set to %.4f (tree-size scaling exponent)\n", ml_eta);
				}
			}
		if(strcmp(parsed_command[i], "bsweight") == 0)
			{
			if(criterion == 2 || criterion == 3)
				{
				bsweight = atoi(parsed_command[i+1]);
				if(bsweight < 0 || bsweight > 1)
					{
					printf2("Error: bsweight must be 0 or 1\n");
					bsweight = 0; error = TRUE;
					}
				else
					printf2("Bootstrap split weights: %s\n", bsweight ? "on" : "off");
				}
			else
				printf2("Warning: bsweight only applies to sfit and qfit criteria\n");
			}
		if(strcmp(parsed_command[i], "visitedtrees") == 0)
			{
			strncpy(g_landscape_file, parsed_command[i+1], sizeof(g_landscape_file) - 1);
			g_landscape_file[sizeof(g_landscape_file) - 1] = '\0';
			}
		if(strcmp(parsed_command[i], "clusterlandscape") == 0)
			{
			if(strcmp(parsed_command[i+1], "yes") == 0)
				g_cluster_enabled = 1;
			else if(strcmp(parsed_command[i+1], "no") == 0)
				g_cluster_enabled = 0;
			else
				{ printf2("Error: clusterlandscape must be yes or no\n"); error = TRUE; }
			}
		if(strcmp(parsed_command[i], "clusteroutput") == 0)
			{
			strncpy(g_cluster_output, parsed_command[i+1], sizeof(g_cluster_output) - 1);
			g_cluster_output[sizeof(g_cluster_output) - 1] = '\0';
			}
		if(strcmp(parsed_command[i], "clusterthreshold") == 0)
			{
			g_cluster_threshold = tofloat(parsed_command[i+1]);
			if(g_cluster_threshold < 0.0f || g_cluster_threshold > 1.0f)
				{ printf2("Error: clusterthreshold must be between 0.0 and 1.0 (normalized RF distance)\n"); error = TRUE; }
			}
		if(strcmp(parsed_command[i], "clusterorderby") == 0)
			{
			if(strcmp(parsed_command[i+1], "score") == 0)
				g_cluster_orderby = 0;
			else if(strcmp(parsed_command[i+1], "visits") == 0)
				g_cluster_orderby = 1;
			else
				{ printf2("Error: clusterorderby must be score or visits\n"); error = TRUE; }
			}


        }
    for(i=0; i<num_commands; i++)   /* This has to be outside the normal for loop because if the assignment of hsreps (or nreps) is made after the assignment of smaple, then an error can be reported */
        {
 
        if(strcmp(parsed_command[i], "sample") == 0)  /* the user can assign the number of hs reps from the boot command */
            {
            sample = toint(parsed_command[i+1]);
            if(sample < nreps)
                {
                printf2("Error: '%s' is an invalid integer number to be assigned to sample\n It must be equal to or larger than", parsed_command[i+1]);
                if(user)printf2(" nreps (%d)\n", nreps);
                if(!user)printf2(" hsreps (%d)\n", nreps);
                error = TRUE;
                }
            }
        }

#ifdef _OPENMP
    /* Re-apply nested-parallelism guard: if called from inside a bootstrap
     * parallel region the user-specified nthreads= is ignored and we fall
     * back to single-threaded to prevent g_hs_*_snap races. */
    if(omp_in_parallel()) nthreads = 1;
#endif

    if(!error)
        {

		if(sample > sup) sample=sup;
		if(nreps > sup) nreps=sup;

        /* Auto-scale maxskips to 2*N² when no explicit value was given.
         * This makes the stopping criterion proportional to tree size:
         * N=10 → 200, N=20 → 800, N=30 → 1800, N=50 → 5000. */
        if(hs_maxskips_is_auto)
            hs_maxskips = 2 * number_of_taxa * number_of_taxa;

        /* Initialise global landscape map if visitedtrees= was specified */
        if(g_landscape_file[0])
            g_landscape_map = lm_create(8192);
		
        if(hsprint==TRUE)
            {
            /***** Print out summary of settings for this run */
            printf2("\n\nHeuristic Search settings:\n");
            printf2("\tCriterion = ");
            switch(criterion)
                {
                case 1:
                    printf2("Matrix Representation Using Parsimony (MRP)\n");
                    break;
                case 0:
                    printf2("Most Similar Supertree (dfit)\n");
                    break;
                case 2:
                    printf2("Maximum split fit (SFIT)\n");
                    break;
                case 4:
                    printf2("Average consensus (AVCON)\n");
                    break;
				case 5:
					printf2("Reconstruction of duplications and losses (RECON)\n");
					break;
                case 6:
                    printf2("Robinson-Foulds distance (RF)\n");
                    break;
                case 7:
                    printf2("Maximum likelihood supertree (beta=%.4f, eta=%.4f, scale=%s)\n",
                            ml_beta, ml_eta,
                            ml_scale == 1 ? "lust (Akanni et al. 2014)" : ml_scale == 2 ? "lnl" : "Steel & Rodrigo 2008");
                    break;
                default:
                    printf2("Maximum quartet fit (QFIT)\n");
                    break;
                    
                }
            if(criterion != 1 && criterion != 4)
                {
                printf2("\tHeuristic search algorithm = ");
                if(method == 1) printf2("Nearest Neighbor Interchange (NNI)\n");
                if(method == 2) printf2("Sub-tree Pruning and Regrafting (SPR)\n");
				if(method == 3) printf2("Tree Bisection and Reconnection (TBR)\n");
                printf2("\tMaximum Number of Steps (nsteps) = %d\n", number_of_steps);
                printf2("\tMaximum Number of Swaps (maxswaps) = %d\n", numswaps);
                printf2("\tNumber of repetitions of Heuristic search = %d\n", nreps);
#ifdef _OPENMP
                printf2("\tOpenMP threads = %d (of %d available)\n", nthreads, omp_get_num_procs());
#endif
                if(hs_maxskips > 0)
                    {
                    if(hs_maxskips_is_auto)
                        printf2("\tMax consecutive visited skips (maxskips) = %d (auto: 2×N²=2×%d²)\n", hs_maxskips, number_of_taxa);
                    else
                        printf2("\tMax consecutive visited skips (maxskips) = %d\n", hs_maxskips);
                    }
                else
                    printf2("\tMax consecutive visited skips (maxskips) = disabled\n");
#ifdef _OPENMP
                if(hs_progress_interval == 0)
                    printf2("\tParallel progress reporting (progress) = every improvement\n");
                else
                    printf2("\tParallel progress reporting (progress) = every %d seconds\n", hs_progress_interval);
                if(hs_droprep > 0.0f)
                    printf2("\tDrop lagging reps (droprep) = %.4f (abandon rep if >%.1f%% above global best)\n",
                            hs_droprep, hs_droprep * 100.0f);
                if(hs_thread_report_interval > 0)
                    printf2("\tPer-thread status reporting (threadreport) = every %d seconds\n", hs_thread_report_interval);
#endif
                if(criterion != 5)
					printf2("\tWeighting Scheme = ");
                if(criterion==0)
                    {
                    if(dweight == 0) printf2("equal\n");
                    if(dweight == 1) printf2("comparisons\n");
                    }
                if(criterion == 2)
                    {
                    if(splits_weight == 1) printf2("equal\n");
                    if(splits_weight == 2) printf2("splits\n");
                    }
                if(criterion == 3)
                    {
                    if(quartet_normalising == 1) printf2("equal\n");
                    if(quartet_normalising == 2) printf2("taxa\n");
                    if(quartet_normalising == 3) printf2("quartets\n");
                    }
                if(criterion == 2 || criterion == 3)
                    printf2("\tBS split weights (bsweight) = %s\n", bsweight ? "on" : "off");
                printf2("\tStarting trees = ");
                if(start == 0)
                    printf2("Top %d random trees chosen from %d random samples\n",nreps, sample);
                if(start == 1)
                    printf2("User defined from file named %s\n", userfilename);
				if(start == 2)
					{
					printf2("neighbor-joining tree from Average consensus distances\n\tMissing data estimated using ");
					if(missing_method == 1) printf2(" 4 point condition distances\n");
					if(missing_method == 0) printf2(" ultrametric distances\n");
					}
				if(start == 3)
					printf2("Best tree currently in memory (start=memory)\n");
				if(user) printf2("\tOutput file = %s\n", useroutfile);
                }
            }
		if(criterion == 0 && do_histogram == TRUE)
			{
			printf2("\tSource tree scores to be plotted against best supertree(s)\n");
			printf2("\tNumber of bins to summarise the source tree scores = %d\n", bins);
			printf2("\tHistogram to be written to file '%s'\n", histogramfile_name);
			}
		if(criterion==5 && hsprint == TRUE)
			{
			printf2("\n\tDuplication weight = %f\n\tLosses weight = %f\n", dup_weight, loss_weight);
			printf2("\tNumber of species tree rootings = ");
			if(numspectries == -1) printf2("all possible\n");
			else printf2("%d\n", numspectries);
			printf2("\tNumber of gene tree rootings = ");
			if(numgenetries == -1) printf2("all possible\n");
			else printf2("%d\n", numgenetries);
				
			}
		
#ifdef _OPENMP
        omp_set_num_threads(nthreads);
#endif
        /* criterion=recon scores all trees (single-copy and multicopy); skip the filter */
        if(criterion != 5)
            saved_tags = apply_singlecopy_filter();

        for(i=0; i<Total_fund_trees; i++) sourcetree_scores[i] = -1;
        /* Reset per-search state so each heuristic_search call starts fresh.
         * In sequential mode hs_alloc_thread_state is never called, so stale
         * values from a previous run would otherwise persist as the acceptance
         * threshold, distorting the second run's search and reported scores. */
        for(k=0; k<number_retained_supers; k++)
            {
            retained_supers[k][0]      = '\0';
            scores_retained_supers[k]  = -1;
            }
        BESTSCORE = -1;
        NUMSWAPS  = 0;
        if(criterion == 1 || criterion == 4)
            {
            if(criterion==1)  /* DO MRP */
                {
				BR_file = fopen("coding.nex", "w");
				error = coding(0, 1, 0);
				fclose(BR_file);
				if(error == FALSE)
					{
					if(print_log == FALSE) 
						strcpy(system_call, "paup coding.nex");
					else
						{
						strcpy(system_call, "paup coding.nex");
						strcat(system_call, " | tee -a ");
						strcat(system_call, logfile_name);
						}

					if(system(system_call) != 0) printf2("Error calling paup, please execute the file coding.nex in paup to perform the parsimony step\n");

					}
				if(error == FALSE)
					{
					remove("clanntmp.chr");
					remove("clanntree.chr");
					/*pars("coding.nex", "clanntmp.chr"); */
					printf2("\n");
					}

                }
            if(criterion == 4)   /* DO AVERAGE CONSENSUS */
                {
					
				if(!got_weights)
					{
					printf2("Warning: There were no weights included in the input trees.\nThe analysis will assume that all branch lengths are set to unity (1)\n");
					}
				paupfile = fopen("average_consensus.nex", "w");
				error = average_consensus(0, missing_method,  useroutfile, paupfile);
				fclose(paupfile);

				if(print_log == FALSE) 
					strcpy(system_call, "paup average_consensus.nex");
				else
					{
					strcpy(system_call, "paup average_consensus.nex");
					strcat(system_call, " | tee -a ");
					strcat(system_call, logfile_name);
					}

				if(system(system_call) != 0) printf2("Error calling PAUP*\n\tPlease execute the file average_consensus.nex in PAUP to complete the analysis\n");
				           
                }
            }  /***** FINISH AVERAGE CONSENSUS */
        else
            {
            
            psfile = fopen("supertree.ps", "w");
            
            tree = malloc(TREE_LENGTH*sizeof(char));
            if(!tree) memory_error(53);
            tree[0] = '\0';

            
            if(criterion == 2 || criterion == 3 || criterion == 6 || criterion == 7)
                rf_precompute_fund_biparts();
            
            if(outfile == NULL && print)
                {
                if(outfile == NULL)
                    {
                    if((outfile = fopen("Heuristic_result.txt", "w")) == NULL)
                        {
                        printf2("Error opening file named Heuristic_result.txt\n");
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
			
			
        

        /* Initialise thread_seed for sequential mode (overridden per-thread inside parallel regions). */
        thread_seed = (unsigned int)seed;

            if(start == 0)  /* if we are to use a random starting tree */
                {
                swaps = 0;
                i=0;
                tried_regrafts = 0;
				NUMSWAPS=0;
                interval1 = time(NULL);
                if(print)
                    {
                    printf2("\nRandom sampling progress indicator:");
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
						printf2("An error occurred while setting a signal handler\n");
						}
                for(GC=0; GC<sample; GC++)
                    {
					i++;
                    interval2 = time(NULL);
                    if(difftime(interval2, interval1) > 5) /* every 5 seconds print a dot to the screen */
                        {
                        printf2("*");
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
                    { int _to = 0; tree_build(1, tree, tree_top, TRUE, -1, &_to); }
                    tree_top = temp_top;
                    temp_top = NULL;
            
                    strcpy(best_tree, "");
                    print_named_tree(tree_top, best_tree);
                    while(unroottree(best_tree));
                    strcat(best_tree, ";");
                    /* evaluate the random tree */
                    if(criterion == 0) distance = compare_trees(FALSE);  /* calculate the distance from the super tree to all the fundamental trees */
                    if(criterion == 2) distance = compare_trees_sfit(FALSE);
                    if(criterion == 3) distance = compare_trees_qfit(FALSE);
					if(criterion==5)
						{
						strcpy(temptree, "");
                        print_tree(tree_top, temptree);
                        strcat(temptree, ";");
                        unroottree(temptree);
						distance = get_recon_score(temptree, numspectries, numgenetries);
						}
				if(criterion == 6)
					{
					if(fund_bipart_sets == NULL) rf_precompute_fund_biparts();
					distance = compare_trees_rf(FALSE);
					}
				if(criterion == 7)
					{
					if(fund_bipart_sets == NULL) rf_precompute_fund_biparts();
					distance = compare_trees_ml(FALSE);
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
				
                if(signal(SIGINT, controlc2) == SIG_ERR)
					{
					printf2("An error occurred while setting a signal handler\n");
					}
				for(i=0; i<Total_fund_trees; i++)sourcetree_scores[i] = -1;
				i=0;
				 if(print) printf2("\nHeuristic search progress indicator:\n");
                    fflush(stdout);
#ifdef _OPENMP
                if(nthreads > 1)
                    {
                    /* ---- Parallel start==0 nreps loop ---- */
                    int    prl_i;
                    char **par_retained = NULL;
                    float *par_scores   = NULL;
                    int    par_n        = 10;
                    float  par_best     = -1;
                    int    par_NUMSWAPS = 0;
                    VisitedSet *par_visited_set = vs_create(1024);

                    par_retained = malloc(par_n * sizeof(char *));
                    par_scores   = malloc(par_n * sizeof(float));
                    for(k=0; k<par_n; k++)
                        {
                        par_retained[k] = malloc(TREE_LENGTH * sizeof(char));
                        par_retained[k][0] = '\0';
                        par_scores[k] = -1;
                        }

                    if(criterion == 2 || criterion == 3 || criterion == 6 || criterion == 7) rf_precompute_fund_biparts();
                    par_progress_best    = -1.0f;
                    par_last_print_score = -1.0f;
                    par_last_progress_time = 0;
                    par_search_start     = time(NULL);
                    /* Snapshot master's input-data pointers for worker threads */
                    g_hs_fundamentals_snap = fundamentals;
                    g_hs_presence_snap     = presence_of_taxa;
                    g_hs_fund_scores_snap  = fund_scores;
                    g_hs_num_comp_snap     = number_of_comparisons;
                    g_hs_fund_bipart_snap  = fund_bipart_sets;
                    #pragma omp parallel num_threads(nthreads) default(shared) \
                            private(prl_i)
                        {
                        int thr_swaps = 0;

                        hs_alloc_thread_state();

                        #pragma omp for schedule(dynamic, 1)
                        for(prl_i=0; prl_i<nreps; prl_i++)
                            {
                            if(user_break) continue;
                            hs_par_rep = prl_i + 1;
                            thr_swaps += do_search(starths[prl_i], TRUE, FALSE, numswaps, outfile, numspectries, numgenetries);

                            /* --- per-rep completion report (all threads) --- */
                            #pragma omp critical (par_progress_print)
                                {
                                int new_best = 0;
                                if(sprscore >= 0.0f && (par_progress_best < 0.0f || sprscore < par_progress_best))
                                    {
                                    par_progress_best = sprscore;
                                    new_best = 1;
                                    }
                                double hs_el = difftime(time(NULL), par_search_start);
                                int    hs_em = (int)(hs_el / 60), hs_es = (int)hs_el % 60;
                                const char *stop_reason;
                                if(user_break)                                           stop_reason = "interrupted";
                                else if(rep_abandon)                                     stop_reason = "suboptimal";
                                else if(hs_maxskips > 0 && skip_streak >= hs_maxskips)  stop_reason = "maxskips";
                                else if(tried_regrafts >= numswaps)                      stop_reason = "maxswaps";
                                else                                                     stop_reason = "converged";
                                if(hsprint)
                                    {
                                    if(sprscore >= 0.0f)
                                        printf2("  [%2d:%02d]  rep %d/%d  %s= %s%.6f  (%d spr/tbr + %d nni swaps, %s)\n",
                                                hs_em, hs_es, prl_i+1, nreps,
                                                ml_score_label(), new_best ? "* " : "  ",
                                                ml_display_score(sprscore),
                                                tried_regrafts, nni_swaps, stop_reason);
                                    else
                                        printf2("  [%2d:%02d]  rep %d/%d  (no new topologies)  (%d spr/tbr + %d nni swaps, %s)\n",
                                                hs_em, hs_es, prl_i+1, nreps,
                                                tried_regrafts, nni_swaps, stop_reason);
                                    fflush(stdout);
                                    }
                                }
                            }

                        #pragma omp critical (hs_merge)
                            {
                            hs_merge_results(&par_retained, &par_scores, &par_n, &par_best, &par_NUMSWAPS, par_visited_set);
                            swaps += thr_swaps;
                            }

                        hs_free_thread_state();
                        }
                    /* ---- end parallel region ---- */

                    /* Install merged results into master globals. */
                    for(k=0; k<number_retained_supers; k++)
                        { if(retained_supers[k]) free(retained_supers[k]); }
                    free(retained_supers);
                    free(scores_retained_supers);
                    retained_supers        = par_retained;
                    scores_retained_supers = par_scores;
                    number_retained_supers = par_n;
                    BESTSCORE              = par_best;
                    NUMSWAPS               = par_NUMSWAPS;
                    swaps = (int)par_visited_set->count + (par_visited_set->zero_present ? 1 : 0);
                    vs_free(par_visited_set);

                    if(user_break && print)
                        printf2("All threads stopped — writing best tree to output file.\n");

                    /* Re-allocate super_scores for post-loop use (freed by hs_free_thread_state). */
                    super_scores = malloc(number_of_taxa * sizeof(int *));
                    for(k=0; k<number_of_taxa; k++)
                        super_scores[k] = malloc(number_of_taxa * sizeof(int));
                    }
                else
#endif
                    {
                    /* ---- Sequential while loop ---- */
                    if(criterion == 2 || criterion == 3 || criterion == 6 || criterion == 7) rf_precompute_fund_biparts();
                    par_progress_best = -1.0f;
                    par_last_progress_time = 0;
                    par_search_start  = time(NULL);
                    while(i != nreps && !user_break) /* Start searching tree space */
                        {
                        swaps += do_search(starths[i], TRUE, print, numswaps, outfile, numspectries, numgenetries);
                        /* --- per-rep completion report --- */
                        {
                        int new_best = 0;
                        if(sprscore >= 0.0f && (par_progress_best < 0.0f || sprscore < par_progress_best))
                            {
                            par_progress_best = sprscore;
                            new_best = 1;
                            }
                        double hs_el = difftime(time(NULL), par_search_start);
                        int    hs_em = (int)(hs_el / 60), hs_es = (int)hs_el % 60;
                        const char *stop_reason;
                        if(user_break)                                           stop_reason = "interrupted";
                        else if(rep_abandon)                                     stop_reason = "suboptimal";
                        else if(hs_maxskips > 0 && skip_streak >= hs_maxskips)  stop_reason = "maxskips";
                        else if(tried_regrafts >= numswaps)                      stop_reason = "maxswaps";
                        else                                                     stop_reason = "converged";
                        if(hsprint)
                            {
                            if(sprscore >= 0.0f)
                                printf2("  [%2d:%02d]  rep %d/%d  %s= %s%.6f  (%d spr/tbr + %d nni swaps, %s)\n",
                                        hs_em, hs_es, i+1, nreps,
                                        ml_score_label(), new_best ? "* " : "  ",
                                        ml_display_score(sprscore),
                                        tried_regrafts, nni_swaps, stop_reason);
                            else
                                printf2("  [%2d:%02d]  rep %d/%d  (no new topologies)  (%d spr/tbr + %d nni swaps, %s)\n",
                                        hs_em, hs_es, i+1, nreps,
                                        tried_regrafts, nni_swaps, stop_reason);
                            fflush(stdout);
                            }
                        }
                        i++;
                        }
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
						printf2("An error occurred while setting a signal handler\n");
						}
					for(i=0; i<Total_fund_trees; i++)sourcetree_scores[i] = -1;
					 if(print) printf2("\nHeuristic search progress indicator:\n");
#ifdef _OPENMP
					if(nthreads > 1)
						{
						/* ---- Parallel nreps loop (start==2, NJ+perturbation mode) ---- */
						int    prl_i;
						char **par_retained = NULL;
						float *par_scores   = NULL;
						int    par_n        = 10;
						float  par_best     = -1;
						int    par_NUMSWAPS = 0;
						VisitedSet *par_visited_set = vs_create(1024);

						/* Pre-compute nreps starting trees.
						   Rep 0    = NJ tree from all source trees (best single starting point).
						   Reps 1+  = NJ tree from a bootstrap resample of source trees — each
						              resample produces a different distance matrix and therefore
						              a genuinely different NJ topology, giving independent basin
						              sampling without relying on SPR perturbation.
						   Fallback = if a resample cannot cover all taxa after BOOT_NJ_MAX_ATTEMPTS
						              tries, fall back to NJ + random SPR moves for that rep. */
						ensure_stored_source_trees();
						char **start_trees = malloc(nreps * sizeof(char *));
						char *nj_tree_save = malloc(TREE_LENGTH * sizeof(char));
						strcpy(nj_tree_save, temptree);   /* save pristine NJ tree */
						for(k=0; k<nreps; k++)
							{
							start_trees[k] = malloc(TREE_LENGTH * sizeof(char));
							if(k == 0)
								{
								strcpy(start_trees[k], nj_tree_save);
								}
							else
								{
								if(!build_bootstrap_nj(start_trees[k], missing_method))
									{
									/* Fallback: SPR perturbation of original NJ */
									strcpy(temptree, nj_tree_save);
									random_num = 10 + (int)fmod(rand(), 21);   /* 10-30 moves */
									for(j=0; j<random_num; j++)
										string_SPR(temptree);
									strcpy(start_trees[k], temptree);
									}
								}
							}
						free(nj_tree_save);

						/* Initialise shared merge buffers. */
						par_retained = malloc(par_n * sizeof(char *));
						par_scores   = malloc(par_n * sizeof(float));
						for(k=0; k<par_n; k++)
							{
							par_retained[k] = malloc(TREE_LENGTH * sizeof(char));
							par_retained[k][0] = '\0';
							par_scores[k] = -1;
							}

						if(criterion == 2 || criterion == 3 || criterion == 6 || criterion == 7) rf_precompute_fund_biparts();
						par_progress_best    = -1.0f;
						par_last_print_score = -1.0f;
						par_last_progress_time = 0;
						par_search_start     = time(NULL);
                    /* Snapshot master's input-data pointers for worker threads */
                    g_hs_fundamentals_snap = fundamentals;
                    g_hs_presence_snap     = presence_of_taxa;
                    g_hs_fund_scores_snap  = fund_scores;
                    g_hs_num_comp_snap     = number_of_comparisons;
                    g_hs_fund_bipart_snap  = fund_bipart_sets;
						#pragma omp parallel num_threads(nthreads) default(shared) \
						        private(prl_i)
							{
							char *thr_tree = malloc(TREE_LENGTH * sizeof(char));
							int  thr_swaps = 0;

							hs_alloc_thread_state();

							#pragma omp for schedule(dynamic, 1)
							for(prl_i=0; prl_i<nreps; prl_i++)
								{
								if(user_break) continue;
								hs_par_rep = prl_i + 1;
								strcpy(thr_tree, start_trees[prl_i]);
								returntree(thr_tree);
								thr_swaps += do_search(thr_tree, TRUE, FALSE, numswaps, outfile, numspectries, numgenetries);

								/* --- per-rep completion report (all threads) --- */
								#pragma omp critical (par_progress_print)
									{
									int new_best = 0;
									if(sprscore >= 0.0f && (par_progress_best < 0.0f || sprscore < par_progress_best))
										{
										par_progress_best = sprscore;
										new_best = 1;
										}
									double hs_el = difftime(time(NULL), par_search_start);
									int    hs_em = (int)(hs_el / 60), hs_es = (int)hs_el % 60;
									const char *stop_reason;
									if(user_break)                                           stop_reason = "interrupted";
									else if(rep_abandon)                                     stop_reason = "suboptimal";
									else if(hs_maxskips > 0 && skip_streak >= hs_maxskips)  stop_reason = "maxskips";
									else if(tried_regrafts >= numswaps)                      stop_reason = "maxswaps";
									else                                                     stop_reason = "converged";
									if(hsprint)
										{
										if(sprscore >= 0.0f)
											printf2("  [%2d:%02d]  rep %d/%d  %s= %s%.6f  (%d spr/tbr + %d nni swaps, %s)\n",
												hs_em, hs_es, prl_i+1, nreps,
												ml_score_label(), new_best ? "* " : "  ",
												ml_display_score(sprscore),
												tried_regrafts, nni_swaps, stop_reason);
										else
											printf2("  [%2d:%02d]  rep %d/%d  (no new topologies)  (%d spr/tbr + %d nni swaps, %s)\n",
												hs_em, hs_es, prl_i+1, nreps,
												tried_regrafts, nni_swaps, stop_reason);
										fflush(stdout);
										}
									}
								}

							#pragma omp critical (hs_merge)
								{
								hs_merge_results(&par_retained, &par_scores, &par_n, &par_best, &par_NUMSWAPS, par_visited_set);
								swaps += thr_swaps;
								}

							hs_free_thread_state();
							free(thr_tree);
							}
						/* ---- end parallel region ---- */

						/* Install merged results into master globals. */
						for(k=0; k<number_retained_supers; k++)
							{ if(retained_supers[k]) free(retained_supers[k]); }
						free(retained_supers);
						free(scores_retained_supers);
						retained_supers        = par_retained;
						scores_retained_supers = par_scores;
						number_retained_supers = par_n;
						BESTSCORE              = par_best;
						NUMSWAPS               = par_NUMSWAPS;
						swaps = (int)par_visited_set->count + (par_visited_set->zero_present ? 1 : 0);
						vs_free(par_visited_set);

						if(user_break && print)
							printf2("All threads stopped — writing best tree to output file.\n");

						/* Re-allocate super_scores for post-loop use (freed by hs_free_thread_state). */
						super_scores = malloc(number_of_taxa * sizeof(int *));
						for(k=0; k<number_of_taxa; k++)
							super_scores[k] = malloc(number_of_taxa * sizeof(int));

						for(k=0; k<nreps; k++) free(start_trees[k]);
						free(start_trees);
						}
					else
#endif
						{
						/* ---- Sequential nreps loop ---- */
						/* Rep 0 = NJ from all source trees; reps 1+ = bootstrap-resampled NJ.
						   Fallback to SPR perturbation if resample cannot cover all taxa. */
						{
						ensure_stored_source_trees();
						char *nj_tree_save = malloc(TREE_LENGTH * sizeof(char));
						strcpy(nj_tree_save, temptree);
						if(criterion == 2 || criterion == 3 || criterion == 6 || criterion == 7) rf_precompute_fund_biparts();
						par_progress_best = -1.0f;
						par_last_progress_time = 0;
						par_search_start  = time(NULL);
						for(i=0; i<nreps; i++)
							{
							if(i != 0)
								{
								if(!build_bootstrap_nj(temptree, missing_method))
									{
									/* Fallback: SPR perturbation of original NJ */
									strcpy(temptree, nj_tree_save);
									random_num = 10 + (int)fmod(rand(), 21);   /* 10-30 moves */
									for(j=0; j<random_num; j++)
										string_SPR(temptree);
									}
								}
							strcpy(tree, temptree);
							returntree(tree);
							swaps += do_search(tree, TRUE, print, numswaps, outfile, numspectries,numgenetries);
							/* --- per-rep completion report --- */
							{
							int new_best = 0;
							if(sprscore >= 0.0f && (par_progress_best < 0.0f || sprscore < par_progress_best))
								{
								par_progress_best = sprscore;
								new_best = 1;
								}
							double hs_el = difftime(time(NULL), par_search_start);
							int    hs_em = (int)(hs_el / 60), hs_es = (int)hs_el % 60;
							const char *stop_reason;
							if(user_break)                                           stop_reason = "interrupted";
							else if(rep_abandon)                                     stop_reason = "suboptimal";
							else if(hs_maxskips > 0 && skip_streak >= hs_maxskips)  stop_reason = "maxskips";
							else if(tried_regrafts >= numswaps)                      stop_reason = "maxswaps";
							else                                                     stop_reason = "converged";
							if(hsprint)
								{
								if(sprscore >= 0.0f)
									printf2("  [%2d:%02d]  rep %d/%d  %s= %s%.6f  (%d spr/tbr + %d nni swaps, %s)\n",
										hs_em, hs_es, i+1, nreps,
										ml_score_label(), new_best ? "* " : "  ",
										ml_display_score(sprscore),
										tried_regrafts, nni_swaps, stop_reason);
								else
									printf2("  [%2d:%02d]  rep %d/%d  (no new topologies)  (%d spr/tbr + %d nni swaps, %s)\n",
										hs_em, hs_es, i+1, nreps,
										tried_regrafts, nni_swaps, stop_reason);
								fflush(stdout);
								}
							}
							}
						free(nj_tree_save);
						}
						}
					}
				else if(start == 3)  /* use best tree currently in memory as starting point */
					{
					if(criterion == 2 || criterion == 3 || criterion == 6 || criterion == 7) rf_precompute_fund_biparts();
					if(signal(SIGINT, controlc2) == SIG_ERR)
						printf2("An error occurred while setting a signal handler\n");
					for(i = 0; i < Total_fund_trees; i++) sourcetree_scores[i] = -1;
					if(print) printf2("\nHeuristic search progress indicator:\n");
					/* Pre-compute nreps starting trees in numbered format (mirrors NJ start path).
					 * Rep 0    = memory tree converted to numbered format (exact in-memory topology).
					 * Reps 1+  = same base tree with 10-30 random string_SPR moves applied, giving
					 *             independent starting points around the known-best neighbourhood.
					 * returntree() is called inside the parallel/sequential loops (as for start=nj)
					 * to convert numbered trees to named format immediately before do_search(). */
					{
					char **start_trees   = malloc(nreps * sizeof(char *));
					char *mem_numbered   = malloc(TREE_LENGTH * sizeof(char));

					/* Convert named memory tree → numbered Newick once.
					 * Use memory_start_tree (saved at parse time) — retained_supers[0]
					 * has been cleared by the reset loop before we reach this point. */
					if(tree_top != NULL) { dismantle_tree(tree_top); tree_top = NULL; }
					temp_top = NULL;
					{ int _to = 0; tree_build(1, memory_start_tree, tree_top, TRUE, -1, &_to); }
					tree_top = temp_top;
					temp_top = NULL;
					number_tree(tree_top, 0);
					strcpy(mem_numbered, "");
					print_tree(tree_top, mem_numbered);
					strcat(mem_numbered, ";");

					for(k = 0; k < nreps; k++)
						{
						start_trees[k] = malloc(TREE_LENGTH * sizeof(char));
						strcpy(start_trees[k], mem_numbered);
						if(k > 0)
							{
							/* Perturb with 10-30 random SPR moves for diversity */
							random_num = 10 + (int)fmod(rand(), 21);
							for(j = 0; j < random_num; j++)
								string_SPR(start_trees[k]);
							}
						}
					free(mem_numbered);

#ifdef _OPENMP
					if(nthreads > 1)
						{
						/* ---- Parallel start==3 nreps loop (mirrors start==2 parallel path) ---- */
						int    prl_i;
						char **par_retained = NULL;
						float *par_scores   = NULL;
						int    par_n        = 10;
						float  par_best     = -1;
						int    par_NUMSWAPS = 0;
						VisitedSet *par_visited_set = vs_create(1024);

						par_retained = malloc(par_n * sizeof(char *));
						par_scores   = malloc(par_n * sizeof(float));
						for(k = 0; k < par_n; k++)
							{
							par_retained[k] = malloc(TREE_LENGTH * sizeof(char));
							par_retained[k][0] = '\0';
							par_scores[k] = -1;
							}

						par_progress_best    = -1.0f;
						par_last_print_score = -1.0f;
						par_last_progress_time = 0;
						par_search_start     = time(NULL);
						g_hs_fundamentals_snap = fundamentals;
						g_hs_presence_snap     = presence_of_taxa;
						g_hs_fund_scores_snap  = fund_scores;
						g_hs_num_comp_snap     = number_of_comparisons;
						g_hs_fund_bipart_snap  = fund_bipart_sets;
						#pragma omp parallel num_threads(nthreads) default(shared) \
						        private(prl_i)
							{
							char *thr_tree = malloc(TREE_LENGTH * sizeof(char));
							int  thr_swaps = 0;

							hs_alloc_thread_state();

							#pragma omp for schedule(dynamic, 1)
							for(prl_i = 0; prl_i < nreps; prl_i++)
								{
								if(user_break) continue;
								hs_par_rep = prl_i + 1;
								strcpy(thr_tree, start_trees[prl_i]);
								returntree(thr_tree);
								thr_swaps += do_search(thr_tree, TRUE, FALSE, numswaps, outfile, numspectries, numgenetries);

								#pragma omp critical (par_progress_print)
									{
									int new_best = 0;
									if(sprscore >= 0.0f && (par_progress_best < 0.0f || sprscore < par_progress_best))
										{ par_progress_best = sprscore; new_best = 1; }
									double hs_el = difftime(time(NULL), par_search_start);
									int    hs_em = (int)(hs_el / 60), hs_es = (int)hs_el % 60;
									const char *stop_reason;
									if(user_break)                                           stop_reason = "interrupted";
									else if(rep_abandon)                                     stop_reason = "suboptimal";
									else if(hs_maxskips > 0 && skip_streak >= hs_maxskips)  stop_reason = "maxskips";
									else if(tried_regrafts >= numswaps)                      stop_reason = "maxswaps";
									else                                                     stop_reason = "converged";
									if(hsprint)
										{
										if(sprscore >= 0.0f)
											printf2("  [%2d:%02d]  rep %d/%d  %s= %s%.6f  (%d spr/tbr + %d nni swaps, %s)\n",
												hs_em, hs_es, prl_i+1, nreps,
												ml_score_label(), new_best ? "* " : "  ",
												ml_display_score(sprscore),
												tried_regrafts, nni_swaps, stop_reason);
										else
											printf2("  [%2d:%02d]  rep %d/%d  (no new topologies)  (%d spr/tbr + %d nni swaps, %s)\n",
												hs_em, hs_es, prl_i+1, nreps,
												tried_regrafts, nni_swaps, stop_reason);
										fflush(stdout);
										}
									}
								}

							#pragma omp critical (hs_merge)
								{
								hs_merge_results(&par_retained, &par_scores, &par_n, &par_best, &par_NUMSWAPS, par_visited_set);
								swaps += thr_swaps;
								}

							hs_free_thread_state();
							free(thr_tree);
							}
						/* ---- end parallel region ---- */

						for(k = 0; k < number_retained_supers; k++)
							{ if(retained_supers[k]) free(retained_supers[k]); }
						free(retained_supers);
						free(scores_retained_supers);
						retained_supers        = par_retained;
						scores_retained_supers = par_scores;
						number_retained_supers = par_n;
						BESTSCORE              = par_best;
						NUMSWAPS               = par_NUMSWAPS;
						swaps = (int)par_visited_set->count + (par_visited_set->zero_present ? 1 : 0);
						vs_free(par_visited_set);

						if(user_break && print)
							printf2("All threads stopped — writing best tree to output file.\n");

						/* Re-allocate super_scores for post-loop use (freed by hs_free_thread_state). */
						super_scores = malloc(number_of_taxa * sizeof(int *));
						for(k = 0; k < number_of_taxa; k++)
							super_scores[k] = malloc(number_of_taxa * sizeof(int));
						}
					else
#endif
						{
						/* ---- Sequential nreps loop ---- */
						par_progress_best = -1.0f;
						par_last_progress_time = 0;
						par_search_start  = time(NULL);
						for(i = 0; i < nreps && !user_break; i++)
							{
							strcpy(tree, start_trees[i]);
							returntree(tree);
							swaps += do_search(tree, TRUE, print, numswaps, outfile, numspectries, numgenetries);
							{
							int new_best = 0;
							if(sprscore >= 0.0f && (par_progress_best < 0.0f || sprscore < par_progress_best))
								{ par_progress_best = sprscore; new_best = 1; }
							double hs_el = difftime(time(NULL), par_search_start);
							int    hs_em = (int)(hs_el / 60), hs_es = (int)hs_el % 60;
							const char *stop_reason;
							if(user_break)                                           stop_reason = "interrupted";
							else if(rep_abandon)                                     stop_reason = "suboptimal";
							else if(hs_maxskips > 0 && skip_streak >= hs_maxskips)  stop_reason = "maxskips";
							else if(tried_regrafts >= numswaps)                      stop_reason = "maxswaps";
							else                                                     stop_reason = "converged";
							if(hsprint)
								{
								if(sprscore >= 0.0f)
									printf2("  [%2d:%02d]  rep %d/%d  %s= %s%.6f  (%d spr/tbr + %d nni swaps, %s)\n",
										hs_em, hs_es, i+1, nreps,
										ml_score_label(), new_best ? "* " : "  ",
										ml_display_score(sprscore),
										tried_regrafts, nni_swaps, stop_reason);
								else
									printf2("  [%2d:%02d]  rep %d/%d  (no new topologies)  (%d spr/tbr + %d nni swaps, %s)\n",
										hs_em, hs_es, i+1, nreps,
										tried_regrafts, nni_swaps, stop_reason);
								fflush(stdout);
								}
							}
							}
						}

					for(k = 0; k < nreps; k++) free(start_trees[k]);
					free(start_trees);
					free(memory_start_tree);
					memory_start_tree = NULL;
					}
					}
				else/* if we are to use a tree (or trees) from a file as the starting tree */
					{
					if(criterion == 2 || criterion == 3 || criterion == 6 || criterion == 7) rf_precompute_fund_biparts();
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
						if(print)printf2("\t\nusertree %d progress indicator: =", k+1);
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
						{ int _to = 0; tree_build(1, tree, tree_top, user, -1, &_to); }
						tree_top = temp_top;
						temp_top = NULL;
						
									
						strcpy(best_tree, "");
						print_named_tree(tree_top, best_tree);
						strcat(best_tree, ";");
						while(unroottree(best_tree)==TRUE);
						/*  printf2("!:assigned tree: %s\n", best_tree); */
								
								
						if(criterion == 0) distance = compare_trees(FALSE);  /* calculate the distance from the super tree to all the fundamental trees */
						if(criterion == 2) distance = compare_trees_sfit(FALSE);
						if(criterion == 3) distance = compare_trees_qfit(FALSE);
						if(criterion==5)
							{
							strcpy(temptree, "");
							print_tree(tree_top, temptree);
							strcat(temptree, ";");

							while(unroottree(temptree)==TRUE);
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
								while(unroottree(best_tree)== TRUE);
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
									while(unroottree(best_tree)== TRUE);
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
                if(print)printf2("\n");
                if(print)printf2("Number of unique topologies scored: %d\n",swaps);

                /* In sequential mode landscape_map is not NULL and holds all accumulated
                 * topology records.  Merge into g_landscape_map now (parallel mode already
                 * merged per-thread into g_landscape_map inside hs_merge_results). */
                if(g_landscape_file[0] && landscape_map)
                    {
                    lm_merge(g_landscape_map, landscape_map);
                    lm_free(landscape_map);
                    landscape_map = NULL;
                    }


                /**** Rescore retained supertrees to fix stale-cache scores from SPR search ****/
                /* presenceof_SPRtaxa[] is never set to TRUE/FALSE, so compare_trees(TRUE)
                 * always uses cached scores; scores_retained_supers[] may therefore reflect
                 * a stale score from a different tree's context. Re-evaluate each retained
                 * tree fresh (spr=FALSE) to get the correct score before displaying.
                 * The landscape map entries for these trees are also updated so the visitedtrees
                 * file reports the same score as the HS output. */
                for(i = 0; i < number_retained_supers && scores_retained_supers[i] != -1; i++)
                    {
                    if(criterion == 0 || criterion == 2 || criterion == 3 || criterion == 6 || criterion == 7)
                        {
                        if(tree_top != NULL) { dismantle_tree(tree_top); tree_top = NULL; }
                        temp_top = NULL;
                        { int _to = 0; tree_build(1, retained_supers[i], tree_top, TRUE, -1, &_to); }
                        tree_top = temp_top;
                        temp_top = NULL;
                        for(k = 0; k < Total_fund_trees; k++) sourcetree_scores[k] = -1;
                        if(criterion == 0)       scores_retained_supers[i] = compare_trees(FALSE);
                        else if(criterion == 2)  scores_retained_supers[i] = compare_trees_sfit(FALSE);
                        else if(criterion == 3)  scores_retained_supers[i] = compare_trees_qfit(FALSE);
                        else if(criterion == 6)  scores_retained_supers[i] = compare_trees_rf(FALSE);
                        else                     scores_retained_supers[i] = compare_trees_ml(FALSE);
                        /* Update the landscape entry so the visitedtrees file score matches the
                         * final reported score rather than the stale-cache value from during search. */
                        if(g_landscape_file[0] && g_landscape_map && tree_top != NULL)
                            {
                            uint64_t rs_hash = tree_topo_hash(tree_top);
                            lm_update_score(g_landscape_map, rs_hash,
                                            ml_display_score(scores_retained_supers[i]));
                            }
                        }
                    }

                /**** Filter retained supertrees: keep only those equal to the best rescored score ****/
                /* Trees added as "equal-score" during NNI may have been scored with stale cached
                 * source-tree scores; after the fresh rescore above some will have higher scores.
                 * Compact the retained_supers[] array, keeping only those at the minimum score. */
                {
                float best_rs = -1.0f;
                int out = 0;
                for(i = 0; i < number_retained_supers && scores_retained_supers[i] != -1; i++)
                    {
                    if(best_rs < 0.0f || scores_retained_supers[i] < best_rs)
                        best_rs = scores_retained_supers[i];
                    }
                for(i = 0; i < number_retained_supers && scores_retained_supers[i] != -1; i++)
                    {
                    if(scores_retained_supers[i] <= best_rs * 1.000001f &&
                       scores_retained_supers[i] >= best_rs * 0.999999f)
                        {
                        if(out != i)
                            {
                            char *tmp_swap = retained_supers[out];
                            retained_supers[out] = retained_supers[i];
                            retained_supers[i] = tmp_swap;
                            scores_retained_supers[out] = scores_retained_supers[i];
                            }
                        out++;
                        }
                    }
                for(; out < number_retained_supers && scores_retained_supers[out] != -1; out++)
                    {
                    retained_supers[out][0] = '\0';
                    scores_retained_supers[out] = -1;
                    }
                }

                /**** Write visited-topology landscape file if requested ******/
                if(g_landscape_file[0] && g_landscape_map)
                    {
                    lm_write(g_landscape_map, g_landscape_file);
                    printf2("Visited topology landscape written to: %s\n", g_landscape_file);
                    printf2("  Unique topologies recorded: %zu\n", g_landscape_map->count);
                    if(g_cluster_enabled)
                        lm_cluster(g_landscape_map, g_cluster_output,
                                   g_cluster_threshold, g_cluster_orderby);
                    lm_free(g_landscape_map);
                    g_landscape_map = NULL;
                    }
                else if(g_cluster_enabled)
                    printf2("Warning: clusterlandscape=yes ignored (visitedtrees= not set)\n");

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
                        if(outfile != NULL) fprintf(outfile, "%s\t[%f]\n", retained_supers[i], ml_display_score(scores_retained_supers[i]) );
                        tree_coordinates(retained_supers[i], FALSE, TRUE, FALSE, -1);
                        printf2("\nSupertree %d of %d %s = %f\n", i+1, j, ml_score_label(), ml_display_score(scores_retained_supers[i]) );
						


						if(do_histogram && criterion == 0)  /* If we are going to draw a histogram of the scores of the source trees when compared to the best supertree */
							{		/* right now we can only do this for mssa (dfit) */
							if(tree_top != NULL)
								{
								dismantle_tree(tree_top);
								tree_top = NULL;
								}
							temp_top = NULL;
							
							{ int _to = 0; tree_build(1, retained_supers[i], tree_top, TRUE, -1, &_to); }
							tree_top = temp_top;
							temp_top = NULL;
							
							/**** evaluate its fit to the source trees in memory *****/
							temptree[0] = '\0';
							for(j=0; j<Total_fund_trees; j++) sourcetree_scores[j] = -1;
							compare_trees(FALSE);  /* calculate the distance from the super tree to all the fundamental trees */
							printf2("\nPlot of source tree scores against supertree number %d", i+1);
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
	
	
    restore_singlecopy_filter(saved_tags);
    free(temptree);
    free(best_tree);
    }




int do_search(char *tree, int user, int print, int maxswaps, FILE *outfile, int numspectries, int numgenetries)
    {
    int swaps = 0, i=0, better_score = TRUE;
	char *temporary_tree = malloc(TREE_LENGTH * sizeof(char));
	char *rep_spr_result = NULL;  /* best topology Newick for this rep after SPR/TBR converges */

	temporary_tree[0] = '\0';
    while(unroottree(tree));  /* fully unroot — one call may not suffice for bifurcating root */
        /****** We now need to build the Supertree in memory *******/
        if(tree_top != NULL)
            {
            dismantle_tree(tree_top);
            tree_top = NULL;
            }
        temp_top = NULL;
        { int _to = 0; tree_build(1, tree, tree_top, user, -1, &_to); }
        tree_top = temp_top;
        temp_top = NULL;
        /* Allocate per-replicate visited-topology hash set */
        visited_set = vs_create(1u << 20);
        /* Initialise periodic progress-line state */
        rep_start_time    = time(NULL);
        hs_do_print       = print;
        last_status_score = -1.0f;
        skip_streak       = 0;
        rep_abandon       = 0;
/*		print_tree(tree_top, temporary_tree);
		printf2("built tree = %s\n", temporary_tree);
  */      if(method == 1)
            {
            /* if NNI is to be carried out */       
            swaps = branchswap(maxswaps,-1, numspectries, numgenetries);
            }
        if(method == 2 || method == 3)
            {
            sprscore = -1;
            tried_regrafts = 0;
            for(i=0; i<Total_fund_trees; i++) sourcetree_scores[i] = -1;  /* invalidate per-rep cache */
            /* if SPR is to be carried out */
/*            while(better_score == TRUE && tried_regrafts < maxswaps && !user_break)
                {
                do
                    {
                    branchpointer = NULL;
                    better_score = spr(tree_top, maxswaps, numspectries, numgenetries);
                    }while(better_score == FALSE && remaining_spr(tree_top) > 0 && tried_regrafts < maxswaps && !user_break);
					printf2("better_score = %d\t, remaining_spr = %d\ttried_regrafts = %d\tuserbreak = %d\n", better_score, remaining_spr(tree_top), tried_regrafts, user_break);
                reset_spr(tree_top);
                }
            swaps = tried_regrafts;
 */
 

			/* Allocate a buffer to capture the best topology Newick for this rep.
			 * spr_new() leaves tree_top in a partially-pruned state after its final
			 * (non-improving) pass — the last branch pruned in its inner loop is never
			 * restored.  We save the valid topology after each pass that finds an
			 * improvement so we can rebuild tree_top cleanly before the NNI pass. */
			rep_spr_result = malloc(TREE_LENGTH * sizeof(char));
			rep_spr_result[0] = '\0';
			print_named_tree(tree_top, rep_spr_result);
			strcat(rep_spr_result, ";");
			while(unroottree(rep_spr_result));

			do/* for a single rep */
				{
				/* Ensure this thread sees the latest user_break value (ARM64 / Apple
				 * Silicon requires an explicit OMP flush for cross-core visibility;
				 * volatile alone only prevents compiler register-caching). */
#ifdef _OPENMP
				#pragma omp flush(user_break)
#endif
				if(user_break) break;   /* bail out before the expensive spr_new() call */
				branchpointer = NULL;
				if(method == 3)
					better_score = tbr_new2(tree_top, maxswaps, numspectries, numgenetries);
				else
					better_score = spr_new3(tree_top, maxswaps, numspectries, numgenetries);
				if(better_score)
					{
					/* tree_top is valid here — save before next pass may corrupt it */
					rep_spr_result[0] = '\0';
					print_named_tree(tree_top, rep_spr_result);
					strcat(rep_spr_result, ";");
					while(unroottree(rep_spr_result));
					}
                /* droprep: only when this SPR pass found no improvement — the rep has
                 * converged to its local optimum and we can meaningfully compare it to
                 * the global best.  Checking mid-climb would drop brand-new reps whose
                 * starting tree is legitimately far from any optimum. */
#ifdef _OPENMP
                if(!better_score && hs_droprep > 0.0f && omp_get_num_threads() > 1 &&
                   sprscore != -1.0f)
                    {
                    float _pbest;
                    #pragma omp atomic read
                    _pbest = par_progress_best;
                    if(_pbest != -1.0f)
                        {
                        float _denom = _pbest > 0.0f ?  _pbest
                                     : _pbest < 0.0f ? -_pbest : 1.0f;
                        if((sprscore - _pbest) / _denom > hs_droprep)
                            rep_abandon = 1;
                        }
                    }
#endif
				}while(better_score == TRUE && tried_regrafts < maxswaps && !user_break && !rep_abandon
			       && (hs_maxskips == 0 || skip_streak < hs_maxskips));

            swaps = tried_regrafts;

            /* Post-SPR/TBR NNI refinement pass: polish the local optimum found
             * by SPR/TBR with fast NNI swaps until no further improvement.
             * NNI is O(N) per round and can tighten hill-climbing after the
             * coarser SPR/TBR search, improving inter-rep convergence. */
            nni_swaps = 0;
            if(!user_break && !rep_abandon)
                {
                /* Rebuild tree_top from the last known-good SPR/TBR topology.
                 * spr_new() leaves tree_top in a partially-pruned state after its
                 * final non-improving pass.  Rebuild from the saved Newick so NNI
                 * operates on the full tree.  Also fixes parent pointers broken
                 * by regraft() (which sets position->parent = NULL). */
                if(tree_top != NULL) { dismantle_tree(tree_top); tree_top = NULL; }
                temp_top = NULL;
                { int _to = 0; tree_build(1, rep_spr_result, tree_top, TRUE, -1, &_to); }
                tree_top = temp_top;
                temp_top = NULL;
                fix_parent_pointers(tree_top, NULL);
                branchpointer = NULL;
                /* branchswap() (NNI) uses its own internal counter; capture it directly */
                nni_swaps = branchswap(maxswaps, -1, numspectries, numgenetries);
                swaps += nni_swaps;
                /* If NNI improved the score, update sprscore so the per-rep
                 * completion report reflects the final post-NNI score. */
                if(scores_retained_supers[0] >= 0.0f &&
                   (sprscore < 0.0f || scores_retained_supers[0] < sprscore))
                    sprscore = scores_retained_supers[0];
                }
			free(rep_spr_result);
			rep_spr_result = NULL;
			
            }
    /* Accumulate this rep's unique topologies into the thread-level accumulator,
     * then free the per-replicate set. */
    if(thread_visited_acc && visited_set) vs_merge(thread_visited_acc, visited_set);
    vs_free(visited_set);
    visited_set = NULL;
    free(temporary_tree);
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


        while(i < number_of_swaps && !user_break
              && (hs_maxskips == 0 || skip_streak < hs_maxskips))
            {
            last_number = number;
            i+= find_swaps(&number, tree_top, number_of_swaps, numspectries, numgenetries);
            last = i;
#ifdef _OPENMP
            if(hs_thread_report_interval > 0 && omp_get_num_threads() > 1)
                {
                time_t _tr_now = time(NULL);
                if(difftime(_tr_now, thread_report_last) >= hs_thread_report_interval)
                    {
                    thread_report_last = _tr_now;
                    double tr_el = difftime(_tr_now, par_search_start);
                    int    tr_em = (int)(tr_el / 60), tr_es = (int)tr_el % 60;
                    float  disp  = (number >= 0.0f) ? number : sprscore;
                    float _pbest4;
                    #pragma omp atomic read
                    _pbest4 = par_progress_best;
                    if(disp < 0.0f)
                        printf2("  [%2d:%02d]  thread %2d  rep %d  score=(searching)  nni_tried=%d  [nni]\n",
                                tr_em, tr_es, omp_get_thread_num(), hs_par_rep, i);
                    else
                        printf2("  [%2d:%02d]  thread %2d  rep %d  %s=  %s%.6f  nni_tried=%d  [nni]\n",
                                tr_em, tr_es, omp_get_thread_num(), hs_par_rep,
                                ml_score_label(),
                                ml_display_score(disp) < ml_display_score(_pbest4) ? "* " : "  ",
                                ml_display_score(disp), i);
                    fflush(stdout);
                    }
                }
#endif
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
	while(count < i && swaps < number_of_swaps && !better_score && !user_break
	         && (hs_maxskips == 0 || skip_streak < hs_maxskips))
		{
                /* choose the pointer sibling to travel down */
                j =(int)fmod(rand_r(&thread_seed), i);  /* the jth pointer sibling is the chosen one */
                while(siblings[j])
                    j =(int)fmod(rand_r(&thread_seed), i);  /* the jth pointer sibling is the chosen one */

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
		 		
		 		
/* Repair parent pointers after SPR/TBR.  spr_new's regraft() sets
   position->parent = NULL for the grafted node (the scoring path
   converts to Newick and doesn't need parent pointers).  Before NNI
   runs directly on tree_top we must fix those so is_ancestor_of()
   works correctly.  Call as fix_parent_pointers(tree_top, NULL). */
static void fix_parent_pointers(struct taxon *pos, struct taxon *parent)
    {
    if(!pos) return;
    while(pos != NULL)
        {
        pos->parent = parent;
        if(pos->daughter)
            fix_parent_pointers(pos->daughter, pos);
        pos = pos->next_sibling;
        }
    }

/* Helper: returns TRUE if ancestor is in the parent-chain of node.
   Used to prevent do_swap from placing a node inside its own subtree,
   which would create a cycle in the tree structure. */
static int is_ancestor_of(struct taxon *ancestor, struct taxon *node)
	{
	struct taxon *pos = node ? node->parent : NULL;
	int _depth = 0, _max = number_of_taxa * 4 + 16;
	while(pos != NULL && _depth++ < _max)
		{
		if(pos == ancestor) return TRUE;
		pos = pos->parent;
		}
	/* If _depth hit _max the parent chain is cyclic — return FALSE (safe: skip the swap/graft) */
	return FALSE;
	}

/* This function actually does the swapping */
/* Part 3 */
void do_swap(struct taxon * first, struct taxon * second)
	{
	struct taxon *f_parent, *f_prev, *f_next;
	struct taxon *s_parent, *s_prev, *s_next;

        interval2 = time(NULL);
        if(difftime(interval2, interval1) > 5) /* every 5 seconds print a dot to the screen */
            {
           /* printf2("="); */
            fflush(stdout);
            interval1 = time(NULL);
            }

	/* Save ALL original pointer values before any modifications.
	   This prevents order-of-operations bugs when the two nodes are adjacent
	   siblings (first->next_sibling == second or second->next_sibling == first),
	   where naive in-place updates create pointer self-loops. */
	f_parent = first->parent;
	f_prev   = first->prev_sibling;
	f_next   = first->next_sibling;
	s_parent = second->parent;
	s_prev   = second->prev_sibling;
	s_next   = second->next_sibling;

	/* === Update EXTERNAL pointers that reference first or second === */
	/* Use saved values only.  Skip updates where the "neighbour" IS the
	   other node (adjacent-sibling case) — those are fixed up below. */

	/* Parent's first-child (daughter) pointer */
	if(f_parent != NULL && f_prev == NULL)
		f_parent->daughter = second;
	if(s_parent != NULL && s_prev == NULL &&
	   !(s_parent == f_parent && f_prev == NULL))   /* avoid double-update when same parent */
		s_parent->daughter = first;

	/* Previous sibling's next pointer */
	if(f_prev != NULL && f_prev != second)
		f_prev->next_sibling = second;
	if(s_prev != NULL && s_prev != first)
		s_prev->next_sibling = first;

	/* Next sibling's prev pointer */
	if(f_next != NULL && f_next != second)
		f_next->prev_sibling = second;
	if(s_next != NULL && s_next != first)
		s_next->prev_sibling = first;

	/* === Update first's own pointers to second's old position === */
	first->parent       = s_parent;
	first->prev_sibling = s_prev;
	first->next_sibling = s_next;

	/* === Update second's own pointers to first's old position === */
	second->parent       = f_parent;
	second->prev_sibling = f_prev;
	second->next_sibling = f_next;

	/* === Fix up cross-references when the two nodes were adjacent siblings ===
	   After the own-pointer swap above, adjacent nodes point at each other
	   using the OLD saved value.  Correct those now. */
	if(f_next == second)
		{
		/* first was immediately before second.
		   After swap: second occupies first's old slot, first occupies second's.
		   They are still adjacent but in the opposite order. */
		second->next_sibling = first;
		first->prev_sibling  = second;
		}
	else if(s_next == first)
		{
		/* second was immediately before first — symmetric case. */
		first->next_sibling  = second;
		second->prev_sibling = first;
		}

	/* === Handle tree_top change === */
	if(tree_top == first)
		tree_top = second;
	else if(tree_top == second)
		tree_top = first;
	}
	


/* this is the function which will allow the swapping to take place x number of steps away, this will make the branch swapping that much more robust. The function from its start point will travel UP from its present position and then will travel up or down the tree (but never back to its start point) looking for viable swaps x number of steps away */

int swapper(struct taxon * position,struct taxon * prev_pos, int stepstaken, struct taxon * first_swap, struct taxon * second_swap, float * number, int * swaps, int number_of_swaps, int numspectries, int numgenetries)
	{
	
  	float  distance = *number;
	int better_score = FALSE, i=0, j=0, different = TRUE;
	struct taxon *start2 = NULL, *start1 = NULL;
        uint64_t topo_h_nni = 0;
        char *best_tree = NULL, *temptree = NULL;
        
        temptree = malloc(TREE_LENGTH*sizeof(char));
        if(!temptree) memory_error(66);
        temptree[0] = '\0';
        best_tree = malloc(TREE_LENGTH*sizeof(char));
        if(!best_tree) memory_error(52);
        best_tree[0] = '\0';
	
	{ int _g = 0, _gmax = number_of_taxa * 4 + 16;
	  while(position->prev_sibling != NULL && _g++ < _gmax) position = position->prev_sibling; /* rewind to start of this level */
	  if(_g >= _gmax) rep_abandon = 1; }
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
			
		if(second_swap == NULL) printf2("second_swap is null!\n");
		{ int _g = 0, _gmax = number_of_taxa * 4 + 16;
		  while(second_swap->prev_sibling != NULL && _g++ < _gmax) second_swap = second_swap->prev_sibling; /* rewinding */
		  if(_g >= _gmax) rep_abandon = 1; }
		
		start2 = second_swap; /* start2 now points to the start levle to be swapped against */
	
                /* now first_swap is at the first sibling at this level, and second_swap is at the first sibling at the other level */
                /* The next step is to swap every sibling at this level with every sibling at the level above */
		while(first_swap != NULL && *swaps < number_of_swaps && !better_score && !user_break && !rep_abandon
		          && (hs_maxskips == 0 || skip_streak < hs_maxskips))
			{
			while(second_swap != NULL && *swaps < number_of_swaps && !better_score && !user_break && !rep_abandon
			          && (hs_maxskips == 0 || skip_streak < hs_maxskips))
				{
				if(second_swap->daughter != prev_pos && !is_ancestor_of(second_swap, first_swap) && !is_ancestor_of(first_swap, second_swap))   /* Do not swap if either node is an ancestor of the other (would create a cycle) */
					{
					do_swap(first_swap, second_swap);

                                        if(check_taxa(tree_top) == number_of_taxa)
                                            {
                                            topo_h_nni = tree_topo_hash(tree_top);
                                            if(visited_set != NULL && vs_contains(visited_set, topo_h_nni))
                                                {
                                                /* Already-visited NNI topology: count skip, record landscape, revert swap. */
                                                skip_streak++;
                                                if(g_landscape_file[0])
                                                    lm_record(landscape_map ? landscape_map : g_landscape_map,
                                                              topo_h_nni, 0.0f, NULL);
                                                do_swap(first_swap, second_swap);
                                                }
                                            else
                                                {
                                                if(visited_set != NULL) { vs_insert(visited_set, topo_h_nni); skip_streak = 0; }
											(*swaps)++;
                                           /* *swaps = *swaps + 1; */ /* count the number of swaps done */
											NUMSWAPS++;
                                            strcpy(best_tree, "");
                                            print_named_tree(tree_top, best_tree);
                                            strcat(best_tree, ";");
                                            while(unroottree(best_tree));
                                            /*  printf2("!:assigned tree: %s\n", best_tree);  */
                                            
                                            
                                            if(criterion == 0) distance = compare_trees(FALSE);  /* calculate the distance from the super tree to all the fundamental trees */
                                            if(criterion == 2) distance = compare_trees_sfit(FALSE);
                                            if(criterion == 3) distance = compare_trees_qfit(FALSE);

											if(criterion==5)
												{
												strcpy(temptree, "");
												print_tree(tree_top, temptree);
												strcat(temptree, ";");
												while(unroottree(temptree));
												distance = get_recon_score(temptree, numspectries, numgenetries);
												}

											if(criterion==6)
												{
												distance = compare_trees_rf(FALSE);
												}

											if(criterion==7)
												{
												distance = compare_trees_ml(FALSE);
												}

                                            /* Landscape recording: NNI-scored topology */
                                            if(g_landscape_file[0])
                                                {
                                                lm_record(landscape_map ? landscape_map : g_landscape_map,
                                                          topo_h_nni, ml_display_score(distance), best_tree);
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
#ifdef _OPENMP
                                                        /* NNI parallel progress: only when THIS rep set the new
                                                         * global best — the rep number then matches the rep that
                                                         * actually achieved the displayed score. */
                                                        if(omp_get_num_threads() > 1)
                                                            {
                                                            time_t _now = time(NULL);
                                                            int new_global_best = 0;
                                                            #pragma omp critical (par_progress)
                                                                {
                                                                if(par_progress_best < 0.0f || distance < par_progress_best)
                                                                    {
                                                                    par_progress_best = distance;
                                                                    new_global_best = 1;
                                                                    }
                                                                }
                                                            if(new_global_best &&
                                                               (hs_progress_interval == 0 ||
                                                                difftime(_now, par_last_progress_time) >= hs_progress_interval))
                                                                {
                                                                double hs_el = difftime(_now, par_search_start);
                                                                int    hs_em = (int)(hs_el / 60), hs_es = (int)hs_el % 60;
                                                                int    hs_imp = (par_last_print_score < 0.0f ||
                                                                                 par_progress_best < par_last_print_score);
                                                                if(hsprint)
                                                                    {
                                                                    printf2("  [%2d:%02d]  best so far = %s%.6f\tnni\trep %d\n",
                                                                            hs_em, hs_es,
                                                                            hs_imp ? "* " : "  ", ml_display_score(par_progress_best), hs_par_rep);
                                                                    par_last_print_score = par_progress_best;
                                                                    par_last_progress_time = _now;
                                                                    fflush(stdout);
                                                                    }
                                                                }
                                                            }
#endif
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
                                                }  /* end else (new topology) */
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
    char filename[1000];
    char *best_tree = malloc(TREE_LENGTH * sizeof(char));
    FILE *yaptpfile = NULL;

    if(!best_tree) { printf2("Error: out of memory in yaptp_search\n"); return; }
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
                printf2("Error: '%s' is an invalid number of repetitions\n", parsed_command[i+1]);
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
			printf2("Using the citerion recon it is only possible to use the markovian method for randomisation\n");
			yaptp_method = 2;
			}
		}
            else
                {
                if(strcmp(parsed_command[i+1], "markovian") == 0)
                    yaptp_method =2;
                else
                    {
                    printf2("Error: method option '%s' unknown\n", parsed_command[i+1]);
                        error = TRUE;
                    }
                }
            }
            
        if(strcmp(parsed_command[i], "search") == 0)
            {
            if(strcmp(parsed_command[i+1], "all") == 0) search = 0;
            else
                {
                if(strcmp(parsed_command[i+1], "nni") == 0) { search = 1; method = 1; }
                else
                    {
                    if(strcmp(parsed_command[i+1], "spr") == 0) { search = 1; method = 2; }
                    else
                        {
                        printf2("Error: search option '%s' unknown\n", parsed_command[i+1]);
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
                printf2("Error: Unable to open file named '%s'\n", parsed_command[i+1]);
                error = TRUE;
                strcpy(filename, "yaptp.txt");
                }
            }
        }
        
    if(!error)
        {
		
        printf2("\n\nYAPTP Settings:\n\tNumber of YAPTP repetitions: %d\n\tSearching Supertree-space: ", Nreps);
        if(search == 0)printf2("Exhaustive Search of Supertree-spce\n");
        if(search == 1)printf2("Heuristic Search\n");
        printf2("\tRandomisation method (method) = ");
        if(yaptp_method==1) printf2("equiprobable\n");
        if(yaptp_method==2) printf2("markovian\n");
        printf2("\tYAPTP output file: yaptp.txt\n\n");
        
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

			if(print_log == FALSE) 
				strcpy(system_call, "paup coding.nex");
			else
				{
				strcpy(system_call, "paup coding.nex");
				strcat(system_call, " | tee -a ");
				strcat(system_call, logfile_name);
				}

			if(system(system_call) != 0) printf2("Error calling paup, please execute the file coding.nex in paup to perform the parsimony step\n");

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
            printf2("YAPTP progress indicator: "); 
			GC = 0;
            for(i=0; i<Nreps; i++)  /* for all reps of the yaptp algorithm */
                {
				if(!user_break)
					{
					if(i>0)hsprint=FALSE;
					printf2("\n\trepetition %d:=", i+1);
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
							{ int _to = 0; tree_build(1, retained_supers[k], tree_top, FALSE, -1, &_to); }
							tree_top = temp_top;
							temp_top = NULL;
						
							strcpy(best_tree, "");
							print_named_tree(tree_top, best_tree);
							/*printf("%s;\t[%f]\n", best_tree, scores_retained_supers[k]);  */
							fprintf(yaptpfile, "%s\t[%f %d]\n", best_tree, scores_retained_supers[k], i);
							}
						else
							{
						/* printf2("%s;\t[%f]\n", retained_supers[k], scores_retained_supers[k]);  */
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
    free(best_tree);
    hsprint=TRUE;
    }


void randomise_tree(char *tree)
    {
    int i=0, j=0, k=0, l=0, x=0, y=0, treecount = 0, random=0, supers = 0, actual_num = 0;
    char **array = NULL, *temptree = malloc(TREE_LENGTH * sizeof(char)), *newtree = NULL, *tmp;
    if(!temptree) { printf2("Error: out of memory in randomise_tree\n"); return; }
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
    free(temptree);
    free(tmp);
    }

void randomise_taxa(char *tree)
    {
    int i=0, j=0, k=0, l=0, x=0, y=0, treecount = 0, random=0, tottax;
    char **array = NULL, *temptree = malloc(TREE_LENGTH * sizeof(char));
    if(!temptree) { printf2("Error: out of memory in randomise_taxa\n"); return; }

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
    free(temptree);
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
    char *tree1 = malloc(TREE_LENGTH * sizeof(char));
    char *tree2 = malloc(TREE_LENGTH * sizeof(char));
    char *name = NULL;

    if(!tree1 || !tree2) { free(tree1); free(tree2); memory_error(62); }

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
    free(scores1);
    free(scores2);
    free(name);
    free(tree1); free(tree2);
    return(different);
    }


/* Next is the code needed for the Baum/ragan coding scheme */




 

/* this function calculates the score of a supertree to a set of source trees using the criterion of Matrix representation using Compatibility */
/* This is analogous to a clique analysis */
                         


/* this function calculates the score of a supertree to a set of source trees using the criterion of Matrix representation using Compatibility */
/* This is analogous to a clique analysis */









/* -------------------------------------------------------------------------
 * TBR helpers — string-based rerooting
 *
 * These three static functions implement the new TBR rerooting path.
 * The key invariant: every tree passed to regraft() is built by tree_build()
 * from a Newick string, guaranteeing it is cycle-free.  reroot_tree() is
 * never called.
 * ----------------------------------------------------------------------- */

/* true_parent: return the actual structural parent of a node.
 *
 * Clann tree invariant: only the FIRST sibling at each level has its
 * parent pointer set; all subsequent siblings have parent == NULL.
 * Walk back via prev_sibling to the first sibling, then return its parent.
 * For genuine root-level nodes (first sibling at root also has parent==NULL)
 * this correctly returns NULL. */
static struct taxon *true_parent(struct taxon *node)
	{
	while(node->prev_sibling != NULL) node = node->prev_sibling;
	return node->parent;
	}


/* serialise_complement: write the Newick fragment for every node in the tree
 * EXCEPT `node` and its descendants.  `root_node` is the first sibling at
 * the top level of the freshly-built working tree.
 * Returns 0 on success, -1 if the depth guard fires (sets rep_abandon). */
static int serialise_complement(struct taxon *node, struct taxon *root_node,
                                char *buf, int depth, int maxdepth)
	{
	struct taxon *child = NULL, *sib = NULL;
	int first = 1;
	struct taxon *parent = NULL;

	if(depth > maxdepth) { rep_abandon = 1; return -1; }

	/* Use true_parent() — clann invariant: only first siblings carry parent */
	parent = true_parent(node);

	if(parent == NULL)
		{
		/* node is at root level — emit all root siblings except node */
		sib = root_node;
		while(sib != NULL)
			{
			if(sib != node)
				{
				if(!first) strcat(buf, ",");
				first = 0;
				print_single_subtree(sib, buf);
				}
			sib = sib->next_sibling;
			}
		return 0;
		}

	/* General case: wrap children-of-parent (minus path node) + complement above */
	strcat(buf, "(");

	child = parent->daughter;
	while(child != NULL)
		{
		if(child != node)
			{
			if(!first) strcat(buf, ",");
			first = 0;
			print_single_subtree(child, buf);
			}
		child = child->next_sibling;
		}

	if(true_parent(parent) != NULL)
		{
		/* recurse upward for complement of parent */
		if(!first) strcat(buf, ",");
		if(serialise_complement(parent, root_node, buf, depth + 1, maxdepth) != 0)
			{
			strcat(buf, ")");
			return -1;
			}
		}
	else
		{
		/* parent is at root level — emit parent's root-level siblings */
		sib = root_node;
		while(sib != NULL)
			{
			if(sib != parent)
				{
				if(!first) strcat(buf, ",");
				first = 0;
				print_single_subtree(sib, buf);
				}
			sib = sib->next_sibling;
			}
		}

	strcat(buf, ")");
	return 0;
	}


/* newick_reroot_at_tag: given an unrooted Newick string and node tag q
 * (assigned by number_tree()), produce a rerooted Newick with the edge to
 * node q as the root edge.  Uses an independent working copy built by
 * tree_build() — no pointer surgery on the caller's tree.
 * Returns 0 on success, -1 on failure (caller should set rep_abandon). */
static int newick_reroot_at_tag(const char *base_nwk, int q, char *out_nwk)
	{
	struct taxon *root = NULL, *target = NULL;
	char *working = NULL;
	int status;

	working = malloc(TREE_LENGTH * sizeof(char));
	if(!working) return -1;
	strcpy(working, base_nwk);

	/* Build an independent, cycle-free working copy */
	temp_top = NULL;
	{ int _to = 0; tree_build(1, working, NULL, FALSE, -1, &_to); }
	root = temp_top;
	temp_top = NULL;
	free(working);

	if(root == NULL) return -1;

	number_tree(root, 0);

	target = get_branch(root, q);
	if(target == NULL) { dismantle_tree(root); return -1; }

	/* Produce rerooted Newick: (target_subtree, complement) */
	out_nwk[0] = '\0';
	strcat(out_nwk, "(");
	print_single_subtree(target, out_nwk);
	strcat(out_nwk, ",");

	status = serialise_complement(target, root, out_nwk, 0, number_of_taxa + 4);

	strcat(out_nwk, ");");

	dismantle_tree(root);

	if(status != 0) return -1;

	/* Normalise to unrooted form */
	while(unroottree(out_nwk));

	return 0;
	}


/* bisect_newick: split tree at branch tag x, producing two valid unrooted
 * Newick strings without any pointer surgery on the input tree.
 *
 * subtree_nwk  — the pruned subtree rooted at the node with tag x
 * remaining_nwk — everything else (complement), normalised to unrooted form
 *
 * Returns:  0 = success
 *          -1 = skip (x is the root node — no meaningful bisection)
 *          -2 = serialise_complement depth overflow (sets rep_abandon)
 *
 * Both output buffers must be at least TREE_LENGTH bytes. */
static int bisect_newick(struct taxon *tree_root, int x,
                         char *subtree_nwk, char *remaining_nwk)
	{
	struct taxon *position = NULL;
	int rc = 0, wrap = 0;

	position = get_branch(tree_root, x);
	if(position == NULL || position == tree_root) return -1;

	/* --- pruned subtree --- */
	subtree_nwk[0] = '\0';
	print_single_subtree(position, subtree_nwk);
	strcat(subtree_nwk, ";");
	/* Unroot the subtree Newick.  Use a change-tracking loop instead of
	 * while(unroottree(...)) to avoid infinite loops when unroottree() returns
	 * TRUE but makes no actual change (e.g. for a 2-leaf subtree "(A,B);"). */
	if(strchr(subtree_nwk, '(') != NULL)
		{
		char *_prev = malloc(TREE_LENGTH * sizeof(char));
		if(_prev)
			{
			do { strcpy(_prev, subtree_nwk); unroottree(subtree_nwk); }
			while(strcmp(_prev, subtree_nwk) != 0);
			free(_prev);
			}
		}

	/* --- remaining (complement) ---
	 * serialise_complement() output format depends on position's level:
	 *   non-root (parent != NULL): already wrapped in "(...)" → just add ";"
	 *   root-level (true_parent()==NULL): comma-separated list → wrap in "(...);" */
	remaining_nwk[0] = '\0';
	wrap = (true_parent(position) == NULL);
	if(wrap) strcat(remaining_nwk, "(");
	rc = serialise_complement(position, tree_root, remaining_nwk, 0, number_of_taxa + 4);
	if(rc != 0) { rep_abandon = 1; return -2; }
	if(wrap) strcat(remaining_nwk, ")");
	strcat(remaining_nwk, ";");
	if(strchr(remaining_nwk, '(') != NULL)
		{
		char *_prev = malloc(TREE_LENGTH * sizeof(char));
		if(_prev)
			{
			do { strcpy(_prev, remaining_nwk); unroottree(remaining_nwk); }
			while(strcmp(_prev, remaining_nwk) != 0);
			free(_prev);
			}
		}

	return 0;
	}


/* tbr_new: new TBR search driver.
 *
 * Uses bisect_newick() + tree_build() for both the pruned subtree and the
 * remaining tree — zero pointer surgery, guaranteed cycle-free trees into
 * regraft().  The inner rerooting loop uses newick_reroot_at_tag().
 *
 * Thread safety: all variables are stack-local or threadprivate globals.
 * ctrl+C: user_break checked at every outer- and inner-loop iteration.
 * rep_abandon: propagated from bisect_newick(), newick_reroot_at_tag(),
 * and regraft(). */
static int tbr_new(struct taxon *master, int maxswaps, int numspectries, int numgenetries)
	{
	struct taxon *newbie = NULL, *remaining_tree = NULL;
	struct taxon *master_copy = NULL, *position = NULL, *temper = NULL;
	int better_score = FALSE;
	int x = 0, y = 0, q = 0, r = 0, i = 0;
	int brc = 0, reroot_status = 0;
	char *base_nwk     = NULL;   /* pruned subtree Newick (bisect output / canonical TBR base) */
	char *remaining_nwk = NULL;  /* remaining tree Newick (bisect output) */
	char *rerooted_nwk  = NULL;  /* rerooted subtree Newick (inner loop) */

	base_nwk      = malloc(TREE_LENGTH * sizeof(char));
	remaining_nwk = malloc(TREE_LENGTH * sizeof(char));
	rerooted_nwk  = malloc(TREE_LENGTH * sizeof(char));
	if(!base_nwk || !remaining_nwk || !rerooted_nwk)
		{ free(base_nwk); free(remaining_nwk); free(rerooted_nwk); return FALSE; }

	/* Make a master copy to restore from each outer iteration */
	temp_top = NULL;
	duplicate_tree(master, NULL);
	master_copy = temp_top;
	temp_top = NULL;

	y = number_tree(master, 0);

	for(x = 0; x < y; x++)
		{
#ifdef _OPENMP
		if(hs_thread_report_interval > 0 && omp_get_num_threads() > 1)
			{
			time_t _tr_now = time(NULL);
			if(difftime(_tr_now, thread_report_last) >= hs_thread_report_interval)
				{
				thread_report_last = _tr_now;
				double tr_el = difftime(_tr_now, par_search_start);
				int    tr_em = (int)(tr_el / 60), tr_es = (int)tr_el % 60;
				float  _pbest_tbr;
				#pragma omp atomic read
				_pbest_tbr = par_progress_best;
				if(sprscore < 0.0f)
					printf2("  [%2d:%02d]  thread %2d  rep %d  score=(searching)  tried=%d  skips=%d  [tbr_new branch %d/%d]\n",
					        tr_em, tr_es, omp_get_thread_num(), hs_par_rep,
					        tried_regrafts, skip_streak, x, y);
				else
					printf2("  [%2d:%02d]  thread %2d  rep %d  %s=  %s%.6f  tried=%d  skips=%d  [tbr_new branch %d/%d]\n",
					        tr_em, tr_es, omp_get_thread_num(), hs_par_rep,
					        ml_score_label(),
					        ml_display_score(sprscore) < ml_display_score(_pbest_tbr) ? "* " : "  ",
					        ml_display_score(sprscore), tried_regrafts, skip_streak, x, y);
				fflush(stdout);
				}
			}
#endif

		/* Restore fresh master from master_copy */
		if(master != NULL) dismantle_tree(master);
		master = NULL;
		temp_top = NULL;
		duplicate_tree(master_copy, NULL);
		master = temp_top;
		tree_top = master;
		temp_top = NULL;
		number_tree(master, 0);

		if(rep_abandon) { x = y; continue; }

		/* --- String-based bisect: no pointer surgery on master --- */
		position = get_branch(master, x);
		brc = bisect_newick(master, x, base_nwk, remaining_nwk);

		/* Dismantle master — Newick strings have everything we need */
		dismantle_tree(master); master = NULL; tree_top = NULL;

		if(brc == -1) continue;          /* root node — skip */
		if(brc == -2) { x = y; continue; } /* depth overflow — rep_abandon already set */

		/* Build fresh, cycle-free remaining tree */
		temp_top = NULL;
		{ int _to = 0; tree_build(1, remaining_nwk, NULL, FALSE, -1, &_to); }
		remaining_tree = temp_top; temp_top = NULL;
		if(remaining_tree == NULL) { rep_abandon = 1; x = y; continue; }
		tree_top = remaining_tree;

		/* Build initial newbie from subtree Newick.
		 * For a bare leaf ("N;"), tree_build(1,...) starts at ';' and returns
		 * name=-1.  Detect and create the leaf directly in that case. */
		if(strchr(base_nwk, '(') == NULL)
			{
			newbie = make_taxon();
			newbie->name = atoi(base_nwk);
			}
		else
			{
			temp_top = NULL;
			{ int _to = 0; tree_build(1, base_nwk, NULL, FALSE, -1, &_to); }
			newbie = temp_top; temp_top = NULL;
			}
		if(newbie == NULL) { rep_abandon = 1; goto tbr_cleanup; }

		r = number_tree(newbie, 0);

		for(i = 0; i < number_of_taxa; i++) presenceof_SPRtaxa[i] = FALSE;
		identify_taxa(newbie, presenceof_SPRtaxa);

		if(r <= 4)
			{
			/* Subtree too small to meaningfully reroot — single SPR-style graft */
			temper = make_taxon();
			temper->daughter = newbie;
			newbie->parent = temper;
			newbie = temper;
			temper = NULL;
			newbie->spr = TRUE;
			better_score = regraft(remaining_tree, newbie, NULL, 1, maxswaps, numspectries, numgenetries);
			}
		else
			{
			/* Re-serialise newbie after tree_build() so base_nwk has canonical
			 * tag assignment consistent with newick_reroot_at_tag(). */
			r = number_tree(newbie, 0);
			base_nwk[0] = '\0';
			print_single_subtree(newbie, base_nwk);
			strcat(base_nwk, ";");
			while(unroottree(base_nwk));

			/* Discard the build-for-count newbie — inner loop builds fresh
			 * per-rerooting copies from rerooted_nwk. */
			dismantle_tree(newbie); newbie = NULL;

			/* --- Inner loop: try each rerooting q against remaining_tree --- */
			for(q = 0; q < r; q++)
				{
				if(user_break || rep_abandon) break;

				reroot_status = newick_reroot_at_tag(base_nwk, q, rerooted_nwk);
				if(reroot_status != 0) { rep_abandon = 1; break; }

				temp_top = NULL;
				{ int _to = 0; tree_build(1, rerooted_nwk, NULL, FALSE, -1, &_to); }
				newbie = temp_top; temp_top = NULL;
				if(newbie == NULL) { rep_abandon = 1; break; }

				temper = make_taxon();
				temper->daughter = newbie;
				newbie->parent = temper;
				newbie = temper;
				temper = NULL;
				newbie->spr = TRUE;

				/* regraft() restores remaining_tree on FALSE return — safe to reuse */
				better_score = regraft(remaining_tree, newbie, NULL, 1, maxswaps, numspectries, numgenetries);

				if(better_score || user_break || rep_abandon) break;

				dismantle_tree(newbie);
				newbie = NULL;
				}
			}

tbr_cleanup:
		if(!better_score && newbie != NULL) { dismantle_tree(newbie); newbie = NULL; }
		if(!better_score && remaining_tree != NULL)
			{ dismantle_tree(remaining_tree); remaining_tree = NULL; tree_top = NULL; }

		if(better_score || (position == branchpointer) || tried_regrafts >= maxswaps
		   || user_break || rep_abandon
		   || (hs_maxskips > 0 && skip_streak >= hs_maxskips))
			x = y;
		}

	free(base_nwk);
	free(remaining_nwk);
	free(rerooted_nwk);
	if(master_copy != NULL) { dismantle_tree(master_copy); master_copy = NULL; }
	return better_score;
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

/* spr_new2: SPR search driver.
 *
 * Uses bisect_newick() + tree_build() for both the pruned subtree and the
 * remaining tree — zero pointer surgery, guaranteed cycle-free trees into
 * regraft(), and correct parent pointers throughout the remaining tree so
 * regraft() can traverse all internal branches.
 *
 * All additions preserved: OpenMP thread reporting, user_break, rep_abandon,
 * skip_streak/hs_maxskips convergence guard, branchpointer early-exit,
 * per-iteration master restore from master_copy, identify_taxa() optimisation.
 *
 * spr_new() is retained as dead code; do_search() calls spr_new2() for SPR. */
static int spr_new2(struct taxon *master, int maxswaps, int numspectries, int numgenetries)
	{
	struct taxon *newbie = NULL, *remaining_tree = NULL;
	struct taxon *temper = NULL, *master_copy = NULL, *position = NULL;
	int better_score = FALSE;
	int x = 0, y = 0, i = 0, brc = 0;
	char *subtree_nwk  = NULL;
	char *remaining_nwk = NULL;

	subtree_nwk  = malloc(TREE_LENGTH * sizeof(char));
	remaining_nwk = malloc(TREE_LENGTH * sizeof(char));
	if(!subtree_nwk || !remaining_nwk)
		{ free(subtree_nwk); free(remaining_nwk); return FALSE; }

	/* Make a master copy to restore from each outer iteration */
	temp_top = NULL;
	duplicate_tree(master, NULL);
	master_copy = temp_top;
	temp_top = NULL;

	y = number_tree(master, 0);

	for(x = 0; x < y; x++)
		{
#ifdef _OPENMP
		if(hs_thread_report_interval > 0 && omp_get_num_threads() > 1)
			{
			time_t _tr_now = time(NULL);
			if(difftime(_tr_now, thread_report_last) >= hs_thread_report_interval)
				{
				thread_report_last = _tr_now;
				double tr_el = difftime(_tr_now, par_search_start);
				int    tr_em = (int)(tr_el / 60), tr_es = (int)tr_el % 60;
				float  _pbest_spr2;
				#pragma omp atomic read
				_pbest_spr2 = par_progress_best;
				if(sprscore < 0.0f)
					printf2("  [%2d:%02d]  thread %2d  rep %d  score=(searching)  tried=%d  skips=%d  [spr_new2 branch %d/%d]\n",
					        tr_em, tr_es, omp_get_thread_num(), hs_par_rep,
					        tried_regrafts, skip_streak, x, y);
				else
					printf2("  [%2d:%02d]  thread %2d  rep %d  %s=  %s%.6f  tried=%d  skips=%d  [spr_new2 branch %d/%d]\n",
					        tr_em, tr_es, omp_get_thread_num(), hs_par_rep,
					        ml_score_label(),
					        ml_display_score(sprscore) < ml_display_score(_pbest_spr2) ? "* " : "  ",
					        ml_display_score(sprscore), tried_regrafts, skip_streak, x, y);
				fflush(stdout);
				}
			}
#endif

		/* Restore fresh master from master_copy */
		if(master != NULL) dismantle_tree(master);
		master = NULL;
		temp_top = NULL;
		duplicate_tree(master_copy, NULL);
		master = temp_top;
		tree_top = master;
		temp_top = NULL;
		number_tree(master, 0);

		if(rep_abandon) { x = y; continue; }

		/* --- String-based bisect: no pointer surgery on master --- */
		position = get_branch(master, x);
		brc = bisect_newick(master, x, subtree_nwk, remaining_nwk);

		/* Dismantle master — Newick strings have everything we need */
		dismantle_tree(master); master = NULL; tree_top = NULL;

		if(brc == -1) continue;           /* root node — skip */
		if(brc == -2) { x = y; continue; } /* depth overflow — rep_abandon already set */

		/* Build fresh, cycle-free remaining tree */
		temp_top = NULL;
		{ int _to = 0; tree_build(1, remaining_nwk, NULL, FALSE, -1, &_to); }
		remaining_tree = temp_top; temp_top = NULL;
		if(remaining_tree == NULL) { rep_abandon = 1; x = y; continue; }
		tree_top = remaining_tree;

		/* Build fresh, cycle-free newbie.
		 * tree_build(1,...) skips the leading '(' — for a bare leaf string like
		 * "3;" there is no '(', so it starts at ';' and produces name=-1.
		 * Detect this case and create the leaf node directly. */
		if(strchr(subtree_nwk, '(') == NULL)
			{
			newbie = make_taxon();
			newbie->name = atoi(subtree_nwk);  /* "3;" → 3 */
			}
		else
			{
			temp_top = NULL;
			{ int _to = 0; tree_build(1, subtree_nwk, NULL, FALSE, -1, &_to); }
			newbie = temp_top; temp_top = NULL;
			}
		if(newbie == NULL) { rep_abandon = 1; goto spr2_cleanup; }

		for(i = 0; i < number_of_taxa; i++) presenceof_SPRtaxa[i] = FALSE;
		identify_taxa(newbie, presenceof_SPRtaxa);

		/* Wrap newbie in SPR sentinel and regraft */
		temper = make_taxon();
		temper->daughter = newbie;
		newbie->parent = temper;
		newbie = temper;
		temper = NULL;
		newbie->spr = TRUE;

		better_score = regraft(remaining_tree, newbie, NULL, 1, maxswaps, numspectries, numgenetries);

spr2_cleanup:
		if(!better_score && newbie != NULL) { dismantle_tree(newbie); newbie = NULL; }
		if(!better_score && remaining_tree != NULL)
			{ dismantle_tree(remaining_tree); remaining_tree = NULL; tree_top = NULL; }

		if(better_score || (position == branchpointer) || tried_regrafts >= maxswaps
		   || user_break || rep_abandon
		   || (hs_maxskips > 0 && skip_streak >= hs_maxskips))
			x = y;
		}

	free(subtree_nwk);
	free(remaining_nwk);
	if(master_copy != NULL) { dismantle_tree(master_copy); master_copy = NULL; }
	return better_score;
	}


/* =======================================================================
 * evaluate_candidate: score a candidate tree topology.
 *
 * Builds tree_top from candidate_nwk (integer-label Newick produced by
 * spr_graft()), scores it, and updates all bookkeeping (BESTSCORE,
 * sprscore, scores_retained_supers[], retained_supers[], visited_set,
 * OMP parallel-progress).
 *
 * tmp_fund_scores must be a caller-allocated float[Total_fund_trees]
 * scratch buffer (so we don't malloc/free on every evaluation call).
 *
 * On return:
 *   TRUE  — better score found; tree_top valid (caller owns it).
 *   FALSE — not an improvement; tree_top dismantled and set to NULL.
 *
 * Thread safety: all state is threadprivate or stack-local; no shared
 * mutable state accessed outside the OMP critical(par_progress) guard.
 * ======================================================================= */
/* =======================================================================
 * probe_candidate: score a candidate tree without committing to it.
 *
 * Used by the best-improvement strategy to survey the full SPR/TBR
 * neighbourhood before selecting the single best move.
 *
 * Differences from evaluate_candidate:
 *   - Does NOT insert the topology into visited_set (deferred to commit).
 *   - Always restores sourcetree_scores[] on return.
 *   - Does NOT update sprscore, BESTSCORE, retained_supers, or landscape.
 *   - Does NOT leave tree_top set on return.
 *
 * Returns the candidate's score, or -1.0f if it should be skipped
 * (already visited, build failure, wrong taxon count, user break).
 * ======================================================================= */
static float probe_candidate(const char *candidate_nwk,
                              float      *tmp_fund_scores,
                              int         numspectries,
                              int         numgenetries)
	{
	float    tmpscore = -1.0f;
	uint64_t topo_h;
	char    *temptree = NULL;
	int      i;

	if(tree_top != NULL) { dismantle_tree(tree_top); tree_top = NULL; }
	temp_top = NULL;
	{ int _to = 0; tree_build(1, (char *)candidate_nwk, NULL, FALSE, -1, &_to); }
	tree_top = temp_top;
	temp_top = NULL;
	if(tree_top == NULL) return -1.0f;

	if(check_taxa(tree_top) != number_of_taxa || user_break)
		{ dismantle_tree(tree_top); tree_top = NULL; return -1.0f; }

	/* Check visited-set but do NOT insert — insertion is deferred to the
	 * evaluate_candidate commit step so the winner can be re-processed. */
	topo_h = tree_topo_hash(tree_top);
	if(visited_set != NULL && vs_contains(visited_set, topo_h))
		{
		skip_streak++;
		dismantle_tree(tree_top); tree_top = NULL;
		return -1.0f;
		}

	/* New topology: count it and reset streak */
	skip_streak = 0;
	tried_regrafts++;
	NUMSWAPS++;

	/* Save sourcetree_scores — always restored before return */
	for(i = 0; i < Total_fund_trees; i++) tmp_fund_scores[i] = sourcetree_scores[i];

	/* Score by criterion (mirrors evaluate_candidate) */
	if(criterion == 0)
		tmpscore = compare_trees(FALSE);
	else if(criterion == 2)
		tmpscore = compare_trees_sfit(FALSE);
	else if(criterion == 3)
		tmpscore = compare_trees_qfit(FALSE);
	else if(criterion == 5)
		{
		temptree = malloc(TREE_LENGTH * sizeof(char));
		if(temptree)
			{
			print_tree(tree_top, temptree); strcat(temptree, ";");
			while(unroottree(temptree));
			tmpscore = get_recon_score(temptree, numspectries, numgenetries);
			free(temptree);
			}
		}
	else if(criterion == 6)
		tmpscore = compare_trees_rf(TRUE);
	else if(criterion == 7)
		tmpscore = compare_trees_ml(TRUE);

	/* Always restore sourcetree_scores — no permanent state change */
	for(i = 0; i < Total_fund_trees; i++) sourcetree_scores[i] = tmp_fund_scores[i];

	dismantle_tree(tree_top); tree_top = NULL;
	return tmpscore;
	}


static int evaluate_candidate(const char *candidate_nwk,
                               float      *tmp_fund_scores,
                               int         numspectries,
                               int         numgenetries)
	{
	int   better_score = FALSE;
	int   i, j, different;
	float tmpscore = 0.0f;
	char *best_tree = NULL;
	char *temptree  = NULL;
	uint64_t topo_h;

	best_tree = malloc(TREE_LENGTH * sizeof(char));
	temptree  = malloc(TREE_LENGTH * sizeof(char));
	if(!best_tree || !temptree)
		{ free(best_tree); free(temptree); return FALSE; }
	best_tree[0] = '\0';
	temptree[0]  = '\0';

	/* Build candidate tree_top from integer-label Newick */
	if(tree_top != NULL) { dismantle_tree(tree_top); tree_top = NULL; }
	temp_top = NULL;
	{ int _to = 0; tree_build(1, (char *)candidate_nwk, NULL, FALSE, -1, &_to); }
	tree_top = temp_top;
	temp_top = NULL;
	if(tree_top == NULL) { free(best_tree); free(temptree); return FALSE; }

	/* Sanity: correct taxon count and no ctrl+C */
	if(check_taxa(tree_top) != number_of_taxa || user_break)
		{
		dismantle_tree(tree_top); tree_top = NULL;
		free(best_tree); free(temptree);
		return FALSE;
		}

	/* Deduplication: skip already-visited topologies this replicate */
	topo_h = tree_topo_hash(tree_top);
	if(visited_set != NULL && vs_contains(visited_set, topo_h))
		{
		skip_streak++;
		if(g_landscape_file[0])
			lm_record(landscape_map ? landscape_map : g_landscape_map,
			          topo_h, 0.0f, NULL);
		dismantle_tree(tree_top); tree_top = NULL;
		free(best_tree); free(temptree);
		return FALSE;
		}
	if(visited_set != NULL) vs_insert(visited_set, topo_h);
	skip_streak = 0;
	tried_regrafts++;
	NUMSWAPS++;

	/* Save sourcetree_scores[] — restore on non-improvement */
	for(i = 0; i < Total_fund_trees; i++) tmp_fund_scores[i] = sourcetree_scores[i];

	/* Named-tree string (for retained_supers and some scoring functions) */
	print_named_tree(tree_top, best_tree);
	strcat(best_tree, ";");
	while(unroottree(best_tree));

	/* Score by criterion */
	if(criterion == 0)
		tmpscore = compare_trees(FALSE);
	else if(criterion == 2)
		tmpscore = compare_trees_sfit(TRUE);
	else if(criterion == 3)
		tmpscore = compare_trees_qfit(TRUE);
	else if(criterion == 5)
		{
		print_tree(tree_top, temptree);
		strcat(temptree, ";");
		while(unroottree(temptree));
		tmpscore = get_recon_score(temptree, numspectries, numgenetries);
		}
	else if(criterion == 6)
		tmpscore = compare_trees_rf(TRUE);
	else if(criterion == 7)
		tmpscore = compare_trees_ml(TRUE);

	/* Landscape recording — store display-scale score (lnL for ML, raw for others)
	 * so the TSV column is on the same scale as all other reported values. */
	if(g_landscape_file[0])
		lm_record(landscape_map ? landscape_map : g_landscape_map,
		          topo_h, ml_display_score(tmpscore), best_tree);

	/* --- Bookkeeping --- */
	/* Always initialise sprscore at the start of a new rep */
	if(sprscore == -1)
		sprscore = tmpscore;

	if(BESTSCORE == -1) BESTSCORE = tmpscore;

	if(scores_retained_supers[0] == -1)
		{
		/* Very first tree ever across all reps: bootstrap retained_supers */
		retained_supers[0] = realloc(retained_supers[0],
		                            (strlen(best_tree) + 10) * sizeof(char));
		strcpy(retained_supers[0], best_tree);
		scores_retained_supers[0] = tmpscore;
		better_score = TRUE;
		}
	else
		{
		if(tmpscore < scores_retained_supers[0] || tmpscore < sprscore)
			{
			if(tmpscore < BESTSCORE) BESTSCORE = tmpscore;
			if(tmpscore < sprscore)
				{ better_score = TRUE; sprscore = tmpscore; }
			if(tmpscore < scores_retained_supers[0])
				{
				better_score = TRUE;
				retained_supers[0] = realloc(retained_supers[0],
				                            (strlen(best_tree) + 10) * sizeof(char));
				strcpy(retained_supers[0], best_tree);
				scores_retained_supers[0] = tmpscore;
#ifdef _OPENMP
				if(omp_get_num_threads() > 1)
					{
					time_t _now = time(NULL);
					int    new_global_best = 0;
					#pragma omp critical (par_progress)
						{
						if(par_progress_best < 0.0f || tmpscore < par_progress_best)
							{ par_progress_best = tmpscore; new_global_best = 1; }
						}
					if(new_global_best &&
					   (hs_progress_interval == 0 ||
					    difftime(_now, par_last_progress_time) >= hs_progress_interval))
						{
						double hs_el = difftime(_now, par_search_start);
						int hs_em = (int)(hs_el / 60), hs_es = (int)hs_el % 60;
						int hs_imp = (par_last_print_score < 0.0f ||
						              par_progress_best < par_last_print_score);
						if(hsprint)
							{
							printf2("  [%2d:%02d]  best so far = %s%.6f\tspr/tbr\trep %d\n",
							        hs_em, hs_es,
							        hs_imp ? "* " : "  ",
							        ml_display_score(par_progress_best), hs_par_rep);
							par_last_print_score   = par_progress_best;
							par_last_progress_time = _now;
							fflush(stdout);
							}
						}
					}
#endif
				j = 1;
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
			/* Not strictly better — check for equal-score distinct topology */
			if(tmpscore == scores_retained_supers[0])
				{
				j = 0;
				different = TRUE;
				if(!check_if_diff_tree(best_tree)) different = FALSE;
				if(different)
					{
					while(scores_retained_supers[j] != -1)
						{
						j++;
						if(j + 1 == number_retained_supers) reallocate_retained_supers();
						}
					retained_supers[j] = realloc(retained_supers[j],
					                            (strlen(best_tree) + 10) * sizeof(char));
					strcpy(retained_supers[j], best_tree);
					scores_retained_supers[j] = tmpscore;
					}
				}
			/* Restore sourcetree_scores to last-known-good state */
			for(i = 0; i < Total_fund_trees; i++)
				sourcetree_scores[i] = tmp_fund_scores[i];
			/* Dismantle non-improving tree_top */
			dismantle_tree(tree_top); tree_top = NULL;
			}
		}

	free(best_tree);
	free(temptree);
	return better_score;
	}


/* =======================================================================
 * spr_new3: string-native SPR search driver.
 *
 * Uses spr_tree.c (spr_parse / spr_bisect / spr_graft) for all tree
 * manipulation — zero pointer surgery, no clann tree invariant issues.
 * Candidate topologies are produced as Newick strings and scored via
 * evaluate_candidate(), which builds a fresh tree_top for each candidate.
 *
 * Algorithm per call:
 *   1. Serialise master → integer-label Newick.
 *   2. Parse → spr_node tree (once per call; dismantle master after parse).
 *   3. Enumerate all edges (prune positions).
 *   4. For each prune edge x:
 *      a. Bisect → sub_nwk + rem_nwk.
 *      b. Parse rem_nwk → rem_tree.
 *      c. Enumerate graft edges.
 *      d. For each graft edge g: spr_graft() → candidate_nwk.
 *      e. evaluate_candidate() — stops the inner loop on first improvement.
 *   5. Return TRUE on first improvement (do_search outer loop will re-call).
 *
 * Thread safety: all state is stack-local or threadprivate.
 * ctrl+C: user_break checked at every outer- and inner-loop iteration.
 * ======================================================================= */
static int spr_new3(struct taxon *master, int maxswaps, int numspectries, int numgenetries)
	{
	int   better_score = FALSE;
	int   x, g, n_prune = 0, n_graft = 0;
	int   max_nodes;
	float best_score;   /* best-improvement: best score seen in survey */

	char *master_nwk    = NULL;
	char *sub_nwk       = NULL;
	char *rem_nwk       = NULL;
	char *candidate_nwk = NULL;
	char *best_nwk      = NULL;   /* best-improvement: Newick of best candidate */
	float *tmp_fund_scores = NULL;

	struct spr_node  *spr_tree  = NULL;
	struct spr_node  *rem_tree  = NULL;
	struct spr_node **prune_edges = NULL;
	struct spr_node **graft_edges = NULL;

	/* Allocate Newick string buffers */
	master_nwk    = malloc(TREE_LENGTH * sizeof(char));
	sub_nwk       = malloc(TREE_LENGTH * sizeof(char));
	rem_nwk       = malloc(TREE_LENGTH * sizeof(char));
	candidate_nwk = malloc(TREE_LENGTH * sizeof(char));
	best_nwk      = malloc(TREE_LENGTH * sizeof(char));
	tmp_fund_scores = malloc(Total_fund_trees * sizeof(float));
	if(!master_nwk || !sub_nwk || !rem_nwk || !candidate_nwk || !best_nwk || !tmp_fund_scores)
		goto spr3_cleanup;
	best_nwk[0] = '\0';

	/* Edge arrays: an unrooted tree with N taxa has at most 2N-2 nodes */
	max_nodes   = 2 * number_of_taxa + 8;
	prune_edges = malloc(max_nodes * sizeof(struct spr_node *));
	graft_edges = malloc(max_nodes * sizeof(struct spr_node *));
	if(!prune_edges || !graft_edges) goto spr3_cleanup;

	/* --- Step 1: Serialise master to integer-label Newick --- */
	master_nwk[0] = '\0';
	print_tree(master, master_nwk);
	strcat(master_nwk, ";");

	/* --- Step 2: Parse into spr_node tree then dismantle master ---
	 * Once parsed, we have everything we need in spr_tree; master (which
	 * is tree_top) can be freed so tree_top is clean for evaluate_candidate. */
	spr_tree = spr_parse(master_nwk);
	if(spr_tree == NULL) goto spr3_cleanup;

	dismantle_tree(master); master = NULL; tree_top = NULL;

	/* --- Step 3: Enumerate prune edges --- */
	spr_get_edges(spr_tree, prune_edges, &n_prune);

	/* --- Steps 4–5: For each prune edge, bisect and try all grafts --- */
	for(x = 0; x < n_prune; x++)
		{
		if(better_score) break;   /* a commit was made — stop */
#ifdef _OPENMP
		if(hs_thread_report_interval > 0 && omp_get_num_threads() > 1)
			{
			time_t _tr_now = time(NULL);
			if(difftime(_tr_now, thread_report_last) >= hs_thread_report_interval)
				{
				thread_report_last = _tr_now;
				double tr_el = difftime(_tr_now, par_search_start);
				int    tr_em = (int)(tr_el / 60), tr_es = (int)tr_el % 60;
				float  _pbest_s3;
				#pragma omp atomic read
				_pbest_s3 = par_progress_best;
				if(sprscore < 0.0f)
					printf2("  [%2d:%02d]  thread %2d  rep %d  score=(searching)  tried=%d  skips=%d  [spr_new3 branch %d/%d]\n",
					        tr_em, tr_es, omp_get_thread_num(), hs_par_rep,
					        tried_regrafts, skip_streak, x, n_prune);
				else
					printf2("  [%2d:%02d]  thread %2d  rep %d  %s=  %s%.6f  tried=%d  skips=%d  [spr_new3 branch %d/%d]\n",
					        tr_em, tr_es, omp_get_thread_num(), hs_par_rep,
					        ml_score_label(),
					        ml_display_score(sprscore) < ml_display_score(_pbest_s3) ? "* " : "  ",
					        ml_display_score(sprscore), tried_regrafts, skip_streak, x, n_prune);
				fflush(stdout);
				}
			}
#endif
		if(user_break || rep_abandon) break;

		/* Bisect at prune_edges[x] */
		if(spr_bisect(spr_tree, prune_edges[x], sub_nwk, rem_nwk,
		              TREE_LENGTH) != 0)
			continue;

		/* Parse remaining tree */
		rem_tree = spr_parse(rem_nwk);
		if(rem_tree == NULL) { rep_abandon = 1; break; }

		/* Enumerate graft positions in remaining tree */
		spr_get_edges(rem_tree, graft_edges, &n_graft);

		/* Per-cut best-improvement: reset tracker for this prune edge */
		best_score  = 1e30f;
		best_nwk[0] = '\0';

		/* Try each graft position */
		for(g = 0; g < n_graft; g++)
			{
			if(hs_strategy == 0 && better_score) break;
			if(user_break || rep_abandon) break;
			if(tried_regrafts >= maxswaps) break;
			if(hs_maxskips > 0 && skip_streak >= hs_maxskips) break;

			spr_graft(sub_nwk, rem_tree, graft_edges[g],
			          candidate_nwk, TREE_LENGTH);

			if(hs_strategy == 0)
				{
				/* First-improvement: commit immediately on improvement */
				if(evaluate_candidate(candidate_nwk, tmp_fund_scores,
				                      numspectries, numgenetries))
					better_score = TRUE;
				}
			else
				{
				/* Best-improvement: probe only, track best */
				float s = probe_candidate(candidate_nwk, tmp_fund_scores,
				                          numspectries, numgenetries);
				if(s >= 0.0f && s < best_score)
					{ best_score = s; strcpy(best_nwk, candidate_nwk); }
				}
			}

		spr_free(rem_tree); rem_tree = NULL;

		/* Per-cut best-improvement: commit the best graft found for this cut */
		if(hs_strategy == 1 && best_nwk[0] != '\0')
			{
			float threshold = (sprscore >= 0.0f) ? sprscore : 1e30f;
			if(best_score < threshold)
				better_score = evaluate_candidate(best_nwk, tmp_fund_scores,
				                                  numspectries, numgenetries);
			}

		if(tried_regrafts >= maxswaps) break;
		if(hs_strategy == 0 && hs_maxskips > 0 && skip_streak >= hs_maxskips) break;
		}

spr3_cleanup:
	spr_free(spr_tree);
	spr_free(rem_tree);
	free(master_nwk);
	free(sub_nwk);
	free(rem_nwk);
	free(candidate_nwk);
	free(best_nwk);
	free(tmp_fund_scores);
	free(prune_edges);
	free(graft_edges);

	/* On FALSE return tree_top is already NULL (evaluate_candidate dismantled it).
	 * On TRUE  return tree_top points to the new best tree. */
	if(!better_score) tree_top = NULL;
	return better_score;
	}


/* =======================================================================
 * tbr_new2: string-native TBR search driver.
 *
 * Same outer structure as spr_new3() but adds an inner rerooting loop:
 * for each pruned subtree with more than 4 nodes, every internal edge
 * of the subtree is tried as a rerooting point before attempting the graft.
 * This enumerates the full TBR neighbourhood.
 *
 * For pruned subtrees with <= 4 nodes (0 or 1 internal edges), a single
 * SPR-style graft is performed (equivalent to TBR for small subtrees).
 * ======================================================================= */
static int tbr_new2(struct taxon *master, int maxswaps, int numspectries, int numgenetries)
	{
	int   better_score = FALSE;
	int   x, g, q, n_prune = 0, n_graft = 0, n_internals = 0;
	int   max_nodes;
	float best_score;   /* best-improvement: best score seen in survey */

	char *master_nwk    = NULL;
	char *sub_nwk       = NULL;
	char *rem_nwk       = NULL;
	char *rerooted_nwk  = NULL;
	char *candidate_nwk = NULL;
	char *best_nwk      = NULL;   /* best-improvement: Newick of best candidate */
	float *tmp_fund_scores = NULL;

	struct spr_node  *spr_tree    = NULL;
	struct spr_node  *sub_tree    = NULL;
	struct spr_node  *rem_tree    = NULL;
	struct spr_node **prune_edges = NULL;
	struct spr_node **graft_edges = NULL;
	struct spr_node **sub_internals = NULL;

	/* Allocate Newick string buffers */
	master_nwk    = malloc(TREE_LENGTH * sizeof(char));
	sub_nwk       = malloc(TREE_LENGTH * sizeof(char));
	rem_nwk       = malloc(TREE_LENGTH * sizeof(char));
	rerooted_nwk  = malloc(TREE_LENGTH * sizeof(char));
	candidate_nwk = malloc(TREE_LENGTH * sizeof(char));
	best_nwk      = malloc(TREE_LENGTH * sizeof(char));
	tmp_fund_scores = malloc(Total_fund_trees * sizeof(float));
	if(!master_nwk || !sub_nwk || !rem_nwk || !rerooted_nwk ||
	   !candidate_nwk || !best_nwk || !tmp_fund_scores)
		goto tbr2_cleanup;
	best_nwk[0] = '\0';

	/* Edge / internal-node arrays */
	max_nodes     = 2 * number_of_taxa + 8;
	prune_edges   = malloc(max_nodes * sizeof(struct spr_node *));
	graft_edges   = malloc(max_nodes * sizeof(struct spr_node *));
	sub_internals = malloc(max_nodes * sizeof(struct spr_node *));
	if(!prune_edges || !graft_edges || !sub_internals) goto tbr2_cleanup;

	/* Serialise master → integer-label Newick */
	master_nwk[0] = '\0';
	print_tree(master, master_nwk);
	strcat(master_nwk, ";");

	/* Parse once, then dismantle master */
	spr_tree = spr_parse(master_nwk);
	if(spr_tree == NULL) goto tbr2_cleanup;

	dismantle_tree(master); master = NULL; tree_top = NULL;

	spr_get_edges(spr_tree, prune_edges, &n_prune);

	for(x = 0; x < n_prune; x++)
		{
		if(better_score) break;   /* a commit was made — stop */
#ifdef _OPENMP
		if(hs_thread_report_interval > 0 && omp_get_num_threads() > 1)
			{
			time_t _tr_now = time(NULL);
			if(difftime(_tr_now, thread_report_last) >= hs_thread_report_interval)
				{
				thread_report_last = _tr_now;
				double tr_el = difftime(_tr_now, par_search_start);
				int    tr_em = (int)(tr_el / 60), tr_es = (int)tr_el % 60;
				float  _pbest_t2;
				#pragma omp atomic read
				_pbest_t2 = par_progress_best;
				if(sprscore < 0.0f)
					printf2("  [%2d:%02d]  thread %2d  rep %d  score=(searching)  tried=%d  skips=%d  [tbr_new2 branch %d/%d]\n",
					        tr_em, tr_es, omp_get_thread_num(), hs_par_rep,
					        tried_regrafts, skip_streak, x, n_prune);
				else
					printf2("  [%2d:%02d]  thread %2d  rep %d  %s=  %s%.6f  tried=%d  skips=%d  [tbr_new2 branch %d/%d]\n",
					        tr_em, tr_es, omp_get_thread_num(), hs_par_rep,
					        ml_score_label(),
					        ml_display_score(sprscore) < ml_display_score(_pbest_t2) ? "* " : "  ",
					        ml_display_score(sprscore), tried_regrafts, skip_streak, x, n_prune);
				fflush(stdout);
				}
			}
#endif
		if(user_break || rep_abandon) break;

		if(spr_bisect(spr_tree, prune_edges[x], sub_nwk, rem_nwk,
		              TREE_LENGTH) != 0)
			continue;

		rem_tree = spr_parse(rem_nwk);
		if(rem_tree == NULL) { rep_abandon = 1; break; }

		spr_get_edges(rem_tree, graft_edges, &n_graft);

		/* Parse subtree to count internal nodes */
		sub_tree = spr_parse(sub_nwk);
		if(sub_tree == NULL) { spr_free(rem_tree); rem_tree = NULL; rep_abandon = 1; break; }

		if(sub_tree->first_child == NULL)
			{
			/* Single-leaf subtree: no rerooting needed — treat as SPR */
			n_internals = 0;
			}
		else
			{
			spr_get_internals(sub_tree, sub_internals, &n_internals);
			}

		/* Per-cut best-improvement: reset tracker for this prune edge
		 * (includes all rerootings of the subtree) */
		best_score  = 1e30f;
		best_nwk[0] = '\0';

		if(n_internals <= 1)
			{
			/* Small subtree (<=4 nodes): single SPR-style pass */
			for(g = 0; g < n_graft; g++)
				{
				if(hs_strategy == 0 && better_score) break;
				if(user_break || rep_abandon) break;
				if(tried_regrafts >= maxswaps) break;
				if(hs_maxskips > 0 && skip_streak >= hs_maxskips) break;

				spr_graft(sub_nwk, rem_tree, graft_edges[g],
				          candidate_nwk, TREE_LENGTH);

				if(hs_strategy == 0)
					{
					if(evaluate_candidate(candidate_nwk, tmp_fund_scores,
					                      numspectries, numgenetries))
						better_score = TRUE;
					}
				else
					{
					float s = probe_candidate(candidate_nwk, tmp_fund_scores,
					                          numspectries, numgenetries);
					if(s >= 0.0f && s < best_score)
						{ best_score = s; strcpy(best_nwk, candidate_nwk); }
					}
				}
			}
		else
			{
			/* Larger subtree: for each rerooting of the subtree, try all
			 * graft positions in the remaining tree. */
			for(q = 0; q < n_internals; q++)
				{
				if(hs_strategy == 0 && better_score) break;
				if(user_break || rep_abandon) break;

				/* Produce rerooted subtree Newick */
				spr_write_rerooted(sub_tree, sub_internals[q],
				                   rerooted_nwk, TREE_LENGTH);

				for(g = 0; g < n_graft; g++)
					{
					if(hs_strategy == 0 && better_score) break;
					if(user_break || rep_abandon) break;
					if(tried_regrafts >= maxswaps) break;
					if(hs_maxskips > 0 && skip_streak >= hs_maxskips) break;

					spr_graft(rerooted_nwk, rem_tree, graft_edges[g],
					          candidate_nwk, TREE_LENGTH);

					if(hs_strategy == 0)
						{
						if(evaluate_candidate(candidate_nwk, tmp_fund_scores,
						                      numspectries, numgenetries))
							better_score = TRUE;
						}
					else
						{
						float s = probe_candidate(candidate_nwk, tmp_fund_scores,
						                          numspectries, numgenetries);
						if(s >= 0.0f && s < best_score)
							{ best_score = s; strcpy(best_nwk, candidate_nwk); }
						}
					}

				if(tried_regrafts >= maxswaps) break;
				if(hs_strategy == 0 && hs_maxskips > 0 && skip_streak >= hs_maxskips) break;
				}
			}

		/* Per-cut best-improvement: commit the best graft found for this cut */
		if(hs_strategy == 1 && best_nwk[0] != '\0')
			{
			float threshold = (sprscore >= 0.0f) ? sprscore : 1e30f;
			if(best_score < threshold)
				better_score = evaluate_candidate(best_nwk, tmp_fund_scores,
				                                  numspectries, numgenetries);
			}

		spr_free(sub_tree); sub_tree = NULL;
		spr_free(rem_tree); rem_tree = NULL;

		if(tried_regrafts >= maxswaps) break;
		if(hs_strategy == 0 && hs_maxskips > 0 && skip_streak >= hs_maxskips) break;
		}

tbr2_cleanup:
	spr_free(spr_tree);
	spr_free(sub_tree);
	spr_free(rem_tree);
	free(master_nwk);
	free(sub_nwk);
	free(rem_nwk);
	free(rerooted_nwk);
	free(candidate_nwk);
	free(best_nwk);
	free(tmp_fund_scores);
	free(prune_edges);
	free(graft_edges);
	free(sub_internals);

	if(!better_score) tree_top = NULL;
	return better_score;
	}


int spr_new(struct taxon * master, int maxswaps, int numspectries, int numgenetries)
	{
	struct taxon * start = NULL, * newbie = NULL, *newbie_backup = NULL, * latest = NULL, * tmp = NULL, *copy = NULL, *temper = NULL, *temper1 = NULL, *master_copy = NULL, *position = NULL;
	int better_score = FALSE, lastinline = FALSE, numofsiblings = 0, donenextlevel = FALSE, i=0, j=0, x=0, y=0, q=0, r=0;
	char *debugtree = NULL;
	
	debugtree = malloc(TREE_LENGTH*sizeof(char));
	debugtree[0] = '\0';
    

	start = master;



	/* As opposed to travelling down the tree in a recursive manner we shall implement a methid which replaces the tree used after each regrafting */
	
	/* 1) First make a copy of the master tree passed to the application, called master_copy */
	
	temp_top = NULL;
	duplicate_tree(master, NULL); /* make a copy of the tree */
	master_copy = temp_top;
	temp_top = NULL;
	/* end 1) */
	/* 2) count the number of nodes that it is possible to break this tree (internal and external branches) */
	
	y = number_tree(master, 0); /* count the parts of newbie */

	
	/* end 2) */
	for(x=0; x<y; x++) /* for each branch of the master tree, create an instance where it is bisected at this point */
		{
#ifdef _OPENMP
		/* Per-thread diagnostic: report from spr_new outer loop so silent threads are
		 * visible even when stuck in dismantle/duplicate between regraft() calls. */
		if(hs_thread_report_interval > 0 && omp_get_num_threads() > 1)
			{
			time_t _tr_now = time(NULL);
			if(difftime(_tr_now, thread_report_last) >= hs_thread_report_interval)
				{
				thread_report_last = _tr_now;
				double tr_el = difftime(_tr_now, par_search_start);
				int    tr_em = (int)(tr_el / 60), tr_es = (int)tr_el % 60;
				float  _pbest2;
				#pragma omp atomic read
				_pbest2 = par_progress_best;
				if(sprscore < 0.0f)
					printf2("  [%2d:%02d]  thread %2d  rep %d  score=(searching)  tried=%d  skips=%d  [spr_new branch %d/%d]\n",
					        tr_em, tr_es, omp_get_thread_num(), hs_par_rep,
					        tried_regrafts, skip_streak, x, y);
				else
					printf2("  [%2d:%02d]  thread %2d  rep %d  %s=  %s%.6f  tried=%d  skips=%d  [spr_new branch %d/%d]\n",
					        tr_em, tr_es, omp_get_thread_num(), hs_par_rep,
					        ml_score_label(),
					        ml_display_score(sprscore) < ml_display_score(_pbest2) ? "* " : "  ",
					        ml_display_score(sprscore), tried_regrafts, skip_streak, x, y);
				fflush(stdout);
				}
			}
#endif
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
		{ int _g = 0, _gmax = number_of_taxa * 4 + 16;
		  while(tmp->prev_sibling != NULL && _g++ < _gmax) tmp = tmp->prev_sibling;
		  if(_g >= _gmax) rep_abandon = 1; /* sibling chain is cyclic — corrupt tree, abandon rep */ }
		/* count number of sibligs */
		numofsiblings = 0;
		{ int _g2 = 0, _gmax2 = number_of_taxa * 4 + 16;
		  while(tmp != NULL && _g2++ < _gmax2) { numofsiblings++; tmp = tmp->next_sibling; }
		  if(_g2 >= _gmax2) rep_abandon = 1; }

		if(position != master && !rep_abandon) /* if we are at the top of the tree, don't do anything */
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
			/* Only update parent->daughter when position was the first child.
			 * When position has a prev_sibling it is not the first child and
			 * parent->daughter already points elsewhere.  The redundant
			 * next_sibling->parent assignment is removed: tree_build() already
			 * set those pointers, and it crashes when next_sibling is NULL. */
			if(position->parent != NULL && position->prev_sibling == NULL)
				(position->parent)->daughter = position->next_sibling;
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
			
			tmp = latest;
			{ int _g = 0, _gmax = number_of_taxa * 4 + 16;
			  while(tmp->prev_sibling != NULL && _g++ < _gmax) tmp = tmp->prev_sibling; /* rewinding */
			  if(_g >= _gmax) rep_abandon = 1; }
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
				
				
				better_score = regraft(latest, newbie, NULL, 1, maxswaps, numspectries, numgenetries);
				
				}
			if(method == 3) /* tbr — handled by tbr_new(), called from do_search() */
				{
				/* This branch is dead code when method==3 because do_search()
				 * dispatches to tbr_new() instead of spr_new() for TBR. */
				if(newbie != NULL) { dismantle_tree(newbie); newbie = NULL; }
				}
			if(!better_score && newbie != NULL) dismantle_tree(newbie);
			newbie = NULL;
            }
		if(better_score == TRUE || (position == branchpointer) || donenextlevel || tried_regrafts >= maxswaps || user_break || rep_abandon
		   || (hs_maxskips > 0 && skip_streak >= hs_maxskips))
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
	*/
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
									
									dismantle_tree(newbie);
									newbie = NULL;
									temp_top = NULL;
									{ int _to = 0; tree_build(1, debugtree, newbie, FALSE, -1, &_to); }
									newbie = temp_top;
									temp_top = NULL;
									
									debugtree[0] = '\0';
									
									i = number_tree(newbie, 0); /* count the parts of newbie */
									temp_top = NULL;
									duplicate_tree(newbie, NULL); /* make a copy of the subtree */
									copy = temp_top;
									temp_top = NULL;
									
									
								
									for(j=0; j<i; j++) /* reroot newbie at each of the "i" positions on the subtree */
										{
										temper1 = get_branch(newbie, j);
										
										
										temp_top = newbie;
										reroot_tree(temper1);
										newbie = temp_top;

										
										temper = make_taxon();
										temper->daughter = newbie; 
										newbie->parent = temper;
										newbie = temper;
										temper= NULL;
										newbie->spr = TRUE;


										better_score = regraft(tmp, newbie, NULL, 1, maxswaps, numspectries, numgenetries);
										

										if(better_score || user_break) j = i;
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
									}
								else
									{
									temper = make_taxon();
									temper->daughter = newbie; 
									newbie->parent = temper;
									newbie = temper;
									temper= NULL;
									newbie->spr = TRUE;
									
									
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
            
        
	{ int _g = 0, _gmax = number_of_taxa * 4 + 16;
	  while(position->prev_sibling != NULL && !user_break && _g++ < _gmax) position = position->prev_sibling; /* rewinding */
	  if(_g >= _gmax) rep_abandon = 1; }
	start = position;
	
	while(position != NULL && !better_score && tried_regrafts < maxswaps && !user_break && !rep_abandon
	      && (hs_maxskips == 0 || skip_streak < hs_maxskips))
		{
		if(steps < number_of_steps && position->parent != NULL && position->parent != last && !user_break) better_score = regraft(position->parent, newbie, position, steps+1, maxswaps, numspectries, numgenetries);  /* go up the tree */
		position = position->next_sibling;
		}

	position = start;
	if(!better_score && !user_break && !rep_abandon)
		{
		while(position != NULL && !better_score && tried_regrafts < maxswaps && !user_break && !rep_abandon
		      && (hs_maxskips == 0 || skip_streak < hs_maxskips))
			{
			if(steps < number_of_steps && !user_break)
				if(position->daughter != NULL && position->daughter != last && !user_break) better_score = regraft(position->daughter, newbie, position, steps+1, maxswaps, numspectries, numgenetries);

			position = position->next_sibling;
			}
		position = start;

		while(position != NULL && !better_score && tried_regrafts < maxswaps && !user_break && !rep_abandon
		      && (hs_maxskips == 0 || skip_streak < hs_maxskips))
			{
			if(!better_score && steps <= number_of_steps && position != tree_top && !user_break
			   && !is_ancestor_of(newbie, position))  /* skip if position is inside newbie's subtree — would create a cycle */
				{
				/* regraft newbie here -- skip position==tree_top: inserting there creates a
				 * degree-2 internal node (invalid in a binary unrooted tree). The 3 edges
				 * adjacent to tree_top are covered by grafting at its daughters instead. */
                                
                                
			/*	if((newbie->daughter)->name != -1)printf2("newbie daughter = %d\n", newbie->daughter->name);
				else printf2("newbie daughter = pointer\n");
				if(position->name != -1) printf2("new sibling daughter = %d\n", position->name);
				else printf2("pointer to %d\n", position->daughter->name);
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
#ifdef _OPENMP
                                /* Per-thread diagnostic report: each thread independently prints
                                 * its own current state every hs_thread_report_interval seconds.
                                 * thread_report_last is threadprivate so each thread has its own
                                 * timer and fires independently of all other threads. */
                                if(hs_thread_report_interval > 0 && omp_get_num_threads() > 1 &&
                                   difftime(interval2, thread_report_last) >= hs_thread_report_interval)
                                    {
                                    thread_report_last = interval2;
                                    double tr_el = difftime(interval2, par_search_start);
                                    int    tr_em = (int)(tr_el / 60), tr_es = (int)tr_el % 60;
                                    float _pbest3;
                                    #pragma omp atomic read
                                    _pbest3 = par_progress_best;
                                    if(sprscore < 0.0f)
                                        printf2("  [%2d:%02d]  thread %2d  rep %d  score=(searching)  tried=%d  skips=%d\n",
                                                tr_em, tr_es, omp_get_thread_num(), hs_par_rep,
                                                tried_regrafts, skip_streak);
                                    else
                                        printf2("  [%2d:%02d]  thread %2d  rep %d  %s=  %s%.6f  tried=%d  skips=%d\n",
                                                tr_em, tr_es, omp_get_thread_num(), hs_par_rep,
                                                ml_score_label(), ml_display_score(sprscore) < ml_display_score(_pbest3) ? "* " : "  ",
                                                ml_display_score(sprscore), tried_regrafts, skip_streak);
                                    fflush(stdout);
                                    }
#endif
                                if(hs_progress_interval > 0 && difftime(interval2, interval1) >= hs_progress_interval)
                                    {
                                    if(hs_do_print)
                                        {
                                        /* Sequential: per-replicate status line */
                                        double hs_el = difftime(interval2, rep_start_time);
                                        int    hs_em = (int)(hs_el / 60);
                                        int    hs_es = (int)hs_el % 60;
                                        if(sprscore < 0.0f)
                                            printf2("  [%2d:%02d]  tried=%6d  %s=(searching...)\n",
                                                    hs_em, hs_es, tried_regrafts, ml_score_label());
                                        else
                                            {
                                            int hs_imp = (last_status_score < 0.0f || sprscore < last_status_score);
                                            printf2("  [%2d:%02d]  tried=%6d  %s=%s%.6f\n",
                                                    hs_em, hs_es, tried_regrafts,
                                                    ml_score_label(), hs_imp ? "* " : "  ", ml_display_score(sprscore));
                                            if(hs_imp) last_status_score = sprscore;
                                            }
                                        fflush(stdout);
                                        }
#ifdef _OPENMP
                                    else if(omp_get_num_threads() > 1)
                                        {
                                        /* Parallel heartbeat: update shared best; the first thread to
                                         * exceed the interval window claims the print slot by updating
                                         * par_last_progress_time inside the critical section, so only one
                                         * line is emitted per interval regardless of which thread fires.
                                         * No rep attribution — the global best may have been found by any rep. */
                                        int do_heartbeat = 0;
                                        float heartbeat_score = -1.0f;
                                        int   heartbeat_imp   = 0;
                                        #pragma omp critical (par_progress)
                                            {
                                            if(sprscore >= 0.0f &&
                                               (par_progress_best < 0.0f || sprscore < par_progress_best))
                                                par_progress_best = sprscore;
                                            if(par_progress_best >= 0.0f &&
                                               difftime(interval2, par_last_progress_time) >= hs_progress_interval)
                                                {
                                                do_heartbeat   = 1;
                                                heartbeat_score = par_progress_best;
                                                heartbeat_imp   = (par_last_print_score < 0.0f ||
                                                                   par_progress_best < par_last_print_score);
                                                par_last_print_score   = par_progress_best;
                                                par_last_progress_time = interval2;
                                                }
                                            }
                                        if(do_heartbeat)
                                            {
                                            double hs_el = difftime(interval2, par_search_start);
                                            int    hs_em = (int)(hs_el / 60), hs_es = (int)hs_el % 60;
                                            if(hsprint)
                                                {
                                                printf2("  [%2d:%02d]  best so far = %s%.6f\tspr/tbr\n",
                                                        hs_em, hs_es,
                                                        heartbeat_imp ? "* " : "  ", ml_display_score(heartbeat_score));
                                                fflush(stdout);
                                                }
                                            }
                                        }
#endif
                                    interval1 = time(NULL);
                                    }
                        
                                
                                if(check_taxa(tree_top) == number_of_taxa && !user_break)
                                    {
                                    /* Skip trees already scored in this replicate */
                                    uint64_t topo_h = tree_topo_hash(tree_top);
                                    if(visited_set == NULL || !vs_contains(visited_set, topo_h))
                                        {
                                        if(visited_set != NULL) vs_insert(visited_set, topo_h);
                                        skip_streak = 0;  /* new topology found — reset convergence counter */
                                        tried_regrafts++;
									NUMSWAPS++;
                                    
                                    strcpy(best_tree, "");
                                    print_named_tree(tree_top, best_tree);
                                    strcat(best_tree, ";");
                                    while(unroottree(best_tree));  /* ensure unrooted even if graft was at root */
                                    if(criterion == 0) tmpscore = compare_trees(TRUE);  /* calculate the distance from the super tree to all the fundamental trees */
                                    if(criterion == 2) tmpscore = compare_trees_sfit(TRUE);
                                    if(criterion == 3) tmpscore = compare_trees_qfit(TRUE);
														
									if(criterion==5)
										{
										strcpy(temptree, "");
										print_tree(tree_top, temptree);
										strcat(temptree, ";");
										while(unroottree(temptree));
										tmpscore = get_recon_score(temptree, numspectries, numgenetries);
										}
									if(criterion == 6)
										tmpscore = compare_trees_rf(TRUE);
									if(criterion == 7)
										tmpscore = compare_trees_ml(TRUE);

                                    /* Record this topology in the landscape map (first visit: store score+newick) */
                                    if(g_landscape_file[0])
                                        lm_record(landscape_map ? landscape_map : g_landscape_map,
                                                  topo_h, ml_display_score(tmpscore), best_tree);

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
                                                while(unroottree(best_tree));  /* ensure unrooted */
												retained_supers[0] = realloc(retained_supers[0], (strlen(best_tree)+10)*sizeof(char));
                                                strcpy(retained_supers[0], best_tree);
                                                scores_retained_supers[0] = tmpscore;
#ifdef _OPENMP
                                                /* SPR/TBR parallel progress: only when THIS rep set the new
                                                 * global best — the rep number in the message then matches
                                                 * the rep that actually achieved the displayed score. */
                                                if(omp_get_num_threads() > 1)
                                                    {
                                                    time_t _now = time(NULL);
                                                    int new_global_best = 0;
                                                    #pragma omp critical (par_progress)
                                                        {
                                                        if(par_progress_best < 0.0f || tmpscore < par_progress_best)
                                                            {
                                                            par_progress_best = tmpscore;
                                                            new_global_best = 1;
                                                            }
                                                        }
                                                    if(new_global_best &&
                                                       (hs_progress_interval == 0 ||
                                                        difftime(_now, par_last_progress_time) >= hs_progress_interval))
                                                        {
                                                        double hs_el = difftime(_now, par_search_start);
                                                        int    hs_em = (int)(hs_el / 60), hs_es = (int)hs_el % 60;
                                                        int    hs_imp = (par_last_print_score < 0.0f ||
                                                                         par_progress_best < par_last_print_score);
                                                        if(hsprint)
                                                            {
                                                            printf2("  [%2d:%02d]  best so far = %s%.6f\tspr/tbr\trep %d\n",
                                                                    hs_em, hs_es,
                                                                    hs_imp ? "* " : "  ", ml_display_score(par_progress_best), hs_par_rep);
                                                            par_last_print_score = par_progress_best;
                                                            par_last_progress_time = _now;
                                                            fflush(stdout);
                                                            }
                                                        }
                                                    }
#endif
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
                                                while(unroottree(best_tree));  /* ensure unrooted */
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
                                        }  /* end !already_visited */
                                    else /* topology already scored this replicate — undo graft and skip */
                                        {
                                        skip_streak++;  /* one more consecutive already-visited hit */
                                        /* Increment visit count in landscape map for this topology */
                                        if(g_landscape_file[0])
                                            lm_record(landscape_map ? landscape_map : g_landscape_map,
                                                      topo_h, 0.0f, NULL);
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
                                        }
                                    } /* if check_taxa */
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
			printf2("%f\n", position->length);
            if(position->length > largest_length) largest_length = position->length;
            }
        position = position->next_sibling;
        }
    }   
    
       
       
          

/*** The next three functions is designed to produce x  co-ordinates for each node in the tree for graphical representation ***/

void generatetrees(void)
	{
	int i, j, k, ntrees = 100, error = FALSE, gen_method = 1, random = TRUE, n = 20, data = 1, tree_rand_method = 1, super = 1, print_all_scores = FALSE, saveideal = FALSE;
	float *results = NULL;
	char *temptree = NULL, *rand_tree = NULL, filename[100], *pruned_tree = NULL, *tmp = malloc(TREE_LENGTH * sizeof(char)), superfilename[1000], c;
	if(!tmp) { printf2("Error: out of memory in generatetrees\n"); return; }
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
					printf2("Error: %s is not a valid method\n\n", parsed_command[i+1]);
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
					printf2("%d is not a valid value for ntrees\n", ntrees);
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
				printf2("%d is not a valid value for nbins\n", ntrees);
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
						printf2("Error: %s not valid sourcedata option\n", parsed_command[i+1]);
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
					printf2("Error: %s not valid savesourcetrees option\n", parsed_command[i+1]);
					error = TRUE;
					}
				}
			}
		}
		

	if(criterion == 4 || criterion == 1)
		{
		printf2("Error: This command is not currently implemented for MRP or Average consensus\n\tplease select another criterion to continue\n");
		error = TRUE;
		}
	if(data == 3 && super == 1 && trees_in_memory == 0)
		{
		printf2("\n\nError: there are no trees in memory, a supertree is required to create ideal data\n\tYou have two choices:\n\t1) choose to run a Heuristic search to create a supertree\n\t2) specify a file containing the supertree to use\n\n\n");
		error = TRUE;
		}
	if(!error)
		{
		printf2("\n\nGeneratetrees Settings:\n");
		printf2("\tGenerating");
		if(random) printf2(" %d random trees\n", ntrees);
		else printf2(" all %.0f possible trees\n", sup);
		if(random)
			{
			printf2("\tRandom supertree generator used: ");
			if(gen_method == 1) printf2("Equiprobable\n");
			else printf2("Markovian\n");
			}
		printf2("\tNumber of bins = %d\n", n);
		printf2("\tCalculating the tree scores using ");
		if(data == 1) printf2("the real source trees\n");
		if(data == 2) printf2("randomised versions of the source trees\n");
		if(data == 3) 
			{
			printf2("idealised versions of the source trees\n");
			printf2("\tIdeal data created from first supertree in ");
			if(super == 1)printf2("memory\n");
			if(super == 2)printf2("file named %s\n", superfilename);
			if(saveideal) printf2("\tIdeal data to be saved to the file 'idealtrees.ph'\n");
			}
		printf2("\tCriterion used to calculate supertree scores = ");
		if(criterion == 0) printf2("dfit\n");
		if(criterion == 2) printf2("sfit\n");
		if(criterion == 3) printf2("qfit\n");
		printf2("\tOutput file = %s\n", filename);
		if(print_all_scores)printf2("\tSaving all %d supertree scores to file 'allscores.txt'\n", ntrees);
		
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
			printf2("An error occurred while setting a signal handler\n");
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
					printf2("finished printing out trees\n");
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
						{ int _to = 0; tree_build(1, retained_supers[0], tree_top, TRUE, -1, &_to); }
					else
						{
						if((superfile = fopen(superfilename, "r")) == NULL)
							{
							printf2("Error opening supertree file %s\n", superfilename);
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
							{ int _to = 0; tree_build(1, tmp, tree_top, TRUE, -1, &_to); }
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
							if(print_pruned_tree(tree_top, 0, pruned_tree, FALSE, 0) >1)
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
		
			if(criterion == 2 || criterion == 3 || criterion == 6 || criterion == 7)
				rf_precompute_fund_biparts();
				
			interval1 = time(NULL);
			printf2("\nProgress Indicator:\n");
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
						/* printf2("="); */
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
					
					if(gen_method == 1){ int _to = 0; tree_build(1, rand_tree, tree_top, FALSE, -1, &_to); }
					if(gen_method == 2){ int _to = 0; tree_build(1, rand_tree, tree_top, TRUE, -1, &_to); }
					tree_top = temp_top;
					temp_top = NULL;
			/**** evaluate its fit to the source trees in memory *****/
					temptree[0] = '\0';
					if(criterion == 0) results[i] = compare_trees(FALSE);  /* calculate the distance from the super tree to all the fundamental trees */
					if(criterion == 2) results[i] = compare_trees_sfit(FALSE);
					if(criterion == 3) results[i] = compare_trees_qfit(FALSE);
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
	free(tmp);
	}




/* xposition1, middle_number, xposition2, yposition0, yposition1, yposition2,
 * printcolour, print_coordinates, tree_coordinates, draw_histogram:

    
/*** consensus is a function to carry out a majority-rule consensus on the results of a bootstrap analysis ****/











void sourcetree_dists(void)
	{
	int i=0, j=0, k=0, l=0, m=0, y=0, ***trees_coding = NULL, **included = NULL, *tracking = NULL, count = 0, x, total, *tree1 = NULL, *tree2 = NULL, *t1tag = NULL, *t2tag = NULL, *t1score = NULL, *t2score = NULL, same, same1, same2, same3;
	char number[100];
	char *string = malloc(TREE_LENGTH * sizeof(char));
	char RFfilename[100];
	int p1 = 0, p2 = 0, counter = 0, remaining_taxa = 0, num_falsed = 0, error = FALSE, r = 0, **shared_taxa = NULL, output_format = 0, missing_method = 2, here = TRUE, found = TRUE;
	float **results = NULL;
	FILE *RFfile = NULL;

	if(!string) { printf2("Error: out of memory in sourcetree_dists\n"); return; }
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
					printf2("Error: '%s' not valid option for output\n", parsed_command[i+1]);
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
						printf2("Error: '%s' not valid option for missing\n", parsed_command[i+1]);
						error = TRUE;
						}
					}
				}
			}
		}
	
	if((RFfile = fopen(RFfilename, "w")) == NULL)		/* check to see if the file is there */
			{								/* Open the source tree file */
			printf2("Cannot open file %s\n", RFfilename);
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

		printf2("Calculating the Robinson-Foulds distance between the trees .... \n\n");
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
			printf2("Estimating the missing data using ");
			if(missing_method == 0) printf2("ultrametric distances\n");
			else printf2("4 point condition distances\n");
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
				printf2("\n\nERROR: the overlap in the data is too sparse to calculate missing cells using ");
				if(missing_method == 0) printf2("an ultrametric estimate\n");
				if(missing_method == 1) printf2("a 4 point condition estimate\n");
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
			
		printf2("\nRobinson-Foulds distances of the source trees have been writen to file named %s\n\n", RFfilename);
		
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
	free(string);
	}


	/* Prune tree: This is a recursive function that is called for every node position of the supertree
	it then checks to see if any of the siblings on this node are not contained in the fundamental tree, these siblings are then turned off.
	This only turns off taxa, pointer siblings will have to be turned off using a separate program */


void spr_dist(void)
	{
	float real_score = 0, sprscore = 0, bestreal = 0,  amountspr, previous, totalnow, bestfake = 0, *results = NULL;
	int i=0, j=0, k=0, l=0, x=0, y=0, error = FALSE, diff, *originaldiff = NULL, numbersprs =0, nreps=100, bestnumSPR = 0, minsprs = 0, best = FALSE, now = 0, bestscore = -1, ***scores_original = NULL, **scores_changed = NULL;
	char *pruned_tree = NULL;
	char *tmp = malloc(TREE_LENGTH * sizeof(char));
	char *ideal = malloc(TREE_LENGTH * sizeof(char));
	char userinfile[1000], c, outputfile[1000];
	char *inputtree = malloc(TREE_LENGTH * sizeof(char));
	if(!tmp || !ideal || !inputtree) { free(tmp); free(ideal); free(inputtree); printf2("Error: out of memory in spr_dist\n"); return; }
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
						printf2("Error there are no trees in memory\n");
						error = TRUE;
						}
					}
				else
					{
					if((infile = fopen(parsed_command[i+1], "r")) == NULL)
						{
						printf2("Error opening file named %s\n", parsed_command[i+1]);
						error = TRUE;
						strcpy(userinfile, parsed_command[i+1]);
						}
					else
						{
						printf2("opened input file %s\n", parsed_command[i+1]);
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
					printf2("Error: %s not valid option for dorandomisation\n", parsed_command[i+1]);
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
		printf2("Error opening file named %s\n", outputfile);
		error = TRUE;
		}
	else
		{
		printf2("opened output file %s\n", outputfile);
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
			printf2("\nbestreal = %f\n", bestreal);
			hsprint = FALSE;
			if(signal(SIGINT, controlc4) == SIG_ERR)
				printf2("An error occurred while setting a signal handler\n");
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
			

			if(starting_super != 2) { int _to = 0; tree_build(1, retained_supers[0], tree_top, TRUE, -1, &_to); }
			else
				{ int _to = 0; tree_build(1, inputtree, tree_top, TRUE, -1, &_to); }

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
					if(print_pruned_tree(tree_top, 0, pruned_tree, FALSE, 0) >1)
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
			if(y==0)printf2("amount SPR per tree required to make Ideal as incongruent as real data = %f\n", totalnow/(Total_fund_trees-num_excluded_trees));
			else printf2("amount SPR per tree required to make Ideal as incongruent as random data = %f\n", totalnow/(Total_fund_trees-num_excluded_trees));
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

	free(tmp);
	free(ideal);
	free(inputtree);
	}


int string_SPR(char * string)
	{
	int i=0, j=0, k=0, l=0, components = 0, random_num = 0, done = FALSE, found = FALSE, **scores_original = NULL, **scores_changed = NULL, attempts = 0;
	char *extracted = malloc(TREE_LENGTH * sizeof(char));
	char *temptree = NULL, *tmp = NULL;
	char *original = malloc(TREE_LENGTH * sizeof(char));
	char *string1 = NULL;

	if(!extracted || !original) { free(extracted); free(original); printf2("Error: out of memory in string_SPR\n"); return 0; }
	
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
		 random_num = (int)fmod(rand_r(&thread_seed), components)+1;
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
		
		{ int _to = 0; tree_build(1, string1, tree_top, FALSE, -1, &_to); }

		tree_top = temp_top;
		temp_top = NULL;

		/**** shrink the tree ***/
		shrink_tree(tree_top);
		temptree[0] = '\0';
		/** printf the result **/
		if(print_pruned_tree(tree_top, 0, temptree, FALSE, 0) >1)
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
		 random_num = (int)fmod(rand_r(&thread_seed), components)+1;
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
		for(j=l; j >= i + (int)strlen(extracted); j--) {
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
	free(extracted);
	free(original);

	return(attempts);
	}





void exhaustive_SPR(char * string)
	{
	int i, j, k, l, x, y, q, r=-1, labelonly = FALSE, exnum, components =0, pruned_components =0, num, *component_index = NULL, *pruned_component_index = NULL, **scores_original = NULL, **scores_changed = NULL;
	char *labeledtree = malloc(TREE_LENGTH * sizeof(char));
	char *tmp_labeledtree = malloc(TREE_LENGTH * sizeof(char));
	char *extractedpart = malloc(TREE_LENGTH * sizeof(char));
	char *tmp_tree = malloc(TREE_LENGTH * sizeof(char));
	char *pruned_tree = malloc(TREE_LENGTH * sizeof(char));
	char *tmp = malloc(TREE_LENGTH * sizeof(char));
	char taxaname[100], filename[100], filename1[100], pasted_name[10000], cut_name[100];
	FILE *sproutfile = NULL, *sprdescriptor = NULL, *labeledtreefile = NULL;

	if(!labeledtree || !tmp_labeledtree || !extractedpart || !tmp_tree || !pruned_tree || !tmp)
		{
		free(labeledtree); free(tmp_labeledtree); free(extractedpart);
		free(tmp_tree); free(pruned_tree); free(tmp);
		printf2("Error: out of memory in exhaustive_SPR\n"); return;
		}

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
	printf2("written labelledtree to file\n");
	
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
			{ int _to = 0; tree_build(1, tmp_tree, tree_top, FALSE, -1, &_to); }

			tree_top = temp_top;
			temp_top = NULL;
			/**** shrink the tree ***/
			shrink_tree(tree_top);
			pruned_tree[0] = '\0';
			
			/** print the result **/
			if(print_pruned_tree(tree_top, 0, pruned_tree, FALSE, 0) >1)
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
	free(labeledtree);
	free(tmp_labeledtree);
	free(extractedpart);
	free(tmp_tree);
	free(pruned_tree);
	free(tmp);

	}

void neighbor_joining(int brlens, char *tree, int names)
	{
	int num_nodes = number_of_taxa, *deleted = NULL, i, j, smallest_i = -1, smallest_j = -1;
	float *transformed = NULL, smallest = 0, vi, vj, vtmp;
	char **tree_structure = NULL;
	char *string = malloc(TREE_LENGTH * sizeof(char));
	char *tmp = malloc(TREE_LENGTH * sizeof(char));
	char c;

	if(!string || !tmp) { free(string); free(tmp); printf2("Error: out of memory in neighbor_joining\n"); return; }
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
	free(string);
	free(tmp);
	}


    


void tips(int num)
	{
	switch(num)
		{
		case(0):
			printf2("\n\t1. Clann can be used to transform nexus formatted tree files into newick formatted files.\n\tThis is done by executing the nexus file as normal and then using the command: \"showtrees savetrees=yes\"\n\tIt is also possible to set the name of the file to which the trees are saved, and to stop clann from displaying a graphical representation of each source tree while this is done\n\n");
			break;
		
		case(1):
			printf2("\n\t2. Clann can be told only to read the first few characters of each taxa name whenreading the source trees into memory\n\tThis is useful when it is necessary to have unique identifiers (for instance gene IDs) on the source trees\n\tThe option \'maxnamelen\' in the \"exe\" command sets this value\n\tIf the names are not fixed widths, the option \'maxnamelen=delimited\' tells Clann to look for the fist dot \".\" which will specifying the end of the taxon ID in the trees\n\n\t\tFor instance using \"exe maxnamelen=delimited\" on the following tree:\n\t\t(apple.00121,(orange.1435,lemon.3421),pear.1032);\n\t\tResults in clann ignoring the numbers after the dots in the taxa names\n");
			break;

		case(2):
			printf2("\n\t3. The equals sign (=), hyphen (-) and space ( ) are special characters in Clann and by default cannot be used in filenames to be read by clann\n\tIf a filename contain some of these characters Clann can only read the name of the file properly by putting the name in inverted commas.\n\t\tFor example: exe \"my-file.txt\"\n");
			break;

		case(3):
			printf2("\n\t4. The first command that you should run if you don’t know what to do is \"help\".\n\tThis will display the list of the commands that are available.\n\tCalling any of the commands followed by a question mark (for instance \"hs ?\"), will display the options and defaults associated with that command\n");
			break;
		
		case(4):
			printf2("\n\t5. The command \"!\" runs a shell terminal on Unix and Mac operating systems allowing system commands can be run without having to quit Clann\n");
			break;

		case(5):
			printf2("\n\t6. Clann can assess supertrees created using other programs\n\tUsing the \"usertrees\" command, clann will read in the file specified and assess all the trees it contains\n\tThe best supertree found in the file is displayed along with its score\n");
			break;

		case(6):
			printf2("\n\t7. All commands in Clann should be written completely in lowercase, typing the command \"boot\" is not the same as \"Boot\" and only the first will be recognised as a valid command\n");
			break;

		case(7):
			printf2("\n\t8. Heuristic and exhaustive searches of supertree space can be interrupted using the key combination \"control-c\"\n\tThis dispays the score if the best tree found so far and give the user the option to stop the search now or continue.\n\tIf this is done during the random sampling phase of a heuristic search, it will allow the user to move straight to the heuristic search without completing the random sampling\n");
			break;

		case(8):
			printf2("\n\t9. Users can assess different configurations of their data by excluding (or including) certain source trees from subsequent commands using the \"excludetrees\" and \"includetrees\" commands\n\tSource trees can be selected based on their name,the taxa they contain, their size (number of taxa they contain) or their score when compared to a supertree\n");
			break;

		case(9):
			printf2("\n\t10. Individual (or multiple) taxa can be pruned from the source trees using the command \"deletetaxa\"\n\tBranch lengths are adjusted to take the deletion of thetaxa into account\n\tIf the deletion of taxa from a source tree means that there are less than 4 taxa remaining, that source tree is removed from the analysis\n\tClann will display the names of the source trees removed if this occurs\n");
			break;

		default:
			printf2("\n\t1. Clann can be used to transform nexus formatted tree files into newick formatted files.\n\tThis is done by executing the nexus file as normal and then using the command: \"showtrees savetrees=yes\"\n\tIt is also possible to set the name of the file to which the trees are saved, and to stop clann from displaying a graphical representation of each source tree while this is done\n\n");
			printf2("\n\t2. Clann can be told only to read the first few characters of each taxa name when reading the source trees into memory\n\tThis is useful when it is necessary to have unique identifiers (for instance gene IDs) on the source trees\n\tThe option \'maxnamelen\' in the \"exe\" command sets this value\n\tIf the names are not fixed widths, the option \'maxnamelen=delimited\' tells Clann to look for the fist dot \".\" which will specify the end of the taxon ID in the trees\n\n\t\tFor instance using \"exe maxnamelen=delimited\" on the following tree:\n\t\t(apple.00121,(orange.1435,lemon.3421),pear.1032);\n\t\tResults in clann ignoring the numbers after the dots in the taxa names\n");
			printf2("\n\t3. The equals sign (=), hyphen (-) and space ( ) are special characters in Clann and by default cannot be used in filenames to be read by clann\n\tIf a filename contain some of these characters Clann can only read the name of the file properly by putting the name in inverted commas.\n\t\tFor example: exe \"my-file.txt\"\n");
			printf2("\n\t4. The first command that you should run if you don’t know what to do is \"help\".\n\tThis will display the list of the commands that are available.\n\tCalling any of the commands followed by a question mark (for instance \"hs ?\"), will display the options and defaults associated with that command\n");
			printf2("\n\t5. The command \"!\" runs a shell terminal on Unix and Mac operating systems allowing system commands can be run without having to quit Clann\n");
			printf2("\n\t6. Clann can assess supertrees created using other programs\n\tUsing the \"usertrees\" command, clann will read in the file specified and assess all the trees it contains\n\tThe best supertree found in the file is displayed along with its score\n");
			printf2("\n\t7. All commands in Clann should be written completely in lowercase, typing the command \"boot\" is not the same as \"Boot\" and only the first will be recognised as a valid command\n");
			printf2("\n\t8. Heuristic and exhaustive searches of supertree space can be interrupted using the key combination \"control-c\"\n\tThis dispays the score if the best tree found so far and give the user the option to stop the search now or continue.\n\tIf this is done during the random sampling phase of a heuristic search, it will allow the user to move straight to the heuristic search without completing the random sampling\n");
			printf2("\n\t9. Users can assess different configurations of their data by deleting certain source trees using the \"deletetrees\" command\n\tSource trees can be selected based on their name,the taxa they contain, their size (number of taxa they contain) or their score when compared to a supertree\n");
			printf2("\n\t10. Individual (or multiple) taxa can be pruned from the source trees using the command \"deletetaxa\"\n\tBranch lengths are adjusted to take the deletion of thetaxa into account\n\tIf the deletion of taxa from a source tree means that there are less than 4 taxa remaining, that source tree is removed from the analysis\n\tClann will display the names of the source trees removed if this occurs\n");
			break;

		}

	}


/* This function controls the redirection of output from CLann to a log file (off by default), allow the use of the overwritten "printf" function (below) */


/* -----------------------------------------------------------------------
 * execute_recluster
 *
 * Standalone command: read a landscape TSV previously written by
 * 'hs visitedtrees=<file>' and cluster its topologies at a user-supplied
 * RF threshold.  No gene trees need to be loaded.
 *
 * Syntax (REPL):
 *   recluster <landscapefile> [clusterthreshold=<f>] [clusteroutput=<file>]
 *              [clusterorderby=score|visits]
 *
 * CLI:
 *   clann recluster landscape.tsv clusterthreshold=0.1 clusteroutput=out.tsv
 * ----------------------------------------------------------------------- */
void execute_recluster(void)
    {
    int   i, error = FALSE;
    char  landscape_file[4096] = "";
    char  cluster_output[4096] = "treeclusters.tsv";
    float cluster_threshold    = 0.2f;
    int   cluster_orderby      = 0;
    LandscapeMap *lm;

    /* parsed_command[0] == "recluster"
     * parsed_command[1] is the landscape file (positional) unless it looks
     * like an option key.  Option keys never contain '/' or '.' as first
     * char, but filenames can — so we treat parsed_command[1] as the file
     * if it contains no '=' AND is non-empty.                              */

    if(num_commands >= 2 && parsed_command[1][0] != '\0' &&
       strchr(parsed_command[1], '=') == NULL &&
       strcmp(parsed_command[1], "clusterthreshold") != 0 &&
       strcmp(parsed_command[1], "clusteroutput")    != 0 &&
       strcmp(parsed_command[1], "clusterorderby")   != 0)
        {
        strncpy(landscape_file, parsed_command[1], sizeof(landscape_file) - 1);
        landscape_file[sizeof(landscape_file) - 1] = '\0';
        }

    /* Parse key=value options from parsed_command[1..] */
    for(i = 1; i < num_commands; i++)
        {
        if(strcmp(parsed_command[i], "clusterthreshold") == 0)
            {
            cluster_threshold = tofloat(parsed_command[i+1]);
            if(cluster_threshold < 0.0f || cluster_threshold > 1.0f)
                { printf2("Error: clusterthreshold must be between 0.0 and 1.0\n"); error = TRUE; }
            }
        else if(strcmp(parsed_command[i], "clusteroutput") == 0)
            {
            strncpy(cluster_output, parsed_command[i+1], sizeof(cluster_output) - 1);
            cluster_output[sizeof(cluster_output) - 1] = '\0';
            }
        else if(strcmp(parsed_command[i], "clusterorderby") == 0)
            {
            if(strcmp(parsed_command[i+1], "score") == 0)
                cluster_orderby = 0;
            else if(strcmp(parsed_command[i+1], "visits") == 0)
                cluster_orderby = 1;
            else
                { printf2("Error: clusterorderby must be score or visits\n"); error = TRUE; }
            }
        }

    if(error) return;

    if(landscape_file[0] == '\0')
        {
        printf2("Error: recluster requires a landscape file as its first argument.\n");
        printf2("  Usage: recluster <landscapefile> [clusterthreshold=<f>]\n");
        printf2("                   [clusteroutput=<file>] [clusterorderby=score|visits]\n");
        return;
        }

    printf2("Reading landscape file: %s\n", landscape_file);
    lm = lm_read(landscape_file);
    if(!lm)
        { printf2("Error: failed to read landscape file '%s'\n", landscape_file); return; }

    printf2("  Unique topologies loaded: %zu\n", lm->count);

    /* ---- Setup taxa if no gene trees are loaded ----
     * lm_cluster uses the global taxa_names[] / taxon_hash_vals[] to parse
     * Newick bipartitions.  If no gene trees have been loaded yet, derive
     * the taxon list from the Newick strings in the landscape map so that
     * RF distances can be computed correctly.
     */
    {
    int   save_ntaxa   = number_of_taxa;
    char **save_names  = taxa_names;
    uint64_t *save_hv  = taxon_hash_vals;
    int   rc_ntaxa     = 0;
    char **rc_taxa     = NULL;
    uint64_t *rc_hv    = NULL;

    if(number_of_taxa == 0)
        {
        /* Collect unique leaf names from all Newick strings in the map */
        size_t  j;
        int     cap = 256;
        rc_taxa = malloc((size_t)cap * sizeof(char *));
        if(rc_taxa)
            {
            for(j = 0; j < lm->capacity; j++)
                {
                const char *nwk = lm->slots[j].newick;
                int k;
                if(lm->slots[j].hash == 0 || !nwk) continue;
                /* Walk the Newick string, collecting leaf tokens */
                k = 0;
                while(nwk[k] && nwk[k] != ';')
                    {
                    if(nwk[k] == '(' || nwk[k] == ')' || nwk[k] == ',')
                        { k++; continue; }
                    if(nwk[k] == ':')
                        { while(nwk[k] && nwk[k] != ',' && nwk[k] != ')' && nwk[k] != ';') k++; continue; }
                    /* Read a name token */
                    {
                    char name[NAME_LENGTH + 1];
                    int  nlen = 0, found, t;
                    while(nwk[k] && nwk[k] != '(' && nwk[k] != ')' &&
                          nwk[k] != ',' && nwk[k] != ':' && nwk[k] != ';')
                        { if(nlen < NAME_LENGTH) name[nlen++] = nwk[k]; k++; }
                    name[nlen] = '\0';
                    if(nlen == 0) continue;
                    /* Skip if this is an internal-node label (after a ')') */
                    /* Already checked: we only get here from a non-'(' non-')' char */
                    found = 0;
                    for(t = 0; t < rc_ntaxa; t++)
                        if(rc_taxa[t] && strcmp(rc_taxa[t], name) == 0) { found = 1; break; }
                    if(!found)
                        {
                        if(rc_ntaxa == cap)
                            {
                            int newcap = cap * 2;
                            char **tmp = realloc(rc_taxa, (size_t)newcap * sizeof(char *));
                            if(!tmp) break;
                            rc_taxa = tmp; cap = newcap;
                            }
                        rc_taxa[rc_ntaxa++] = strdup(name);
                        }
                    }
                    }
                }

            if(rc_ntaxa > 0)
                {
                /* Assign splitmix64 weights (same formula as execute_command) */
                int t;
                rc_hv = malloc((size_t)rc_ntaxa * sizeof(uint64_t));
                if(rc_hv)
                    {
                    for(t = 0; t < rc_ntaxa; t++)
                        {
                        uint64_t x = (uint64_t)(t + 1);
                        x += 0x9e3779b97f4a7c15ULL;
                        x  = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
                        x  = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
                        rc_hv[t] = x ^ (x >> 31);
                        }
                    number_of_taxa  = rc_ntaxa;
                    taxa_names      = rc_taxa;
                    taxon_hash_vals = rc_hv;
                    printf2("  Derived %d taxa from landscape Newick strings\n", rc_ntaxa);
                    }
                }
            }
        }

    lm_cluster(lm, cluster_output, cluster_threshold, cluster_orderby);

    /* Restore global taxa state if we replaced it */
    if(save_ntaxa == 0 && rc_ntaxa > 0)
        {
        int t;
        number_of_taxa  = save_ntaxa;
        taxa_names      = save_names;
        taxon_hash_vals = save_hv;
        free(rc_hv);
        for(t = 0; t < rc_ntaxa; t++) free(rc_taxa[t]);
        free(rc_taxa);
        }
    }

    lm_free(lm);
    }


/* printf2: moved to utils.c */
