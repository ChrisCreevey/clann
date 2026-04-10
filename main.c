/*
 *  main.c — Clann v5.0.0
 *  Entry point, command dispatch, and CLI interface.
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

#include "clann.h"
#include "utils.h"
#include "tree_io.h"
#include "viz.h"
#include "main.h"


/* ===== CLI direct-command interface ===================================
 *
 * Allows CLANN to be called as:
 *   clann <command> <treefile> [key=value ...]
 *   clann usertrees <source> <candidates> [key=value ...]
 *   clann --help
 *   clann hs --help
 *
 * Options with '=' are parsed as key=value pairs.
 * Named file flags: --source=<f>/--trees=<f> set the gene-tree file;
 *   --topologies=<f>/--candidates=<f> set the usertrees candidates file.
 * Global options (criterion, nthreads, mlbeta, mlscale, seed) become
 * "set X Y" commands; all others are appended to the command string.
 * Leading '--' is stripped so GNU-style --criterion=ml also works.
 * ===================================================================== */

static const char *cli_known_commands[] = {
    "hs", "hsearch", "alltrees", "usertrees", "consensus", "nj", "boot", "bootstrap", NULL
};

static int is_cli_command(const char *s)
    {
    int i;
    for(i = 0; cli_known_commands[i]; i++)
        if(strcmp(s, cli_known_commands[i]) == 0) return TRUE;
    return FALSE;
    }

/* Print usage.  cmd==NULL: general help.  cmd!=NULL: per-command help
 * (synthesises "<cmd> ?" into stored_commands[] for the REPL to handle). */
static void cli_usage(const char *cmd)
    {
    if(cmd == NULL)
        {
        printf("\nUsage:\n");
        printf("  clann <command> <treefile> [options]\n");
        printf("  clann usertrees --source=<gene-trees> --topologies=<candidates> [options]\n\n");
        printf("Commands:\n");
        printf("  hs / hsearch   Heuristic supertree search\n");
        printf("  alltrees       Exhaustive supertree search (small datasets)\n");
        printf("  usertrees      Score / test user-supplied topologies (requires two files)\n");
        printf("  consensus      Consensus tree from source trees\n");
        printf("  nj             Neighbour-joining supertree\n\n");
        printf("File options (all commands):\n");
        printf("  --source=<f> / --trees=<f>        Source gene-tree file (positional arg also accepted)\n");
        printf("  --topologies=<f> / --candidates=<f>  Candidate-topology file (usertrees only)\n\n");
        printf("Global options (applied before the command):\n");
        printf("  criterion=<c>   dfit (default), ml, rf, sfit, qfit, avcon\n");
        printf("  nthreads=<n>    threads for parallel search (default: all CPUs)\n");
        printf("  mlbeta=<f>      beta for ML criterion (default 1.0)\n");
        printf("  mlscale=<s>     lnl (default), paper, lust\n");
        printf("  seed=<n>        random seed\n\n");
        printf("Examples:\n");
        printf("  clann hs gene_trees.ph criterion=ml nthreads=4\n");
        printf("  clann usertrees --source=gene_trees.ph --topologies=candidates.ph criterion=ml tests=yes\n\n");
        printf("Run 'clann <command> --help' for command-specific options.\n");
        printf("Run 'clann' with no arguments to start the interactive session.\n\n");
        }
    else
        {
        /* Let the existing REPL help system handle per-command help */
        snprintf(stored_commands[parts], 10000, "%s ?", cmd);
        parts++;
        }
    }

/* Populate stored_commands[] with synthesised CLANN commands derived
 * from argv[start..argc-1].
 *
 * Returns  TRUE  : CLI mode; caller should set command_line / doexecute_command
 *          -1    : general --help printed; caller should clean_exit(0)
 *          -2    : per-command help queued; run REPL without loading a tree
 *
 * out_treefiles[0] = source-trees file (for exe / execute_command)
 * out_treefiles[1] = candidate-topologies file (usertrees only)
 */
static int cli_dispatch(int argc, char *argv[], int start,
                        char out_treefiles[2][1000])
    {
    static const char *global_keys[] = {
        "criterion", "mlbeta", "mlscale", "seed", NULL
    };

    const char *cmdname = argv[start];
    char main_cmd[10000];
    int file_count = 0, k, g, is_usertrees;

    is_usertrees = (strcmp(cmdname, "usertrees") == 0);

    /* Scan for --help / -h */
    for(k = start; k < argc; k++)
        {
        const char *a = argv[k];
        while(*a == '-') a++;
        if(strcmp(a, "help") == 0 || strcmp(a, "h") == 0)
            {
            /* Per-command help: queue "<cmd> ?" for the REPL to print */
            snprintf(stored_commands[parts], 10000, "%s ?", cmdname);
            parts++;
            return -2;  /* run REPL, no tree needed */
            }
        }

    snprintf(main_cmd, sizeof(main_cmd), "%s", cmdname);

    for(k = start + 1; k < argc; k++)
        {
        char stripped[10000];
        char key[512], val[512];
        char *p, *eq;

        strncpy(stripped, argv[k], sizeof(stripped) - 1);
        stripped[sizeof(stripped) - 1] = '\0';

        /* Strip leading -- or - */
        p = stripped;
        if(p[0] == '-' && p[1] == '-')      p += 2;
        else if(p[0] == '-' && p[1] != '\0') p += 1;

        /* Split on first '=' */
        eq = strchr(p, '=');
        if(eq)
            {
            int klen = (int)(eq - p);
            strncpy(key, p, klen); key[klen] = '\0';
            strncpy(val, eq + 1, sizeof(val) - 1); val[sizeof(val) - 1] = '\0';
            }
        else
            {
            strncpy(key, p, sizeof(key) - 1); key[sizeof(key) - 1] = '\0';
            val[0] = '\0';
            }

        /* No '=' → file argument */
        if(val[0] == '\0')
            {
            if(file_count == 0)
                strncpy(out_treefiles[0], key, 999);
            else if(file_count == 1 && is_usertrees)
                {
                /* usertrees needs the candidates file as its first argument */
                snprintf(main_cmd, sizeof(main_cmd), "%s %s", cmdname, key);
                }
            file_count++;
            continue;
            }

        /* Named file flags: --source=<f> / --trees=<f> → source-trees file
         *                   --topologies=<f> / --candidates=<f> → usertrees candidates */
        if(strcasecmp(key, "source") == 0 || strcasecmp(key, "trees") == 0)
            { strncpy(out_treefiles[0], val, 999); continue; }
        if(is_usertrees &&
           (strcasecmp(key, "topologies") == 0 || strcasecmp(key, "candidates") == 0))
            { snprintf(main_cmd, sizeof(main_cmd), "%s %s", cmdname, val); continue; }

        /* Is it a global set option? */
        int is_global = FALSE;
        for(g = 0; global_keys[g]; g++)
            if(strcasecmp(key, global_keys[g]) == 0) { is_global = TRUE; break; }

        if(is_global)
            {
            snprintf(stored_commands[parts], 10000, "set %s %s", key, val);
            parts++;
            }
        else
            {
            /* Command-specific: append to main_cmd */
            snprintf(main_cmd + strlen(main_cmd),
                     sizeof(main_cmd) - strlen(main_cmd),
                     " %s %s", key, val);
            }
        }

    /* Enqueue the main command */
    snprintf(stored_commands[parts], 10000, "%s", main_cmd);
    parts++;
    return TRUE;
    }

/* ===== end CLI helpers =============================================== */


int main(int argc, char *argv[])
    {
	
    int i = 0, j=0, k=0, l=0, m=0, error=FALSE, x, doexecute_command = FALSE, command_line = FALSE, tipnum=0, autocommands=FALSE;
    char c, s, *command = NULL, *prev_command = NULL, HOME[1000], PATH[1000], exefilename[1000], string[10000];
    time_t time1, time2;
    double diff=0;    
    FILE *tmpclann = NULL;
	
	exefilename[0] = '\0';
	inputfilename[0] = '\0';
    saved_supertree[0] = '\0';
    logfile_name[0] = '\0';
	strcpy(logfile_name, "clann.log");
	system_call[0] ='\0';

	test_array = malloc(TREE_LENGTH*sizeof(int));
	test_array[0] = '\0';

	
	/********** START OF READLINE STUFF ***********/
	#ifdef HAVE_READLINE    
	if(command_line == FALSE)
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
    seed=(int)((time(NULL)/2)+getpid());
    srand((unsigned) (seed));
    #else
    seed=(int)(time(NULL)/2);
    srand((unsigned) (seed));
    #endif

    command = malloc(1000000*sizeof(char));
    if(!command) memory_error(74);
    command[0] = '\0';

    prev_command = malloc(1000000*sizeof(char));
    if(!prev_command) memory_error(74);
    prev_command[0] = '\0';

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
    for(i=0; i<10000; i++)
        {
        parsed_command[i] = malloc(10000*sizeof(char));
        parsed_command[i][0] = '\0';
        }


	/* Pre-getopt: handle --help before BSD getopt tries to parse double-dash */
	if(argc >= 2 && strcmp(argv[1], "--help") == 0)
		{
		cli_usage(NULL);
		clean_exit(0);
		}

	if(argc > 1)
		{
		while ((c = getopt(argc, argv, "nlhc:")) != -1)
		    {   
			switch (c) 
			      {
			      case 'n':
				        command_line = TRUE;
						printf("\nNon-Internactive mode entered - commands must be provided in a nexus \'clann block\' or with \'-c\'\n");
				        break;
			      case 'l':
				        printf2("opening logfile %s\n", logfile_name);	
			        	if((logfile = fopen(logfile_name, "w")) == NULL)
			        		{	
			        		printf2("Error opening log file named %s\n", logfile_name);
			        		}
			        	else
			        		{
			        		printf2("starting logging to logfile %s\n", logfile_name );
			        		print_log = TRUE;
			        		}
				        break;
				  case 'c':
				  		commands_filename = optarg;
				  		printf2("opening commands file %s\n", commands_filename);

			        	if((commands_file = fopen(commands_filename, "r")) == NULL)
			        		{	
			        		printf2("Error opening commands file named %s\n", commands_filename);
			        		}
			        	else
			        		{
			        		s = getc(commands_file);
							while(!feof(commands_file) && !error)
								{
								i=0;
								while((s == ' ' || s == '\n' || s == '\r' || s == '\t' || s == ';') && !feof(commands_file)) s = getc(commands_file);
								if(s == '#') /* skip over lines that start with a '#' -- this allows the user to add comments */
									{
									while(!feof(commands_file) && s != '\n' && s != '\r') s = getc(commands_file);
									}	
								while(s != ';' && s != '\n' && !feof(commands_file)) /* commands can finish with a ';' or a newline character */
									{
									string[i] = s;
									i++;
									s = getc(commands_file);
									}
								string[i] = ';';
								string[i+1] = '\0';
								strcpy(stored_commands[parts], string);
								parts++;
								}
							autocommands=TRUE; /* this is used indicate that we need to print the commands both to screen and log, as otherwise they will not appear */
							fclose(commands_file);
			        		}
			        	break;
			      case 'h':
			      		print_splash();
			      		

			      		printf("\nUsage: \"clann -lnh [-c commands file] [tree file]\"\n");
			      		printf("\n\tWhere [tree file] is an optional Nexus or Phylip formatted file of phylogenetic trees\n");
			      		printf("\t-l turn on logging of screen output to file \"clann.log\"\n");
			      		printf("\t-n turns off interactive mode - requires commands to be provided in a nexus \'clann block\' or with \'-c\'\n");
			      		printf("\t-c <file name> specifies a file with commands to be executed (one command per line)\n");
			      		printf("\t-h prints this message\n\n");

				        print_commands(0);

				        clean_exit(0);
				        break;
			      }
		 	}
		for (i = optind; i < argc; i++)
			{
			doexecute_command = TRUE;
			strcpy(exefilename, argv[i]);
			}

		/* ---- CLI direct-command mode ------------------------------------ *
		 * If the first non-flag argument is a known command name, synthesise *
		 * the stored_commands queue and run non-interactively.              *
		 * ------------------------------------------------------------------ */
		if(optind < argc && is_cli_command(argv[optind]))
			{
			char cli_treefiles[2][1000];
			int cli_r;
			cli_treefiles[0][0] = '\0';
			cli_treefiles[1][0] = '\0';
			/* Reset any state set by the earlier non-flag argv loop */
			doexecute_command = FALSE;
			exefilename[0] = '\0';
			cli_r = cli_dispatch(argc, argv, optind, cli_treefiles);
			if(cli_r == -1) clean_exit(0);   /* general --help printed */
			if(cli_r == -2)
				{
				/* Per-command help: run REPL with "cmd ?" queued, no tree */
				command_line = TRUE;
				}
			if(cli_r == TRUE)
				{
				command_line = TRUE;
				if(cli_treefiles[0][0] != '\0')
					{
					doexecute_command = TRUE;
					strcpy(exefilename, cli_treefiles[0]);
					}
				}
			/* prevent the earlier non-flag loop from double-loading */
			optind = argc;
			}
		}
	if(command_line == TRUE)
		{
		strcpy(stored_commands[parts], "quit");
		parts++;
		}

    print_splash(); /* Print the Clann splash screen */

     
    if(doexecute_command) /* if the user has specified an input file at the command line */
        {
        execute_command(exefilename, TRUE);
        }
    

    i=0;
    while(strcmp(parsed_command[0], "quit") != 0)  /* while the user has not chosen to quit */
        {
        
        if(num_commands > 0)
            {
            
            if(print_log && strcmp(command, "") != 0)
            	{ 
            	fprintf(logfile, "clann> %s\n", command);
            	if(autocommands==TRUE) printf("clann> %s\n", command);
            	}
            if(strcmp(parsed_command[0], "execute") == 0 || strcmp(parsed_command[0], "exe") == 0)
                {
                if(num_commands == 2 && parsed_command[1][0] == '?')
                    print_commands(1);
                else
                    {
                
                    if(num_commands > 1)
                        {
                        delimiter = TRUE;   /* default: on; user can override with maxnamelen=full */
                        execute_command(parsed_command[1], TRUE);
							
						/*spr_dist(); */ 
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
                                if(criterion == 0 || criterion == 2 || criterion == 3 || criterion == 5 || criterion == 6 || criterion == 7)
                                    alltrees_search(TRUE);
                                else
                                    {
                                    BR_file = fopen("coding.nex", "w");
                                    x = coding(0, 0, 0);
                                    fclose(BR_file);
									if(x == FALSE)
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
                                printf2("Error: You need to load source trees before using this command\n");
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
                                        printf2("Error: You cannont do a user tree search in MRP\n");
                                    else
                                        usertrees_search();
                                    }
                                else
                                    {
                                    printf2("Error: You need to load source trees before using this command\n");
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
                                        printf2("Error: You need to load source trees before using this command\n");
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
                                            printf2("Error: You need to load source trees before using this command\n");
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
                                                printf2("Error: You need to load source trees before using this command\n");
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
                                                    printf2("Error: You need to load source trees before using this command\n");
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
                                                    if(strcmp(parsed_command[0], "excludetrees") == 0)
                                                        {
                                                        if(num_commands == 2 && parsed_command[1][0] == '?')
                                                            print_commands(17);
                                                        else
															{
															if(number_of_taxa > 0)
																exclude(TRUE);
															else
																printf2("Error: You need to load source trees before using this command\n");
															}
														}
                                                    else
                                                        {
                                                        if(strcmp(parsed_command[0], "includetrees") == 0)
                                                            {
                                                            if(num_commands == 2 && parsed_command[1][0] == '?')
                                                                print_commands(18);
                                                            else
																{
																if(number_of_taxa > 0)
																		include(TRUE);
																else
																	printf2("Error: You need to load source trees before using this command\n");
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
                                                                    printf2("\n\tType exit to return to Clann\n\n");
                                                                    system("bash");
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
																			printf2("Error: You need to load source trees before using this command\n");
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
																					showtrees(FALSE);
																			else
																				printf2("Error: You need to load source trees before using this command\n");
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
																					printf2("Error: You need to load source trees before using this command\n");
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
																						printf2("Error: You need to load source trees before using this command\n");
																					}
																				}

																			else
																				{
																				if(strcmp(parsed_command[0], "restoretaxa") == 0)
																					{
																					if(number_of_taxa > 0)
																						{
																						if(restoretaxa_available)
																							restoretaxa(TRUE);
																						else
																							printf2("Error: no deletetaxa snapshot available to restore\n");
																						}
																					else
																						printf2("Error: You need to load source trees before using this command\n");
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
																							printf2("Error: You need to load source trees before using this command\n");
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
																								printf2("Error: You need to load source trees before using this command\n");
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
																									printf2("Error: You need to load source trees before using this command\n");	
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
																										printf2("Error: You need to load source trees before using this command\n");	
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
																											printf2("Error: You need to load source trees before using this command\n");	
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
																												printf2("Error: You need to load source trees before using this command\n");	
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
																													printf2("Error: You need to load source trees before using this command\n");	
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
																														printf2("Error: You need to load source trees before using this command\n");	
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
																															printf2("The source trees have now been randomised\n");
																															}
																														else
																															printf2("Error: You need to load source trees before using this command\n");	
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
																																printf2("Error: You need to load source trees before using this command\n");	
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
                                                                                                                                    printf2("Error: You need to load source trees before using this command\n"); 
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
	                                                                                                                        	{
	                                                                                                                        	if(strcmp(parsed_command[0], "log") == 0)
	                                                                                                                            	{                                     
	                                                                                                                            	if(num_commands == 2 && parsed_command[1][0] == '?')
		                                                                                                                                print_commands(29);
		                                                                                                                            else
		                                                                                                                            	{
		                                                                                                                            	do_log();
		                                                                                                                            	}
		                                                                                                                        	}
		                                                                                                                        else
		                                                                                                                            {
	                                                                                                                        		if(strcmp(parsed_command[0], "savetrees") == 0)
	                                                                                                                            		{                                     
	                                                                                                                            		if(num_commands == 2 && parsed_command[1][0] == '?')
		                                                                                                                                	print_commands(30);
		                                                                                                                            	else
		                                                                                                                            		{
		                                                                                                                            		if(number_of_taxa > 0)
																																				showtrees(TRUE);
																																			else
																																				{
																																				printf2("Error: You need to load source trees before using this command\n");
																																				}
		                                                                                                                            		}
		                                                                                                                        		}
		                                                                                                                        	else
	                                                                                                                        	    {
	                                                                                                                        	    if(strcmp(parsed_command[0], "mlscores") == 0)
	                                                                                                                        	        {
	                                                                                                                        	        if(num_commands == 2 && parsed_command[1][0] == '?')
	                                                                                                                        	            print_commands(31);
	                                                                                                                        	        else
	                                                                                                                        	            mlscores();
	                                                                                                                        	        }
	                                                                                                                        	    else
	                                                                                                                        			printf2("Error: command not known.\n\tType help at the prompt to get a list of available commands.\n");
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
                            }
                        }
                    }
                }
            } /* end if*/

        time2 = time(NULL);
        diff = difftime(time2, time1);
        if(diff > 0)
            {
            printf2("\nTime taken:"); 
            if(diff > 60)
                {
                if(diff > 3600)
                    {
                    if(diff > 86400)
                        {
                        printf2(" %0.0lf Day", diff/86400);
                        if(diff/86400 > 1) printf2("s");
                        diff = diff- (86400 * ((int)diff/86400));
                        }
                    printf2(" %0.0lf Hour", diff/3600);
                    if(diff/3600 > 1) printf2("s");
                    diff = diff - (3600 * ((int)diff/3600));
                    }
                printf2(" %0.0lf Minute", diff/60);
                if(diff/60 > 1) printf2("s");
                diff = diff - (60 * ((int)diff/60));
                }
            printf2(" %0.0lf second", diff);
            if(diff > 1) printf2("s");
            printf2("\n\n");
            }

        if(logfile!=NULL) fflush(logfile);

        if(i == parts)
            {
            parts =0;
			#ifdef HAVE_READLINE
/****/      free(command);
			#endif
			user_break = FALSE;
			if(signal(SIGINT, controlc4) == SIG_ERR)
				{
				printf2("An error occurred while setting a signal handler\n");
				}
		   #ifdef HAVE_READLINE
/*****/    command = readline("clann> "); 
/*****/    command = realloc(command, 10000*sizeof(char));
		   #else
 	        printf2("clann> ");
 	        fflush(stdout);
 	        command = xgets(command);
		   #endif
          	autocommands=FALSE; /* indicate that from this point on any commands have been entered by the user and have not come from a commands file */
           
            /* if the commmand doesn't end in a ";" then put one on */
            k=0;
            while(command[k] != '\0') k++;
            if(command[k-1] != ';')
                strcat(command, ";");  

            
            
			#ifdef HAVE_READLINE
            if(strcmp(command, ";") != 0 && command_line == FALSE)   /* if the command is not just a blank line */
                {
            	if(strcmp(prev_command, command) != 0)  /* This gets rid of multiple instances of the exact same command in a row being recorded */
            		{
/******/            add_history(command);
            		strcpy(prev_command, command);  
            		}
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
    free(prev_command);
    free(tempsuper);
	#ifdef HAVE_READLINE
/****/ if(command_line == FALSE) write_history(PATH);  /* write the history of commands to file */
	#endif
/*	printf2("malloc_check = %d x (%d)\n", malloc_check, sizeof(taxon_type));
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

void print_splash(void)
	{
	printf2("\n\n\t*********************************************************************");
    printf2("\n\t*                                                                   *");
    printf2("\n\t*                           Clann  v%s                           *", VERSION);
    printf2("\n\t* Investigating phylogenetic information through supertree analyses *");
    printf2("\n\t*                                                                   *");
    printf2("\n\t*                  lab: http://www.creeveylab.org                   *");
    printf2("\n\t*                  email: %s                   *", PACKAGE_BUGREPORT);
    printf2("\n\t*                                                                   *");
    printf2("\n\t*                 Copyright Chris Creevey 2003-2020                 *");
    printf2("\n\t*                                                                   *");
    printf2("\n\t*          HINT: Type \"help\" to see all available commands          *");
    printf2("\n\t*********************************************************************\n\n");

	}
 
    
void print_commands(int num)
    {
    
    if(num == 0)
        {
        printf2("\nAvailable Commands:\n\n");
        printf2("\nThe following commands are always available:\n\n");
        printf2("\texecute (exe)\t- Read in a file of source trees\n");
        printf2("\thelp\t\t- Display this message\n");
        printf2("\tquit\t\t- Quit Clann\n");
        printf2("\tset\t\t- Set global parameters such as optimality criterion for carrying reconstructing a supertree\n");
        printf2("\t!\t\t- Run a shell session, while preserving the current Clann session (type \'exit\' to return)\n");
        printf2("\ttips\t\t- Show tips and hints for better use of Clann\n");
        printf2("\tlog\t\t- Control logging of screen output to a log file\n");

        printf2("\nThe following commands are only available when there are source trees in memory:\n");

        printf2("\nSupertree reconstruction:\n");
        printf2("\ths\t\t- Carry out a heuristic search for the best supertree usign the criterion selected\n");
        printf2("\tbootstrap\t- Carry out a bootstrap supertree analysis using the criterion selected\n");
        printf2("\tnj\t\t- Construct a neighbour-joining supertree\n");
        printf2("\talltrees\t- Exhaustively search all possible supertrees\n");
        printf2("\tusertrees\t- Assess user-defined supertrees (from seperate file), to find the best scoring\n");
        printf2("\tconsensus\t- Calculate a consensus tree of all trees containing all taxa\n");

        printf2("\nSource tree selection and modification:\n");
        printf2("\tsavetrees\t- Save source trees to file in phylip format (subsets of trees can be chosen based on a number of criteria)\n");
        printf2("\tshowtrees\t- Visualise selected source trees in ASCII format (also can save selected trees to file)\n");
       	printf2("\texcludetrees\t- Exclude source trees from analyses (based on a variety of criteria)\n");
       	printf2("\tincludetrees\t- Restore previously excluded source trees (based on a variety of criteria)\n");
        printf2("\tdeletetaxa\t- Specify taxa to delete from all source trees in memory (i.e. prune from the trees while preserving branch lengths)\n");
        printf2("\trestoretaxa\t- Restore the original trees before the last deletetaxa operation\n");
        printf2("\trandomisetrees\t- Randomises the source trees in memory, while preserving taxa composition in each tree\n");

        printf2("\nMiscellaneous calculations:\n");
        printf2("\trfdists\t\t- Calculate Robinson-Foulds distances between all source trees\n");
        printf2("\tgeneratetrees\t- Generate random supertrees & assess  against source trees in memory\n");
        printf2("\tyaptp\t\t- \"Yet another permutation-tail-probability\" test - performs a randomisation test\n");



        printf2("\nExperimental Options:\n");
        printf2("\treconstruct\t- Carry out a gene-tree reconciliation (source trees against a species tree)\n");
        printf2("\tprunemonophylies - Prunes clades which consist of multiple sequences from the same species, to a single representative\n");
        printf2("\tsprdists\t- Carry out estimation of SPR distances of real data versus ideal and randomised versions of the data\n");

        printf2("\n\n\nType a command followed by '?' in interactive mode to get information on the options available i.e.: \"exe ?\"\n");
        printf2("Full descriptions of the commands are available in the manual\n\n\n");


        }
   
        

    if(num == 1)
        {
        printf2("\nexecute (exe)\t<filename>\n");
        printf2("\nexe\t<filename>\n\n");
		 printf2("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf2("\t===========================================================\n");
		printf2("\n\tmaxnamelen\t<integer> | delimited | full\t*delimited");
		printf2("\n\t  (delimiter mode is ON by default: species names extracted before '.'");
		printf2("\n\t   use 'maxnamelen=full' to disable and use full taxon names)");
		printf2("\n\tdelimiter_char\t<character>\t\t\t'.'");
		printf2("\n\tsummary\t\tshort | long\t\t\t*long");
		printf2("\n\tautoprunemono\tyes | no\t\t\t*no");
		printf2("\n\t  (prune monophyletic same-species clades from multicopy trees at load time;");
		printf2("\n\t   trees that become single-copy after pruning join the supertree pool;");
		printf2("\n\t   original unpruned trees are preserved for 'reconstruct')");
		printf2("\n\tautoweight\tclan | splitviol\t\t*off");
		printf2("\n\t  clan:      weights = compatible_clans / testable_clans  (per-clan metric)");
		printf2("\n\t  splitviol: weights = 1 - violated_splits / total_splits  (per-split metric)");
		printf2("\n\t  Both modes report split score range alongside weights.");
		printf2("\n\t  A split is 'violated' only when non-clan taxa appear on both sides");
		printf2("\n\t  (internal clan splits are not counted as violations).");
		printf2("\n\t  Requires clanfile=<file>.");
		printf2("\n\tclanfile\t<filename>\t\t\t*none");
		printf2("\n\t  (clan file: one clan per line, taxon names separated by spaces or commas,");
		printf2("\n\t   lines beginning with '#' are comments)");
        }
    if(num == 2)
        {
        printf2("\nalltrees [options]\n\n");
        printf2("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf2("\t===========================================================\n");
        if(delimiter) printf2("\n\t(Note: delimiter mode is ON -- multicopy gene trees will be automatically\n\texcluded from the search; use 'maxnamelen=full' to disable)\n");
        if(criterion == 0 || criterion == 2 || criterion == 3 || criterion == 6 || criterion == 7)
            {
            printf2("\n\trange\t\t<treenumber> - <treenumber>\t*all\n\tsavetrees\t<filename>\t\t\ttop_alltrees.txt\n\tcreate\t\tyes | no\t\t\t*no\n");
            if(criterion == 0)
                {
                printf2("\tweight\t\tequal | comparisons\t\t");
                if(dweight == 1) printf2("comparisons\n");else printf2("euqal\n");
                }
            if(criterion == 2)
                {
                printf2("\tweight\t\tequal | splits\t\t\t");
                if(splits_weight == 1) printf2("equal\n");if(splits_weight == 2) printf2("splits\n");
                }
            if(criterion == 3)
                {
                printf2("\tweight\t\tequal | taxa | quartets\t\t");
                if(quartet_normalising == 1) printf2("equal\n");if(quartet_normalising == 2) printf2("taxa\n");if(quartet_normalising == 3) printf2("quartets\n");
                }
            }
        
        }
    if(num == 3)
        {
        printf2("\nusertrees <candidates-file> [options]\n");
        printf2("  (source gene trees must be loaded first with 'exe <file>')\n");
        printf2("  CLI form: clann usertrees --source=<gene-trees> --topologies=<candidates> [options]\n");
        printf2("        or: clann usertrees <gene-trees> <candidates> [options]\n\n");
        printf2("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf2("\t===========================================================\n");
        printf2("\n\toutfile\t\t<filename>\t\t\tUsertrees_result.txt");
        if(criterion == 0)
            {
            printf2("\n\tweight\t\tequal | comparisons\t\t");
            if(dweight == 1) printf2("comparisons\n");else printf2("splits\n");
            }
        if(criterion == 2)
            {
            printf2("\n\tweight\t\tequal | splits\t\t\t");
            if(splits_weight == 1) printf2("equal\n");if(splits_weight == 2) printf2("splits\n");
            }
        if(criterion == 3)
            {
            printf2("\n\tweight\t\tequal | taxa | quartets\t\t");
            if(quartet_normalising == 1) printf2("equal\n");if(quartet_normalising == 2) printf2("taxa\n");if(quartet_normalising == 3) printf2("quartets\n");
            }
		printf2("\n\tprintsourcescores\tyes | no\t\t\t*no");
	if(criterion == 7)
		{
		printf2("\n\ttests\t\tyes | no\t\t\t*no (ML topology tests; criterion=ml only)");
		printf2("\n\tnboot\t\t<integer>\t\t\t*1000 (SH bootstrap replicates)");
#ifdef _OPENMP
		printf2("\n\tnthreads\t<integer>\t\t\t*%-3d (parallel SH bootstrap; default=all CPUs)", omp_get_num_procs());
#endif
		printf2("\n\ttestsfile\t<filename>\t\t\t*mltest_results.txt");
		printf2("\n\tnormcorrect\t(flag)\t\t\t\t*off");
		printf2("\n\t  Apply Bryant & Steel (2008) normalising constant correction to lnL.");
		printf2("\n\t  Subtracts log(Z_{T|X_i}) per source tree, where Z_T = sum_m b_m(T)*e^{-beta*m}.");
		printf2("\n\t  Uses truncated large-beta expansion (b_0+b_2+b_4 terms); accurate for beta>~1.5.");
		printf2("\n\t  Improves absolute lnL accuracy; minor effect on tree rankings (see Bryant & Steel 2008).");
		}
	printf2("\n");
        }
    if(num == 4)
        {
        printf2("\nhs (or hsearch) [options]  \n");
        printf2("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf2("\t===========================================================\n");
        if(delimiter) printf2("\n\t(Note: delimiter mode is ON -- multicopy gene trees will be automatically\n\texcluded from the search; use 'maxnamelen=full' to disable)\n");
        if(criterion == 0 || criterion == 2 || criterion == 3 || criterion==5 || criterion==6 || criterion==7)
            {
            printf2("\n\tsample\t\t<integer number>\t\t*10,000\n\tnreps\t\t<integer number>\t\t*10");
            printf2("\n\tswap\t\tnni | spr | tbr\t\t\t");
            if(method == 1) printf2("nni");
            if(method == 2) printf2("spr"); 
            if(method == 3) printf2("tbr");
            printf2("\n\tnsteps\t\t<integer number>\t\t%d", number_of_steps);
			printf2("\n\tstart\t\tnj | random | memory | <filename>\tnj");
			printf2("\n\t  memory  = use best tree currently in memory as starting point");

			printf2("\n\tmaxswaps\t<integer number>\t\t*1,000,000\n\tsavetrees\t<filename>\t\t\tHeuristic_result.txt");
#ifdef _OPENMP
            printf2("\n\tnthreads\t<integer number>\t\t*%-3d (OpenMP threads; default=all CPUs; not available for criterion=recon)", omp_get_num_procs());
#endif
            printf2("\n\tmaxskips\t<integer number>\t\t*auto=2N² (stop replicate after this many consecutive already-visited moves; 0=disabled)");
#ifdef _OPENMP
            if(hs_progress_interval == 0)
                printf2("\n\tprogress\t<integer seconds | 0>\t\t0 (report every improvement)");
            else
                printf2("\n\tprogress\t<integer seconds | 0>\t\t%d (parallel best-so-far update interval; 0=every improvement)", hs_progress_interval);
            if(hs_droprep > 0.0f)
                printf2("\n\tdroprep\t\t<float 0-1>\t\t\t%.4f (abandon rep if score >X%% above global best; 0=disabled)", hs_droprep);
            else
                printf2("\n\tdroprep\t\t<float 0-1>\t\t\t0 (disabled; abandon rep if score >X%% above global best; parallel only)");
#endif
            printf2("\n\tautoprunemono\t(set via exe) prune monophyletic same-species clades from multicopy trees at load time");
            printf2("\n\tvisitedtrees\t<filename>\t\t\t(disabled) Record all visited topologies (TSV: newick, score, visit_count)");
            if(criterion == 0)
                {
                printf2("\n\tweight\t\tequal | comparisons\t\t");
                if(dweight == 1) printf2("comparisons");else printf2("splits");
                }
            if(criterion == 2)
                {
                printf2("\n\tweight\t\tequal | splits\t\t\t");
                if(splits_weight == 1) printf2("equal");if(splits_weight == 2) printf2("splits");
                }
            if(criterion == 3)
                {
                printf2("\n\tweight\t\tequal | taxa | quartets\t\t");
                if(quartet_normalising == 1) printf2("equal");if(quartet_normalising == 2) printf2("taxa");if(quartet_normalising == 3) printf2("quartets");
                }
			if(criterion == 0)
				{
				printf2("\n\tdrawhistogram\tyes | no\t\t\t*no\n\tnbins\t\t<integer number>\t\t*20\n\thistogramfile\t<filename>\t\t\t*Heuristic_histogram.txt");
				}
            
            }
        if(criterion == 1)
            {
            printf2("\n\tanalysis\tparsimony | nj\t\t\t*parsimony\n");
            printf2("\n\tParsimony options:\n\tweighted\tyes | no\t\t\t*no\n\tswap\t\tnni | spr | tbr\t\t\t*tbr\n\taddseq\t\tsimple | closest | asis |\n\t\t\trandom | furthest\t\t*random\n\tnreps\t\t<integer number>\t\t*10\n\n\tGeneral Options:\n\tsavetrees\t<filename>\t\t\tMRP.tree\n");
            }
		if(criterion == 4)
			{
			printf2("\n\tmissing\t4point | ultrametric\t\t\t*4point\n");
			}
		if(criterion == 5)
			{
			printf2("\n\tduplications\t<value>\t\t\t\t*1.0\n\tlosses\t\t<value>\t\t\t\t*1.0");
			printf2("\n\tnumspeciesrootings\t<value> | all\t\t*2\n\tnumgenerootings\t\t<value> | all\t\t*2\n");
			}
		if(criterion == 7)
			{
			printf2("\n\tmlbeta\t\t<float > 0>\t\t\t%.4f", ml_beta);
			printf2("\n\tmlscale\t\tpaper | lust | lnl\t\t%s", ml_scale==1?"lust":ml_scale==2?"lnl":"paper");
			printf2("\n\t  paper = Steel & Rodrigo (2008) formula: minimise beta*RF directly");
			printf2("\n\t  lust  = L.U.st (Akanni et al. 2014) log10 scaling (beta*d*log10(e))");
			printf2("\n\t  lnl   = report as lnL = -beta*RF  [default; matches ML tool conventions]");
			printf2("\n\tmleta\t\t<float >= 0>\t\t\t%.4f  [experimental]", ml_eta);
			printf2("\n\t  0     = Steel & Rodrigo (2008) raw RF (default)");
			printf2("\n\t  1     = divide RF by source-tree split count (full normalisation)");
			printf2("\n\t  >1    = actively down-weight large trees beyond normalisation\n");
			printf2("\n\t  Use 'mlscores eta=auto' to estimate optimal eta from data\n");
			}
        }
    if(num == 5)
        {
        printf2("\nbootstrap (or boot) [options]  \n\n\n");
        printf2("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf2("\t===========================================================\n");
        if(criterion != 1)
            {
            printf2("\n\tnreps\t\t<integer number>\t\t*100\n\thsreps\t\t<integer number>\t\t*1\n\tsample\t\t<integer number>\t\t*10,000\n\tswap\t\tnni | spr | tbr | all\t\t");
            if(method == 1) printf2("nni");
            if(method == 2) printf2("spr");
            printf2("\n\tstart\t\trandom | <filename>\t\trandom");
            
            printf2("\n\tnsteps\t\t<integer number>\t\t%d\n\ttreefile\t<output treefile name>\t\tbootstrap.txt\n\tmaxswaps\t<integer number>\t\t*1,000,000\n", number_of_steps);
            if(criterion == 0)
                {
                printf2("\tweight\t\tequal | comparisons\t\t");
                if(dweight == 1) printf2("comparisons");else printf2("splits");
                }
            if(criterion == 2)
                {
                printf2("\tweight\t\tequal | splits\t\t\t");
                if(splits_weight == 1) printf2("equal");if(splits_weight == 2) printf2("splits");
                }
            if(criterion == 3)
                {
                printf2("\tweight\t\tequal | taxa | quartets\t\t");
                if(quartet_normalising == 1) printf2("equal");if(quartet_normalising == 2) printf2("taxa");if(quartet_normalising == 3) printf2("quartets");
                }
			if(criterion == 5)
				{
				printf2("\n\tduplications\t<value>\t\t\t\t*1.0\n\tlosses\t\t<value>\t\t\t\t*1.0");
				printf2("\n\tnumspeciesrootings\t<value> | all\t\t*2\n\tnumgenerootings\t\t<value> | all\t\t*2\n");
				}
			if(criterion == 7)
				{
				printf2("\n\tmlbeta\t\t<float > 0>\t\t\t%.4f", ml_beta);
				printf2("\n\tmlscale\t\tpaper | lust | lnl\t\t%s", ml_scale==1?"lust":ml_scale==2?"lnl":"paper");
				printf2("\n\tmleta\t\t<float >= 0>\t\t\t%.4f  [experimental]", ml_eta);
				}
#ifdef _OPENMP
			if(criterion == 0 || criterion == 2 || criterion == 3)
				printf2("\n\tnthreads\t<integer number>\t\t*%-3d (parallel bootstrap replicates; default=all CPUs)", omp_get_num_procs());
			else if(criterion == 6 || criterion == 7)
				printf2("\n\tnthreads\t(single-threaded for RF/ML; parallelism not yet implemented)");
			if(hs_progress_interval == 0)
				printf2("\n\tprogress\t<integer seconds | 0>\t\t0 (report every improvement)");
			else
				printf2("\n\tprogress\t<integer seconds | 0>\t\t%d (parallel best-so-far update interval; 0=every improvement)", hs_progress_interval);
			if(hs_droprep > 0.0f)
				printf2("\n\tdroprep\t\t<float 0-1>\t\t\t%.4f (abandon rep if score >X%% above global best; 0=disabled)", hs_droprep);
			else
				printf2("\n\tdroprep\t\t<float 0-1>\t\t\t0 (disabled; abandon rep if score >X%% above global best; parallel only)");
#endif
			printf2("\n\tautoprunemono\t(set via exe) prune monophyletic same-species clades from multicopy trees at load time");
			printf2("\n\tconsensus\tstrict | majrule | minor | <proportion>\t*majrule");
			printf2("\n\tconsensusfile\t<filename>\t\t\tconsensus.ph\n");
            printf2("\n");
            }
        if(criterion == 1)
            printf2("\n\tnreps\t\t<integer number>\t\t*100\n\tswap\t\tnni | spr | tbr | all\t\t*tbr\n\taddseq\t\tsimple | closest | asis |\n\t\t\trandom | furthest\t\t*random\n\ttreefile\t<output treefile name>\t\tMRP.tree\n");
        }
    if(num == 6)
        {
        printf2("\nyaptp [options]  \n\n\n");
        printf2("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf2("\t===========================================================\n");
        if(criterion != 1)
            {
            printf2("\n\tmethod\t\tequiprobable | markovian\t*equiprobable");
            printf2("\n\tnreps\t\t<integer number>\t\t*100\n\thsreps\t\t<integer number>\t\t*1\n\tsample\t\t<integer number>\t\t*10,000\n\tsearch\t\tnni | spr | all\t\t\t");
            if(method == 1) printf2("nni");
            if(method == 2) printf2("spr");
            printf2("\n\tnsteps\t\t<integer number>\t\t%d\n\ttreefile\t<output treefile name>\t\tyaptp.ph\n\tmaxswaps\t<integer number>\t\t*1,000,000\n", number_of_steps);
            if(criterion == 0)
                {
                printf2("\tweight\t\tequal | comparisons\t\t");
                if(dweight == 1) printf2("comparisons");else printf2("splits");
                }
            if(criterion == 2)
                {
                printf2("\tweight\t\tequal | splits\t\t\t");
                if(splits_weight == 1) printf2("equal");if(splits_weight == 2) printf2("splits");
                }
            if(criterion == 3)
                {
                printf2("\tweight\t\tequal | taxa | quartets\t\t");
                if(quartet_normalising == 1) printf2("equal");if(quartet_normalising == 2) printf2("taxa");if(quartet_normalising == 3) printf2("quartets");
                }
	if(criterion == 5)
		{
		printf2("\n\tduplications\t<value>\t\t\t\t*1.0\n\tlosses\t\t<value>\t\t\t\t*1.0");
		printf2("\n\tnumspeciesrootings\t<value> | all\t\t*2\n\tnumgenerootings\t\t<value> | all\t\t*2\n");
		}

            printf2("\n");
            }
        else
            {
            printf2("\n\tnreps\t\t<integer number>\t\t*100\n");
            }
            
        }
    if(num == 7)
        {
        
        }
    if(num == 8)
        {
        printf2("\nset [options]  \n\n\n");
        printf2("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf2("\t===========================================================\n");
        printf2("\n\tcriterion\tdfit | sfit | qfit | mrp | avcon | rf | ml\t");
        printf2("\n\tmlbeta\t\t<float > 0>\t\t\t\t%.4f", ml_beta);
        printf2("\n\tmlscale\t\tpaper | lust | lnl\t\t\t%s", ml_scale == 1 ? "lust (Akanni et al. 2014)" : ml_scale == 2 ? "lnl" : "Steel & Rodrigo 2008");
        printf2("\n\tmleta\t\t<float >= 0>\t\t\t\t%.4f  [experimental]", ml_eta);
        if(criterion == 0) printf2("dfit");
        if(criterion == 1) printf2("mrp");
        if(criterion == 2) printf2("sfit");
        if(criterion == 3) printf2("qfit");
        if(criterion == 4) printf2("avcon");
        if(criterion == 6) printf2("rf");
        if(criterion == 7) printf2("ml");
        printf2("\n");
        printf2("\n\t  dfit  = Distance Fit (path-length distances; default)");
        printf2("\n\t  sfit  = Splits Fit (bipartition compatibility)");
        printf2("\n\t  qfit  = Quartet Fit (quartet topology compatibility)");
        printf2("\n\t  mrp   = Matrix Representation Parsimony (requires PAUP*)");
        printf2("\n\t  avcon = Average Consensus distances (requires PAUP*)");
        printf2("\n\t  recon = Duplication/Loss Reconciliation");
        printf2("\n\t  rf    = Robinson-Foulds distance (normalised, sum across gene trees)");
        printf2("\n\t  ml    = Maximum Likelihood supertree (Steel & Rodrigo 2008; exponential RF model)\n");
        printf2("\n\tseed\t\t<integer number>\t\t\t%d", seed);
/*        printf2("\n\n\t\t\tdfit = best Distance Fit\n\t\t\tsfit = maximum Splits Fit\n\t\t\tqfit = maximum Quartet Fit\n\t\t\tmrp = Matrix representation using parsimony\n");
  */      }
    if(num == 9)
        printf2("\nquit \n\tquit from the program and return to the operating system\n");
    if(num == 10)
        {
        printf2("\nexclude [options]  \nNOT IMPLEMENTED YET\n\n");
        printf2("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf2("\t===========================================================\n");

        printf2("\n\tall\t\tyes | no\t\t\tno\n\tlessthan\t<integer number>  \n\tgreaterthan\t<integer number>\n\tequalto\t\t<integer number>\n\tcontaining\t<taxanumbers>\n");
        }
    if(num == 11)
        {
        printf2("\ninclude [options]  \nNOT IMPLEMENTED YET\n\n");
        printf2("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf2("\t===========================================================\n");

        printf2("\n\tall\t\tyes | no\t\t\tno\n\tlessthan\t<integer number>  \n\tgreaterthan\t<integer number>\n\tequalto\t\t<integer number>\n\tcontaining\t<taxanumbers>\n");
        }
     if(num == 12)
        {
        printf2("\n!  \n\n\n");
        printf2("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf2("\t===========================================================\n");

        printf2("\n\tTemporarily returns the user to the operating system by running a C shell\n\tType 'exit' to return to clann");
        }
    if(num == 13)
        printf2("\nquit\n\tThis command quits Clann\n");
	
	if(num == 14)
		{
        printf2("\ngeneratetrees [options]  \n\n\n");
        printf2("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf2("\t===========================================================\n");

        printf2("\n\tmethod\t\tequiprobable | markovian\t*equiprobable\n\tntrees\t\tall | <integer number>\t\t*100\n\tnbins\t\t<integer number>\t\t*20\n\toutfile\t\t<output file name>\t\t*histogram.txt\n\tsourcedata\treal | randomised | ideal\t*real\n\tsavescores\tyes | no\t\t\t*no\n\tsupertree\tmemory | <supertree file name>\t*memory\n\tsavesourcetrees\tyes | no\t\t\t*no\n\n\tThe option 'supertree' is only used when idealised data are being created\n\n");
		
		
		}
     
	 if(num == 15)
		{
		printf2("\nconsensus [options]\n\n\n");
		printf2("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf2("\t===========================================================\n");
		
		printf2("\n\tmethod\tstrict | majrule | minor | <value>\t*minor\n\tguidetree\t<guide tree file name>\t\t*<none>\n\tfilename\t<output file name>\t\t*consensus.ph\n");
		}
	 
	 if(num == 16)
		{
		printf2("\nshowtrees\trange | size | namecontains | containstaxa | score \n\n");
		printf2("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf2("\t===========================================================\n");
		
		printf2("\n\trange\t\t<integer value> - <integer value> \t*all\n\tsize\t\tequalto <integer value>\n\t\t\tlessthan <integer value>\n\t\t\tgreaterthan <integer value>\t\t*none\n\tnamecontians\t<character string>\t\t\t*none\n\tcontainstaxa\t<character string>\t\t\t*none\n\tscore\t\t<min score> - <max score>\t\t*none\n\tsavetrees\tyes | no\t\t\t\t*no\n\tfilename\t<output file name>\t\t\t*showtrees.txt\n\tdisplay\t\tyes | no\t\t\t\t*yes\n");
		}

	 if(num == 17)
		{
		printf2("\texcludetrees\trange | size | namecontains | containstaxa | score \n\n");
		printf2("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf2("\t===========================================================\n");
		
		printf2("\n\tsinglecopy\tN/A\n\tmulticopy\tN/A\n\trange\t\t<integer value> - <integer value> \t*all\n\tsize\t\tequalto <integer value>\n\t\t\tlessthan <integer value>\n\t\t\tgreaterthan <integer value>\t\t*none\n\tnamecontains\t<character string>\t\t\t*none\n\tcontainstaxa\t<character string>\t\t\t*none\n\tscore\t\t<min score> - <max score>\t\t*none\n");
		}

	 if(num == 18)
		{
		printf2("\nincludetrees\trange | size | namecontains | containstaxa | score \n\n");
		printf2("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf2("\t===========================================================\n");
		
		printf2("\n\trange\t\t<integer value> - <integer value> \t*all\n\tsize\t\tequalto <integer value>\n\t\t\tlessthan <integer value>\n\t\t\tgreaterthan <integer value>\t\t*none\n\tnamecontians\t<character string>\t\t\t*none\n\tcontainstaxa\t<character string>\t\t\t*none\n\tscore\t\t<min score> - <max score>\t\t*none\n");
		}
	
	 if(num == 19)
		{
		printf2("\nrfdists [options]\t\n\n");
		printf2("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf2("\t===========================================================\n");
		printf2("\n\tfilename\t<output file name>\t\t*robinson_foulds.txt\n\toutput\t\tmatrix | vector\t\t*matrix\n\tmissing\t\tnone | 4point | ultrametric\t*none\n");
		}
	if(num == 20)
		{
		printf2("\ndeletetaxa\t <taxa name> <taxa name> etc...\n\n");
		printf2("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf2("\t===========================================================\n");
		}
	if(num == 21)
		{
		printf2("\nnj [options]\t\n\n");
		printf2("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf2("\t===========================================================\n");
		if(delimiter) printf2("\n\t(Note: delimiter mode is ON -- multicopy gene trees will be automatically\n\texcluded from the search; use 'maxnamelen=full' to disable)\n");
		printf2("\n\tmissing\t\t4point | ultrametric\t\t*4point\n\tsavetrees\t<file name>\t\t\t*NJtree.ph\n");
		}
	 if(num == 22)
		{
		printf2("\nsprdists [options]\t\n\n");
		printf2("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf2("\t===========================================================\n");
		printf2("\n\tsupertree\tcreate | memory | <file name>\t*create\n\tdorandomisation\tyes | no\t\t\t*yes\n\toutfile\t\t<file name>\t\t\t*SPRdistances.txt\n");
		}
	 if(num == 23)
		{
		printf2("\nmapunknowns \t\n\n");
		printf2("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf2("\t===========================================================\n");
		printf2("\n\t\n");
		}
	 if(num == 24)
		{
		printf2("\nreconstruct [options]\t\n\n");
		printf2("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf2("\t===========================================================\n");
		printf2("\n\tduplications\t<value>\t\t\t\t*1.0\n\tlosses\t\t<value>\t\t\t\t*1.0");
		printf2("\n\tshowrecon\tyes | no\t\t\t*no\n\tbasescore\t<value>\t\t\t\t*1.0\n\tprintfiles\tyes | no\t\t\t*yes\n\tspeciestree\tmemory | first | <file>\t\t*memory");
		}
     if(num == 25)
        {
        printf2("\nprunemonophylies [options]\t\n\n");
        printf2("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf2("\t===========================================================\n");
        printf2("\n\tfilename\t<output file name>\t\t*prunedtrees.txt");
        printf2("\n\tselection\trandom | length\t\t\t*random\n\n\tIf \"length\" is chosen, then name MUST have a number directly following the name of the species\n\t representing the sequence length. i.e.: \"Speces.length.XXXXXX\n\n" );
        }
    if(num == 26)
        {
        printf2("\ntips [options]\t\n\n");
        printf2("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf2("\t===========================================================\n");
        printf2("\n\tnumber\t<value between 1 and 10>\t\trandom");
        }
    if(num == 27)
        {
        printf2("\nrandomisetrees \t\n\n");
        printf2("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf2("\t===========================================================\n");
        printf2("\n\tThis command randomises all source trees in memory\n");
        }
    if(num == 28)
        {
        printf2("\nrandomprune \t\n\n");
        printf2("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf2("\t===========================================================\n");
        printf2("\n\tThis command randomly prunes the source trees in memory\n");
        }
   if(num == 29)
        {
        printf2("\nlog [options]\t\n\n");
        printf2("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf2("\t===========================================================\n");
        printf2("\n\tstatus\t\ton|off\t\t\t\t");
        if(print_log == TRUE) printf2("on");
        else printf2("off");
        printf2("\n\tfile\t\t<output file name>\t\t%s", logfile_name);
        }
	 if(num == 31)
		{
		printf2("\nmlscores\t[outfile=<file>] [scan=<n>] [scanmin=<f>] [scanmax=<f>]\n");
		printf2("\t\t[eta=auto] [escan=<n>] [etamax=<f>] [fixbeta]\n\n");
		printf2("  Estimates the Steel & Rodrigo (2008) ML beta parameter for the current\n");
		printf2("  supertree by closed-form MLE: beta = W / WD, where W is the sum of\n");
		printf2("  source-tree weights and WD is the weighted sum of (scaled) RF distances.\n");
		printf2("  Updates ml_beta so subsequent hs/boot runs use the estimated value.\n\n");
		printf2("  eta=auto [experimental]: jointly estimates the tree-size scaling\n");
		printf2("  exponent eta via 1-D grid search. Model:\n");
		printf2("  log L = W*log(beta) - eta*sum(log k_i) - beta*sum(w_i*d_i/k_i^eta)\n");
		printf2("  eta=0: Steel 2008 (default); eta=1: normalised by split count;\n");
		printf2("  eta>1: down-weights large trees. Updates ml_eta after estimation.\n\n");
		printf2("\tOptions\t\t\tSettings\t\t\tCurrent\n");
        printf2("\t===========================================================\n");
		printf2("\n\toutfile\t\t\t<filename>\t\t\t*none\n");
		printf2("\t\t\t\tWrite beta (and eta) log-likelihood profile\n");
		printf2("\tscan\t\t\t<integer>\t\t\t*100\n");
		printf2("\t\t\t\tNumber of points in beta profile\n");
		printf2("\tscanmin\t\t\t<float>\t\t\t\t*beta/100\n");
		printf2("\t\t\t\tLower bound of beta scan\n");
		printf2("\tscanmax\t\t\t<float>\t\t\t\t*beta*10\n");
		printf2("\t\t\t\tUpper bound of beta scan\n");
		printf2("\teta\t\t\tauto\t\t\t\t*off\n");
		printf2("\t\t\t\tEstimate optimal eta from data\n");
		printf2("\tescan\t\t\t<integer>\t\t\t*50\n");
		printf2("\t\t\t\tNumber of points in eta grid\n");
		printf2("\tetamax\t\t\t<float>\t\t\t\t*3.0\n");
		printf2("\t\t\t\tUpper bound of eta grid\n");
		printf2("\tfixbeta\t\t\t(no value)\t\t\t*off\n");
		printf2("\t\t\t\tHold beta fixed at current ml_beta; do not recompute.\n");
		printf2("\t\t\t\tWith eta=auto: grid-search eta with beta fixed.\n");
		printf2("\tsourcescores\t\t<filename>\t\t\t*none\n");
		printf2("\t\t\t\tWrite per-source-tree lnL contributions to a TSV file.\n");
		printf2("\t\t\t\tColumns: name, weight, lnL  (one row per active source tree).\n");
		}

	 if(num == 30)
		{
		printf2("\nsavetrees\trange | size | namecontains | containstaxa | score \n\n");
		printf2("\tOptions\t\tSettings\t\t\tCurrent\n");
        printf2("\t===========================================================\n");
		
		printf2("\n\trange\t\t<integer value> - <integer value> \t*all\n\tsize\t\tequalto <integer value>\n\t\t\tlessthan <integer value>\n\t\t\tgreaterthan <integer value>\t\t*none\n\tnamecontians\t<character string>\t\t\t*none\n\tcontainstaxa\t<character string>\t\t\t*none\n\tscore\t\t<min score> - <max score>\t\t*none\n\tfilename\t<output file name>\t\t\t*savedtrees.txt\n");
		}
	
	
		
     if(num != 0)
        {
        printf2("\n\n\t\t\t\t\t\t\t*Option is nonpersistent\n");
        printf2("\t===========================================================\n\n");    
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


        
/* --- autoprunemono helpers --- */

/* Recursively walk in-memory tree and recount presence_of_taxa[treenum][j] from surviving leaves.
 * Checks pos->tag: pruned leaves have tag==FALSE and are skipped (same logic as print_pruned_tree).
 * Uses pos->name which holds the taxa_names[] index set by find_taxon() at tree-build time. */
static void recount_from_tree(struct taxon *pos, int treenum)
	{
	while(pos != NULL)
		{
		if(pos->daughter != NULL)
			{
			/* internal node: always recurse into daughters */
			recount_from_tree(pos->daughter, treenum);
			}
		else
			{
			/* leaf node: only count if not pruned (tag != FALSE) */
			if(pos->tag != FALSE && pos->name >= 0 && pos->name < number_of_taxa)
				presence_of_taxa[treenum][pos->name]++;
			}
		pos = pos->next_sibling;
		}
	}

/* After loading trees, prune monophyletic same-species clades from multicopy trees.
 * Promoted (become single-copy) trees join the supertree pool.
 * Originals are saved in original_fundamentals[] for reconstruct. */
static void autoprunemono_apply(void)
	{
	int i, j, n_multicopy = 0, n_promoted = 0, n_still_multi = 0, n_pruned = 0;
	int numt = 0, clannID = 0, num_nodes = 0, leaf_count = 0;
	int *taxa_fate = NULL;
	char *temptree = malloc(TREE_LENGTH * sizeof(char));
	char *pruned_tree = NULL, *tmp = NULL;

	if(!temptree) { printf2("Error: out of memory for autoprunemono\n"); return; }

	/* Allocate original_fundamentals if not already done */
	if(original_fundamentals == NULL)
		{
		original_fundamentals = (char **)calloc(fundamental_assignments * FUNDAMENTAL_NUM, sizeof(char *));
		if(!original_fundamentals)
			{
			printf2("Error: out of memory for autoprunemono\n");
			free(temptree);
			return;
			}
		}

	pruned_tree = malloc(TREE_LENGTH * sizeof(char));
	if(!pruned_tree) { free(temptree); printf2("Error: out of memory for autoprunemono\n"); return; }
	tmp = malloc(TREE_LENGTH * sizeof(char));
	if(!tmp) { free(temptree); free(pruned_tree); printf2("Error: out of memory for autoprunemono\n"); return; }

	for(i = 0; i < Total_fund_trees; i++)
		{
		if(!sourcetreetag[i]) continue;

		/* Check if multicopy */
		int multicopy = FALSE;
		for(j = 0; j < number_of_taxa; j++)
			if(presence_of_taxa[i][j] > 1) { multicopy = TRUE; break; }
		if(!multicopy) continue;
		n_multicopy++;

		/* Dismantle any tree already in memory */
		if(tree_top != NULL) { dismantle_tree(tree_top); tree_top = NULL; }
		temp_top = NULL;

		/* Build tree from stored representation using full names */
		temptree[0] = '\0';
		strcpy(temptree, fundamentals[i]);
		returntree_fullnames(temptree, i);
		basic_tree_build(1, temptree, tree_top, TRUE);
		tree_top = temp_top;
		temp_top = NULL;
		reset_tree(tree_top);

		/* Number all taxa (sets tag2 on leaves) */
		numt = 0; clannID = 0;
		numt = number_tree2(tree_top, numt);

		/* Mark species-specific clades */
		taxa_fate = malloc(numt * sizeof(int));
		if(!taxa_fate)
			{
			printf2("Error: out of memory for autoprunemono\n");
			free(temptree); free(pruned_tree); free(tmp);
			if(tree_top != NULL) { dismantle_tree(tree_top); tree_top = NULL; }
			return;
			}
		for(j = 0; j < numt; j++) taxa_fate[j] = 0;
		clannID = identify_species_specific_clades(tree_top, numt, taxa_fate, clannID);
		free(taxa_fate); taxa_fate = NULL;

		/* Shrink and build pruned Newick (numeric IDs = fundamentals[] format) */
		shrink_tree(tree_top);
		pruned_tree[0] = '\0';
		leaf_count = print_pruned_tree(tree_top, 0, pruned_tree, FALSE, i);

		/* Count commas to determine remaining taxon count */
		num_nodes = 0;
		for(j = 0; j < (int)strlen(pruned_tree); j++)
			if(pruned_tree[j] == ',') num_nodes++;

		/* Wrap with outer parens if multi-node (same as prune_monophylies) */
		if(leaf_count > 1)
			{
			tmp[0] = '\0';
			strcpy(tmp, "(");
			strcat(tmp, pruned_tree);
			strcat(tmp, ")");
			strcpy(pruned_tree, tmp);
			}
		strcat(pruned_tree, ";");

		if(num_nodes > 2)  /* tree has 4+ taxa: worth keeping */
			{
			n_pruned++;

			/* Save original fundamentals[i] */
			original_fundamentals[i] = malloc((strlen(fundamentals[i]) + 2) * sizeof(char));
			if(!original_fundamentals[i])
				{
				printf2("Error: out of memory for autoprunemono\n");
				free(temptree); free(pruned_tree); free(tmp);
				if(tree_top != NULL) { dismantle_tree(tree_top); tree_top = NULL; }
				return;
				}
			strcpy(original_fundamentals[i], fundamentals[i]);

			/* Replace fundamentals[i] with pruned version */
			free(fundamentals[i]);
			fundamentals[i] = malloc((strlen(pruned_tree) + 100) * sizeof(char));
			if(!fundamentals[i]) memory_error(220);
			strcpy(fundamentals[i], pruned_tree);

			/* Recount presence_of_taxa[i][*] from in-memory pruned tree */
			for(j = 0; j < number_of_taxa; j++)
				presence_of_taxa[i][j] = 0;
			recount_from_tree(tree_top, i);

			/* Check if now single-copy */
			int still_multicopy = FALSE;
			for(j = 0; j < number_of_taxa; j++)
				if(presence_of_taxa[i][j] > 1) { still_multicopy = TRUE; break; }
			if(still_multicopy)
				n_still_multi++;
			else
				n_promoted++;
			}

		/* Free in-memory tree for next iteration */
		if(tree_top != NULL) { dismantle_tree(tree_top); tree_top = NULL; }
		temp_top = NULL;
		}

	free(temptree);
	free(pruned_tree);
	free(tmp);

	if(n_multicopy > 0)
		{
		printf2("\nAutoprunemono: pruned monophyletic same-species clades in %d multicopy trees.\n", n_pruned);
		printf2("  Promoted to single-copy pool:   %d trees\n", n_promoted);
		printf2("  Still multicopy after pruning:  %d trees (retained for reconstruct)\n", n_still_multi);
		if(n_pruned > 0)
			printf2("  Original (unpruned) trees stored for reconstruct.\n");
		}
	}

/* -----------------------------------------------------------------------
 * is_monophyletic_subtree: recursive helper.
 * Returns number of marked leaves in the subtree rooted at pos.
 * Sets *mono=FALSE if the marked leaves are not monophyletic.
 * total_marked: total number of marked leaves in the whole tree.
 * ----------------------------------------------------------------------- */
static int is_monophyletic_subtree(struct taxon *pos, const int *marked,
                                   int total_marked, int *mono)
    {
    if(pos == NULL) return 0;
    if(pos->daughter == NULL)   /* leaf */
        return (marked[pos->name] ? 1 : 0);

    int count = 0;
    struct taxon *child = pos->daughter;
    while(child != NULL)
        {
        count += is_monophyletic_subtree(child, marked, total_marked, mono);
        child = child->next_sibling;
        }
    /* If this subtree contains some but not all marked leaves, not monophyletic */
    if(count > 0 && count < total_marked) *mono = FALSE;
    return count;
    }

/* Returns TRUE if the taxa marked in marked[] form a monophyletic group
 * in the tree rooted at root.  total_taxa = number_of_taxa. */
static int is_monophyletic(struct taxon *root, const int *marked, int total_taxa)
    {
    (void)total_taxa;
    /* count total marked */
    int total_marked = 0, t;
    for(t = 0; t < total_taxa; t++) if(marked[t]) total_marked++;
    if(total_marked < 2) return TRUE;
    int mono = TRUE;
    is_monophyletic_subtree(root, marked, total_marked, &mono);
    return mono;
    }

/* -----------------------------------------------------------------------
 * count_clan_and_total_in_subtree: recursive helper used by
 * count_violated_splits.  Returns the number of clan members in the
 * subtree rooted at pos, and stores the total leaf count (all taxa,
 * clan or not) in *total_out.
 * ----------------------------------------------------------------------- */
static int count_clan_and_total_in_subtree(struct taxon *pos,
                                            const int *clan_ids, int clan_size,
                                            int *total_out)
    {
    if(pos == NULL) { *total_out = 0; return 0; }
    if(pos->daughter == NULL)   /* leaf */
        {
        *total_out = 1;
        int j;
        for(j = 0; j < clan_size; j++)
            if(clan_ids[j] == pos->name) return 1;
        return 0;
        }
    int clan_cnt = 0, tot = 0;
    struct taxon *ch = pos->daughter;
    while(ch != NULL)
        {
        int sub_tot = 0;
        clan_cnt += count_clan_and_total_in_subtree(ch, clan_ids, clan_size, &sub_tot);
        tot += sub_tot;
        ch = ch->next_sibling;
        }
    *total_out = tot;
    return clan_cnt;
    }

/* -----------------------------------------------------------------------
 * count_violated_splits: recursively count internal nodes (splits) in the
 * tree rooted at pos and determine how many are violated by at least one
 * clan.
 *
 * A split {L | R} violates clan C only when non-clan taxa appear on BOTH
 * sides alongside clan members — i.e. the split genuinely crosses the clan
 * boundary in both directions.  Internal splits within a monophyletic clan
 * (where one side consists entirely of clan members) are NOT counted as
 * violations.
 *
 * Formal condition: split {L|R} violates clan C iff all four hold:
 *   clan_in  >= 1  (clan member(s) in subtree L)
 *   clan_out >= 1  (clan member(s) in R)
 *   non_clan_in  >= 1  (non-clan taxon in L)
 *   non_clan_out >= 1  (non-clan taxon in R)
 *
 * n_present[c] = number of clan c members present in this source tree.
 * n_taxa_tree  = total taxa present in this source tree.
 * violated and total are incremented in-place.
 * ----------------------------------------------------------------------- */
static void count_violated_splits(struct taxon *pos,
                                   int **clan_ids, int *clan_sizes,
                                   int *n_present, int nclans,
                                   int n_taxa_tree,
                                   int *violated, int *total)
    {
    if(pos == NULL || pos->daughter == NULL) return;  /* skip leaves */
    (*total)++;
    int c, is_viol = 0;
    for(c = 0; c < nclans && !is_viol; c++)
        {
        if(n_present[c] < 2) continue;
        int total_in = 0;
        int clan_in  = count_clan_and_total_in_subtree(pos, clan_ids[c],
                                                        clan_sizes[c], &total_in);
        int clan_out     = n_present[c] - clan_in;
        int non_clan_in  = total_in - clan_in;
        int non_clan_out = (n_taxa_tree - total_in) - clan_out;
        if(clan_in >= 1 && clan_out >= 1 && non_clan_in >= 1 && non_clan_out >= 1)
            is_viol = 1;
        }
    if(is_viol) (*violated)++;
    struct taxon *ch = pos->daughter;
    while(ch != NULL)
        {
        count_violated_splits(ch, clan_ids, clan_sizes, n_present, nclans,
                               n_taxa_tree, violated, total);
        ch = ch->next_sibling;
        }
    }

/* -----------------------------------------------------------------------
 * compute_autoweights_clan: assign per-tree weights based on compatibility
 * with a set of user-defined clans (irrefutable monophyletic groups).
 *
 * For each source tree i, w_i = (number of clans compatible with T_i) /
 * (number of clans testable in T_i), where a clan C is testable in T_i if
 * at least 2 members of C are present.  Trees with no testable clans receive
 * weight 1.0.  Weights are normalised so the mean across tagged trees is 1.0.
 *
 * Clan file format: one clan per line, taxon names separated by whitespace
 * or commas.  Lines beginning with '#' are comments.  Taxon names must match
 * those in the source trees (after any delimiter processing).
 * ----------------------------------------------------------------------- */
/* mode 0 = clan-compatibility weight (compatible/testable clans)
 * mode 1 = split-violation weight (1 - violated_splits/total_splits)
 * Both modes always compute and report the split-violation metric; the
 * chosen mode determines which value is stored in tree_weights[]. */
static void compute_autoweights_clan(const char *clanfile, int mode)
    {
    int i, j, c, t;
    FILE *fp = fopen(clanfile, "r");
    if(!fp)
        { printf2("autoweight=clan: cannot open clan file '%s'\n", clanfile); return; }

    /* ---- read clans from file ---- */
    int    max_clans = 256, nclans = 0;
    int  **clan_ids   = malloc(max_clans * sizeof(int *));
    int   *clan_sizes = malloc(max_clans * sizeof(int));
    if(!clan_ids || !clan_sizes)
        { printf2("autoweight=clan: out of memory\n"); fclose(fp); return; }

    char line[100000];
    while(fgets(line, sizeof(line), fp))
        {
        /* strip newline and skip comments / blank lines */
        int len = (int)strlen(line);
        while(len > 0 && (line[len-1] == '\n' || line[len-1] == '\r')) line[--len] = '\0';
        if(line[0] == '#' || len == 0) continue;

        /* tokenise on whitespace and commas */
        int  clan_buf[10000], clan_n = 0;
        char *tok = strtok(line, " \t,");
        while(tok != NULL)
            {
            /* look up taxon name in taxa_names[] */
            int found = -1;
            for(t = 0; t < number_of_taxa; t++)
                if(taxa_names[t] != NULL && strcmp(taxa_names[t], tok) == 0)
                    { found = t; break; }
            if(found >= 0 && clan_n < 10000)
                clan_buf[clan_n++] = found;
            tok = strtok(NULL, " \t,");
            }

        if(clan_n < 2) continue;  /* need at least 2 members to be testable */

        if(nclans >= max_clans)
            {
            max_clans *= 2;
            clan_ids   = realloc(clan_ids,   max_clans * sizeof(int *));
            clan_sizes = realloc(clan_sizes, max_clans * sizeof(int));
            if(!clan_ids || !clan_sizes)
                { printf2("autoweight=clan: out of memory\n"); fclose(fp); return; }
            }
        clan_ids[nclans]   = malloc(clan_n * sizeof(int));
        if(!clan_ids[nclans]) { printf2("autoweight=clan: out of memory\n"); fclose(fp); return; }
        for(t = 0; t < clan_n; t++) clan_ids[nclans][t] = clan_buf[t];
        clan_sizes[nclans] = clan_n;
        nclans++;
        }
    fclose(fp);

    if(nclans == 0)
        { printf2("autoweight=clan: no valid clans found in '%s'\n", clanfile); goto cleanup; }

    printf2("  Clan file: %s  (%d clans loaded)\n", clanfile, nclans);

    /* ---- score each source tree ---- */
    /* For each tree i and clan c:
     *   testable if >= 2 clan members present in T_i
     *   compatible if all clan members present in T_i form a clade —
     *   i.e. none of the non-clan taxa present in T_i appear *between*
     *   clan members.  We test this via presence_of_taxa: a clan is
     *   incompatible if any taxon NOT in the clan is present in T_i and
     *   the clan members do not form a contiguous subtree.
     *
     *   Simple testable compatibility check using only presence_of_taxa:
     *   build the tree in memory and check actual monophyly.
     */
    char *temptree = malloc(TREE_LENGTH * sizeof(char));
    if(!temptree) { printf2("autoweight=clan: out of memory\n"); goto cleanup; }

    /* Save current tree_top */
    struct taxon *saved_top = tree_top;
    tree_top = NULL;

    double weight_sum = 0.0;
    int    n_tagged   = 0;

    /* per-tree split scores for summary (allocated once, reused) */
    double *split_scores = malloc(Total_fund_trees * sizeof(double));
    if(!split_scores) { printf2("autoweight=clan: out of memory\n"); free(temptree); tree_top = saved_top; goto cleanup; }
    for(i = 0; i < Total_fund_trees; i++) split_scores[i] = 1.0;

    /* pre-compute n_present[c] for all clans per tree (reused array) */
    int *np = malloc(nclans * sizeof(int));
    if(!np) { printf2("autoweight=clan: out of memory\n"); free(split_scores); free(temptree); tree_top = saved_top; goto cleanup; }

    for(i = 0; i < Total_fund_trees; i++)
        {
        if(!sourcetreetag[i]) { continue; }

        /* Build source tree into tree_top */
        temptree[0] = '\0';
        strcpy(temptree, fundamentals[i]);
        returntree_fullnames(temptree, i);
        basic_tree_build(1, temptree, tree_top, TRUE);
        tree_top = temp_top;
        temp_top = NULL;
        reset_tree(tree_top);

        /* pre-compute clan member counts for this tree */
        for(c = 0; c < nclans; c++)
            {
            np[c] = 0;
            for(j = 0; j < clan_sizes[c]; j++)
                if(presence_of_taxa[i][clan_ids[c][j]] > 0) np[c]++;
            }

        int compatible = 0, testable = 0;
        for(c = 0; c < nclans; c++)
            {
            if(np[c] < 2) continue;  /* not testable */
            testable++;

            /* count non-clan taxa present */
            int n_other = 0;
            for(t = 0; t < number_of_taxa; t++)
                {
                if(presence_of_taxa[i][t] == 0) continue;
                int in_clan = FALSE;
                for(j = 0; j < clan_sizes[c]; j++)
                    if(clan_ids[c][j] == t) { in_clan = TRUE; break; }
                if(!in_clan) n_other++;
                }

            /* if no other taxa, clan is trivially compatible */
            if(n_other == 0) { compatible++; continue; }

            /* mark clan members and check monophyly */
            int *clan_mark = calloc(number_of_taxa, sizeof(int));
            if(!clan_mark) continue;
            for(j = 0; j < clan_sizes[c]; j++)
                clan_mark[clan_ids[c][j]] = 1;

            int mono = is_monophyletic(tree_top, clan_mark, number_of_taxa);
            if(mono) compatible++;

            free(clan_mark);
            }

        /* compute split-violation metric for this tree */
        int n_taxa_tree = 0;
        for(t = 0; t < number_of_taxa; t++)
            if(presence_of_taxa[i][t] > 0) n_taxa_tree++;
        int vsplit = 0, tsplit = 0;
        count_violated_splits(tree_top, clan_ids, clan_sizes, np, nclans,
                               n_taxa_tree, &vsplit, &tsplit);
        split_scores[i] = (tsplit > 0) ? 1.0 - (double)vsplit / (double)tsplit : 1.0;

        if(tree_top != NULL) { dismantle_tree(tree_top); tree_top = NULL; }

        double clan_w  = (testable > 0) ? (double)compatible / (double)testable : 1.0;
        double split_w = split_scores[i];
        tree_weights[i] = (float)(mode == 1 ? split_w : clan_w);
        weight_sum += tree_weights[i];
        n_tagged++;
        }

    free(np);
    free(temptree);
    tree_top = saved_top;

    if(mode == 1)
        printf2("  Split-violation weights assigned to %d trees (range 0–1: fraction of splits not violated by any clan)\n", n_tagged);
    else
        printf2("  Clan weights assigned to %d trees (range 0–1: fraction of clans supported)\n", n_tagged);

    /* print weight range and split score range */
    {
    float wmin = 1e30f, wmax = -1e30f;
    double smin = 1e30, smax = -1e30;
    for(i = 0; i < Total_fund_trees; i++)
        if(sourcetreetag[i])
            {
            if(tree_weights[i] < wmin) wmin = tree_weights[i];
            if(tree_weights[i] > wmax) wmax = tree_weights[i];
            if(split_scores[i] < smin) smin = split_scores[i];
            if(split_scores[i] > smax) smax = split_scores[i];
            }
    printf2("  Weight range: %.4f – %.4f\n", wmin, wmax);
    if(mode == 0)
        printf2("  Split score range: %.4f – %.4f  (fraction of splits not violated by any clan)\n", (float)smin, (float)smax);
    }

    free(split_scores);

cleanup:
    for(c = 0; c < nclans; c++) free(clan_ids[c]);
    free(clan_ids);
    free(clan_sizes);
    }

void execute_command(char *commandline, int do_all)
    {
    int i = 0, j=0, k=0, printfundscores = FALSE, error = FALSE, do_autoprunemono = FALSE;
    int do_clan_weights = 0;  /* 0=off, 1=clan-compatibility, 2=split-violation */
    char c = '\0', temp[NAME_LENGTH], filename[10000], *newbietree, string_num[1000];
    char clanfile_name[10000];
    int newbietree_alloc = 0;  /* tracks allocated size of newbietree for overflow detection */
	float num = 0;
    clanfile_name[0] = '\0';

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
				else if(strcmp(parsed_command[i+1], "full") == 0)
					{
					max_name_length = NAME_LENGTH;
					delimiter = FALSE;
					printf2("Full taxon names in use (delimiter mode disabled).\n");
					}
				else
					{
					error = TRUE;
					printf2("Error: value %s not valid for maxnamelen\n\n", parsed_command[i+1]);
					}
				}
            }
        if(strcmp(parsed_command[i], "delimiter_char") == 0)
            {
            max_name_length = toint(parsed_command[i+1]);
            if(max_name_length == 0)
                {
                if(parsed_command[i+1][0] != '\0')
					{
					delimiter_char=parsed_command[i+1][0];
					}					
				else
					{
					error = TRUE;
					printf2("Error: A character must be provided as a delimiter\n");
					}
				}
            }
        if(strcmp(parsed_command[i], "autoprunemono") == 0)
            {
            if(strcmp(parsed_command[i+1], "yes") == 0)
                do_autoprunemono = TRUE;
            else if(strcmp(parsed_command[i+1], "no") == 0)
                do_autoprunemono = FALSE;
            else
                {
                printf2("Error: value '%s' not valid for autoprunemono (use yes or no)\n", parsed_command[i+1]);
                error = TRUE;
                }
            }
        if(strcmp(parsed_command[i], "autoweight") == 0)
            {
            if(strcmp(parsed_command[i+1], "clan") == 0)
                do_clan_weights = 1;
            else if(strcmp(parsed_command[i+1], "splitviol") == 0)
                do_clan_weights = 2;
            else
                printf2("Warning: autoweight value '%s' not recognised (use 'clan' or 'splitviol')\n", parsed_command[i+1]);
            }
        if(strcmp(parsed_command[i], "clanfile") == 0)
            {
            strncpy(clanfile_name, parsed_command[i+1], 9999);
            clanfile_name[9999] = '\0';
            }
        }

    if(delimiter == TRUE) printf2("\nDelimiter mode active: species names extracted before '%c' (use maxnamelen=full to disable)\n", delimiter_char);

    /********** open the input file ************/
    filename[0] = '\0';
 
    if(commandline == NULL)
        {
        strcpy(filename, parsed_command[1]);
        }
    else
        strcpy(filename, commandline);
    
    if((infile = fopen(filename, "r")) == NULL)		/* check to see if the file is there */
        {								/* Open the source tree file */
        printf2("Cannot open file %s\n", filename);
        error = TRUE;
        }
	strcpy(inputfilename, filename);
    if(!error)
        {
		largest_tree = 0;
		if(yaptp_results != NULL)
			{
			free(yaptp_results);
			yaptp_results = NULL;
			}
		num_excluded_trees = 0;
		num_excluded_taxa = 0;
		trees_in_memory = 0;
            /************************ Assign the dynamic arrays *************************/
        newbietree_alloc = TREE_LENGTH;
        newbietree = malloc(newbietree_alloc*sizeof(char));
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
        /* Reset autoprunemono state when new trees are loaded */
        if(original_fundamentals != NULL)
            {
            for(i=0; i<fundamental_assignments*FUNDAMENTAL_NUM; i++)
                {
                if(original_fundamentals[i] != NULL)
                    {
                    free(original_fundamentals[i]);
                    original_fundamentals[i] = NULL;
                    }
                }
            free(original_fundamentals);
            original_fundamentals = NULL;
            }
        autoprunemono_active = 0;

        
		
		
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
            printf2("\nReading Nexus format source tree file \n");


			if(nexusparser(infile) == TRUE)
				{
				printf2("error reading nexus file\n");
				}


            
            
            
            }  /* End reading NEXUS format tree */
        else
            {
			printf2("\nReading Newhampshire (Phylip) format source trees\n");
            
			if(c != '(')
                {
                printf2("Error: Treefile not in correct format\n");
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
							/* Grow buffer if approaching the limit (handles trees longer than TREE_LENGTH) */
							if(i >= newbietree_alloc - 3)
								{
								newbietree_alloc *= 2;
								newbietree = realloc(newbietree, newbietree_alloc * sizeof(char));
								if(newbietree == NULL) memory_error(111);
								}
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
		/* Initialise per-taxon weights for bipartition tree hashing (splitmix64 finaliser) */
		taxon_hash_vals = realloc(taxon_hash_vals, number_of_taxa * sizeof(uint64_t));
		{
		uint64_t thv_x;
		for(i=0; i<number_of_taxa; i++)
			{
			thv_x = (uint64_t)(i + 1);
			thv_x += 0x9e3779b97f4a7c15ULL;
			thv_x = (thv_x ^ (thv_x >> 30)) * 0xbf58476d1ce4e5b9ULL;
			thv_x = (thv_x ^ (thv_x >> 27)) * 0x94d049bb133111ebULL;
			taxon_hash_vals[i] = thv_x ^ (thv_x >> 31);
			}
		}
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

		/* Apply autoprunemono if requested */
		if(do_autoprunemono && Total_fund_trees > 0)
			{
			autoprunemono_active = 1;
			autoprunemono_apply();
			}

        /* Apply clan-based autoweights if requested */
        if(do_clan_weights > 0 && Total_fund_trees > 0)
            {
            if(clanfile_name[0] == '\0')
                printf2("autoweight=clan/splitviol: please specify clanfile=<filename>\n");
            else
                compute_autoweights_clan(clanfile_name, do_clan_weights - 1);
            }

		free(newbietree);

        }
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
    if(partition_number != NULL) free(partition_number);
    if(from_tree != NULL) free(from_tree);

    if(sourcetree_scores != NULL) free(sourcetree_scores);
	if(presenceof_SPRtaxa) free(presenceof_SPRtaxa);
	
	if(sourcetreetag != NULL) free(sourcetreetag);
    /* Flush all stdio buffers before _exit so no output is lost.
     * We use _exit() rather than exit() to bypass libgomp's atexit cleanup
     * handler, which calls pthread_join on idle worker threads parked in
     * Mach semaphore_wait() (an uninterruptible kernel wait on macOS).
     * If the wake signal races past the thread, the join blocks indefinitely,
     * leaving the process in the UE (uninterruptible+exiting) state.
     * _exit() lets the OS tear down all threads atomically on process exit. */
    fflush(NULL);
    _exit(error ? 1 : 0);

	if(logfile!= NULL)
		fclose(logfile);
    }


void controlc1(int signal)
	{
	char *c = NULL;
	
	printf2("\n\nCompleted %d random samples before CTRL-C was caught\nDo you want to stop the random sampling and start the heuristuc searches now? (Y/N): \n", GC);
	c = malloc(10000*sizeof(char));
	xgets(c);
	while(strcmp(c, "Y") != 0 && strcmp(c, "y") != 0 && strcmp(c, "N") != 0 && strcmp(c, "n") != 0)
		{
		printf2("Error '%s' is not a valid response, please respond Y or N\n", c);
		xgets(c);
		}
	
	if(strcmp(c, "Y") == 0 || strcmp(c, "y") == 0)
		{
		GC = 100000000;
		printf2("Starting Heuristic searches\n");
		}

	if(strcmp(c, "N") == 0 || strcmp(c, "n") == 0)
			printf2("Continuing random samples.....\n");
	
	free(c);

	}

void controlc2(int sig)
	{
	(void)sig;
	/* async-signal-safe: only write() and _exit() are used here.
	 * printf/fflush/signal are NOT async-signal-safe and can deadlock
	 * if the signal interrupts a thread holding the stdio mutex. */
	static volatile sig_atomic_t ctrl_c_count = 0;
	ctrl_c_count++;
	if(ctrl_c_count == 1)
		{
		user_break = TRUE;
		static const char msg1[] =
			"\nSearch interrupted \xe2\x80\x94 threads will stop after their current swap; "
			"best tree will be written to output.\n"
			"Press Ctrl+C again to force-quit immediately (output will NOT be saved).\n";
		write(STDOUT_FILENO, msg1, sizeof(msg1) - 1);
		/* Handler persists on macOS/Linux (BSD signal semantics); no re-arm needed. */
		}
	else
		{
		static const char msg2[] = "\nForce quit \xe2\x80\x94 output not saved.\n";
		write(STDOUT_FILENO, msg2, sizeof(msg2) - 1);
		_exit(1);
		}
	}

void controlc3(int signal)
	{
	char *c = NULL;
	
	c = malloc(10000*sizeof(char));
	
	printf2("\n\nDo you want to stop the generation of trees now? (Y/N): ");
	xgets(c);
	while(strcmp(c, "Y") != 0 && strcmp(c, "y") != 0 && strcmp(c, "N") != 0 && strcmp(c, "n") != 0)
		{
		printf2("Error '%s' is not a valid response, please respond Y or N\n", c);
		xgets(c);
		}
	
	if(strcmp(c, "Y") == 0 || strcmp(c, "y") == 0)
		{
		user_break = TRUE;
		printf2("Stopping the generation of treesn\n");
		}


	if(strcmp(c, "N") == 0 || strcmp(c, "n") == 0)
			printf2("Continuing the generation of trees.....\n");
	
	free(c);
	}

void controlc4(int signal)
	{
	char *c = NULL;
	
	c = malloc(10000*sizeof(char));
	
	printf2("\n\nDo you really wish to quit? (Y/N): ");
	xgets(c);
	while(strcmp(c, "Y") != 0 && strcmp(c, "y") != 0 && strcmp(c, "N") != 0 && strcmp(c, "n") != 0)
		{
		printf2("Error '%s' is not a valid response, please respond Y or N\n", c);
		xgets(c);
		}
	
	if(strcmp(c, "Y") == 0 || strcmp(c, "y") == 0)
		{
		user_break = TRUE;
		printf2("Bye\n");
		clean_exit(0);
		}
	
	free(c);
	}

void controlc5(int signal)
	{
	char *c = NULL;
	
	c = malloc(10000*sizeof(char));
	
	printf2("\n\nDo you want to stop the exhaustive searches now? (Y/N): ");
	xgets(c);
	while(strcmp(c, "Y") != 0 && strcmp(c, "y") != 0 && strcmp(c, "N") != 0 && strcmp(c, "n") != 0)
		{
		printf2("Error '%s' is not a valid response, please respond Y or N\n", c);
		xgets(c);
		}
	if(strcmp(c, "Y") == 0 || strcmp(c, "y") == 0)
		{
		user_break = TRUE;
		printf2("Stopping exhaustive search. The resulting tree may be sub-optimal\n\n");
		}


	if(strcmp(c, "N") == 0 || strcmp(c, "n") == 0)
			printf2("Continuing heuristic searches.....\n");

	free(c);
	}


void set_parameters(void)
    {
    int i=0, j=0, isdigit=TRUE, this=FALSE;
    for(i=0; i<num_commands; i++)
        {
        if(strcmp(parsed_command[i], "criterion") == 0)
            {
            if(strcmp(parsed_command[i+1], "mrp") == 0 || strcmp(parsed_command[i+1], "MRP") == 0)
                {
                if(criterion != 1) printf2("Scoring criterion set to matrix representation using parsimony (MRP)\n");
                else printf2("Scoring criterion is already set to MRP\n");
                criterion = 1;
                }
            else
                {
                if(strcmp(parsed_command[i+1], "dfit") == 0 || strcmp(parsed_command[i+1], "DFIT") == 0)
                    {
                    if(criterion != 0) printf2("Scoring criterion set to best distance fit (DFIT)\n");
                    else printf2("Scoring criterion is already set to DFIT\n");
                    criterion = 0;
                    }
                else
                    {
                    if(strcmp(parsed_command[i+1], "sfit") == 0 || strcmp(parsed_command[i+1], "SFIT") == 0)
                        {
                        if(criterion != 2) printf2("Scoring criterion set to maximum split fit (SFIT)\n");
                        else printf2("Scoring criterion is already set to SFIT\n");
                        criterion = 2;
                        }
                    else
                        {
                        if(strcmp(parsed_command[i+1], "qfit") == 0 || strcmp(parsed_command[i+1], "QFIT") == 0)
                            {
                            if(criterion != 3) printf2("Scoring criterion set to maximum quartet fit (QFIT)\n");
                            else printf2("Scoring criterion is already set to QFIT\n");
                            criterion = 3;
                            }
                        else
                            {
                            if(strcmp(parsed_command[i+1], "avcon") == 0 || strcmp(parsed_command[i+1], "AVCON") == 0)
                                {
								if(criterion != 4) printf2("Scoring criterion set to average consensus (AVCON)\n");
								else printf2("Scoring criterion is already set to AVCON\n");
								criterion = 4;
                                }
                            else
								{
								if(strcmp(parsed_command[i+1], "recon") == 0 || strcmp(parsed_command[i+1], "RECON") == 0)
									{
									if(criterion != 5) printf2("Scoring criterion set to duplication and loss reconstruction (RECON)\n");
									else printf2("Scoring criterion is already set to RECON\n");
									criterion = 5;
									}
								else if(strcmp(parsed_command[i+1], "rf") == 0 || strcmp(parsed_command[i+1], "RF") == 0)
									{
									if(criterion != 6) printf2("Scoring criterion set to Robinson-Foulds distance (RF)\n");
									else printf2("Scoring criterion is already set to RF\n");
									criterion = 6;
									}
								else if(strcmp(parsed_command[i+1], "ml") == 0 || strcmp(parsed_command[i+1], "ML") == 0
									     || strcmp(parsed_command[i+1], "lust") == 0)
									{
									if(criterion != 7) printf2("Scoring criterion set to ML supertree likelihood (Steel & Rodrigo 2008)\n");
									else printf2("Scoring criterion is already set to ML\n");
									criterion = 7;
									}
								else
									printf2("Error: %s not known as criterion tpye\n", parsed_command[i+1]);
								}
							}
                        }
                    }
                }
            }
        if(strcmp(parsed_command[i], "seed") == 0)
        	{
        	isdigit=TRUE;
        	for(j=0; j<strlen(parsed_command[i+1]); j++)
        		{
        		this=FALSE;
        		if( parsed_command[i+1][j] == '0' || parsed_command[i+1][j] == '1' || parsed_command[i+1][j] == '2' || parsed_command[i+1][j] == '3' || parsed_command[i+1][j] == '4' || parsed_command[i+1][j] == '5' || parsed_command[i+1][j] == '6' || parsed_command[i+1][j] == '7' || parsed_command[i+1][j] == '8' || parsed_command[i+1][j] == '9' ) this=TRUE;
        		if(this == FALSE) isdigit=FALSE;
        		}
      		if(isdigit==FALSE)
  				{
      			printf2("ERROR: the value you entered (%s) is not an integer\n", parsed_command[i+1] );
  				}
  			else
  				{
  				seed=atoi(parsed_command[i+1]);
  				srand((unsigned) (seed));
  				printf2("The seed value for the random number generator has been set to %d\n", seed);
  				}
        	}
        if(strcmp(parsed_command[i], "mlbeta") == 0)
            {
            ml_beta = atof(parsed_command[i+1]);
            if(ml_beta <= 0)
                {
                printf2("Error: mlbeta must be > 0\n");
                ml_beta = 1.0f;
                }
            else
                printf2("ML beta parameter set to %.4f\n", ml_beta);
            }
        if(strcmp(parsed_command[i], "mlscale") == 0)
            {
            if(strcmp(parsed_command[i+1], "lnl") == 0 || strcmp(parsed_command[i+1], "LNL") == 0)
                { ml_scale = 2; printf2("ML scale set to lnl\n"); }
            else if(strcmp(parsed_command[i+1], "lust") == 0 || strcmp(parsed_command[i+1], "LUST") == 0)
                { ml_scale = 1; printf2("ML scale set to lust\n"); }
            else if(strcmp(parsed_command[i+1], "paper") == 0 || strcmp(parsed_command[i+1], "PAPER") == 0)
                { ml_scale = 0; printf2("ML scale set to paper\n"); }
            else
                printf2("Error: mlscale must be 'paper', 'lust', or 'lnl'\n");
            }
        if(strcmp(parsed_command[i], "mleta") == 0)
            {
            double a = atof(parsed_command[i+1]);
            if(a < 0.0)
                printf2("Error: mleta must be >= 0 (0=Steel 2008, 1=normalised, >1=downweight large trees)\n");
            else
                { ml_eta = a; printf2("[Experimental] ml_eta set to %.4f\n", ml_eta); }
            }
        }

    }


void do_log(void)
	{
	int error = FALSE, newlogfile=FALSE, start=FALSE, stop=FALSE, i=0;
	

    for(i=0; i<num_commands; i++)
        {

        if(strcmp(parsed_command[i], "file") == 0)
             {
             newlogfile=TRUE;
             if(logfile!= NULL)
		     	{
		     	printf2("closing log file %s \n", logfile_name);
		     	fclose(logfile);
		     	}
		     
	     	strcpy(logfile_name,parsed_command[i+1]);
	     	
	     	if((logfile = fopen(logfile_name, "w")) == NULL)
				{
				printf2("Error opening log file named %s\n", logfile_name);
				error = TRUE;
				}
			else
				printf2("opening log file %s \n", logfile_name);
		     
         	 }
        }
    for(i=0; i<num_commands; i++)
        {
        if(strcmp(parsed_command[i], "status") == 0)
        	{	
	        if(strcmp(parsed_command[i+1], "on") == 0)
		        {
		        if(print_log == FALSE)
		        	{	
			        if(logfile == NULL)
			        	{
			        	printf2("opening logfile %s\n", logfile_name);	
			        	if((logfile = fopen(logfile_name, "w")) == NULL)
			        		printf2("Error opening log file named %s\n", logfile_name);
			        	}
			        printf2("starting logging to logfile %s\n", logfile_name );
			        print_log = TRUE;
			        start=TRUE;
			    	}
			    else
			    	{
			    	printf2("Already logging output, no changes have been made\n");
			    	}
			    }
			else
				{
			    if(strcmp(parsed_command[i+1], "off") == 0)
			        {
			        if(print_log == TRUE)
			        	{	
			        	printf2("stopping logging to file %s\n", logfile_name);
			        	print_log = FALSE;
			        	stop=TRUE;
			        	}
			        else
			        	{	
			        	printf2("Cannot stop logging as there is no logging currently happenning\n");
			        	}
			        }
			    else
			    	{
			    	printf2("Error: \"%s\" is not a valid modifer for option \"status\", please choose \"on\" or \"off\"", parsed_command[i+1]);
			    	}	
			    }
		    }

	    }
	if(newlogfile==FALSE && start==FALSE && stop == FALSE)
		{
		printf2("\nLog file status not changed (");
		if(logfile!= NULL) 
			{
			printf2("log file \"%s\" is open ", logfile_name);
			if(print_log == TRUE)
				printf2("and logging is currently on)\n");
			else
				printf2("and logging is currently off)\n");

			}
		else printf2("log file not open)\n");
		}	


	}


