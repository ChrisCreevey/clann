/*
 *  clann_api.c — Clann v5.0.0  (library-mode API implementation)
 *
 *  Implements the public C API declared in clann_api.h.  This file is only
 *  compiled when building the shared library target (Makefile rule
 *  "libclann.so"), which passes -DCLANN_LIBRARY_MODE to every compilation
 *  unit.
 *
 *  Copyright (C) 2003-2026 Chris Creevey <chris.creevey@gmail.com>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */

#include "clann_api.h"

#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "globals.h"
#include "main.h"
#include "tree_io.h"
#include "utils.h"

/* -----------------------------------------------------------------------
 * Error-recovery globals (definitions; declarations are in clann_api.h).
 * ----------------------------------------------------------------------- */
jmp_buf clann_jmp_exit;
int     clann_exit_code = 0;

/* -----------------------------------------------------------------------
 * Internal: library command dispatcher.
 *
 * Called after parse_command() has filled parsed_command[] and num_commands.
 * Mirrors the if-else dispatch in main()'s REPL loop, covering the most
 * commonly used analysis commands.  Returns 0 on success, 1 on error.
 * ----------------------------------------------------------------------- */
static int clann_library_dispatch(void)
    {
    const char *cmd = parsed_command[0];

    if(strcmp(cmd, "quit") == 0 || strcmp(cmd, "exit") == 0)
        return 0;

    /* ---------- tree loading ----------------------------------------- */
    if(strcmp(cmd, "execute") == 0 || strcmp(cmd, "exe") == 0)
        {
        if(num_commands > 1)
            execute_command(parsed_command[1], TRUE);
        return 0;
        }

    /* ---------- global settings -------------------------------------- */
    if(strcmp(cmd, "set") == 0)
        {
        set_parameters();
        return 0;
        }

    /* All remaining commands require loaded trees. */
    if(number_of_taxa <= 0)
        {
        printf2("Error: You need to load source trees before using this command\n");
        return 1;
        }

    /* ---------- heuristic search ------------------------------------- */
    if(strcmp(cmd, "hs") == 0 || strcmp(cmd, "hsearch") == 0)
        {
        /* nreps is read from parsed_command[] inside heuristic_search().
         * Pass 10 as the default; heuristic_search() overrides it. */
        heuristic_search(TRUE, TRUE, 10000, 10);
        return 0;
        }

    /* ---------- bootstrap -------------------------------------------- */
    if(strcmp(cmd, "bootstrap") == 0 || strcmp(cmd, "boot") == 0)
        {
        bootstrap_search();
        return 0;
        }

    /* ---------- neighbour-joining tree ------------------------------- */
    if(strcmp(cmd, "nj") == 0)
        {
        nj();
        return 0;
        }

    /* ---------- consensus tree --------------------------------------- */
    if(strcmp(cmd, "consensus") == 0)
        {
        do_consensus();
        return 0;
        }

    /* ---------- user-tree scoring ------------------------------------ */
    if(strcmp(cmd, "usertrees") == 0)
        {
        if(criterion == 1)
            {
            printf2("Error: You cannot do a user-tree search in MRP\n");
            return 1;
            }
        usertrees_search();
        return 0;
        }

    /* ---------- exhaustive topology search --------------------------- */
    if(strcmp(cmd, "alltrees") == 0)
        {
        if(criterion == 0 || criterion == 2 || criterion == 3 ||
           criterion == 5 || criterion == 6 || criterion == 7)
            alltrees_search(TRUE);
        else
            {
            printf2("Error: alltrees is not supported for the current criterion in library mode\n");
            return 1;
            }
        return 0;
        }

    /* ---------- show / save trees ------------------------------------ */
    if(strcmp(cmd, "showtrees") == 0)
        {
        showtrees(FALSE);
        return 0;
        }
    if(strcmp(cmd, "savetrees") == 0)
        {
        showtrees(TRUE);
        return 0;
        }

    /* ---------- pairwise distance analyses --------------------------- */
    if(strcmp(cmd, "rfdists") == 0)
        {
        sourcetree_dists();
        return 0;
        }
    if(strcmp(cmd, "sprdists") == 0)
        {
        spr_dist();
        return 0;
        }

    /* ---------- ML scoring ------------------------------------------- */
    if(strcmp(cmd, "mlscores") == 0)
        {
        mlscores();
        return 0;
        }

    /* ---------- reconciliation --------------------------------------- */
    if(strcmp(cmd, "reconstruct") == 0)
        {
        reconstruct(1);
        return 0;
        }

    /* ---------- recluster -------------------------------------------- */
    if(strcmp(cmd, "recluster") == 0)
        {
        execute_recluster();
        return 0;
        }

    /* ---------- taxon management ------------------------------------- */
    if(strcmp(cmd, "excludetrees") == 0)
        {
        exclude(TRUE);
        return 0;
        }
    if(strcmp(cmd, "includetrees") == 0)
        {
        include(TRUE);
        return 0;
        }
    if(strcmp(cmd, "deletetaxa") == 0)
        {
        exclude_taxa(TRUE);
        return 0;
        }
    if(strcmp(cmd, "restoretaxa") == 0)
        {
        if(restoretaxa_available)
            restoretaxa(TRUE);
        else
            printf2("Error: no deletetaxa snapshot available to restore\n");
        return 0;
        }

    /* ---------- misc ------------------------------------------------- */
    if(strcmp(cmd, "generatetrees") == 0)
        {
        generatetrees();
        return 0;
        }

    printf2("Error: command '%s' not supported in library mode\n", cmd);
    return 1;
    }

/* -----------------------------------------------------------------------
 * clann_init
 * ----------------------------------------------------------------------- */
int clann_init(void)
    {
    int i;

    /* Initialise static string buffers that are accessed before trees are
     * loaded.  These are global char arrays defined in treecompare2.c. */
    inputfilename[0]    = '\0';
    saved_supertree[0]  = '\0';
    logfile_name[0]     = '\0';
    strcpy(logfile_name, "clann.log");
    system_call[0]      = '\0';

    /* Seed the random-number generator. */
    seed = (int)(time(NULL) / 2);
    srand((unsigned)seed);

    /* test_array: integer workspace used by topology operations. */
    if(test_array == NULL)
        {
        test_array = malloc(TREE_LENGTH * sizeof(int));
        if(!test_array) return -1;
        test_array[0] = '\0';
        }

    /* tempsuper: temporary Newick buffer. */
    if(tempsuper == NULL)
        {
        tempsuper = malloc(TREE_LENGTH * sizeof(char));
        if(!tempsuper) return -1;
        tempsuper[0] = '\0';
        }

    /* stored_commands: queue for commands pre-loaded from scripts or CLI. */
    if(stored_commands == NULL)
        {
        stored_commands = malloc(100 * sizeof(char *));
        if(!stored_commands) return -1;
        for(i = 0; i < 100; i++)
            {
            stored_commands[i] = malloc(10000 * sizeof(char));
            if(!stored_commands[i]) return -1;
            stored_commands[i][0] = '\0';
            }
        }

    /* retained_supers / scores_retained_supers: best-topology storage. */
    if(retained_supers == NULL)
        {
        retained_supers = malloc(number_retained_supers * sizeof(char *));
        if(!retained_supers) return -1;
        for(i = 0; i < number_retained_supers; i++)
            {
            retained_supers[i] = malloc(TREE_LENGTH * sizeof(char));
            if(!retained_supers[i]) return -1;
            retained_supers[i][0] = '\0';
            }
        }
    if(scores_retained_supers == NULL)
        {
        scores_retained_supers = malloc(number_retained_supers * sizeof(float));
        if(!scores_retained_supers) return -1;
        for(i = 0; i < number_retained_supers; i++)
            scores_retained_supers[i] = -1;
        }

    /* parsed_command: tokenised command array read by analysis functions. */
    if(parsed_command == NULL)
        {
        parsed_command = malloc(10000 * sizeof(char *));
        if(!parsed_command) return -1;
        for(i = 0; i < 10000; i++)
            {
            parsed_command[i] = malloc(10000 * sizeof(char));
            if(!parsed_command[i]) return -1;
            parsed_command[i][0] = '\0';
            }
        }

    return 0;
    }

/* -----------------------------------------------------------------------
 * clann_set_output_fn
 * ----------------------------------------------------------------------- */
void clann_set_output_fn(clann_output_fn_t fn, void *userdata)
    {
    clann_output_fn   = fn;
    clann_output_data = userdata;
    }

/* -----------------------------------------------------------------------
 * clann_load_trees
 * ----------------------------------------------------------------------- */
int clann_load_trees(const char *filename, const char *parse_opts)
    {
    char cmdline[20200];

    if(setjmp(clann_jmp_exit) != 0)
        return clann_exit_code ? clann_exit_code : -1;

    if(parse_opts && parse_opts[0] != '\0')
        snprintf(cmdline, sizeof(cmdline), "exe %s %s", filename, parse_opts);
    else
        snprintf(cmdline, sizeof(cmdline), "exe %s", filename);

    num_commands = parse_command(cmdline);
    if(num_commands > 1)
        execute_command(parsed_command[1], TRUE);

    return 0;
    }

/* -----------------------------------------------------------------------
 * clann_run_command
 * ----------------------------------------------------------------------- */
int clann_run_command(const char *commandline)
    {
    if(setjmp(clann_jmp_exit) != 0)
        return clann_exit_code ? clann_exit_code : -1;

    num_commands = parse_command((char *)commandline);
    if(num_commands <= 0)
        return 0;

    return clann_library_dispatch();
    }

/* -----------------------------------------------------------------------
 * clann_reset
 *
 * Uses clean_exit(0) to free all analysis state, catching the longjmp it
 * raises, then nulls the pointers that clean_exit() freed without zeroing,
 * and calls clann_init() to re-establish a clean baseline.
 * ----------------------------------------------------------------------- */
void clann_reset(void)
    {
    if(setjmp(clann_jmp_exit) == 0)
        clean_exit(0);   /* frees everything, then longjmps back here */

    clann_exit_code = 0;

    /* clean_exit() frees these but does not set the pointers to NULL.
     * Zero them so clann_init() re-allocates them cleanly. */
    stored_commands        = NULL;
    parsed_command         = NULL;
    retained_supers        = NULL;
    scores_retained_supers = NULL;
    tempsuper              = NULL;
    test_array             = NULL;

    clann_init();
    }
