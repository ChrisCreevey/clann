/*
 *  clann_api.h — Clann v5.0.0  (library-mode public API)
 *
 *  This header is only used when Clann is built as a shared library
 *  (compile with -DCLANN_LIBRARY_MODE).  It declares the public C API
 *  for in-process use by ctypes/CFFI bindings or a future Python C
 *  extension, and exposes the longjmp target that allows clean_exit()
 *  to unwind gracefully without terminating the host process.
 *
 *  Copyright (C) 2003-2026 Chris Creevey <chris.creevey@gmail.com>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */

#ifndef CLANN_API_H
#define CLANN_API_H

#include <setjmp.h>

#ifdef __cplusplus
extern "C" {
#endif

/* -----------------------------------------------------------------------
 * Output callback type (shared with utils.h / utils.c).
 *
 * void fn(const char *msg, void *userdata)
 *   msg      — null-terminated formatted output string.
 *   userdata — opaque pointer supplied to clann_set_output_fn().
 * ----------------------------------------------------------------------- */
typedef void (*clann_output_fn_t)(const char *msg, void *userdata);

/* -----------------------------------------------------------------------
 * Error-recovery globals used by clean_exit() in library mode.
 *
 * clann_run_command() and clann_load_trees() establish a setjmp() frame
 * before invoking any Clann code.  When clean_exit() is called (e.g. on a
 * fatal error), it stores the exit code in clann_exit_code and longjmps
 * to clann_jmp_exit rather than calling _exit(), keeping the host process
 * alive.
 * ----------------------------------------------------------------------- */
extern jmp_buf clann_jmp_exit;
extern int     clann_exit_code;

/* -----------------------------------------------------------------------
 * Public library API
 * ----------------------------------------------------------------------- */

/**
 * Initialise Clann global state.
 *
 * Allocates the buffers (parsed_command, stored_commands, retained_supers,
 * etc.) that Clann's analysis functions depend on.  Seeds the random-number
 * generator from the system clock.
 *
 * Must be called once before any other clann_* function.  Safe to call
 * again after clann_reset() to re-initialise for a fresh analysis.
 *
 * Returns 0 on success, -1 on allocation failure.
 */
int clann_init(void);

/**
 * Register (or clear) the output callback.
 *
 * When fn is non-NULL, every printf2() message is forwarded to
 * fn(msg, userdata) instead of being written to stdout / logfile.
 * Pass fn=NULL to revert to normal stdout output.
 */
void clann_set_output_fn(clann_output_fn_t fn, void *userdata);

/**
 * Load source trees from a file.
 *
 * Equivalent to the interactive "exe <filename> [parse_opts]" command.
 * parse_opts may be NULL or "" to use defaults.
 *
 * Must be called before analysis commands (hs, nj, consensus, …).
 *
 * Returns 0 on success, non-zero on error (including any Clann error that
 * would have called _exit() in the standalone binary).
 */
int clann_load_trees(const char *filename, const char *parse_opts);

/**
 * Execute a single Clann command string.
 *
 * The commandline must be in the format accepted by the interactive REPL,
 * e.g. "hs nreps=10 criterion=ml" or "set criterion dfit".
 *
 * Supported commands: exe/execute, set, hs/hsearch, bootstrap/boot, nj,
 * consensus, usertrees, alltrees, showtrees, savetrees, rfdists, sprdists,
 * mlscores, reconstruct, recluster, excludetrees, includetrees,
 * deletetaxa, restoretaxa, generatetrees, quit (no-op).
 *
 * Returns 0 on success, non-zero on error.
 */
int clann_run_command(const char *commandline);

/**
 * Reset all analysis state.
 *
 * Frees all dynamically allocated analysis data (trees, scores, taxa, …)
 * and re-initialises global state so that a new clann_init() +
 * clann_load_trees() cycle can begin in the same process.
 */
void clann_reset(void);

#ifdef __cplusplus
}
#endif

#endif /* CLANN_API_H */
