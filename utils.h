/*
 *  utils.h — Clann v5.0.0
 *  Utility functions: numeric conversion, I/O helpers, printf2
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

#ifndef CLANN_UTILS_H
#define CLANN_UTILS_H

#include "clann.h"

/* Globals defined in treecompare2.c, used by printf2 */
extern FILE *logfile;
extern int   print_log;

/* -----------------------------------------------------------------------
 * Library-mode output callback.
 *
 * When CLANN_LIBRARY_MODE is defined at compile time, printf2() forwards
 * every formatted message to this callback instead of writing to stdout or
 * logfile.  Set via clann_set_output_fn() (declared in clann_api.h).
 *
 * clann_output_fn_t : pointer to   void fn(const char *msg, void *userdata)
 * clann_output_fn   : current callback (NULL = use normal stdout/log path)
 * clann_output_data : opaque pointer forwarded to every callback invocation
 * ----------------------------------------------------------------------- */
#ifdef CLANN_LIBRARY_MODE
typedef void (*clann_output_fn_t)(const char *msg, void *userdata);
extern clann_output_fn_t clann_output_fn;
extern void             *clann_output_data;
#endif

/* Character/numeric conversion */
int   texttoint(char c);
char  inttotext(int c);
int   toint(char *number);
float tofloat(char *number);

/* I/O helpers */
char *xgets(char *s);

/* Logging printf */
void  printf2(char *format, ...);

/* Error handler — defined in treecompare2.c; declared here for universal access */
void  memory_error(int error_num);

#endif /* CLANN_UTILS_H */
