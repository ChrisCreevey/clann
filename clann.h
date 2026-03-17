/*
 *  clann.h — Clann v5.0.0
 *  Shared types, constants, and structure definitions
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

#ifndef CLANN_H
#define CLANN_H

#include "config.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <stdarg.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <signal.h>
#include <unistd.h> /* For Unix Systems Builds */
/*#include <process.h> */ /*  For Windows systems Builds */

#ifdef HAVE_READLINE /* from config.h */
#include <readline/readline.h>
#include <readline/history.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

/****** Constants ********/

#define FUNDAMENTAL_NUM 1000
#define TREE_LENGTH 1000000
#define TAXA_NUM 1000
#define NAME_LENGTH 1000
#define LOG10E 0.43429448190325182765  /* log10(e): for L.U.st ML likelihood conversion */

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

	int name;				/* This holds number which refers to the name of the taxon at this node if it doesn't point to a daughter node */
	char *fullname;				/* This holds the original name as it was given in the tree file, this allows us to record if there was parts of the name excluded earlier */
	float score;
	struct taxon *daughter;   	/* This points to a daughter node, but is only set if name is null */
	struct taxon *parent;		/* This points to the parent node, however its only set on the first sibling at each level */
	struct taxon *prev_sibling;	/* Points to the previous sibling at this node */
	struct taxon *next_sibling;	/* Points to the next sibling at this node */
	int tag;				/* A tag which is included as something which might be useful at some point */
	int tag2;				/* sometimes one of something just isn't enough ;o) */
	float loss;				/* A tag which will tell whether this node was lost.... only used for tree-mapping procedure */
        float xpos;   /* used to fix the position of the node for graphical representation */
        float ypos;
        int spr;
	char weight[100];  /* this will be used when we have a weight to add to the node like with a Bootstrap proportion */
	float length;
	int *donor; /* This is for the HGT analysis */
};

/* Open-addressed hash set for visited tree topologies */
typedef struct {
    uint64_t *keys;        /* 0 = empty slot */
    size_t    cap;         /* always a power of 2 */
    size_t    count;
    int       zero_present; /* special-case: was key==0 inserted? */
} VisitedSet;

/* Open-addressed hash map for visited-topology landscape recording */
typedef struct {
    uint64_t  hash;        /* topology hash (key; 0 = empty slot)  */
    float     score;       /* score on first visit                  */
    int       visit_count; /* total visits across all reps          */
    char     *newick;      /* heap-allocated named-taxon Newick     */
} LandscapeEntry;

typedef struct {
    LandscapeEntry *slots;
    size_t          capacity;  /* always a power of 2 */
    size_t          count;
} LandscapeMap;

typedef struct {
    uint64_t *hashes;   /* sorted array of canonical bipartition hashes */
    int       count;
} BipartSet;

#endif /* CLANN_H */
