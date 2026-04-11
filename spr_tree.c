/*
 *  spr_tree.c — Clann v5.0.0
 *  String-native SPR/TBR tree representation and operations.
 *
 *  All tree operations work by pure Newick serialisation — the tree is
 *  treated as an immutable value and every SPR/TBR move produces a new
 *  string.  No in-place pointer surgery is ever performed.
 *
 *  Tree invariant maintained by spr_parse():
 *    - EVERY non-root node has parent != NULL.
 *    - ALL siblings in a sibling list have parent pointing to their
 *      common parent (no "first-sibling-only" clann invariant here).
 *
 *  Copyright (C) 2003-2026 Chris Creevey <chris.creevey@gmail.com>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */

#include "spr_tree.h"

#include <stdlib.h>
#include <string.h>
#include <ctype.h>


/* =======================================================================
 * Memory management
 * ======================================================================= */

struct spr_node *spr_node_new(void)
    {
    return calloc(1, sizeof(struct spr_node));
    }

void spr_free(struct spr_node *root)
    {
    struct spr_node *c, *next;
    if(root == NULL) return;
    c = root->first_child;
    while(c != NULL)
        {
        next = c->next_sib;
        spr_free(c);
        c = next;
        }
    free(root->label);
    free(root);
    }


/* =======================================================================
 * Newick parser
 * ======================================================================= */

/* parse_subtree: recursive descent parser.
 * *pos is the current position in the Newick string s.
 * Returns the root of the parsed subtree, or NULL on error. */
static struct spr_node *parse_subtree(const char *s, int *pos)
    {
    struct spr_node *node, *child;
    struct spr_node **tail;
    int start, len;

    node = spr_node_new();
    if(node == NULL) return NULL;

    /* Skip leading whitespace */
    while(s[*pos] == ' ' || s[*pos] == '\t') (*pos)++;

    if(s[*pos] == '(')
        {
        /* Internal node */
        (*pos)++;  /* skip '(' */
        tail = &node->first_child;

        while(1)
            {
            /* Skip whitespace before child */
            while(s[*pos] == ' ' || s[*pos] == '\t') (*pos)++;

            child = parse_subtree(s, pos);
            if(child == NULL) { spr_free(node); return NULL; }
            child->parent = node;
            *tail = child;
            tail  = &child->next_sib;

            /* Skip whitespace after child */
            while(s[*pos] == ' ' || s[*pos] == '\t') (*pos)++;

            if(s[*pos] == ')') { (*pos)++; break; }   /* end of children */
            if(s[*pos] == ',') { (*pos)++; continue; } /* next child */

            /* Unexpected character — malformed Newick */
            spr_free(node);
            return NULL;
            }
        }
    else
        {
        /* Leaf node: read the label (any characters up to ',', ')', ';', ':') */
        start = *pos;
        while(s[*pos] != '\0' && s[*pos] != ',' && s[*pos] != ')' &&
              s[*pos] != ';'  && s[*pos] != ':')
            (*pos)++;
        len = *pos - start;
        node->label = malloc(len + 1);
        if(node->label == NULL) { free(node); return NULL; }
        memcpy(node->label, s + start, len);
        node->label[len] = '\0';
        }

    /* Skip optional branch length ":value" */
    if(s[*pos] == ':')
        {
        (*pos)++;
        while(s[*pos] != '\0' && s[*pos] != ',' &&
              s[*pos] != ')' && s[*pos] != ';')
            (*pos)++;
        }

    return node;
    }

struct spr_node *spr_parse(const char *nwk)
    {
    int pos = 0;
    struct spr_node *root;
    if(nwk == NULL) return NULL;
    /* Skip leading whitespace */
    while(nwk[pos] == ' ' || nwk[pos] == '\t') pos++;
    root = parse_subtree(nwk, &pos);
    return root;
    }


/* =======================================================================
 * Newick serialisation
 * ======================================================================= */

/* write_node: append the Newick fragment for node's subtree to buf. */
static void write_node(const struct spr_node *node, char *buf)
    {
    const struct spr_node *c;
    int first;

    if(node->first_child == NULL)
        {
        /* Leaf */
        strcat(buf, node->label ? node->label : "?");
        }
    else
        {
        strcat(buf, "(");
        first = 1;
        for(c = node->first_child; c != NULL; c = c->next_sib)
            {
            if(!first) strcat(buf, ",");
            first = 0;
            write_node(c, buf);
            }
        strcat(buf, ")");
        }
    }

void spr_write(const struct spr_node *root, char *buf, size_t bufsz)
    {
    (void)bufsz;   /* caller is responsible for sufficient size */
    buf[0] = '\0';
    if(root == NULL) { strcat(buf, ";"); return; }
    write_node(root, buf);
    strcat(buf, ";");
    }


/* =======================================================================
 * Tree query
 * ======================================================================= */

int spr_count_nodes(const struct spr_node *root)
    {
    const struct spr_node *c;
    int n;
    if(root == NULL) return 0;
    n = 1;
    for(c = root->first_child; c != NULL; c = c->next_sib)
        n += spr_count_nodes(c);
    return n;
    }


/* =======================================================================
 * Edge and internal-node enumeration
 * ======================================================================= */

static void collect_edges_rec(const struct spr_node *node,
                               struct spr_node **out, int *count)
    {
    struct spr_node *c;
    for(c = node->first_child; c != NULL; c = c->next_sib)
        {
        out[(*count)++] = c;
        collect_edges_rec(c, out, count);
        }
    }

void spr_get_edges(const struct spr_node *root,
                   struct spr_node **out, int *count)
    {
    *count = 0;
    if(root == NULL) return;
    collect_edges_rec(root, out, count);
    }

static void collect_internals_rec(const struct spr_node *node,
                                   struct spr_node **out, int *count)
    {
    struct spr_node *c;
    if(node->first_child != NULL)
        {
        out[(*count)++] = (struct spr_node *)node;
        for(c = node->first_child; c != NULL; c = c->next_sib)
            collect_internals_rec(c, out, count);
        }
    }

void spr_get_internals(const struct spr_node *root,
                       struct spr_node **out, int *count)
    {
    *count = 0;
    if(root == NULL) return;
    collect_internals_rec(root, out, count);
    }


/* =======================================================================
 * Complement serialisation — the engine behind bisect and rerooting.
 *
 * write_complement(node, exclude, buf):
 *   Append the Newick fragment for all of node's descendants EXCEPT
 *   the subtree rooted at exclude.
 *
 *   Degree-2 contraction: if an internal node would be left with exactly
 *   one child after the exclusion, it is contracted (no parens emitted —
 *   the sole remaining child is written in-place).  This applies to ALL
 *   nodes, including the root: when the root has exactly 2 children and
 *   one is excluded, writing "(remaining_child)" would produce a spurious
 *   degree-2 wrapper in the Newick string.  The correct output is just
 *   "remaining_child_subtree" with no wrapping parens.
 * ======================================================================= */

static void write_complement(const struct spr_node *node,
                              const struct spr_node *exclude,
                              char *buf)
    {
    const struct spr_node *c;
    int n, first;

    if(node == exclude) return;   /* skip the excluded subtree entirely */

    if(node->first_child == NULL)
        {
        /* Leaf (and not excluded) */
        strcat(buf, node->label ? node->label : "?");
        return;
        }

    /* Count children that survive the exclusion */
    n = 0;
    for(c = node->first_child; c != NULL; c = c->next_sib)
        if(c != exclude) n++;

    if(n == 0)
        return;   /* degenerate: excluded subtree was the only child */

    if(n == 1)
        {
        /* Degree-2 contraction: skip this node's parens and recurse into
         * the sole surviving child.  Applies to root AND non-root nodes.
         * Without this, a 2-child root that loses one child would emit
         * "(remaining_child)" — a degree-2 wrapper — instead of just
         * "remaining_child_subtree". */
        for(c = node->first_child; c != NULL; c = c->next_sib)
            if(c != exclude)
                { write_complement(c, exclude, buf); break; }
        return;
        }

    /* General case: emit "(child1, child2, ...)" for all surviving children */
    strcat(buf, "(");
    first = 1;
    for(c = node->first_child; c != NULL; c = c->next_sib)
        {
        if(c == exclude) continue;
        if(!first) strcat(buf, ",");
        first = 0;
        write_complement(c, exclude, buf);
        }
    strcat(buf, ")");
    }


/* =======================================================================
 * Bisection
 * ======================================================================= */

int spr_bisect(const struct spr_node *root,
               const struct spr_node *edge_node,
               char *sub_nwk, char *rem_nwk, size_t bufsz)
    {
    (void)bufsz;

    if(edge_node == NULL || edge_node == root) return -1;

    /* --- Pruned subtree: serialise edge_node's own subtree --- */
    sub_nwk[0] = '\0';
    write_node(edge_node, sub_nwk);
    strcat(sub_nwk, ";");

    /* --- Remaining tree: complement of edge_node in root --- */
    rem_nwk[0] = '\0';
    write_complement(root, edge_node, rem_nwk);
    strcat(rem_nwk, ";");

    return 0;
    }


/* =======================================================================
 * SPR regraft — pure serialisation graft
 *
 * write_graft(node, graft_node, sub_content, buf):
 *   Walk node's subtree.  When graft_node is encountered, instead of
 *   writing its subtree directly, write:
 *       (sub_content, graft_node_subtree)
 *   which inserts a new internal node (the re-attachment point) between
 *   graft_node and its parent.
 * ======================================================================= */

static void write_graft(const struct spr_node *node,
                        const struct spr_node *graft_node,
                        const char            *sub_content,
                        char                  *buf)
    {
    const struct spr_node *c;
    int first;

    if(node == graft_node)
        {
        /* Insert new internal node: (pruned_subtree, graft_node_subtree) */
        strcat(buf, "(");
        strcat(buf, sub_content);
        strcat(buf, ",");
        write_node(node, buf);
        strcat(buf, ")");
        return;
        }

    if(node->first_child == NULL)
        {
        strcat(buf, node->label ? node->label : "?");
        return;
        }

    strcat(buf, "(");
    first = 1;
    for(c = node->first_child; c != NULL; c = c->next_sib)
        {
        if(!first) strcat(buf, ",");
        first = 0;
        write_graft(c, graft_node, sub_content, buf);
        }
    strcat(buf, ")");
    }

void spr_graft(const char            *sub_nwk,
               const struct spr_node *rem_root,
               const struct spr_node *graft_node,
               char *out_nwk, size_t bufsz)
    {
    char *sub_content;
    char *semi;

    (void)bufsz;

    /* Strip trailing ';' from sub_nwk for embedding as a subtree fragment */
    sub_content = malloc(strlen(sub_nwk) + 1);
    if(sub_content == NULL) { out_nwk[0] = '\0'; return; }
    strcpy(sub_content, sub_nwk);
    semi = strrchr(sub_content, ';');
    if(semi) *semi = '\0';

    out_nwk[0] = '\0';
    write_graft(rem_root, graft_node, sub_content, out_nwk);
    strcat(out_nwk, ";");

    free(sub_content);
    }


/* =======================================================================
 * TBR subtree rerooting
 *
 * spr_write_rerooted(sub_root, reroot_node, buf):
 *   Produce a Newick string "(reroot_node_subtree, complement);".
 *
 *   This corresponds to cutting the branch between reroot_node and its
 *   parent inside sub_root, then treating reroot_node's subtree as one
 *   component and the rest as the other.  Used by TBR to enumerate all
 *   internal edges of the pruned subtree as possible re-attachment roots.
 * ======================================================================= */

void spr_write_rerooted(const struct spr_node *sub_root,
                        const struct spr_node *reroot_node,
                        char *buf, size_t bufsz)
    {
    (void)bufsz;

    buf[0] = '\0';

    if(reroot_node == NULL || reroot_node == sub_root)
        {
        /* No rerooting: write sub_root as-is */
        write_node(sub_root, buf);
        strcat(buf, ";");
        return;
        }

    /* "(reroot_node_subtree, complement_at_reroot_node);" */
    strcat(buf, "(");
    write_node(reroot_node, buf);
    strcat(buf, ",");
    write_complement(sub_root, reroot_node, buf);
    strcat(buf, ");");
    }
