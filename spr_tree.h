/*
 *  spr_tree.h — Clann v5.0.0
 *  String-native SPR/TBR tree representation and operations.
 *
 *  A clean phylogenetic tree struct with always-valid parent pointers,
 *  free from the clann-specific invariant that only first siblings carry
 *  their parent pointer.  All operations work by pure Newick serialisation —
 *  no in-place pointer surgery on any tree.
 *
 *  Copyright (C) 2003-2026 Chris Creevey <chris.creevey@gmail.com>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */

#ifndef CLANN_SPR_TREE_H
#define CLANN_SPR_TREE_H

#include <stddef.h>

/* -----------------------------------------------------------------------
 * spr_node: minimal phylogenetic tree node.
 *
 * Invariant: EVERY node has parent != NULL except the root.
 *            ALL siblings have their parent pointer set.
 *            (This differs from the clann taxon struct, where only the
 *            first sibling at each level has parent set.)
 *
 * Leaf nodes:     label != NULL,  first_child == NULL
 * Internal nodes: label == NULL,  first_child != NULL
 * ----------------------------------------------------------------------- */
struct spr_node {
    char            *label;        /* taxon label (leaf) or NULL (internal) */
    struct spr_node *parent;       /* NULL only for root                    */
    struct spr_node *first_child;  /* first child; NULL for leaves          */
    struct spr_node *next_sib;     /* next sibling; NULL for last child     */
};

/* -----------------------------------------------------------------------
 * Memory management
 * ----------------------------------------------------------------------- */

/* Allocate and zero-initialise a single spr_node */
struct spr_node *spr_node_new(void);

/* Recursively free the tree rooted at root */
void spr_free(struct spr_node *root);

/* -----------------------------------------------------------------------
 * Parsing and serialisation
 * ----------------------------------------------------------------------- */

/* Parse a Newick string (integer-label or name-label, with or without
 * trailing ";") into an spr_node tree.
 * Returns the root node, or NULL on allocation failure or parse error. */
struct spr_node *spr_parse(const char *nwk);

/* Write the tree rooted at root to buf as a Newick string with trailing ";".
 * buf must be at least bufsz bytes (use TREE_LENGTH from clann.h). */
void spr_write(const struct spr_node *root, char *buf, size_t bufsz);

/* -----------------------------------------------------------------------
 * Tree query
 * ----------------------------------------------------------------------- */

/* Count all nodes (leaves + internals) in the tree */
int spr_count_nodes(const struct spr_node *root);

/* -----------------------------------------------------------------------
 * Edge enumeration
 *
 * In an unrooted tree stored as a rooted tree, each non-root node
 * represents an edge (the edge between that node and its parent).
 * ----------------------------------------------------------------------- */

/* Enumerate all edges as an array of child-node pointers.
 * out[] must hold at least (spr_count_nodes(root) - 1) entries.
 * *count receives the number of entries written. */
void spr_get_edges(const struct spr_node *root,
                   struct spr_node **out, int *count);

/* Enumerate all internal nodes (for TBR rerooting).
 * out[] must hold at least spr_count_nodes(root) entries.
 * *count receives the number of entries written.
 * The root itself IS included if it is internal. */
void spr_get_internals(const struct spr_node *root,
                       struct spr_node **out, int *count);

/* -----------------------------------------------------------------------
 * Bisection (SPR/TBR outer step)
 * ----------------------------------------------------------------------- */

/* Split the tree at the edge above edge_node (edge_node ↔ its parent).
 *
 * sub_nwk  ← Newick for edge_node's subtree       (includes trailing ";")
 * rem_nwk  ← Newick for the complement             (includes trailing ";")
 *            Degree-2 internal nodes created by the removal are contracted.
 *
 * Both buffers must be at least bufsz bytes.
 * Returns  0 on success.
 * Returns -1 if edge_node is NULL or is the root (no bisection possible). */
int spr_bisect(const struct spr_node *root,
               const struct spr_node *edge_node,
               char *sub_nwk, char *rem_nwk, size_t bufsz);

/* -----------------------------------------------------------------------
 * Regrafting (SPR inner step) — pure serialisation, no tree modification
 * ----------------------------------------------------------------------- */

/* Produce a candidate Newick by grafting sub_nwk onto the edge ABOVE
 * graft_node in the remaining tree rem_root.
 *
 * The graft inserts a new internal node between graft_node and its parent,
 * with the new node having two children: graft_node and the sub_nwk subtree.
 *
 * sub_nwk  — Newick string for the pruned subtree (trailing ";" optional).
 * out_nwk  — output buffer, must be at least bufsz bytes.
 *
 * graft_node must not be the root of rem_root (use spr_get_edges() to
 * enumerate valid graft positions). */
void spr_graft(const char      *sub_nwk,
               const struct spr_node *rem_root,
               const struct spr_node *graft_node,
               char *out_nwk, size_t bufsz);

/* -----------------------------------------------------------------------
 * TBR subtree rerooting
 * ----------------------------------------------------------------------- */

/* Write sub_root's subtree rerooted so that reroot_node is the root of
 * one component, to buf.  The output is a valid unrooted Newick
 * "(reroot_subtree, complement_of_reroot);".
 *
 * Used by TBR to enumerate all internal edges of the pruned subtree
 * and try each as a different re-attachment point.
 *
 * buf must be at least bufsz bytes. */
void spr_write_rerooted(const struct spr_node *sub_root,
                        const struct spr_node *reroot_node,
                        char *buf, size_t bufsz);

#endif /* CLANN_SPR_TREE_H */
