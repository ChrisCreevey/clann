## Makefile.am -- Process this file with automake to produce Makefile.in
bin_PROGRAMS = clann
clann_SOURCES = treecompare2.c utils.c topology.c viz.c tree_io.c main.c \
                tree_ops.c scoring.c consensus.c reconcile.c prune.c spr_tree.c \
                treecluster.c
clann_LDADD = $(LDADD) $(LIBREADLINE) -lm
clann_CFLAGS  = -O2 -Wall -funroll-loops $(OPENMP_CFLAGS)
clann_LDFLAGS = $(OPENMP_CFLAGS)

