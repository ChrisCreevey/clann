## Makefile.am -- Process this file with automake to produce Makefile.in
bin_PROGRAMS = clann
clann_SOURCES = treecompare2.c utils.c topology.c viz.c tree_io.c main.c \
                tree_ops.c scoring.c consensus.c reconcile.c prune.c spr_tree.c \
                treecluster.c
clann_LDADD = $(LDADD) $(LIBREADLINE) -lm
clann_CFLAGS  = -O2 -Wall -funroll-loops $(OPENMP_CFLAGS)
clann_LDFLAGS = $(OPENMP_CFLAGS)

## -----------------------------------------------------------------------
## libclann.so — shared library for Python (ctypes/CFFI) and C extension
## bindings.  Built separately to avoid contaminating the normal binary.
##
##   make libclann.so      — build the shared library
##   make libclann.so.clean — remove it and its object files
##
## The library uses the same sources as the binary plus clann_api.c, but
## compiles them with -DCLANN_LIBRARY_MODE which:
##   • excludes main() from main.c
##   • replaces _exit() with longjmp() in clean_exit()
##   • enables the output-callback path in printf2() (utils.c)
## -----------------------------------------------------------------------

LIBCLANN_SOURCES = treecompare2.c utils.c topology.c viz.c tree_io.c main.c \
                   tree_ops.c scoring.c consensus.c reconcile.c prune.c \
                   spr_tree.c treecluster.c clann_api.c

LIBCLANN_CFLAGS  = -DCLANN_LIBRARY_MODE -fPIC -O2 -Wall -funroll-loops \
                   $(OPENMP_CFLAGS) -I.

LIBCLANN_OBJS    = $(patsubst %.c,libclann_objs/%.o,$(LIBCLANN_SOURCES))

libclann_objs/%.o: %.c
	@mkdir -p libclann_objs
	$(CC) $(LIBCLANN_CFLAGS) -c -o $@ $<

libclann.so: $(LIBCLANN_OBJS)
	$(CC) -shared -o $@ $^ $(OPENMP_CFLAGS) -lm

libclann.so.clean:
	rm -rf libclann_objs libclann.so

.PHONY: libclann.so.clean
