## Makefile.am -- Process this file with automake to produce Makefile.in
bin_PROGRAMS = clann
clann_SOURCES = treecompare2.c
clann_LDADD = $(LDADD) $(LIBREADLINE) -lm
clann_LDFLAGS = -g -O0 -Wall -funroll-loops -fsanitize=address

