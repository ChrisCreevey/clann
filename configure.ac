dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)

AC_INIT([clann], [4.2.5], [chris.creevey@gmail.com], [clann], [https://github.com/ChrisCreevey/clann])

AM_INIT_AUTOMAKE([1.5 foreign])
AM_CONFIG_HEADER(config.h)

dnl do check for readline

AC_SUBST([LIBREADLINE])

AC_CHECK_LIB(readline, readline, 
		       [AC_SUBST([LIBREADLINE], ["-lreadline -lncurses"])
			AC_DEFINE([HAVE_READLINE], [], [READLINE DEF])
             		], [
			echo "No readline found, Clann will use standard input instead"
			echo "Please install readline for full IO functionality" ] )


dnl continue

AC_PROG_CC

AC_CONFIG_FILES([Makefile])
AC_OUTPUT




          
 
