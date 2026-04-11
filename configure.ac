dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)

AC_INIT([clann], [5.0.0], [chris.creevey@gmail.com], [clann], [https://github.com/ChrisCreevey/clann])

AM_INIT_AUTOMAKE([1.5 foreign])
AM_MAINTAINER_MODE([disable])
AC_CONFIG_HEADERS([config.h])

dnl -----------------------------------------------------------------------
dnl Readline detection.
dnl
dnl Prefer GNU readline over macOS libedit:
dnl   1. Check common Homebrew prefixes for GNU readline first.
dnl      If found, prepend to CPPFLAGS/LDFLAGS so AC_CHECK_LIB finds it.
dnl   2. Fall back to the system readline (libedit on macOS, GNU on Linux).
dnl   3. If nothing is found, build without readline (plain fgets prompt).
dnl -----------------------------------------------------------------------

AC_SUBST([LIBREADLINE])

dnl -- Step 1: look for Homebrew GNU readline (macOS only) -----------------
READLINE_PREFIX=""
for candidate in \
        "$(brew --prefix readline 2>/dev/null)" \
        /usr/local/opt/readline \
        /opt/homebrew/opt/readline ; do
    if test -n "$candidate" && test -f "$candidate/include/readline/readline.h"; then
        READLINE_PREFIX="$candidate"
        break
    fi
done

if test -n "$READLINE_PREFIX"; then
    AC_MSG_NOTICE([Found GNU readline at $READLINE_PREFIX — preferring over system libedit])
    CPPFLAGS="-I$READLINE_PREFIX/include $CPPFLAGS"
    LDFLAGS="-L$READLINE_PREFIX/lib $LDFLAGS"
fi

dnl -- Step 2: find terminal library (readline's link-time dependency) ------
dnl On Linux readline typically needs -ltinfo or -lncurses to link.
dnl Detect whichever is present so we can pass it to the readline test.
READLINE_TERMLIB=""
AC_CHECK_LIB([tinfo], [tgetent],
    [READLINE_TERMLIB="-ltinfo"],
    [AC_CHECK_LIB([ncurses], [initscr],
        [READLINE_TERMLIB="-lncurses"],
        [READLINE_TERMLIB=""])])

dnl -- Step 3: check for readline, linking against the terminal library -----
AC_CHECK_LIB([readline], [readline],
    [AC_SUBST([LIBREADLINE], ["-lreadline $READLINE_TERMLIB"])
     AC_DEFINE([HAVE_READLINE], [], [Define if GNU readline or libedit is available])
    ],
    [echo "No readline found — Clann will use standard input instead"
     echo "For full CLI functionality: brew install readline  (macOS)"
     echo "                     or:   apt install libreadline-dev  (Debian/Ubuntu)"],
    [$READLINE_TERMLIB])


dnl continue

AC_PROG_CC

AC_OPENMP

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
