dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)

AC_INIT([clann], [5.0.0], [chris.creevey@gmail.com], [clann], [https://github.com/ChrisCreevey/clann])

AM_INIT_AUTOMAKE([1.5 foreign])
AM_MAINTAINER_MODE([disable])
AC_CONFIG_HEADERS([config.h])

dnl AC_PROG_CC must come before any AC_CHECK_LIB calls so the C compiler
dnl is properly configured for link tests on all platforms.
AC_PROG_CC

dnl -----------------------------------------------------------------------
dnl Readline detection.
dnl
dnl Strategy:
dnl   1. Detect terminal library needed for linking readline on Linux
dnl      (-ltinfo preferred, -lncurses as fallback).
dnl   2. Try Homebrew readline locations (macOS) with saved/restored flags
dnl      so a wrong-arch library does not block the system readline.
dnl      Apple Silicon (/opt/homebrew) is checked before Intel (/usr/local).
dnl   3. Fall back to the system readline (libedit on macOS, GNU on Linux).
dnl   4. If nothing works, build without readline (plain fgets prompt).
dnl -----------------------------------------------------------------------

AC_SUBST([LIBREADLINE])

dnl -- Step 1: detect terminal library (needed on Linux) -------------------
READLINE_TERMLIB=""
AC_CHECK_LIB([tinfo],   [tgetent], [READLINE_TERMLIB="-ltinfo"],
    [AC_CHECK_LIB([ncurses], [initscr], [READLINE_TERMLIB="-lncurses"])])

dnl -- Step 2: try Homebrew readline locations (macOS) ---------------------
dnl Prefer /opt/homebrew (Apple Silicon ARM) over /usr/local (Intel/Rosetta).
dnl If the Homebrew library is the wrong architecture the link test fails;
dnl we restore CPPFLAGS/LDFLAGS and fall through to the system check.
READLINE_FOUND=no
READLINE_TRIED=""
for candidate in \
        "$(brew --prefix readline 2>/dev/null)" \
        /opt/homebrew/opt/readline \
        /usr/local/opt/readline ; do
    dnl skip empty or already-tried paths
    case " $READLINE_TRIED " in *" $candidate "*) continue ;; esac
    READLINE_TRIED="$READLINE_TRIED $candidate"
    if test -n "$candidate" && test -f "$candidate/include/readline/readline.h"; then
        AC_MSG_NOTICE([Trying GNU readline at $candidate])
        SAVED_CPPFLAGS="$CPPFLAGS"
        SAVED_LDFLAGS="$LDFLAGS"
        CPPFLAGS="-I$candidate/include $CPPFLAGS"
        LDFLAGS="-L$candidate/lib $LDFLAGS"
        AS_UNSET([ac_cv_lib_readline_readline])
        AC_CHECK_LIB([readline], [readline],
            [AC_SUBST([LIBREADLINE], ["-lreadline $READLINE_TERMLIB"])
             AC_DEFINE([HAVE_READLINE], [], [Define if GNU readline is available])
             READLINE_FOUND=yes],
            [dnl link failed (e.g. wrong architecture) — restore and try next
             CPPFLAGS="$SAVED_CPPFLAGS"
             LDFLAGS="$SAVED_LDFLAGS"],
            [$READLINE_TERMLIB])
        if test "$READLINE_FOUND" = "yes"; then break; fi
    fi
done

dnl -- Step 3: fall back to system readline --------------------------------
dnl Try plain, then with each terminal lib in turn.  Use AC_LINK_IFELSE
dnl directly to avoid AC_CHECK_LIB cache collisions across iterations.
if test "$READLINE_FOUND" = "no"; then
    for _rl_tl in "" "-ltinfo" "-lncurses" "-lcurses"; do
        _rl_save_LIBS="$LIBS"
        LIBS="-lreadline $_rl_tl $LIBS"
        AC_MSG_CHECKING([for readline with extra libs: '$_rl_tl'])
        AC_LINK_IFELSE(
            [AC_LANG_PROGRAM([[extern char *readline(const char *);]],
                             [[readline(""); return 0;]])],
            [AC_MSG_RESULT([yes])
             AC_SUBST([LIBREADLINE], ["-lreadline $_rl_tl"])
             AC_DEFINE([HAVE_READLINE], [], [Define if GNU readline is available])
             READLINE_FOUND=yes],
            [AC_MSG_RESULT([no])])
        LIBS="$_rl_save_LIBS"
        if test "$READLINE_FOUND" = "yes"; then break; fi
    done
fi
if test "$READLINE_FOUND" = "no"; then
    echo "No readline found — Clann will use standard input instead"
    echo "Install: apt install libreadline-dev  (Debian/Ubuntu)"
    echo "         yum install readline-devel   (RedHat/CentOS/Fedora)"
    echo "         brew install readline        (macOS)"
fi


AC_OPENMP

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
