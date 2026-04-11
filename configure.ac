dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)

AC_INIT([clann], [5.0.0], [chris.creevey@gmail.com], [clann], [https://github.com/ChrisCreevey/clann])

AM_INIT_AUTOMAKE([1.5 foreign])
AM_MAINTAINER_MODE([disable])
AC_CONFIG_HEADERS([config.h])

dnl -----------------------------------------------------------------------
dnl OpenMP compiler selection (macOS only).
dnl
dnl Apple's Clang does not ship an OpenMP runtime, so AC_OPENMP reports
dnl "unsupported" on a stock macOS install.  We inspect the effective
dnl compiler's --version output; if it identifies itself as Apple Clang /
dnl Apple LLVM we search Homebrew prefixes for a real GCC (gcc-15 ... gcc-11)
dnl on both Apple Silicon (/opt/homebrew) and Intel (/usr/local).  When one
dnl is found we set CC and clear the autoconf cache so AC_PROG_CC and then
dnl AC_OPENMP pick it up automatically.
dnl
dnl Reasons this is more robust than checking "${CC+set}":
dnl   * CC may be set to "cc" or "gcc" which are Apple Clang symlinks on
dnl     a stock macOS install — the old guard would silently skip the probe.
dnl   * config.cache from a previous run can hold ac_cv_prog_CC=clang; we
dnl     clear it whenever we override CC so the cache stays consistent.
dnl
dnl On Linux GCC supports OpenMP natively; no special handling is needed.
dnl A user who explicitly supplies a non-Clang CC is left entirely alone.
dnl -----------------------------------------------------------------------
case "$(uname -s 2>/dev/null)" in
Darwin*)
    dnl Determine the effective compiler: honour an explicit CC, else "cc".
    _eff_cc="${CC:-cc}"
    dnl If it reports itself as Apple Clang/LLVM it has no OpenMP runtime.
    if "$_eff_cc" --version 2>&1 | grep -q 'Apple clang\|Apple LLVM'; then
        _brew_gcc=
        for _gcc_ver in 15 14 13 12 11; do
            for _brew_prefix in /opt/homebrew /usr/local; do
                if test -x "$_brew_prefix/bin/gcc-$_gcc_ver"; then
                    _brew_gcc="$_brew_prefix/bin/gcc-$_gcc_ver"
                    break 2
                fi
            done
        done
        if test -n "$_brew_gcc"; then
            CC="$_brew_gcc"
            dnl Clear any cached compiler choice from a previous configure run.
            AS_UNSET([ac_cv_prog_CC])
            AC_MSG_NOTICE([macOS: Apple Clang detected; switching to Homebrew GCC at $CC for OpenMP support])
        else
            AC_MSG_NOTICE([macOS: Apple Clang detected and no Homebrew GCC found; clann will be built single-threaded])
            AC_MSG_NOTICE([  To enable multi-threading: brew install gcc  (then re-run ./configure)])
        fi
    fi
    ;;
esac

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

dnl Advisory when no OpenMP runtime is available.
dnl OPENMP_CFLAGS is empty when the compiler does not support -fopenmp.
if test -z "$OPENMP_CFLAGS"; then
    AC_MSG_WARN([OpenMP is not supported by $CC — clann will run single-threaded.])
    AC_MSG_WARN([On macOS:  brew install gcc   then re-run ./configure])
    AC_MSG_WARN([On Linux:  ensure gcc (not clang) is available, or load the])
    AC_MSG_WARN([           appropriate compiler module on your HPC system.])
fi

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
