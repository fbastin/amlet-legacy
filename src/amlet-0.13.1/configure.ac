dnl Process this file with autoconf to produce a configure script.
AC_INIT(src/amlet.c, 0.13.1, fabian.bastin@cirrelt.ca)

AM_CONFIG_HEADER(config.h)

dnl Allows cross compilation
AC_CANONICAL_SYSTEM
AC_ARG_PROGRAM

AM_INIT_AUTOMAKE(amlet, 0.13.1)

AC_ISC_POSIX
AC_PROG_CC
AC_PROG_CXX
AC_C_CONST
AC_EXEEXT
AC_LIBTOOL_DLOPEN

AM_PROG_LIBTOOL
AC_PROG_INSTALL
AC_PROG_MAKE_SET

dnl We need Lex and Yacc in order to parse driver files
AM_PROG_LEX
AC_PROG_YACC

dnl Checks for programs
AC_HEADER_STDC

dnl Checks C compiler capacities
AC_C_INLINE
AC_C_CONST

AC_CONFIG_MACRO_DIR([m4])

dnl Before configuring libtool check for --enable-all-static option
AC_ARG_ENABLE(all-static,
[ --enable-all-static build completely static (standalone) binaries],
[
 if test "$enableval" = "yes" ; then
 LT_LDFLAGS="-all-static $LT_LDFLAGS"
 fi
])

AC_SUBST(LT_LDFLAGS) 

dnl Check that LaTeX is available for documentation generation
AC_PATH_PROGS(LATEX, latex)
AC_PATH_PROGS(PS2PDF14, ps2pdf14)
AC_PATH_PROGS(PS2PDF, ps2pdf)
AC_PATH_PROG(LATEX2HTML, latex2html)

doc_type=pdf
AC_ARG_ENABLE(doc,
[ --enable-doc[=doc_type] Specify the documentation style:
	pdf [default]: Adobe Portable Document Format
	ps: Postscript
	html: HTML],
[if test "$enableval" = "pdf"; then
   doc_type = pdf
 elif test "$enableval" = "ps" ; then
   doc_type=ps
 elif test "$enableval" = "html" ; then
   doc_type=html
 elif test "$enableval" != "yes"; then
   echo
   echo "Error!"
   echo "Unknown documentation format"
   exit -1
 fi
])
AM_CONDITIONAL(PDF, test x$doc_type = xpdf)
AM_CONDITIONAL(POSTSCRIPT, test x$doc_type = xps)
AM_CONDITIONAL(HTML, test x$doc_type = xhtml)

dnl perl path
AC_ARG_WITH(perl, [  --with-perl=PATH        the Perl executable that Automake should use],
PERL=$withval, [AC_PATH_PROG(PERL, perl)])

AM_CONDITIONAL(PERLEXT, false)
if test -x "$PERL"; then
  AM_CONDITIONAL(PERLEXT, true)
  AC_MSG_RESULT(Search perl libray for integration with C)
dnl We could probably do better things by explicitely checking the modules
  PERLCCOPTS=`perl -MExtUtils::Embed -e ccopts`
  PERLLDOPTS=`perl -MExtUtils::Embed -e ldopts`
else
  AC_MSG_WARN([perl not found (try the --with-perl flag)])
fi
 
dnl Check for header files
AC_HEADER_TIME
AC_CHECK_HEADERS(sys/time.h)
AC_CHECK_HEADERS(signal.h)
AC_CHECK_HEADERS(stdlib.h)

echo $CC $target $prefix

CXXFLAGS="$CXXFLAGS"

dnl C compiler flags
case $CC in
gcc*)
  DEBUG="-g -Wall"
  CFLAGS="$CFLAGS -I${includedir} -Wall -fomit-frame-pointer"
  CXXFLAGS="$CFLAGS -I${includedir} -Wall -fomit-frame-pointer"
  PROFILE="-Wall"
  case $target in
  *-*-linux*)
    CFLAGS="$CFLAGS -O3"
    PROFILE="$PROFILE -W -pg -g -O3"
    case $target-cpu in
      i586-*)
        CFLAGS="$CFLAGS -march=pentium"
        ;;
      i686-*)
        CFLAGS="$CFLAGS -march=pentiumpro"
        ;;
      esac
    ;;
  esac
  ;;
pgcc*)
  DEBUG="-g"
  if test "$prefix"; then
  CFLAGS="$CFLAGS -Minline -c9x -L$prefix/lib -I$prefix/include"
  CXXFLAGS="$CFLAGS -Minline -c9x -L$prefix/lib -I$prefix/include"
  else
  CFLAGS="$CFLAGS -Minline -c9x "
  CXXFLAGS="$CFLAGS -Minline -c9x "
  fi
  PROFILE=""
  ;;
cc)
  DEBUG="-g"
  CFLAGS="$CFLAGS"
  CXXFLAGS="$CFLAGS"
  PROFILE=""
  case $target in
  *-*-osf*)
    DEBUG="-g -Wall"
    CFLAGS="$CFLAGS -O3 -ieee"
    CXXFLAGS="$CFLAGS -O3 -ieee"
    ;;
  esac
  ;;
*)
  DEBUG="-g"
  CFLAGS="$CFLAGS"
  CXXFLAGS="$CFLAGS"
  PROFILE=""
  ;;
esac

dnl Check if MPI is available
AC_ARG_ENABLE(mpicc,
[ --enable-mpicc Specify that mpicc should be used for compilation],
[AC_CHECK_PROG(MPICC,mpicc,yes,no)

if test "$MPICC" = "yes"; then
  CFLAGS="$CFLAGS -DMPI"
  CC="mpicc";
else
  AC_MSG_WARN([MPI not found; there will be no parallel code generation.])
fi
])

dnl Check for libraries
AC_CHECK_LIB(m, main,, AC_MSG_ERROR([Cannot find mathematical library]))

dnl BLAS LIBRARY
AC_CHECK_LIB(gfortran, main, LIBS="$LIBS -lgfortran",)
AC_CHECK_LIB(openblas, main, LIBS="$LIBS -lopenblas -lpthread",
  AC_CHECK_LIB(goto2, main, LIBS="$LIBS -lgoto2 -lpthread",
  AC_CHECK_LIB(cblas, main, LIBS="$LIBS -lcblas",
  AC_CHECK_LIB(gslcblas, main, AC_CHECK_LIB(blas, main, LIBS="$LIBS -lblas -lgslcblas",),
  AC_CHECK_LIB(atlas, main, LIBS=" -latlas $LIBS"
    AC_CHECK_LIB(ptcblas, main, LIBS="-lptcblas $LIBS"
      AC_CHECK_LIB(pthread, main,,)),
        AC_MSG_ERROR([Unable to find a supported BLAS library]))))))

AC_CHECK_LIB(ctools, c_malloc,, AC_MSG_ERROR([Cannot find ctools library]))
AC_CHECK_LIB(dstruct, ds_slist_new,, AC_MSG_ERROR([Cannot find dstruct library]))
AC_CHECK_LIB(ntools, nt_matrix_new,, AC_MSG_ERROR([Cannot find ntools library]))
AC_CHECK_LIB(stat, st_mean,, AC_MSG_ERROR([Cannot find stat stat library]))
AC_CHECK_LIB(ran, ran_random_new,, AC_MSG_ERROR([Cannot find ran library]))
AC_CHECK_LIB(moca, main,, AC_MSG_ERROR([Cannot find moca library]))
AC_CHECK_LIB(gfortran, main,,)
AC_CHECK_LIB(ophelia, op_btr_init,, AC_MSG_ERROR([Cannot find OPHELIA]))

AC_ARG_WITH(verbose,
[ --with-verbose	Display debugging information],[
	CFLAGS="$CFLAGS -DDEBUG"
	CXXFLAGS="$CXXFLAGS -DDEBUG"
	echo "Verbose mode on"
])

AC_ARG_WITH(debug,
[  --with-debug      Compile with debuggage options],[
  CFLAGS="$DEBUG"
  CXXFLAGS="$DEBUG"
  echo "Use debuggage options"
])

AC_ARG_WITH(aggressive,
[ --with-aggressive	Do aggressive optimizations of the code - use it with care!],[
	CFLAGS="$CFLAGS -ffast-math"
	CXXFLAGS="$CXXFLAGS -ffast-math"
	echo "Do aggressive optimizations of the code - use it with care!"
])

AC_ARG_WITH(single,
[ --with-single	Use single precision numbers for sample generation],[
	CFLAGS="$CFLAGS -DUSESINGLE"
	CXXFLAGS="$CXXFLAGS -DUSESINGLE"
	echo "Use single float numbers for sample generation"
])

dnl Make substitutions
AC_SUBST(LIBTOOL_DEPS)
AC_SUBST(DEBUG)
AC_SUBST(CFLAGS)
AC_SUBST(CXXFLAGS)
AC_SUBST(PROFILE)
AC_SUBST(PERLCCOPTS)
AC_SUBST(PERLLDOPTS)

AC_OUTPUT(Makefile amlet/Makefile src/Makefile doc/Makefile test/Makefile)
