AC_INIT(oratio/ran.h, 0.7.5, bastin@iro.umontreal.ca)

AM_CONFIG_HEADER(config.h)

dnl Allows cross compilation
AC_CANONICAL_SYSTEM
AC_ARG_PROGRAM

AM_INIT_AUTOMAKE(oratio, 0.7.5)

AC_ISC_POSIX
AC_PROG_CC
AC_PROG_CXX
AC_C_CONST
AC_EXEEXT
AC_LIBTOOL_DLOPEN

AM_PROG_LIBTOOL
AC_PROG_INSTALL

LT_INIT([win32-dll])

echo $CC $target $prefix

dnl C compiler flags
case $CC in
gcc*)
  DEBUG="-g -Wall"
  CFLAGS="$CFLAGS -Wall -fomit-frame-pointer"
  PROFILE="-Wall"
  case $target in
  *-*-linux*)
    CFLAGS="$CFLAGS -O2"
    PROFILE="$PROFILE -W -pg -g -O3"
  esac
  ;;
pgcc*)
  DEBUG="-g"
  if test "$prefix"; then
  CFLAGS="$CFLAGS -Minline -c9x -L$prefix/lib -I$prefix/include"
  else
  CFLAGS="$CFLAGS -Minline -c9x "
  fi
  PROFILE=""
  ;;
cc)
  DEBUG="-g"
  CFLAGS="$CFLAGS"
  PROFILE=""
  case $target in
  *-*-osf*)
    DEBUG="-g -Wall"
    CFLAGS="$CFLAGS -O3 -ieee"
    ;;
  esac
  ;;
*)
  DEBUG="-g"
  CFLAGS="$CFLAGS"
  PROFILE=""
  ;;
esac

dnl BLAS LIBRARY
AC_CHECK_LIB(gfortran, main, LIBS="$LIBS -lgfortran",)
AC_CHECK_LIB(openblas, main, LIBS="$LIBS -lopenblas -lpthread -lm",
  AC_CHECK_LIB(goto2, main, LIBS="$LIBS -lgoto2 -lpthread -lm",
  AC_CHECK_LIB(cblas, main, LIBS="$LIBS -lcblas -lm",
  AC_CHECK_LIB(gslcblas, main, AC_CHECK_LIB(blas, main, LIBS="$LIBS -lblas -lgslcblas -lm",),
  AC_CHECK_LIB(atlas, main, LIBS=" -latlas $LIBS -lm"
    AC_CHECK_LIB(ptcblas, main, LIBS="-lptcblas $LIBS -lm"
      AC_CHECK_LIB(pthread, main,,)),
        AC_MSG_ERROR([Unable to find a supported BLAS library]))))))

dnl Checks for headers
AC_STDC_HEADERS
AC_HAVE_HEADERS(sys/types.h)

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T
AC_STRUCT_TM

AC_ARG_WITH(verbose,
[ --with-verbose	Display debugging information],[
	CFLAGS="$CFLAGS -DDEBUG"
	echo "Verbose mode on"
])

AC_ARG_WITH(aggressive,
[ --with-aggressive	Do aggressive optimizations of the code - use it with care!],[
	CFLAGS="$CFLAGS -ffast-math"
	echo "Do aggressive optimizations of the code - use it with care!"
])

dnl Make substitutions
AC_SUBST(LIBTOOL_DEPS)
AC_SUBST(DEBUG)
AC_SUBST(CFLAGS)
AC_SUBST(PROFILE)

AC_PATH_PROGS(LATEX, latex)

AC_CHECK_PROGS([DOXYGEN], [doxygen])
if test -z "$DOXYGEN";
   then AC_MSG_WARN([Doxygen not found - continuing without Doxygen support])
fi
AM_CONDITIONAL([HAVE_DOXYGEN], 
[test -n "$DOXYGEN"])AM_COND_IF([HAVE_DOXYGEN], [AC_CONFIG_FILES([doc/Doxyfile])])

AC_OUTPUT(oratio/Makefile
    ctools/Makefile
    dstruct/Makefile
    moca/Makefile
    ntools/Makefile ntools/machar/Makefile
    ran/Makefile
    sim/Makefile
    graph/Makefile
    stat/Makefile
    test/Makefile
    test/ctools/Makefile
    test/dstruct/Makefile
    test/graph/Makefile
    test/moca/Makefile
    test/ntools/Makefile
    test/ran/Makefile
    test/sim/Makefile
    test/stat/Makefile
    doc/Makefile
    examples/Makefile
    Makefile)
