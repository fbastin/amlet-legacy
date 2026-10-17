AC_INIT([ophelia], [0.7.4])

AC_CANONICAL_SYSTEM

AC_CONFIG_SRCDIR([ophelia/nlp.h])
AM_INIT_AUTOMAKE([subdir-objects])

AM_CONFIG_HEADER(config.h)

AC_PROG_CC
AC_PROG_CXX
AC_PROG_F77
AC_ISC_POSIX

dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE

AC_EXEEXT
AC_LIBTOOL_DLOPEN

AM_PROG_LIBTOOL
AC_PROG_INSTALL
AC_PROG_MAKE_SET

AC_PATH_PROGS(LATEX, latex)

dnl Check for hearder files
AC_STDC_HEADERS
AC_CHECK_HEADER(values.h)
AC_CHECK_HEADER(sys/types.h)


echo $CC $target $prefix

dnl C compiler flags
case $CC in
gcc*)
  DEBUG="-g -Wall"
  CFLAGS="$CFLAGS -I${includedir} -Wall -fomit-frame-pointer"
  PROFILE="-Wall"
  case $target in
  *-*-linux*)
    ;;
  ppcMac-darwin)
    AC_CHECK_LIB(dl, dlopen)
    if test "$ac_cv_lib_dl_dlopen" = yes; then
    LD_DYN_FLAGS="-ldl -dynamic"
    SFLAGS="-fpic -DPIC"
SLDFLAGS="-dynamic -twolevel_namespace - weak_reference_mismatches weak -undefined dynamic_lookup"
LD=libtool
ac_lib_dl=yes
fi
;;
ix86Mac-darwin)
AC_CHECK_LIB(dl, dlopen)
if test "$ac_cv_lib_dl_dlopen" = yes; then
LD_DYN_FLAGS="-ldl -dynamic"
SFLAGS="-fpic -DPIC"
SLDFLAGS="-dynamic -twolevel_namespace - weak_reference_mismatches weak -undefined dynamic_lookup"
LD=libtool
ac_lib_dl=yes
fi
;;  esac
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
    CFLAGS="$CFLAGS -ieee"
    ;;
  esac
  ;;
*)
  DEBUG="-g"
  CFLAGS="$CFLAGS"
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
AC_CHECK_LIB(dstruct, ds_list_new,, AC_MSG_ERROR([Cannot find dstruct library]))
AC_CHECK_LIB(ntools, nt_matrix_new,, AC_MSG_ERROR([Cannot find ntools library]))
AC_CHECK_LIB(stat, st_normal_icdf,, AC_MSG_ERROR([Cannot find ran library]))
AC_CHECK_LIB(ran, ran_random_new,, AC_MSG_ERROR([Cannot find ran library]))

AC_ARG_WITH(debug,
[  --with-debug      Compile with debuggage options],[
  CFLAGS="$CFLAGS -g -DDEBUG"
  echo "Use debuggage options"
])

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

AC_OUTPUT(Makefile       \
          ophelia/Makefile test/Makefile doc/Makefile \
	  lapack/Makefile \
	  collection/Makefile collection/HS/Makefile collection/SP/Makefile \
	  collection/leastsquares/Makefile \
	  collection/minpack2/Makefile collection/minpack2/tprobs/Makefile \
	  collection/unconstrained/Makefile \
	  examples/Makefile
         )
