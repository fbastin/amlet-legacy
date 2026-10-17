the_headers = ophelia.h nlp.h \
	common_struct.h wrapper.h \
	errors.h io.h solutions.h utilities.h \
	problem.h op_sp_collection.h \
	constraints.h projection.h isotonic.h \
	btr.h btr_struct.h btr_defaults.h defaults.h \
	rtr.h \
	model.h radius.h stopping.h \
	step.h btr_step.h pcauchy.h pcg.h pcg_monotonicity.h random_step.h \
	bfgs.h bfgs_struct.h linesearch.h \
	gss.h sa.h sa_step.h \
	secant.h sr1.h \
	bounds.h \
	barrier.h data.h \
	outer.h \
	nlp_least_squares.h \
	activesets.h \
	matrix.h \
	clapack.h \
	bsplines.h nonparametric.h \
	finance.h option.h \
	collection.h \
	ophelia.hpp btr.hpp problem.hpp solver.hpp

the_src = problem.c \
	wrapper.c io.c solutions.c utilities.c \
	gss.c \
	btr.c rho.c rtr.c \
	model.c btr_struct.c \
	step.c cauchy.c pcauchy.c pcg.c pcg_monotonicity.c random_step.c \
	stopping.c \
	radius.c \
	constraints.c projection.c isotonic.c \
	bfgs.c bfgs_struct.c dcsrch.c linesearch.c \
	secant.c sr1.c psb.c \
	bounds.c barrier.c \
	outer.c \
	sa.c \
	nlp_least_squares.c \
	activesets.c \
	matrix.c \
	opproblem.cpp solver.cpp \
	op_btr.cpp \
	../lapack/ieeeck.f ../lapack/ilaenv.f ../lapack/iparmq.f \
	../lapack/dgelq2.f ../lapack/dgelqf.f \
	../lapack/dgeqr2.f ../lapack/dgeqrf.f \
	../lapack/dlabad.f ../lapack/dlacn2.f ../lapack/dlange.f \
	../lapack/dlarf.f ../lapack/dlarfb.f ../lapack/dlarfg.f \
	../lapack/dlascl.f ../lapack/dlaset.f \
	../lapack/disnan.f \
	../lapack/dlacpy.f ../lapack/dlaisnan.f ../lapack/dlamch.f \
	../lapack/dlansy.f ../lapack/dlassq.f ../lapack/dlaqsy.f ../lapack/dlapy2.f \
	../lapack/dlarft.f \
	../lapack/dlasyf.f ../lapack/dlatrs.f ../lapack/dlauu2.f \
	../lapack/dlauum.f \
	../lapack/dorm2r.f ../lapack/dorml2.f ../lapack/dormlq.f ../lapack/dormqr.f \
	../lapack/dpocon.f ../lapack/dpoequ.f ../lapack/dporfs.f ../lapack/dposvx.f \
	../lapack/dpotf2.f ../lapack/dpotrf.f ../lapack/dpotri.f ../lapack/dpotrs.f \
	../lapack/drscl.f ../lapack/dsytf2.f ../lapack/dsytrf.f ../lapack/dsytri.f \
	../lapack/dtrti2.f ../lapack/dtrtri.f ../lapack/dtrtrs.f \
	../lapack/dgels.f \
	../lapack/lsame.f

lib_LTLIBRARIES = libophelia.la

libophelia_la_SOURCES = $(the_headers) $(the_src)

libophelia_la_LIBADD =
libophelia_la_LDFLAGS =

opheliaincludedir = $(pkgincludedir)
opheliainclude_HEADERS = $(the_headers)

EXTRA_DIST = $(the_headers)
