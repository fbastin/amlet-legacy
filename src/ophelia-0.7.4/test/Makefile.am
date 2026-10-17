check_PROGRAMS = test_nlp

test_nlp_headers = test_errors.h \
	test_bfgs.h \
	test_btr_unconstrained.h \
	test_btr_constrained_bounds.h test_btr_constrained_monotonicity.h \
	test_bounds.h bounds.h \
	test_outer.h test_ls.h \
	test_siman.h \
	test_isotonic.h \
	test_minpack2.h \
	ophelia_test.h

test_nlp_SOURCES = test_nlp.c \
	test_bfgs.c \
	test_btr_unconstrained.c test_btr_constrained_bounds.c \
	test_btr_constrained_monotonicity.c \
	bounds.c \
	test_outer.c test_ls.c \
	test_gss.c \
	siman.c test_matrix.c \
	test_isotonic.c
# \
#	test_minpack2.c test_tron.c

# Depending of the system, we have to specify gfortran library, or ignore it.
#test_nlp_LDADD = -lophelia -lgfortran -lopcoll -lcminpack2 -lmoca
test_nlp_LDADD = -lophelia -lopcoll -lcminpack2 -lmoca

AM_CPPFLAGS=-I$(top_srcdir)/src

TESTS = $(check_PROGRAMS)

CLEANFILES = barrier.txt

EXTRA_DIST = $(test_nlp_headers) tron.dat
