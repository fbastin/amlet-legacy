check_PROGRAMS = arng_test

arng_test_SOURCES = test_ran.c \
	test_uniform.c test_acceptance_rejection.c \
	cauchy_test.c \
	test_exponential.c test_normal.c test_truncatednormal.c \
	test_mvn.c test_m_unit_sphere.c \
	test_ran_sp.c

arng_test_LDADD = -lran -lstat -lntools -lctools -ldstruct

AM_CPPFLAGS=-I$(top_srcdir)/src

TESTS = $(check_PROGRAMS)

CLEANFILES = toto

EXTRA_DIST = uniforms_test.h \
	cauchy_test.h test_ran.h \
	test_m_unit_sphere.h \
	test_ran_sp.h
