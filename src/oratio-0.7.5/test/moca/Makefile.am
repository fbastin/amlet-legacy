check_PROGRAMS = caffeine

caffeine_SOURCES = caffeine.c test_distributions.c randomization.c \
	test_lattice.c test_sobol.c \
	caffeine.h

caffeine_LDADD = -lmoca -lran -lstat -lntools -ldstruct -lctools

INCLUDES=-I$(top_srcdir)/src

TESTS = $(check_PROGRAMS)

CLEANFILES =
