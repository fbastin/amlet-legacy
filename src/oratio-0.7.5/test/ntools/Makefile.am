check_PROGRAMS = ntools_test

ntools_test_SOURCES = ntools_test.c\
	algebra_test.c bsplines_test.c machar_test.c modular_test.c \
	natural_test.c polynomial.c primes_test.c

ntools_test_LDADD = $(LDADD) $(top_srcdir)/ntools/libntools.la -lctools -ldstruct

INCLUDES = -I$(top_srcdir)/ntools

TESTS = $(check_PROGRAMS)

CLEANFILES =

EXTRA_DIST = ntools_test.h
