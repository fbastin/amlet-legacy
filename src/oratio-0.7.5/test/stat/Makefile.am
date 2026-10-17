check_PROGRAMS = stat_test

stat_test_SOURCES = stat_test.c \
	cauchy_test.c cauchy_test.h \
	descriptive_test.c descriptive_test.h \
	empirical_test.c empirical_test.h \
	lognormal.c lognormal.h poisson.c poisson.h \
	mvn_test.c \
	stat_test.h

stat_test_LDADD = -lran -lstat -lntools -ldstruct -lctools -lm

INCLUDE=-I$(top_srcdir)/oratio

TESTS = $(check_PROGRAMS)
