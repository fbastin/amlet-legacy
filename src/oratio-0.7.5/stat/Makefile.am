the_headers = percentile.h \
	histogram.h truncation.h \
	stat_error.h

the_src = descriptive.c empirical.c estimator.c percentile.c utils.c \
	tally.c tallystore.c \
	uniform.c \
	betasymmetrical.c cauchy.c chi2.c exponential.c fisher.c \
	geometric.c gumbel.c \
	laplace.c  gamma.c lognormal.c logistic.c \
	normal.c pareto.c poisson.c rayleigh.c student.c triangular.c \
	weibull.c \
	mvn.c \
	histogram.c truncation.c \
	gof.c likelihood.c

lib_LTLIBRARIES = libstat.la

libstat_la_SOURCES = $(the_headers) $(the_src)

libstat_la_LIBADD =
libstat_la_LDFLAGS = -no-undefined -L../ctools -lctools -L../dstruct -ldstruct -L../ntools -lntools

statincludedir = $(pkgincludedir)
statinclude_HEADERS = $(the_headers)

EXTRA_DIST = $(the_headers) ChangeLog
