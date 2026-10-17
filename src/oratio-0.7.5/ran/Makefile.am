the_headers = standardminimal.h \
	lecuyer88a.h lecuyer88b.h \
	mrg32k3a.h mt19937.h

the_src = random.c \
	standardminimal.c \
	lecuyer88a.c lecuyer88b.c \
	mrg32k3a.c \
	mt19937.c mt19937ar.c \
	bernouilli.c dirac.c geometric.c intrandom.c \
	arcsine.c cauchy.c exponential.c erlang.c gamma.c gumbel.c \
	johnson.c laplace.c normal.c poisson.c student.c triangular.c \
	truncatednormal.c \
	mvn.c m_unit_sphere.c \
	sp_spline.c spinversion.c \
	bootstrap.c permutations.c ratio.c \
	io.c lcg.c acceptance_rejection.c \
	brownian.c mpar.c \
	copula.c walk.c

lib_LTLIBRARIES = libran.la

libran_la_SOURCES = $(the_headers) $(the_src)

libran_la_LIBADD =
libran_la_LDFLAGS = -no-undefined -L../ctools -lctools -L../ntools -lntools -L../stat -lstat

ranincludedir = $(pkgincludedir)
raninclude_HEADERS = $(the_headers)

EXTRA_DIST = $(the_headers)
