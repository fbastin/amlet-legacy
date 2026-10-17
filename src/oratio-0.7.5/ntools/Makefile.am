SUBDIRS = machar

the_headers = f2c.h cblas.h norm.h log.h

the_src = derive.c \
	matrix.c norm.c vector.c \
	bsplines.c \
	gcd.c io.c utilities.c values.c \
	functions.c gamma.c nt_polynomial.c \
	primes.c simplex.c \
	log.c \
	modular.c

lib_LTLIBRARIES = libntools.la

libntools_la_SOURCES = $(the_headers) $(the_src)

libntools_la_LIBADD = machar/libmachar.la
libntools_la_LDFLAGS = -no-undefined -L../ctools -lctools

ntoolsincludedir = $(pkgincludedir)
ntoolsinclude_HEADERS = $(the_headers)

EXTRA_DIST = $(the_headers)
