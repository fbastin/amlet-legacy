the_headers = faure.h mlhs.h error.h

the_src = moca.c \
	generation.c \
	faure.c halton.c lattice.c mlhs.c sobol.c \
	baker.c randomization.c distributions.c \
	bsplines.c

lib_LTLIBRARIES = libmoca.la

libmoca_la_SOURCES = $(the_headers) $(the_src)

libmoca_la_LIBADD =
libmoca_la_LDFLAGS = -no-undefined -L../ctools -lctools -L../ntools -lntools -L../ran -lran

mocaincludedir = $(pkgincludedir)
mocainclude_HEADERS = $(the_headers)

EXTRA_DIST = $(the_headers)
