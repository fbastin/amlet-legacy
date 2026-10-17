the_headers =

the_src = gsl.c parabola_1.c unconstrained.c \
	cgt1.c rosenbrock.c \
	ackley.c chichinadze.c

noinst_HEADERS =
noinst_LTLIBRARIES = libopcollunconstrained.la

libopcollunconstrained_la_SOURCES = $(the_headers) $(the_src)

EXTRA_DIST = $(the_headers)

# Disable autoheader. 
AUTOHEADER=echo
