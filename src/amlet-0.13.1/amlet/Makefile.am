the_headers = amlet.hpp amlet.h \
	credits.h \
	amlet_config.h amlet_data.h data.h dcl.h model.h \
	lattice.h likelihood.h \
	logit.h mixedlogit.h \
	individual.h observations.h utilities.h \
	bhhh.h \
	optimization.h amlet_btr.h \
	alternative.h nest.h population.h information.h \
	variable.h nonparametric.h \
	bsplines.h \
	amlet_analysis.h validation.h \
	info.h log.h io.h \
	bootstrap.h \
	amlet_stopping.h \
	Train96.h Train99.h \
	external_strategy.h \
	amlet_driver.h driver_parser.h

the_src = model.c credits.c driver.c error.c information.c\
	lattice.c sample.c io.c likelihood.c solve.c variable.c \
	mixedlogit.c \
	individual.c observations.c utilities.c nonparametric.c \
	bhhh.c \
	memory.c \
	logit.c nest.c alternative.c \
	amlet_data.c \
	bsplines.c \
	amlet_analysis.c validation.c \
	bootstrap.c \
	info.c log.c \
	optimization.c amlet_btr.c amlet_stopping.c \
	Train96.c Train99.c \
	external_strategy.c \
	driver_parser.y driver_scanner.l

lib_LTLIBRARIES = libamlet.la

libamlet_la_SOURCES = $(the_headers) $(the_src)

libamlet_la_LIBADD =
libamlet_la_LDFLAGS =

amletincludedir = $(pkgincludedir)
amletinclude_HEADERS = $(the_headers)

EXTRA_DIST = $(the_headers)

AM_CFLAGS=-DDATADIR="$(pkgdatadir)"

AM_YFLAGS = -d -v
AM_LFLAGS = -l

#clean:
#
#	rm -f driver_parser.c driver_parser.h driver_scanner.c
