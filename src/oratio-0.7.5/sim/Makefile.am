the_headers =

the_src = accumulate.cpp esim.cpp event.cpp

lib_LTLIBRARIES = liborsim.la

liborsim_la_SOURCES = $(the_headers) $(the_src)

liborsim_la_LIBADD =
liborsim_la_LDFLAGS = -no-undefined -L../ctools -lctools \
	-L../dstruct -ldstruct -L../ntools -lntools

orsimincludedir = $(pkgincludedir)
orsiminclude_HEADERS = $(the_headers)

EXTRA_DIST = $(the_headers)
