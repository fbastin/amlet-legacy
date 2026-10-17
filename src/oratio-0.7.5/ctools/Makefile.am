the_headers = symbol.h

the_src = ctime.c ctimer.c error.c \
	exception.c file.c c_csv.c \
	memory.c cstring.c symbol.c system.c \
	getopt_long.c

lib_LTLIBRARIES = libctools.la

libctools_la_SOURCES = $(the_headers) $(the_src)

libctools_la_LIBADD =
libctools_la_LDFLAGS = -no-undefined

ctoolsincludedir = $(pkgincludedir)
ctoolsinclude_HEADERS = $(the_headers)

EXTRA_DIST = $(the_headers)
