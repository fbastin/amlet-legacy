SUBDIRS = HS SP leastsquares minpack2 unconstrained

the_headers = leastsquares.h bounds.h unconstrained.h

the_src =

lib_LTLIBRARIES = libopcoll.la

libopcoll_la_SOURCES = $(the_headers) $(the_src)

libopcoll_la_LIBADD = HS/libopcollHS.la SP/libopcollSP.la \
	leastsquares/libopcollleastsquares.la \
	unconstrained/libopcollunconstrained.la
libopcoll_la_LDFLAGS = -release 0.0.1 -version-info 1:0:0

opcollincludedir = $(pkgincludedir)/collection
opcollinclude_HEADERS = $(the_headers)

EXTRA_DIST = $(the_headers)
