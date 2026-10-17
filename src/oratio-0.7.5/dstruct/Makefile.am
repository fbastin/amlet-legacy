the_headers =

the_src = dsarray.c dslist.c dsslist.c dsqueue.c dssort.c dssearch.c \
	dsconvert.c dssarray.c dstree.c dsbtree.c

lib_LTLIBRARIES = libdstruct.la

libdstruct_la_SOURCES = $(the_headers) $(the_src)

libdstruct_la_LIBADD =
libdstruct_la_LDFLAGS = -no-undefined -L../ctools -lctools

dstructincludedir =  $(pkgincludedir)
dstructinclude_HEADERS = $(the_headers)

EXTRA_DIST = $(the_headers) ChangeLog
