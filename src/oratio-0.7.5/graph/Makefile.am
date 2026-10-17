the_headers =

the_src = graph.c dijkstra.c label_correcting.c shortest_path.c \
	edge.c node.c

lib_LTLIBRARIES = liborgraph.la

liborgraph_la_SOURCES = $(the_headers) $(the_src)

liborgraph_la_LIBADD =
liborgraph_la_LDFLAGS = -no-undefined -L../ctools -lctools -L../dstruct -ldstruct -L../ntools -lntools -L../stat -lstat

orgraphincludedir = $(pkgincludedir)
orgraphinclude_HEADERS = $(the_headers)

EXTRA_DIST = $(the_headers)
