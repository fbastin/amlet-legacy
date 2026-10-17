SUBDIRS = tprobs

the_headers = minpack2_collection.h

the_src = minpack2_wrapper.c wrapper.f

lib_LIBRARIES = libcminpack2.a
libcminpack2_a_SOURCES = $(the_headers) $(the_src)

cminpack2includedir = $(pkgincludedir)
cminpack2include_HEADERS = $(the_headers)

EXTRA_DIST = $(the_headers)
