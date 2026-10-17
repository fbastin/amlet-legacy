bin_PROGRAMS = amlet

amlet_SOURCES = amlet.c exec.c

AM_CFLAGS=-DDATADIR="$(pkgdatadir)"

INCLUDES = -I@top_srcdir@/amlet

amlet_LDADD = -L../amlet -lamlet

amlet_LDFLAGS =

AM_YFLAGS = -d -v
AM_LFLAGS = -l

clean:
	rm -f driver_parser.c driver_scanner.c
