SUFFIXES = .dvi

docdir = $(prefix)/share
docname = oratio_manual

EXTRA_DIST = oratio_manual.tex intro.tex \
	ctools.tex dstruct.tex graph.tex moca.tex \
	ntools.tex ran.tex sim.tex stat.tex biblio.bib license.tex \
	delacroix.pdf \
	mathlist.sty eclbkbox.sty mycal.sty mybold.sty url.sty \
	pictex.sty pictex.tex prepictex.tex postpictex.tex

pdf: $(docname).pdf

install-data-local: $(docname).pdf
	$(mkinstalldirs) $(docdir)/oratio
	if [ "$(LATEX)" != "" ]; then $(INSTALL_DATA) $< $(docdir)/oratio/$<; else $(INSTALL_DATA) $(docname).tex $(docdir)/oratio/$(docname).tex; fi

%.pdf: $(docname).tex
	export TEXINPUTS=$(TEXINPUTS):$(top_srcdir)/doc;
	if [ "$(LATEX)" != "" ]; then pdflatex $<; bibtex $*; makeindex $*; pdflatex $<; pdflatex $<; fi

%.dvi: $(docname).tex
	export TEXINPUTS="$(TEXINPUTS)":$(top_srcdir)/doc;
	if [ "$(LATEX)" != "" ]; then $(LATEX) $<; bibtex $*; makeindex $*; $(LATEX) $<; $(LATEX) $<; fi

CLEANFILES = $(wildcard *.dvi *.toc *.log *.bbl *.blg *.aux) $(docname).pdf

if HAVE_DOXYGEN
#directory = $(top_srcdir)/doc/man/man3/

#dist_man_MANS = $(directory)/man_page_1.3 $(directory)/man_page_2.3
#$(directory)/man_page_1.3: doxyfile.stamp
#$(directory)/man_page_2.3: doxyfile.stamp

doxyfile.stamp:
	$(DOXYGEN) Doxyfile
	echo Timestamp > doxyfile.stamp

CLEANFILES = doxyfile.stamp

all-local: doxyfile.stamp
clean-local:
	rm -rf $(top_srcdir)/doc/man
endif

clean:
	rm -f $(CLEANFILES)
