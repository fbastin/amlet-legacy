SUFFIXES = .dvi

docdir = $(prefix)/share
docname = ophelia_manual

EXTRA_DIST = ophelia_manual.tex \
	basics.tex intro.tex analysis.tex constrained.tex \
	metaheuristics.tex tests.tex unconstrained.tex \
	biblio.bib op.pdf osshape.pdf pgpath.pdf \
	url.sty \
	gsl_e1.plot pbr_t1_3.plot gsl_e1.eps pbr_t1_3.eps

CLEANFILES = $(wildcard *.dvi *.ps *.toc *.log *.bbl *.blg *.aux) $(docname).pdf

all: $(docname).pdf

install-data-local: $(docname).pdf
	$(mkinstalldirs) $(docdir)/$(docname)
	if [ "$(LATEX)" != "" ]; then $(INSTALL_DATA) $< $(docdir)/$(docname)/$<; else $(INSTALL_DATA) $(docname).tex $(docdir)/$(docname)/$(docname).tex; fi

%.pdf: $(docname).tex
	export TEXINPUTS=$(TEXINPUTS):$(top_srcdir)/doc;
	if [ "$(LATEX)" != "" ]; then pdflatex $<; bibtex $*; makeindex $*; pdflatex $<; pdflatex $<; fi

%.ps: %.dvi
	if [ "$(LATEX)" != "" ]; then dvips $< -o; fi

%.dvi: $(docname).tex
	export TEXINPUTS=$(TEXINPUTS):$(top_srcdir)/doc;
	if [ "$(LATEX)" != "" ]; then "$(LATEX)" $<; bibtex $*; makeindex $*; "$(LATEX)" $<; "$(LATEX)" $<; fi

clean:
	rm -f $(CLEANFILES)
