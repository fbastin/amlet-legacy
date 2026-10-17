SUFFIXES = .dvi

docdir = $(prefix)/share
docname = amlet_manual

EXTRA_DIST = amlet_manual.tex biblio.bib \
	basics.tex data.tex install.tex intro.tex license.tex models.tex \
	numerical.tex pst.pdf

CLEANFILES = $(wildcard *.dvi *.ps *.toc *.log *.bbl *.blg *.aux)

#if PDF
#  TARGETS = $(docname).pdf
#endif
#if POSTSCRIPT
  TARGETS_PS = $(docname).ps
  TARGETS_PDF = $(docname).pdf
#endif
if HTML
  TARGETS_HTML = $(docname).html
endif

all: $(TARGETS) $(TARGETS_PDF)

install-data-local: $(docname).pdf
	$(mkinstalldirs) $(docdir)/amlet
	if [ -n "$(LATEX)" ]; then $(INSTALL_DATA) $< $(docdir)/amlet/$<; fi

%.html: $(docname).tex
	if [ -n "$(LATEX2HTML)" ]; then latex2html $<; latex2html $<; fi	


%.pdf: $(docname).tex
	export TEXINPUTS=$(TEXINPUTS):$(top_srcdir)/doc;
	if [ "$(LATEX)" != "" ]; then pdflatex -shell-escape $<; bibtex $*; makeindex $*; pdflatex -shell-escape $<; pdflatex -shell-escape $<; fi

%.ps: %.dvi
	if [ -n "$(LATEX)" ]; then dvips -Ppdf -G0 $< -o; fi

%.dvi: $(docname).tex
	if [ -n "$(LATEX)" ]; then export TEXINPUTS="$(TEXINPUTS)":$(top_srcdir)/doc; latex $<; bibtex $*; latex $<; latex $<; fi

clean:
	rm -rf $(CLEANFILES) $(docname).pdf $(docname).ps
