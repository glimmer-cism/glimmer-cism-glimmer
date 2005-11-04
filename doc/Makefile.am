if DO_DOC_PDF
 pdf_doc = glimmer.pdf
endif
if DO_DOC_PS
 ps_doc = glimmer.ps
endif
doc_DATA = $(ps_doc) $(pdf_doc)
docdir = $(prefix)/doc
dist_doc_DATA = $(doc_DATA)
SUBDIRS = dg ug num common

EXTRA_DIST = dg.tex dg_pdf.tex dg_ps.tex glimmer.tex glimmer_pdf.tex glimmer_ps.tex num.tex num_pdf.tex num_ps.tex ug.tex ug_pdf.tex ug_ps.tex 

CLEANFILES = *~ *.aux *.dvi *.log *.out *.toc

code_doc.ps:	code_doc.dvi
		dvips -o code_doc.ps code_doc.dvi

code_doc.tex:	../src/fortran/*.f90
		f90docp.py -t "GLIMMER Code Documentation" -a "Glimmer Developer Team" $^ > $@

.tex.dvi::
	$(LATEX) $<
	$(BIBTEX) $*
	$(LATEX) $<
	$(LATEX) $<

%.ps:	%_ps.dvi
	dvips -o $@ $<

%.pdf:	%_pdf.dvi
	$(DVIPDFT) -o $@ $<