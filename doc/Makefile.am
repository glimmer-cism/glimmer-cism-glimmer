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

EXTRA_DIST = dg.tex dg_pdf.tex dg_ps.tex glimmer.tex glimmer_pdf.tex glimmer_ps.tex num.tex num_pdf.tex num_ps.tex ug.tex ug_pdf.tex ug_ps.tex user_doc.tex

CLEANFILES = *~ *.aux *.dvi *.log *.out *.toc

dg.pdf:	dg_pdf.dvi
	$(DVIPDFT) -o dg.pdf dg_pdf.dvi

dg.ps:	dg_ps.dvi
	dvips -o dg.ps dg_ps.dvi

num.pdf:	num_pdf.dvi
	$(DVIPDFT) -o num.pdf num_pdf.dvi

num.ps:	num_ps.dvi
	dvips -o num.ps num_ps.dvi

ug.pdf:	ug ug_pdf.dvi
	$(DVIPDFT) -o ug.pdf ug_pdf.dvi

ug.ps:	ug ug_ps.dvi
	dvips -o ug.ps ug_ps.dvi

glimmer.pdf:	glimmer_pdf.dvi
	$(DVIPDFT) -o glimmer.pdf glimmer_pdf.dvi

glimmer.ps:	glimmer_ps.dvi
	dvips -o glimmer.ps glimmer_ps.dvi

code_doc.ps:	code_doc.dvi
		dvips -o code_doc.ps code_doc.dvi

code_doc.tex:	../src/fortran/*.f90
		f90docp.py -t "GLIMMER Code Documentation" -a "Glimmer Developer Team" $^ > $@

.tex.dvi::
	latex $<
	latex $<
