if DO_DOC_PDF
 pdf_doc = glimmer.pdf
endif
doc_DATA = glimmer.ps $(pdf_doc)
docdir = $(prefix)/doc
dist_doc_DATA = glimmer.ps $(pdf_doc)
DIRS = ug num

subdirs:	$(DIRS)
$(DIRS):
		$(MAKE) -C $@

if MAINTAINER_MODE
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
endif
clean:
	rm -f *~ *.aux *.dvi *.log *.out *.toc

.PHONY:		clean install subdirs $(DIRS)
