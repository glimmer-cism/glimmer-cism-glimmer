ifndef GLIMMER_PREFIX
GLIMMER_PREFIX=..
endif

SUBDIRS = ug

all:	subdirs glimmer.pdf glimmer.ps num.ps

subdirs:	$(SUBDIRS)
$(SUBDIRS):
		$(MAKE) -C $@

dg.pdf:	dg_pdf.dvi
	dvipdfm -o dg.pdf dg_pdf.dvi

dg.ps:	dg_ps.dvi
	dvips -o dg.ps dg_ps.dvi

num.pdf:	num_pdf.dvi
	dvipdfm -o num.pdf num_pdf.dvi

num.ps:	num_ps.dvi
	dvips -o num.ps num_ps.dvi

ug.pdf:	ug ug_pdf.dvi
	dvipdfm -o ug.pdf ug_pdf.dvi

ug.ps:	ug ug_ps.dvi
	dvips -o ug.ps ug_ps.dvi

glimmer.pdf:	glimmer_pdf.dvi
	dvipdfm -o glimmer.pdf glimmer_pdf.dvi

glimmer.ps:	glimmer_ps.dvi
	dvips -o glimmer.ps glimmer_ps.dvi

code_doc.ps:	code_doc.dvi
		dvips -o code_doc.ps code_doc.dvi

code_doc.tex:	../src/fortran/*.f90
		f90docp.py -t "GLIMMER Code Documentation" -a "Glimmer Developer Team" $^ > $@

install:	glimmer.ps glimmer.pdf
		install -d $(GLIMMER_PREFIX)/share/docs/glimmer
		cp $^ $(GLIMMER_PREFIX)/share/docs/glimmer

.tex.dvi::
	latex $<
	latex $<

clean:
	rm -f *~ *.aux *.dvi *.log *.out *.toc

.PHONY:		clean install subdirs $(SUBDIRS)
