dg.pdf:	dg_pdf.dvi
	dvipdfm -o dg.pdf dg_pdf.dvi

dg.ps:	dg_ps.dvi
	dvips -o dg.ps dg_ps.dvi

ug.pdf:	ug_pdf.dvi
	dvipdfm -o ug.pdf ug_pdf.dvi

ug.ps:	ug_ps.dvi
	dvips -o ug.ps ug_ps.dvi

glimmer.pdf:	glimmer_pdf.dvi
	dvipdfm -o glimmer.pdf glimmer_pdf.dvi

glimmer.ps:	glimmer_ps.dvi
	dvips -o glimmer.ps glimmer_ps.dvi

ug_ps.dvi:	ug_ps.tex
	cd ug; $(MAKE)
	latex $<
	latex $<
ug_pdf.dvi:	ug_pdf.tex
	cd ug; $(MAKE)
	latex $<
	latex $<

.tex.dvi::
	latex $<
	latex $<

clean:
	rm -f *~ *.aux *.dvi *.log *.out *.toc
