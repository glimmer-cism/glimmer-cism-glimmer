EXTRA_DIST = eis_varlist.tex glide_varlist.tex glint_varlist.tex intro.tex paramfiles.tex varlist.tex drivers.tex  example.tex \
             glide.tex glimmer_api.tex installation.tex ug.tex varlist.tex.in

BUILT_SOURCES = glide_varlist.tex eis_varlist.tex glint_varlist.tex
CLEANFILES = $(BUILT_SOURCES)

SUBDIRS = figs

glide_varlist.tex:	../../src/fortran/glide_vars.def varlist.tex.in
			$(PYTHON) ../../src/python/generate_ncvars.py $^
eis_varlist.tex:	../../src/fortran/eis_vars.def varlist.tex.in
			$(PYTHON) ../../src/python/generate_ncvars.py $^
glint_varlist.tex:	../../src/fortran/glint_vars.def varlist.tex.in
			$(PYTHON) ../../src/python/generate_ncvars.py $^
