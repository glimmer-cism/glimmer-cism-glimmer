all:	glide_varlist.tex eis_varlist.tex glint_varlist.tex

glide_varlist.tex:	 ../../src/fortran/glide_vars.def varlist.tex.in
	python ../../src/python/generate_ncvars.py $^
eis_varlist.tex:	 ../../src/fortran/eis_vars.def varlist.tex.in
	python ../../src/python/generate_ncvars.py $^
glint_varlist.tex:	 ../../src/fortran/glint_vars.def varlist.tex.in
	python ../../src/python/generate_ncvars.py $^
