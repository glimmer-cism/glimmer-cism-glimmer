#
THIS_FILE=Makefile
ARC_FILE=../../makefile.arc
#
#
include $(ARC_FILE)
#
GLIMMERFLAGS=$(F77FLAGS) $(GLIMMERPRECOPT)
#
LIBNAME = libglimmer

AUTOGEN = ncdf.f90 ncdf_file.f90 ncdf_infile.f90 ncdf_params.f90

OBJECTS = glimmer_degd.o glimmer_isot.o glimmer_setup.o \
          glimmer_modu.o glimmer_outp.o glimmer_paramets.o glimmer_physcon.o glimmer_restart.o \
          glimmer_temp.o glimmer_thck.o glimmer_velo.o glimmer.o glimmer_global.o \
          glimmer_interp.o glimmer_mbal.o glimmer_object.o glimmer_proj.o \
          glimmer_utils.o gmt.o ncdf.o ncdf_file.o ncdf_params.o ncdf_infile.o glimmer_config.o

SOURCES = glimmer_degd.f90 glimmer_isot.f90 glimmer_setup.f90 \
          glimmer_modu.f90 glimmer_outp.f90 glimmer_paramets.f90 glimmer_physcon.f90 \
          glimmer_restart.f90 glimmer_temp.f90 glimmer_thck.f90 glimmer_velo.f90 \
          glimmer.f90 glimmer_global.f90 glimmer_interp.f90 glimmer_mbal.f90 \
          glimmer_object.f90 glimmer_proj.f90 glimmer_utils.f90 gmt.f90 glimmer_config.f90 $(AUTOGEN)

INCLUDE_FILES=

OTHER_FILES=$(ARC_FILE) $(THIS_FILE)

$(LIBNAME).a: $(OBJECTS) $(OTHER_FILES) $(INCLUDE_FILES)
	if [ -f $(LIBNAME).a ] ; then \rm $(LIBNAME).a ; fi
	ar rcv $(LIBNAME).a $(OBJECTS)
	if [ "$(RANLIB)" != "" ] ; then $(RANLIB) $(LIBNAME).a ; fi
	\cp -f *.mod ../../mod
	\cp -f $(LIBNAME).a ../../lib

example :  $(LIBNAME).a glimmer_example.o
	$(F77) $(GLIMMERFLAGS) -o glimmer_example glimmer_example.o -L. -lglimmer -lslap $(LD_FLAGS)

bin2ncdf: bin2ncdf.o glimmer_global.o  ncdf.o ncdf_file.o glimmer_modu.o glimmer_physcon.o glimmer_paramets.o
	$(F77) -o $@ $^ $(LD_FLAGS) 

input2ncdf:	input2ncdf.o glimmer_global.o  ncdf.o ncdf_file.o glimmer_modu.o
	$(F77) -o $@ $^ $(LD_FLAGS)

test_config:	test_config.o glimmer_config.o
	$(F77) -o $@ $^ $(LD_FLAGS)

clean : 
	\rm -f *.o *.a *.mod f90_dep.mak $(AUTOGEN)

f90_dep.mak: *.f90 $(AUTOGEN)
	python ../python/f90_dependencies.py -o $@ $^

code_doc: $(SOURCES) $(OTHER_FILES)
	f90doc -t "GLIMMER Code Documentation" -a "Ian Rutt" $(SOURCES) > ../../doc/code_doc.tex
	cd ../../doc && latex2html -local_icons -split 4 -toc_depth 3 -link 1 code_doc.tex

#some special rules for automatically generating ncdf src code
ncdf.f90:	ncdf_vars.def ncdf.f90.in
		python ../python/generate_ncvars.py $< $@.in	
ncdf_file.f90:	ncdf_vars.def ncdf_file.f90.in
		python ../python/generate_ncvars.py $< $@.in	
ncdf_infile.f90:	ncdf_vars.def ncdf_infile.f90.in
		python ../python/generate_ncvars.py $< $@.in	
ncdf_params.f90:	ncdf_vars.def ncdf_params.f90.in
		python ../python/generate_ncvars.py $< $@.in	

include f90_dep.mak 

