#
THIS_FILE=Makefile
ARC_FILE=../../makefile.arc
#
#
include $(ARC_FILE)
#
GLIMMERFLAGS=$(F77FLAGS) $(GLIMMERPRECOPT)
#
SUBDIRS = SLAP_library
LIBS = libglimmer.a libglide.a libglint.a libclimate.a
PROGS = glimmer_example simple_glide eis_glide
#relaxed, bin2ncdf, input2ncdf are currently broken...
TESTS = test_config test_setup test_ts
DATA = ncdf_template.f90.in eis_vars.def  glide_vars.def  glint_vars.def

AUTOGEN = glide_io.f90 eis_io.f90 glint_io.f90

SOURCES = *.f90

EIS = eis_forcing.o eis_ela.o eis_temp.o eis_slc.o eis_cony.o glimmer_searchcircle.o eis_io.o eis_types.o

SIMPLE = simple_forcing.o

ISOSTASY = isostasy.o isostasy_el.o isostasy_setup.o isostasy_types.o

INCLUDE_FILES=

all:	$(LIBS) $(TESTS) $(PROGS)


.PHONY:		clean all subdirs $(SUBDIRS)

subdirs:	$(SUBDIRS)

$(SUBDIRS):
		$(MAKE) -C $@


libglimmer.a: 	glimmer_CFproj.o glimmer_config.o glimmer_global.o glimmer_utils.o glimmer_log.o glimmer_ts.o\
		glimmer_ncdf.o glimmer_ncio.o glimmer_ncparms.o glimmer_paramets.o glimmer_physcon.o glimmer_degd.o gmt.o kelvin.o
		ar rcv $@ $^
		if [ "$(RANLIB)" != "" ] ; then $(RANLIB) $@ ; fi

libglide.a:	glide.o glide_setup.o glide_types.o glimmer_temp.o glimmer_thck.o glimmer_velo.o glide_mask.o glide_stop.o glide_io.o glide_nc_custom.o $(ISOSTASY)
		ar rcv $@ $^
		if [ "$(RANLIB)" != "" ] ; then $(RANLIB) $@ ; fi

libglint.a:	glint_climate.o glint_initialise.o glint_io.o glint_timestep.o glint_type.o \
		glimmer.o glimmer_proj.o glimmer_interp.o glimmer_global_grid.o glimmer_mbal.o
		ar rcv $@ $^
		if [ "$(RANLIB)" != "" ] ; then $(RANLIB) $@ ; fi

libclimate.a:	$(SIMPLE) $(EIS)
		ar rcv $@ $^
		if [ "$(RANLIB)" != "" ] ; then $(RANLIB) $@ ; fi

eis_glide:	$(EIS) eis_glide.o libglide.a libglimmer.a SLAP_library
		$(F77) -o $@ $(EIS) eis_glide.o libglide.a libglimmer.a SLAP_library/libslap.a $(LD_FLAGS)

simple_glide:	$(SIMPLE) simple_glide.o libglide.a libglimmer.a SLAP_library
		$(F77) -o $@ $(SIMPLE) simple_glide.o libglide.a libglimmer.a SLAP_library/libslap.a $(LD_FLAGS)

glimmer_example :  glimmer_example.o libglide.a libglimmer.a libglint.a SLAP_library
	$(F77) $(GLIMMERFLAGS) -o $@ $< libglint.a libglide.a libglimmer.a SLAP_library/libslap.a $(LD_FLAGS)

bin2ncdf: bin2ncdf.o libglide.a libglimmer.a
	$(F77) -o $@ $^ $(LD_FLAGS) 

input2ncdf:	input2ncdf.o libglide.a libglimmer.a
	$(F77) -o $@ $^ $(LD_FLAGS)

relaxed: relaxed.o glimmer_global.o ncdf.o ncdf_file.o ncdf_infile.o glide_types.o glide_setup.o glimmer_physcon.o glimmer_paramets.o \
         glimmer_config.o ncdf_params.o glimmer_outp.o glimmer_CFproj.o 
	$(F77) -o $@ $^ $(LD_FLAGS) SLAP_library/libslap.a

test_config:	test_config.o glimmer_config.o glimmer_log.o
	$(F77) -o $@ $^ $(LD_FLAGS)

test_setup:	test_setup.o libglide.a libglimmer.a
		$(F77) -o $@ $^ $(LD_FLAGS) 

test_ts:	glimmer_ts.o glimmer_log.o test_ts.o
		$(F77) -o $@ $^ $(LD_FLAGS)

# install glimmer
install:	install-mod install-libs install-progs install-data
		for dir in $(SUBDIRS); do \
		  $(MAKE) -C $$dir install; \
		done

install-data:	$(DATA)
		install -d $(GLIMMER_PREFIX)/share/glimmer
		cp $(DATA) $(GLIMMER_PREFIX)/share/glimmer

install-mod:	$(OBJECTS)
		install -d $(GLIMMER_PREFIX)/mod
		cp *.mod $(GLIMMER_PREFIX)/mod
install-libs:	$(LIBS)
		install -d $(GLIMMER_PREFIX)/lib
		cp $(LIBS) $(GLIMMER_PREFIX)/lib
install-progs:	$(PROGS)
		install -d $(GLIMMER_PREFIX)/bin
		install $(PROGS) $(GLIMMER_PREFIX)/bin

# clean up
clean : 
	rm -f *.o *.mod f90_dep.mak $(AUTOGEN) $(PROGS) $(LIBS)
	for dir in $(SUBDIRS); do \
	  $(MAKE) -C $$dir clean; \
	done

f90_dep.mak: *.f90 $(AUTOGEN)
	python ../python/f90_dependencies.py -o $@ $^

#some special rules for automatically generating ncdf src code
glide_io.f90:	glide_vars.def ncdf_template.f90.in
		python ../python/generate_ncvars.py $^
eis_io.f90:	eis_vars.def ncdf_template.f90.in
		python ../python/generate_ncvars.py $^
glint_io.f90:	glint_vars.def ncdf_template.f90.in
		python ../python/generate_ncvars.py $^

include f90_dep.mak 
