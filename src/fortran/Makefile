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
LIBS = libglimmer.a
PROGS = glimmer_example bin2ncdf input2ncdf relaxed
DATA = ncdf_vars.def

AUTOGEN = ncdf.f90 ncdf_file.f90 ncdf_infile.f90 ncdf_params.f90

OBJECTS = glimmer_degd.o glimmer_isot.o glimmer_setup.o \
          glimmer_modu.o glimmer_outp.o glimmer_paramets.o glimmer_physcon.o glimmer_restart.o \
          glimmer_temp.o glimmer_thck.o glimmer_velo.o glimmer.o glimmer_global.o \
          glimmer_interp.o glimmer_mbal.o glimmer_object.o glimmer_proj.o \
          glimmer_utils.o gmt.o glimmer_global_grid.o ncdf.o ncdf_file.o ncdf_params.o \
          ncdf_infile.o glimmer_config.o glimmer_CFproj.o glide_messages.o

SOURCES = glimmer_degd.f90 glimmer_isot.f90 glimmer_setup.f90 \
          glimmer_modu.f90 glimmer_outp.f90 glimmer_paramets.f90 glimmer_physcon.f90 \
          glimmer_restart.f90 glimmer_temp.f90 glimmer_thck.f90 glimmer_velo.f90 \
          glimmer.f90 glimmer_global.f90 glimmer_interp.f90 glimmer_mbal.f90 \
          glimmer_object.f90 glimmer_proj.f90 glimmer_utils.f90 gmt.f90 \
          glimmer_global_grid.f90 glimmer_config.f90 glimmer_CFproj.f90 \
          glide_messages.f90 $(AUTOGEN)

INCLUDE_FILES=

OTHER_FILES=$(ARC_FILE) $(THIS_FILE)

all: $(LIBS) $(PROGS)

.PHONY:		clean all subdirs $(SUBDIRS)

subdirs:	$(SUBDIRS)

$(SUBDIRS):
		$(MAKE) -C $@


libglimmer.a: $(OBJECTS) $(OTHER_FILES) $(INCLUDE_FILES)
	ar rcv $@ $(OBJECTS)
	if [ "$(RANLIB)" != "" ] ; then $(RANLIB) $@ ; fi

glimmer_example :  glimmer_example.o SLAP_library libglimmer.a
	$(F77) $(GLIMMERFLAGS) -o $@ $< libglimmer.a SLAP_library/libslap.a $(LD_FLAGS)

bin2ncdf: bin2ncdf.o glimmer_global.o  ncdf.o ncdf_file.o glimmer_modu.o glimmer_physcon.o glimmer_paramets.o glimmer_CFproj.o glide_messages.o
	$(F77) -o $@ $^ $(LD_FLAGS) 

input2ncdf:	input2ncdf.o glimmer_global.o  ncdf.o ncdf_file.o glimmer_modu.o glimmer_paramets.o  glimmer_CFproj.o glide_messages.o
	$(F77) -o $@ $^ $(LD_FLAGS)

relaxed: relaxed.o  SLAP_library libglimmer.a
	$(F77) -o $@ $< libglimmer.a SLAP_library/libslap.a $(LD_FLAGS) 

test_config:	test_config.o glimmer_config.o
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
ncdf.f90:	ncdf_vars.def ncdf.f90.in
		python ../python/generate_ncvars.py $< $@.in	
ncdf_file.f90:	ncdf_vars.def ncdf_file.f90.in
		python ../python/generate_ncvars.py $< $@.in	
ncdf_infile.f90:	ncdf_vars.def ncdf_infile.f90.in
		python ../python/generate_ncvars.py $< $@.in	
ncdf_params.f90:	ncdf_vars.def ncdf_params.f90.in
		python ../python/generate_ncvars.py $< $@.in	

include f90_dep.mak 
