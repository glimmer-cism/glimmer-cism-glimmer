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
LIBS = libglimmer.$(LIB_EXT)
PROGS = glimmer_example bin2ncdf input2ncdf relaxed test_config

OBJECTS = glide_pdd.$(OBJ_EXT) glimmer_isot.$(OBJ_EXT) glimmer_setup.$(OBJ_EXT) \
          glimmer_modu.$(OBJ_EXT) glimmer_outp.$(OBJ_EXT) glimmer_paramets.$(OBJ_EXT) \
          glimmer_physcon.$(OBJ_EXT) glimmer_restart.$(OBJ_EXT) \
          glimmer_temp.$(OBJ_EXT) glimmer_thck.$(OBJ_EXT) glimmer_velo.$(OBJ_EXT) \
          glimmer.$(OBJ_EXT) glimmer_global.$(OBJ_EXT) glimmer_interp.$(OBJ_EXT) \
          glimmer_mbal.$(OBJ_EXT) glimmer_object.$(OBJ_EXT) glimmer_proj.$(OBJ_EXT) \
          glimmer_utils.$(OBJ_EXT) gmt.$(OBJ_EXT) glimmer_global_grid.$(OBJ_EXT) \
          ncdf.$(OBJ_EXT) ncdf_file.$(OBJ_EXT) ncdf_params.$(OBJ_EXT) \
          ncdf_infile.$(OBJ_EXT) glimmer_config.$(OBJ_EXT) glimmer_CFproj.$(OBJ_EXT) \
          glide_messages.$(OBJ_EXT) glide_enmabal.$(OBJ_EXT) glimmer_routing.$(OBJ_EXT)

SOURCES = glide_pdd.f90 glimmer_isot.f90 glimmer_setup.f90 \
          glimmer_modu.f90 glimmer_outp.f90 glimmer_paramets.f90 glimmer_physcon.f90 \
          glimmer_restart.f90 glimmer_temp.f90 glimmer_thck.f90 glimmer_velo.f90 \
          glimmer.f90 glimmer_global.f90 glimmer_interp.f90 glimmer_mbal.f90 \
          glimmer_object.f90 glimmer_proj.f90 glimmer_utils.f90 gmt.f90 \
          glimmer_global_grid.f90 glimmer_config.f90 glimmer_CFproj.f90 \
          glide_messages.f90 ncdf.f90 ncdf_file.f90 ncdf_infile.f90 ncdf_params.f90 \
          glide_enmabal.f90 glimmer_routing.f90

INCLUDE_FILES=

OTHER_FILES=$(ARC_FILE) $(THIS_FILE)

all: $(LIBS) $(PROGS)

.PHONY:		clean all subdirs $(SUBDIRS)

subdirs:	$(SUBDIRS)

$(SUBDIRS):
		$(MAKE) -C $@


libglimmer.$(LIB_EXT): $(OBJECTS) $(OTHER_FILES) $(INCLUDE_FILES)
	ar rcv $@ $(OBJECTS)
	if [ "$(RANLIB)" != "" ] ; then $(RANLIB) $@ ; fi

glimmer_example :  glimmer_example.$(OBJ_EXT) SLAP_library libglimmer.$(LIB_EXT)
	$(F77) $(GLIMMERFLAGS) $(OUT_FLAG)$@ $< libglimmer.$(LIB_EXT) \
        SLAP_library/libslap.$(LIB_EXT) $(LD_FLAGS)

bin2ncdf: bin2ncdf.$(OBJ_EXT) glimmer_global.$(OBJ_EXT)  ncdf.$(OBJ_EXT) ncdf_file.$(OBJ_EXT) \
        glimmer_modu.$(OBJ_EXT) glimmer_physcon.$(OBJ_EXT) glimmer_paramets.$(OBJ_EXT) \
        glimmer_CFproj.$(OBJ_EXT) glide_messages.$(OBJ_EXT)
	$(F77) $(OUT_FLAG)$@ $^ $(LD_FLAGS)

input2ncdf:	input2ncdf.$(OBJ_EXT) glimmer_global.$(OBJ_EXT)  ncdf.$(OBJ_EXT) ncdf_file.$(OBJ_EXT) \
        glimmer_modu.$(OBJ_EXT) glimmer_paramets.$(OBJ_EXT)  glimmer_CFproj.$(OBJ_EXT) glide_messages.$(OBJ_EXT)
	$(F77) $(OUT_FLAG)$@ $^ $(LD_FLAGS) $(GLIMMERFLAGS)

relaxed: relaxed.$(OBJ_EXT)  SLAP_library libglimmer.$(LIB_EXT)
	$(F77) $(OUT_FLAG)$@ $< libglimmer.$(LIB_EXT) SLAP_library/libslap.$(LIB_EXT) \
        $(LD_FLAGS) 

test_config:	test_config.$(OBJ_EXT) glimmer_config.$(OBJ_EXT) glide_messages.$(OBJ_EXT)
	$(F77) $(OUT_FLAG)$@ $^ $(LD_FLAGS)

# install glimmer
install:	install-mod install-libs install-progs
		for dir in $(SUBDIRS); do \
		  $(MAKE) -C $$dir install; \
		done

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
	rm -f *.$(OBJ_EXT) *.mod f90_dep.mak $(PROGS) $(LIBS)
	for dir in $(SUBDIRS); do \
	  $(MAKE) -C $$dir clean; \
	done

f90_dep.mak: *.f90
	python ../python/f90_dependencies.py -o $@ $^

#some special rules for automatically generating ncdf src code	 
#ncdf.f90:       ncdf_vars.def ncdf.f90.in	 
#	python ../python/generate_ncvars.py $< $@.in        	 
#ncdf_file.f90:  ncdf_vars.def ncdf_file.f90.in	 
#	python ../python/generate_ncvars.py $< $@.in        	 
#ncdf_infile.f90:        ncdf_vars.def ncdf_infile.f90.in	 
#	python ../python/generate_ncvars.py $< $@.in        	 
#ncdf_params.f90:        ncdf_vars.def ncdf_params.f90.in	 
#	python ../python/generate_ncvars.py $< $@.in        

include f90_dep.mak 

