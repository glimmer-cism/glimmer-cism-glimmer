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

OBJECTS = glimmer_degd.o glimmer_isot.o glimmer_setup.o \
          glimmer_modu.o glimmer_outp.o glimmer_paramets.o glimmer_physcon.o glimmer_restart.o \
          glimmer_temp.o glimmer_thck.o glimmer_velo.o glimmer.o glimmer_global.o \
          glimmer_interp.o glimmer_mbal.o glimmer_object.o glimmer_proj.o \
          glimmer_utils.o gmt.o

SOURCES = glimmer_degd.f90 glimmer_isot.f90 glimmer_setup.f90 \
          glimmer_modu.f90 glimmer_outp.f90 glimmer_paramets.f90 glimmer_physcon.f90 \
          glimmer_restart.f90 glimmer_temp.f90 glimmer_thck.f90 glimmer_velo.f90 \
          glimmer.f90 glimmer_global.f90 glimmer_interp.f90 glimmer_mbal.f90 \
          glimmer_object.f90 glimmer_proj.f90 glimmer_utils.f90 gmt.f90

ifeq "$(F77)" "ifc"
G_DEGDM=GLIMMER_DEGD.mod
G_ISOTM=GLIMMER_ISOT.mod
G_SETUM=GLIMMER_SETUP.mod
G_MODUM=GLIMMER_TYPES.mod
G_OUTPM=GLIMMER_OUTP.mod
G_PARAM=PARAMETS.mod
G_PHYSM=PHYSCON.mod
G_RESTM=GLIMMER_RESTART.mod
G_TEMPM=GLIMMER_TEMP.mod
G_THCKM=GLIMMER_THCK.mod
G_VELOM=GLIMMER_VELO.mod
GLIMMERM=GLIMMER_MAIN.mod
G_GLOBAM=GLIMMER_GLOBAL.mod
G_INTERM=GLIMMER_INTERP.mod
G_MBALM =GLIMMER_MBAL.mod
G_OBJECM=GLIMMER_OBJECT.mod
G_PROJM =GLIMMER_PROJECT.mod
G_UTILSM=GLIMMER_UTILS.mod
GMTM    =GMT.mod
else
G_DEGDM=glimmer_degd.mod
G_ISOTM=glimmer_isot.mod
G_SETUM=glimmer_setup.mod
G_MODUM=glimmer_types.mod
G_OUTPM=glimmer_outp.mod
G_PARAM=paramets.mod
G_PHYSM=physcon.mod
G_RESTM=glimmer_restart.mod
G_TEMPM=glimmer_temp.mod
G_THCKM=glimmer_thck.mod
G_VELOM=glimmer_velo.mod
GLIMMERM=glimmer_main.mod
G_GLOBAM=glimmer_global.mod
G_INTERM=glimmer_interp.mod
G_MBALM =glimmer_mbal.mod
G_OBJECM=glimmer_object.mod
G_PROJM =glimmer_project.mod
G_UTILSM=glimmer_utils.mod
GMTM    =gmt.mod
endif

INCLUDE_FILES=

OTHER_FILES=$(ARC_FILE) $(THIS_FILE)

$(LIBNAME).a: $(OBJECTS) $(OTHER_FILES) $(INCLUDE_FILES)
	if [ -f $(LIBNAME).a ] ; then \rm $(LIBNAME).a ; fi
	ar rcv $(LIBNAME).a $(OBJECTS)
	if [ "$(RANLIB)" != "" ] ; then $(RANLIB) $(LIBNAME).a ; fi
	\cp -f *.mod ../../mod
	\cp -f $(LIBNAME).a ../../lib

example : glimmer_example.f90 $(LIBNAME).a
	$(F77) -c $(GLIMMERFLAGS) glimmer_example.f90
	$(F77) $(GLIMMERFLAGS) -o glimmer_example glimmer_example.o -L. -lglimmer -lslap

clean : 
	\rm -f *.o *.a *.mod

glimmer_degd.o $(G_DEGDM): glimmer_degd.f90 $(G_MODUM) $(G_GLOBAM) $(G_PHYSM) \
                        $(G_PARAM) $(OTHER_FILES)
	$(F77) -c $(GLIMMERFLAGS) glimmer_degd.f90

glimmer_isot.o $(G_ISOTM): glimmer_isot.f90 $(G_MODUM) $(G_GLOBAM) $(G_PHYSM) \
                        $(G_PARAM) $(OTHER_FILES)
	$(F77) -c $(GLIMMERFLAGS) glimmer_isot.f90

glimmer_setup.o $(G_SETUM): glimmer_setup.f90 $(G_PHYSM) $(G_PARAM) $(G_MODUM) \
                        $(G_PROJM) $(G_OUTPM) $(G_GLOBAM) $(OTHER_FILES) $(G_TEMPM)
	$(F77) -c $(GLIMMERFLAGS) glimmer_setup.f90

glimmer_modu.o $(G_MODUM): glimmer_modu.f90 $(G_PARAM) $(G_GLOBAM) $(OTHER_FILES)
	$(F77) -c $(GLIMMERFLAGS) glimmer_modu.f90

glimmer_outp.o $(G_OUTPM): glimmer_outp.f90 $(G_MODUM) $(G_GLOBAM) $(G_PARAM) $(G_PHYSM) $(OTHER_FILES)
	$(F77) -c $(GLIMMERFLAGS) glimmer_outp.f90

glimmer_paramets.o $(G_PARAM): glimmer_paramets.f90 $(G_GLOBAM) $(G_PHYSM) $(OTHER_FILES)
	$(F77) -c $(GLIMMERFLAGS) glimmer_paramets.f90 

glimmer_physcon.o $(G_PHYSM): glimmer_physcon.f90 $(G_GLOBAM) $(OTHER_FILES)
	$(F77) -c $(GLIMMERFLAGS) glimmer_physcon.f90

glimmer_restart.o $(G_RESTM): glimmer_restart.f90 $(G_OBJECM) $(OTHER_FILES)
	$(F77) -c $(GLIMMERFLAGS) glimmer_restart.f90

glimmer_temp.o $(G_TEMPM): glimmer_temp.f90 $(G_MODUM) $(G_GLOBAM) $(G_PHYSM) $(G_PARAM) \
                        $(G_VELOM) $(G_THCKM) $(OTHER_FILES) $(G_MBALM)
	$(F77) -c $(GLIMMERFLAGS) glimmer_temp.f90

glimmer_thck.o $(G_THCKM): glimmer_thck.f90 $(G_MODUM) $(G_GLOBAM) $(G_PHYSM) $(G_PARAM) $(OTHER_FILES)
	$(F77) -c $(GLIMMERFLAGS) glimmer_thck.f90

glimmer_velo.o $(G_VELOM): glimmer_velo.f90 $(G_MODUM) $(G_GLOBAM) $(G_PHYSM) $(G_PARAM) $(G_UTILSM) $(OTHER_FILES)
	$(F77) -c $(GLIMMERFLAGS) glimmer_velo.f90

glimmer.o $(GLIMMERM): glimmer.f90 $(G_GLOBAM) $(G_OBJECM) $(G_PROJM) $(G_UTILSM) $(G_MBALM) \
                       $(G_INTERM) $(G_RESTM) $(G_PARAM) $(G_PHYSM) $(OTHER_FILES)
	$(F77) -c $(GLIMMERFLAGS) glimmer.f90

glimmer_global.o $(G_GLOBAM): glimmer_global.f90 $(OTHER_FILES)
	$(F77) -c $(GLIMMERFLAGS) glimmer_global.f90

glimmer_interp.o $(G_INTERM): glimmer_interp.f90 $(G_GLOBAM) $(G_UTILSM) $(G_PROJM) $(OTHER_FILES)
	$(F77) -c $(GLIMMERFLAGS) glimmer_interp.f90

glimmer_mbal.o $(G_MBALM): glimmer_mbal.f90 $(G_GLOBAM) $(OTHER_FILES)
	$(F77) -c $(GLIMMERFLAGS) glimmer_mbal.f90

glimmer_object.o $(G_OBJECM): glimmer_object.f90 $(G_GLOBAM) $(G_PROJM) $(G_MODUM) $(G_SETUM) \
                              $(G_TEMPM) $(G_VELOM) $(G_OUTPM) $(G_ISOTM) $(G_THCKM) $(G_INTERM) \
                              $(G_MBALM) $(G_PARAM) $(G_UTILSM) $(GMTM) $(OTHER_FILES)
	$(F77) -c $(GLIMMERFLAGS) glimmer_object.f90

glimmer_proj.o $(G_PROJM): glimmer_proj.f90 $(GMTM) $(G_GLOBAM) $(OTHER_FILES)
	$(F77) -c $(GLIMMERFLAGS) glimmer_proj.f90

glimmer_utils.o $(G_UTILSM): glimmer_utils.f90 $(G_GLOBAM) $(OTHER_FILES)
	$(F77) -c $(GLIMMERFLAGS) glimmer_utils.f90

gmt.o $(GMTM): gmt.f90 $(OTHER_FILES)
	$(F77) -c $(GLIMMERFLAGS) gmt.f90

code_doc: $(SOURCES) $(OTHER_FILES)
	f90doc -t "GLIMMER Code Documentation" -a "Ian Rutt" $(SOURCES) > ../../doc/code_doc.tex
	cd ../../doc && latex2html -local_icons -split 4 -toc_depth 3 -link 1 code_doc.tex
