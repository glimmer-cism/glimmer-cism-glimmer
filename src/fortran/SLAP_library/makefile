#
ARC_FILE=../../../makefile.arc
ifndef GLIMMER_PREFIX
GLIMMER_PREFIX=../../../
endif
include $(ARC_FILE)
#
LIBS = libslap.a
OBJECTS = blas.o dcgs.o dlaputil.o dmset.o dmvops.o mach.o xersla.o

INCLUDE_FILES=

OTHER_FILES=$(ARC_FILE) makefile

libslap.a: $(OBJECTS) $(OTHER_FILES) $(INCLUDE_FILES)
	ar rcv $@ $(OBJECTS)
	if [ "$(RANLIB)" != "" ] ; then $(RANLIB) libslap.a ; fi

install:	$(LIBS)
	install -d $(GLIMMER_PREFIX)/lib
	cp $(LIBS) $(GLIMMER_PREFIX)/lib

clean : 
	rm -f *.o *.a

blas.o: blas.f $(OTHER_FILES)
	$(F77) -c $(GLIMMERFLAGS) blas.f

dcgs.o: dcgs.f $(OTHER_FILES)
	$(F77) -c $(GLIMMERFLAGS) dcgs.f

dlaputil.o: dlaputil.f $(OTHER_FILES)
	$(F77) -c $(GLIMMERFLAGS) dlaputil.f

dmset.o: dmset.f $(OTHER_FILES)
	$(F77) -c $(GLIMMERFLAGS) dmset.f

dmvops.o: dmvops.f $(OTHER_FILES)
	$(F77) -c $(GLIMMERFLAGS) dmvops.f

mach.o: mach.f $(OTHER_FILES)
	$(F77) -c $(GLIMMERFLAGS) mach.f

xersla.o: xersla.f $(OTHER_FILES)
	$(F77) -c $(GLIMMERFLAGS) xersla.f


