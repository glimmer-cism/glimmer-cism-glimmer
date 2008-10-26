# Makefile for building the example climate driver
# Magnus Hagdorn, November 2004

# setup compilation flags: tell the compiler where to look for f90 module files
GLIMMER_FFLAGS = `glimmer-config --fcflags`
# setup linker flags: tell the linker where to look for libraries and which libraries to link
GLIMMER_LIBRARIES = `glimmer-config --libs` -lglide -lglimmer -lslap -lnetcdf

# set variables FCFLAGS and LDFLAGS for extra libraries, e.g. netCDF

FC=gfortran

# how to build the example executable
example:	example.o
		$(FC) -o $@ $< $(LDFLAGS) $(GLIMMER_LIBRARIES)

# implicit rule, telling the make system how to build an object file from a .f90 file
.f90.o:
		$(FC) $(FCFLAGS) $(GLIMMER_FFLAGS) -c $< -o $@


.SUFFIXES:	.f90 .o


# stuff below here is not too interesting. You might as well ignore it...
# build a distribution
DIST_FILES = ChangeLog example.config example.f90 Makefile README
VERSION = 0.6
DIST_NAME = glimmer-example-$(VERSION)
dist::
	mkdir $(DIST_NAME)
	cp $(DIST_FILES) $(DIST_NAME)
	tar cvzf $(DIST_NAME).tar.gz $(DIST_NAME)
	rm -r $(DIST_NAME)

clean:
	rm example example.o
