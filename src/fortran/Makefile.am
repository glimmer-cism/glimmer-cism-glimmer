SUBDIRS = SLAP_library

lib_LIBRARIES = libglimmer.a libglide.a libglint.a libclimate.a
bin_PROGRAMS  =  simple_glide eis_glide glint_example
check_PROGRAMS = test_config test_setup test_ts

variable_DATA = glide_vars.def eis_vars.def glint_vars.def time_vars.def ncdf_template.f90.in
variabledir = $(pkgdatadir)
dist_variable_DATA = $(variable_DATA)

nodist_include_HEADERS = *.mod
EXTRA_DIST = f90_dep.mak
BUILT_SOURCES = glide_io.f90 eis_io.f90 glint_io.f90
CLEANFILES = $(BUILT_SOURCES) f90_dep.mak

libglimmer_a_SOURCES = glimmer_CFproj.f90 glimmer_config.f90 glimmer_global.f90 glimmer_utils.f90 glimmer_log.f90 glimmer_ts.f90 \
                glimmer_ncdf.f90 glimmer_ncio.f90 glimmer_ncparms.f90 glimmer_paramets.f90 glimmer_physcon.f90 kelvin.f90 profile.f90
libglide_a_SOURCES = glide.f90 glide_setup.f90 glide_types.f90 glide_temp.f90 glide_thck.f90 glide_velo.f90 glide_mask.f90 glide_stop.f90 glide_io.f90 \
                glide_nc_custom.f90 isostasy.f90 isostasy_el.f90 isostasy_setup.f90 isostasy_types.f90
libglint_a_SOURCES = glint_climate.f90 glint_initialise.f90 glint_io.f90 glint_timestep.f90 glint_type.f90 glimmer_routing.f90 \
		glint.f90 glint_proj.f90 glint_interp.f90 glint_global_grid.f90 glint_precip_param.f90 glide_enmabal.f90 glide_pdd.f90 \
                glint_mbal.f90 glint_constants.f90 glint_gmt.f90 
libclimate_a_SOURCES = simple_forcing.f90 eis_forcing.f90 eis_ela.f90 eis_temp.f90 eis_slc.f90 eis_cony.f90 glimmer_searchcircle.f90 eis_io.f90 eis_types.f90

simple_glide_SOURCES = simple_glide.f90 
simple_glide_LDADD = libclimate.a libglide.a libglimmer.a SLAP_library/libslap.a
eis_glide_SOURCES = eis_glide.f90
eis_glide_LDADD = libclimate.a libglide.a libglimmer.a SLAP_library/libslap.a
glint_example_SOURCES = glint_example.f90
glint_example_LDADD = libglint.a libglide.a libglimmer.a SLAP_library/libslap.a

test_config_SOURCES = test_config.f90 glimmer_config.f90 glimmer_log.f90
test_setup_SOURCES = test_setup.f90
test_setup_LDADD = libglide.a libglimmer.a
test_ts_SOURCES = glimmer_ts.f90 glimmer_log.f90 test_ts.f90

#some special rules for automatically generating ncdf src code
glide_io.f90:	glide_vars.def ncdf_template.f90.in
		$(PYTHON) ../python/generate_ncvars.py $^
eis_io.f90:	eis_vars.def ncdf_template.f90.in
		$(PYTHON) ../python/generate_ncvars.py $^
glint_io.f90:	glint_vars.def ncdf_template.f90.in
		$(PYTHON) ../python/generate_ncvars.py $^

f90_dep.mak:	*.f90 $(BUILT_SOURCES)
		$(PYTHON) ../python/f90_dependencies.py -o $@ $^

-include f90_dep.mak
