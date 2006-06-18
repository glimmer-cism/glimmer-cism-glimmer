SUBDIRS = SLAP_library

lib_LIBRARIES = libglimmer.a libglide.a libglint.a libclimate.a
bin_PROGRAMS  =  simple_glide eis_glide glint_example glex_ebm eismint3_glide
check_PROGRAMS = test_config test_setup test_ts test_lithot test_integrate

variable_DATA = glide_vars.def eis_vars.def glint_vars.def time_vars.def \
                glint_mbal_vars.def eismint3_vars.def ncdf_template.F90.in
variabledir = $(pkgdatadir)
dist_variable_DATA = $(variable_DATA)

nodist_include_HEADERS := $(wildcard *.mod)
EXTRA_DIST = f90_dep.mak
BUILT_SOURCES = glide_io.F90 eis_io.F90 glint_io.F90 glint_mbal_io.F90 \
                eismint3_io.F90
CLEANFILES = $(BUILT_SOURCES) glimmer_version.F90 f90_dep.mak

libglimmer_a_SOURCES = glimmer_version.F90 glimmer_CFproj.F90 \
                       glimmer_config.F90 glimmer_global.F90 glimmer_utils.F90\
                       glimmer_log.F90 glimmer_ts.F90 glimmer_ncdf.F90 \
                       glimmer_ncio.F90 glimmer_ncparms.F90 \
                       glimmer_paramets.F90 glimmer_physcon.F90 kelvin.F90 \
                       profile.F90 glimmer_pdd.F90 glimmer_routing.F90 \
                       glimmer_daily_pdd.F90 glint_global_interp.F90 \
                       glimmer_sparse.F90 glimmer_coordinates.F90 \
                       glimmer_integrate.F90 ncdf_utils.F90
nodist_libglimmer_a_SOURCES = glimmer_version.F90
libglide_a_SOURCES = glide.F90 glide_setup.F90 glide_types.F90 glide_temp.F90 \
                     glide_thck.F90 glide_velo.F90 glide_mask.F90 \
                     glide_stop.F90 glide_io.F90 glide_nc_custom.F90 \
                     isostasy.F90 isostasy_el.F90 isostasy_setup.F90 \
                     isostasy_types.F90 glide_lithot.F90 glide_lithot3d.F90 \
                     glide_lithot1d.F90 glide_profile.F90
libglint_a_SOURCES = glint_climate.F90 glint_initialise.F90 glint_io.F90 \
                     glint_timestep.F90 glint_type.F90 glint.F90 \
                     glint_proj.F90 glint_interp.F90 glint_global_grid.F90 \
                     glint_precip_param.F90 glint_mbal.F90 glint_constants.F90\
                     glint_gmt.F90 glint_mbal_coupling.F90 glint_mbal_io.F90 \
                     glint_example_clim.F90 glint_smb.F90 glex_ebm_clim.F90
libclimate_a_SOURCES = eismint3_forcing.F90 eismint3_types.F90 \
                       simple_forcing.F90 eis_forcing.F90 eis_ela.F90 \
                       eis_temp.F90 eis_slc.F90 eis_cony.F90 \
                       glimmer_searchcircle.F90 eis_io.F90 eis_types.F90 \
                       eismint3_io.F90

eismint3_glide_SOURCES = eismint3_glide.F90
eismint3_glide_LDADD = libclimate.a libglide.a libglimmer.a \
                       SLAP_library/libslap.a
simple_glide_SOURCES = simple_glide.F90 
simple_glide_LDADD = libclimate.a libglide.a libglimmer.a \
                     SLAP_library/libslap.a
eis_glide_SOURCES = eis_glide.F90
eis_glide_LDADD = libclimate.a libglide.a libglimmer.a SLAP_library/libslap.a
glint_example_SOURCES = glint_example.F90
glint_example_LDADD = libglint.a libglide.a libglimmer.a SLAP_library/libslap.a
glex_ebm_SOURCES = glex_ebm.F90
glex_ebm_LDADD = libglint.a libglide.a libglimmer.a SLAP_library/libslap.a

test_config_SOURCES = test_config.F90 glimmer_config.F90 glimmer_log.F90
test_setup_SOURCES = test_setup.F90
test_setup_LDADD = libglide.a libglimmer.a
test_ts_SOURCES = glimmer_ts.F90 glimmer_log.F90 test_ts.F90
test_lithot_SOURCES = test_lithot.F90
test_lithot_LDADD = libclimate.a libglide.a libglimmer.a SLAP_library/libslap.a
test_integrate_SOURCES = test_integrate.F90 glimmer_global.F90 \
                         glimmer_physcon.F90 glimmer_integrate.F90

#some special rules for automatically generating ncdf src code
glide_io.F90:	glide_vars.def ncdf_template.F90.in
		$(PYTHON) $(top_srcdir)/src/python/generate_ncvars.py $^
eis_io.F90:	eis_vars.def ncdf_template.F90.in
		$(PYTHON) $(top_srcdir)/src/python/generate_ncvars.py $^
eismint3_io.F90: eismint3_vars.def ncdf_template.F90.in
		$(PYTHON) $(top_srcdir)/src/python/generate_ncvars.py $^
glint_io.F90:	glint_vars.def ncdf_template.F90.in
		$(PYTHON) $(top_srcdir)/src/python/generate_ncvars.py $^
glint_mbal_io.F90:	glint_mbal_vars.def ncdf_template.F90.in
		$(PYTHON) $(top_srcdir)/src/python/generate_ncvars.py $^

# handling CVS versioning
if HAVE_CVS
glimmer_version.F90:	$(srcdir)/CVS/Entries glimmer_version.F90.in 
		$(PYTHON) $(top_srcdir)/src/python/cvs_version.py \
                $(srcdir)/CVS/Entries glimmer_version.F90.in $@
else
glimmer_version.F90:	glimmer_version.F90.in
		cp $< $@
endif


f90_dep.mak:	$(srcdir)/*.F90 $(BUILT_SOURCES) glimmer_version.F90
		$(PYTHON) $(top_srcdir)/src/python/f90_dependencies.py -o $@ $^

-include f90_dep.mak
