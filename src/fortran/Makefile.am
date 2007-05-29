SUBDIRS = SLAP_library
RESTART_STEM=glimmer_rst
RESTART_INC_INS=$(RESTART_STEM)_body.inc.in $(RESTART_STEM)_head.inc.in
RESTART_INCS=$(RESTART_INC_INS:.inc.in=.inc)
RESTART_COMMON=glimmer_restart_common.F90 glimmer_restart_pointarr.F90 \
               glimmer_restart_statarr.F90 glimmer_restart_statscal.F90

lib_LIBRARIES = libglimmer.a libglide.a libglint.a libclimate.a
bin_PROGRAMS  =  simple_glide eis_glide glint_example glex_ebm eismint3_glide nc2config
check_PROGRAMS = test_config test_setup test_ts test_lithot test_integrate

variable_DATA = glide_vars.def eis_vars.def glint_vars.def time_vars.def \
                glint_mbal_vars.def eismint3_vars.def ncdf_template.F90.in
variabledir = $(pkgdatadir)
dist_variable_DATA = $(variable_DATA)

nodist_include_HEADERS := $(wildcard *.mod)
EXTRA_DIST = f90_dep.mak $(RESTART_INC_INS)
BUILT_SOURCES = glide_io.F90 eis_io.F90 glint_io.F90 glint_mbal_io.F90 \
                eismint3_io.F90 $(RESTART_INCS)
CLEANFILES = $(BUILT_SOURCES) glimmer_vers.F90 f90_dep.mak
RESTART_CODE = glide_types.F90 glimmer_config.F90 glimmer_coordinates.F90 \
                  glimmer_daily_pdd.F90 glimmer_map_types.F90 glimmer_ncdf.F90 \
                  glimmer_pdd.F90 glimmer_searchcircle.F90 glimmer_sparse.F90 \
                  glimmer_ts.F90 glint_main.F90 glint_global_grid.F90 glint_interp.F90 \
                  glint_mbal_coupling.F90 glint_mbal.F90 glint_smb.F90 glint_type.F90 \
                  isostasy_types.F90 profile.F90 simple_forcing.F90 glimmer_paramets.F90 \
                  glimmer_anomcouple.F90
libglimmer_a_SOURCES = glimmer_vers.F90 \
                       glimmer_config.F90 glimmer_global.F90 glimmer_utils.F90\
                       glimmer_log.F90 glimmer_ts.F90 glimmer_ncdf.F90 \
                       glimmer_ncio.F90 glimmer_ncparams.F90 \
                       glimmer_paramets.F90 glimmer_physcon.F90 kelvin.F90 \
                       profile.F90 glimmer_pdd.F90 glimmer_routing.F90 \
                       glimmer_daily_pdd.F90 glint_global_interp.F90 \
                       glimmer_sparse.F90 glimmer_coordinates.F90 \
                       glimmer_integrate.F90 ncdf_utils.F90 \
                       glimmer_map_CFproj.F90 glimmer_map_proj4.F90 \
                       glimmer_map_types.F90 glimmer_map_init.F90 \
                       glimmer_map_trans.F90 glimmer_restart_common.F90 \
                       glimmer_restart_statscal.F90 glimmer_restart_statarr.F90 \
                       glimmer_restart_pointarr.F90 glimmer_anomcouple.F90 \
                       $(RESTART_INCS) glimmer_restart.F90 glimmer_filenames.F90
nodist_libglimmer_a_SOURCES = glimmer_vers.F90
libglide_a_SOURCES = glide.F90 glide_setup.F90 glide_types.F90 glide_temp.F90 \
                     glide_thck.F90 glide_velo.F90 glide_mask.F90 \
                     glide_stop.F90 glide_io.F90 glide_nc_custom.F90 \
                     isostasy.F90 isostasy_el.F90 isostasy_setup.F90 \
                     isostasy_types.F90 glide_lithot.F90 glide_lithot3d.F90 \
                     glide_lithot1d.F90 glide_profile.F90 glide_diagnostics.F90 \
                     glissade.F90 glissade_remap.F90 glissade_constants.F90
libglint_a_SOURCES = glint_climate.F90 glint_initialise.F90 glint_io.F90 \
                     glint_timestep.F90 glint_type.F90 glint_main.F90 \
                     glint_interp.F90 glint_global_grid.F90 \
                     glint_precip_param.F90 glint_mbal.F90 glint_constants.F90\
                     glint_mbal_coupling.F90 glint_mbal_io.F90 \
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
nc2config_SOURCES = nc2config.F90
nc2config_LDADD = libglimmer.a

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

# Pattern rule for restart include files
if BUILD_RESTARTS
$(RESTART_STEM)_body.inc: $(RESTART_STEM)_body.inc.in Makefile
		cp $< $@
$(RESTART_STEM)_head.inc: $(RESTART_STEM)_head.inc.in Makefile
		cp $< $@
$(RESTART_CODE:.F90=.o): $(RESTART_COMMON:.F90=.o)
else
$(RESTART_STEM)_body.inc $(RESTART_STEM)_head.inc: Makefile
		echo '! Dummy include file - restarts not enabled' > $@
endif

# handling CVS versioning
if HAVE_CVS
glimmer_vers.F90:	$(srcdir)/CVS/Entries glimmer_vers.F90.in 
		$(PYTHON) $(top_srcdir)/src/python/cvs_version.py \
                $(srcdir)/CVS/Entries glimmer_vers.F90.in $@
$(RESTART_STEM)_body.inc.in $(RESTART_STEM)_head.inc.in: $(RESTART_CODE) Makefile
		$(PYTHON) $(top_srcdir)/src/python/make_restarts.py -o $(RESTART_STEM) $(RESTART_CODE)
else
glimmer_vers.F90:	glimmer_vers.F90.in
		cp $< $@
endif


f90_dep.mak:	$(srcdir)/*.F90 $(BUILT_SOURCES) glimmer_vers.F90
		$(PYTHON) $(top_srcdir)/src/python/f90_dependencies.py -o $@ $^

-include f90_dep.mak
