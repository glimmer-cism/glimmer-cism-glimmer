SUBDIRS = SLAP_library
RESTART_STEM=glimmer_rst
RESTART_INC_INS=$(RESTART_STEM)_body.inc.in $(RESTART_STEM)_head.inc.in
RESTART_INCS=$(RESTART_INC_INS:.inc.in=.inc)
RESTART_COMMON=glimmer_restart_common.F90 glimmer_restart_pointarr.F90 \
               glimmer_restart_statarr.F90 glimmer_restart_statscal.F90

lib_LIBRARIES = libglide.a libglimmer.a libglint.a libclimate.a
bin_PROGRAMS  =  simple_glide eis_glide glint_example glex_ebm eismint3_glide nc2config \
                 spin_glide glide_thermal_equib
check_PROGRAMS = test_config test_setup test_ts test_lithot test_integrate test_commandline test_writestats \
                 glimmer_unittest

variable_DATA = glide_vars.def eis_vars.def glint_vars.def time_vars.def \
                glint_mbal_vars.def eismint3_vars.def spin_vars.def ncdf_template.F90.in

variabledir = $(pkgdatadir)
dist_variable_DATA = $(variable_DATA)

nodist_include_HEADERS := $(wildcard *.mod)
EXTRA_DIST = f90_dep.mak $(RESTART_INC_INS)
BUILT_SOURCES = glide_io.F90 eis_io.F90 spin_io.F90 glint_io.F90 glint_mbal_io.F90 \
                eismint3_io.F90 $(RESTART_INCS)
CLEANFILES = $(BUILT_SOURCES) glimmer_vers.F90 f90_dep.mak *.mod
RESTART_CODE = glide_types.F90 glimmer_config.F90 glimmer_coordinates.F90 \
                  glimmer_daily_pdd.F90 glimmer_map_types.F90 glimmer_ncdf.F90 \
                  glimmer_pdd.F90 glimmer_searchcircle.F90 glimmer_sparse.F90 \
                  glimmer_sparse_slap.F90 \
                  glimmer_ts.F90 glint_main.F90 glint_global_grid.F90 glint_interp.F90 \
                  glint_mbal_coupling.F90 glint_mbal.F90 glint_smb.F90 glint_type.F90 \
                  isostasy_types.F90 profile.F90 simple_forcing.F90 glimmer_paramets.F90 \
                  glimmer_anomcouple.F90 glint_mpinterp.F90
libglimmer_a_SOURCES = glimmer_vers.F90 glimmer_commandline.F90 \
                       glimmer_config.F90 glimmer_global.F90 glimmer_utils.F90\
                       glimmer_log.F90 glimmer_ts.F90 glimmer_ncdf.F90 \
                       glimmer_ncio.F90 glimmer_ncparams.F90 \
                       glimmer_paramets.F90 glimmer_physcon.F90 kelvin.F90 \
                       profile.F90 glimmer_pdd.F90 glimmer_routing.F90 \
                       glimmer_daily_pdd.F90 glint_global_interp.F90 \
                       glimmer_sparse_type.F90 glimmer_sparse.F90 \
                       glimmer_sparse_umfpack.F90 \
                       glimmer_sparse_slap.F90 \
		       glimmer_coordinates.F90 \
                       glimmer_integrate.F90 ncdf_utils.F90 \
                       glimmer_map_CFproj.F90 glimmer_map_proj4.F90 \
                       glimmer_map_types.F90 glimmer_map_init.F90 \
                       glimmer_map_trans.F90 glimmer_restart_common.F90 \
                       glimmer_restart_statscal.F90 glimmer_restart_statarr.F90 \
                       glimmer_restart_pointarr.F90 glimmer_anomcouple.F90 \
                       $(RESTART_INCS) glimmer_restart.F90 glimmer_filenames.F90 \
                       glimmer_writestats.F90 writestats.c writestats.h
nodist_libglimmer_a_SOURCES = glimmer_vers.F90

libglide_a_SOURCES = glide.F90 glide_setup.F90 glide_types.F90 glide_temp.F90 \
                     glide_bwater.F90 glide_deriv.F90 xls.F90 ice3d_lib.F90 \
                     glide_velo_higher.F90 glide_thck.F90 glide_velo.F90 \
                     glide_mask.F90 glide_stop.F90 glide_io.F90 \
                     glide_nc_custom.F90 isostasy.F90 isostasy_el.F90\
                     isostasy_setup.F90 isostasy_types.F90 glide_lithot.F90\
                     glide_lithot3d.F90 glide_lithot1d.F90 glide_profile.F90\
                     glide_diagnostics.F90 glissade.F90 glissade_remap.F90\
                     glissade_velo.F90 glissade_constants.F90 glide_vertint.F90\
                     glide_thckmask.F90 glide_nonlin.F90 glide_grids.F90\
                     glam.F90 glam_strs2.F90 glam_thck_ppm.F90\
                     remap_advection.F90 remap_glamutils.F90 glide_ground.F90
# *sfp* Added the last two lines above for addition of Payne/Price HO dynamics and  
# *sfp* LANL inc. remapping for dH/dt calculation

libglint_a_SOURCES = glint_climate.F90 glint_initialise.F90 glint_io.F90 \
                     glint_timestep.F90 glint_type.F90 glint_main.F90 \
                     glint_interp.F90 glint_global_grid.F90 \
                     glint_precip_param.F90 glint_mbal.F90 glint_constants.F90\
                     glint_mbal_coupling.F90 glint_mbal_io.F90 \
                     glint_example_clim.F90 glint_smb.F90 glex_ebm_clim.F90 \
                     glint_mpinterp.F90 glint_commandline.F90
libclimate_a_SOURCES = eismint3_forcing.F90 eismint3_types.F90 \
                       simple_forcing.F90 eis_forcing.F90 eis_ela.F90 \
                       eis_temp.F90 eis_slc.F90 eis_cony.F90 \
                       glimmer_searchcircle.F90 eis_io.F90 eis_types.F90 \
                       eismint3_io.F90 spin_forcing.F90 spin_mb.F90 \
                       spin_temp.F90 spin_slc.F90 spin_io.F90 \
                       testsfg.F90 spin_types.F90 
if PARDISO
RESTART_CODE += glimmer_sparse_pardiso.F90
libglimmer_a_SOURCES += glimmer_sparse_pardiso.F90
endif PARDISO

eismint3_glide_SOURCES = eismint3_glide.F90
eismint3_glide_LDADD = libclimate.a libglide.a libglimmer.a \
                       SLAP_library/libslap.a $(BLAS_LIBS) $(NETCDF_LDFLAGS) $(NETCDF_LIBS)
simple_glide_SOURCES = simple_glide.F90 
simple_glide_LDADD = libclimate.a libglide.a libglimmer.a \
                     SLAP_library/libslap.a $(BLAS_LIBS) $(NETCDF_LDFLAGS) $(NETCDF_LIBS)
eis_glide_SOURCES = eis_glide.F90
eis_glide_LDADD = libclimate.a libglide.a libglimmer.a SLAP_library/libslap.a $(BLAS_LIBS) $(NETCDF_LDFLAGS) $(NETCDF_LIBS)
spin_glide_SOURCES = spin_glide.F90
spin_glide_LDADD = libclimate.a libglide.a libglimmer.a SLAP_library/libslap.a
glint_example_SOURCES = glint_example.F90
glint_example_LDADD = libglint.a libglide.a libglimmer.a SLAP_library/libslap.a $(BLAS_LIBS) $(NETCDF_LDFLAGS) $(NETCDF_LIBS)
glex_ebm_SOURCES = glex_ebm.F90
glex_ebm_LDADD = libglint.a libglide.a libglimmer.a SLAP_library/libslap.a $(BLAS_LIBS) $(NETCDF_LDFLAGS) $(NETCDF_LIBS)
nc2config_SOURCES = nc2config.F90
nc2config_LDADD = libglimmer.a $(NETCDF_LDFLAGS) $(NETCDF_LIBS)

test_config_SOURCES = test_config.F90 glimmer_config.F90 glimmer_log.F90
test_setup_SOURCES = test_setup.F90
test_setup_LDADD = libglide.a libglimmer.a $(NETCDF_LDFLAGS) $(NETCDF_LIBS)
test_ts_SOURCES = glimmer_ts.F90 glimmer_log.F90 test_ts.F90
test_lithot_SOURCES = test_lithot.F90
test_lithot_LDADD = libclimate.a libglide.a libglimmer.a SLAP_library/libslap.a $(BLAS_LIBS) $(NETCDF_LDFLAGS) $(NETCDF_LIBS)
test_integrate_SOURCES = test_integrate.F90 glimmer_global.F90 \
                         glimmer_physcon.F90 glimmer_integrate.F90
test_commandline_SOURCES = test_commandline.F90 glimmer_commandline.F90
test_writestats_SOURCES = test_writestats.F90 glimmer_writestats.F90 writestats.c writestats.h

glimmer_unittest_SOURCES = unittest_framework.F90 unittest.F90 unittest_glide_deriv.F90 \
                           unittest_ice3d_lib.F90 unittest_pbc.F90
glimmer_unittest_LDADD = libglide.a libglimmer.a SLAP_library/libslap.a

glide_thermal_equib_SOURCES = thermal_equib.F90
glide_thermal_equib_LDADD = libclimate.a libglide.a libglimmer.a SLAP_library/libslap.a

ALL_SRC=$(test_integrage_SOURCES) $(test_lithot_SOURCES) $(test_ts_SOURCES) $(test_setup_SOURCES) \
        $(test_config_SOURCES) $(nc2config_SOURCES) $(glex_ebm_SOURCES) $(glint_example_SOURCES) \
	$(eis_glide_SOURCES) $(simple_glide_SOURCES) $(eismint3_glide_SOURCES) $(libglimmer_a_SOURCES) \
	$(libglide_a_SOURCES) $(libglint_a_SOURCES) $(libclimate_a_SOURCES) $(spin_glide_SOURCES) \
        $(glimmer_unittest_SOURCES)

AM_FCFLAGS = $(NETCDF_CPPFLAGS) $(HDF5_CPPFLAGS)

LDADD = $(NETCDF_LDFLAGS) $(NETCDF_LIBS)

AM_CFLAGS=$(CFORTRAN_CFLAGS)

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
spin_io.F90:	spin_vars.def ncdf_template.F90.in
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

# Originally, the makefile dependency builder used all files in $(srcdir)/*.F90
# However, this doesn't work now that mutliple files (e.g. multiple sparse solvers)
# could define the same module.
# So, we give it all of the sources used to build Glimmer
f90_dep.mak:	$(ALL_SRC) $(BUILT_SOURCES) glimmer_vers.F90
		$(PYTHON) $(top_srcdir)/src/python/f90_dependencies.py -o $@ $^

-include f90_dep.mak
