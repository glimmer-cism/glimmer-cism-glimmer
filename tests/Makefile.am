EXTRA_DIST = e2.lin.config  e2.nonlin.config

if RUN_TESTS
  BUILT_SOURCES = e2.lin.nc e2.nonlin.nc
endif

CLEANFILES = $(BUILT_SOURCES)

.config.nc::
	../src/python/glide_launch.py --prefix=.. --src $<
