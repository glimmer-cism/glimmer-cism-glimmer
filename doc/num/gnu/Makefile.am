graphics = w_profile.pslatex wt_sigma.pslatex w_profile.ps wt_sigma.ps

EXTRA_DIST = w_profile.gp wt_sigma.gp velo/calc2.data  velo/e1-mm.1.data  velo/e1-mm2.1.data  velo/temp_cor.data  velo/temp_uncor.data  velo/velo_cor.data  velo/velo_uncor.data $(graphics)
 
if MAINTAINER_MODE
 BUILT_SOURCES = $(graphics)
 CLEANFILES = $(BUILT_SOURCES)
endif
include $(top_srcdir)/doc/make_doc.rules
