graphics = scale.pstex_t grid.pstex_t scale.pstex grid.pstex thick_evo.eps
EXTRA_DIST = grid.fig scale.fig $(graphics)
if MAINTAINER_MODE
 BUILT_SOURCES = $(graphics)
 CLEANFILES = $(BUILT_SOURCES)
endif
include $(top_srcdir)/doc/make_doc.rules
