SUBDIRS = figs gnu

all:	subdirs

subdirs:	$(SUBDIRS)
$(SUBDIRS):
		$(MAKE) -C $@

.PHONY:		subdirs $(SUBDIRS)
