SUBDIRS = src/fortran doc

subdirs:	$(SUBDIRS)

$(SUBDIRS):
		$(MAKE) -C $@

clean:
		for dir in $(SUBDIRS); do \
		  $(MAKE) -C $$dir clean; \
		done

install:
		for dir in $(SUBDIRS); do \
		  $(MAKE) -C $$dir install; \
		done

.PHONY:		clean install subdirs $(SUBDIRS)
