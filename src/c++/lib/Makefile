# define our build subdirectories
SUBDIRS = alignment blt_util common eland_ms kagu

all:
	@for dir in $(SUBDIRS); do \
		echo "- Building in $$dir"; \
		$(MAKE) --no-print-directory -C $$dir; \
		echo ""; \
	done

.PHONY: all
