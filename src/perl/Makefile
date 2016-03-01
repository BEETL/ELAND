# -------------------
# define our includes
# -------------------

# ----------------------------------
# define our source and object files
# ----------------------------------

PROGRAM=ELAND_standalone.pl

all: $(PROGRAM)

.PHONY: all

$(PROGRAM):
	@echo "  * building $@"
	@cp bin/$@ $(BIN_DIR)/$@

clean:
	@echo "Cleaning up."
	@rm -f $(BIN_DIR)/$(PROGRAM)

.PHONY: clean
