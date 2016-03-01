# -------------------
# define our includes
# -------------------

# ----------------------------------
# define our source and object files
# ----------------------------------

PROGRAM=squashGenome
OBJECTS=squashGenome.o SquashGenome.o Exceptions.o Program.o
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))
LIBS=-lz -L$(BOOST_ROOT)/lib -lboost_program_options -lboost_system -lboost_filesystem -lboost_regex -lboost_date_time

all: $(PROGRAM)

.PHONY: all

$(PROGRAM): $(BUILT_OBJECTS)
	@echo "  * linking $(PROGRAM)"
	@$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $(BIN_DIR)/$@ $^ $(LIBS)

clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/* $(BIN_DIR)/*

.PHONY: clean
