# -------------------
# define our includes
# -------------------

# ----------------------------------
# define our source and object files
# ----------------------------------

PROGRAM=FastqConverter
OBJECTS=FastqConverter.o Exceptions.o Timer.o FileConversion.o StringUtilities.o 
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))
LIBS=-lz -L$(BOOST_ROOT)/lib -lboost_program_options -lboost_system -lboost_filesystem -lboost_regex -lboost_date_time -lboost_iostreams

all: $(PROGRAM)

.PHONY: all

$(PROGRAM): $(BUILT_OBJECTS)
	@echo "  * linking $(PROGRAM)"
	@$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $(BIN_DIR)/$@ $^ $(LIBS)

clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/* $(BIN_DIR)/*

.PHONY: clean
