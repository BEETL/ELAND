# -------------------
# define our includes
# -------------------

# ----------------------------------
# define our source and object files
# ----------------------------------

PROGRAM=kagu
OBJECTS=kagu.o AlignmentQuality.o AlignmentResolver.o ExportWriter.o Timer.o AlignmentReader.o AnomalyWriter.o ConfigurationSettings.o XmlTree.o LineReader.o ElandExtendedReader.o StringUtilities.o Exceptions.o FastqReader.o

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
