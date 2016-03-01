# -------------------
# define our includes
# -------------------

INCLUDES = -I. -I$(CASAVA_INCLUDE) -I$(BOOST_ROOT)/include

# ----------------------------------
# define our source and object files
# ----------------------------------

SOURCES=Exceptions.cpp FastqReader.cpp FileConversion.cpp LineReader.cpp Program.cpp Sequence.cpp StreamUtil.cpp StringUtilities.cpp ElandExtendedReader.cpp ExtendedFileReader.cpp
OBJECTS= $(SOURCES:.cpp=.o)
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))

all: $(BUILT_OBJECTS)

.PHONY: all

$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(CXXFLAGS) $(INCLUDES)

clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/* $(BIN_DIR)/*

.PHONY: clean
