# -------------------
# define our includes
# -------------------

INCLUDES = -I. -I$(CASAVA_INCLUDE) -I$(BOOST_ROOT)/include

# ----------------------------------
# define our source and object files
# ----------------------------------

SOURCES=AlignmentQuality.cpp AlignmentResolver.cpp ExportWriter.cpp Timer.cpp AlignmentReader.cpp AnomalyWriter.cpp ConfigurationSettings.cpp XmlTree.cpp

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
