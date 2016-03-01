# -------------------
# define our includes
# -------------------

# ----------------------------------
# define our source and object files
# ----------------------------------

PROGRAM=eland_ms
OBJECTS=eland_ms.o aligner.o BclReader.o GlobalUtilities.o OligoSourceBcl.o OligoSourceFastq.o OligoSourceQseq.o parse_util.o Exceptions.o FastqReader.o LineReader.o Program.o Sequence.o StreamUtil.o ContigNameFinder.o ELAND_options_ms.o Hasher.o MatchTable.o StateMachine.o SuffixScoreTable.o
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
