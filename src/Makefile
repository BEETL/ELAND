# define our object and binary directories
export TOP_DIR = $(PWD)
export OBJ_DIR = $(TOP_DIR)/../obj
export BIN_DIR = $(TOP_DIR)/../bin

# define where samtools lives

# define our common source directories
export CASAVA_INCLUDE=$(TOP_DIR)/c++/include

# define some default flags
FLAGS = -Wall -Wno-unused-local-typedefs -Wno-delete-non-virtual-dtor -ansi -g -O3
#FLAGS = -Wall -ansi -g -O3
export CFLAGS = $(FLAGS)
export CXXFLAGS = $(FLAGS)
#export LDFLAGS 
export CXX ?= g++

SUBDIRS = c++ eland_ms FastqConverter squashGenome kagu orphanAligner samtools perl

all:
	@test -d $(OBJ_DIR) || mkdir $(OBJ_DIR)
	@test -d $(BIN_DIR) || mkdir $(BIN_DIR)
	@for dir in $(SUBDIRS); do \
		echo "- Building in $$dir"; \
		$(MAKE) --no-print-directory -C $$dir; \
		echo ""; \
	done

.PHONY: all

clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/* $(BIN_DIR)/*

.PHONY: clean
