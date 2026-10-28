EXECUTABLE = analyze
OUT_DIR = ..
MKDIR_P = mkdir -p

# Specify the compiler
COMP = g++

# Flags for compilation
FLAGS = -Wall -Wextra -Wpedantic -std=c++20 -O3 -g

# Homebrew dependencies. Any Boost works; prefer boost@1.85, which OpenFPM installs.
BREW_PREFIX ?= $(shell brew --prefix)
BOOST_PREFIX ?= $(firstword $(wildcard $(BREW_PREFIX)/opt/boost@1.85 $(BREW_PREFIX)/opt/boost))
DEP_PREFIXES = $(BOOST_PREFIX) $(BREW_PREFIX)/opt/libomp $(BREW_PREFIX)/opt/gsl

# Include directories
INCS = -I./src/ -I../src/ -I examples/ -I $(OUT_DIR) -I../src/eigen-3.4.0 $(addprefix -I,$(addsuffix /include,$(DEP_PREFIXES))) -I/usr/local/include
# Libraries and library paths
LIBPATH = $(addprefix -L,$(addsuffix /lib,$(DEP_PREFIXES))) -L/usr/local/lib
LIBS = -lgsl -lboost_iostreams -lboost_program_options -lgslcblas -lm

# Homebrew's library directories are not on the Linux loader path (and boost@1.85 is keg-only)
ifeq ($(shell uname -s),Linux)
  comma := ,
  LIBPATH += $(addprefix -Wl$(comma)-rpath$(comma),$(addsuffix /lib,$(DEP_PREFIXES)))
endif

# Header files
HEADERS = ../src/tools.hpp

# Object files
BASE = \
    $(OUT_DIR)/analyze_options.o \
    $(OUT_DIR)/Writer.o \
    $(OUT_DIR)/Vertex.o \
    $(OUT_DIR)/Triangle.o \
    $(OUT_DIR)/Edge.o \
    $(OUT_DIR)/Cell.o \
    $(OUT_DIR)/Triangulation.o \
    $(OUT_DIR)/Tissue.o

.PHONY: all clean cleanall directories debug

# Main targets
all: directories $(EXECUTABLE)

debug: FLAGS += -g -DDEBUG -ggdb3
debug: all

directories: $(OUT_DIR) $(OUT_DIR)/accessories

$(OUT_DIR) $(OUT_DIR)/accessories:
	$(MKDIR_P) $@

# General rule for building objects
$(OUT_DIR)/%.o: src/%.cpp src/%.hpp src/header.hpp.gch
	$(COMP) -c $< -o $@ $(INCS) $(FLAGS) $(LIBPATH) $(LIBS)

$(OUT_DIR)/accessories/%.o: src/accessories/%.cpp src/accessories/%.hpp src/header.hpp.gch
	$(COMP) -c $< -o $@ $(INCS) $(FLAGS) $(LIBPATH) $(LIBS)

# Precompile the header
src/header.hpp.gch: ../src/header.hpp
	$(COMP) -x c++-header $< -o $@ $(INCS) $(FLAGS)

# Linking
$(EXECUTABLE): examples/$(EXECUTABLE).cpp $(BASE) $(HEADERS)
	$(COMP) $< $(BASE) -o $@ $(INCS) $(FLAGS) $(LIBPATH) $(LIBS)

# Cleaning up
clean:
	rm -f $(OUT_DIR)/*.o $(OUT_DIR)/accessories/*.o analyze

cleanall: clean
	rm -rf $(OUT_DIR)
	rm -f src/header.hpp.gch
	rm -f $(EXECUTABLE)
