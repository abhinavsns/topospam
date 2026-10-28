EXECUTABLE = analyze
OUT_DIR = ..
MKDIR_P = mkdir -p

# Specify the compiler
COMP = g++

# Flags for compilation
# gnu++20 rather than c++20: in strict mode, conda's clang does not get INFINITY and NAN from
# the math.h of the macOS 27 SDK
FLAGS = -Wall -Wextra -Wpedantic -std=gnu++20 -O3 -g

comma := ,

# GSL and Boost from a conda environment (pixi.toml) or from Homebrew
ifneq ($(wildcard $(CONDA_PREFIX)/include/gsl),)
  COMP = $(CXX)
  DEP_PREFIXES = $(CONDA_PREFIX)
  # conda's libraries are found through the rpath, on macOS too
  RPATH_DIRS = $(CONDA_PREFIX)/lib
else
  # Any Boost works; prefer boost@1.85, which OpenFPM installs.
  BREW_PREFIX ?= $(shell brew --prefix)
  BOOST_PREFIX ?= $(firstword $(wildcard $(BREW_PREFIX)/opt/boost@1.85 $(BREW_PREFIX)/opt/boost))
  DEP_PREFIXES = $(BOOST_PREFIX) $(BREW_PREFIX)/opt/libomp $(BREW_PREFIX)/opt/gsl
  EXTRA_INCS = -I/usr/local/include
  EXTRA_LIBPATH = -L/usr/local/lib
  # Homebrew's library directories are not on the Linux loader path (and boost@1.85 is keg-only)
  ifeq ($(shell uname -s),Linux)
    RPATH_DIRS = $(addsuffix /lib,$(DEP_PREFIXES))
  endif
endif

# Include directories
INCS = -I./src/ -I../src/ -I examples/ -I $(OUT_DIR) -I../src/eigen-3.4.0 $(addprefix -I,$(addsuffix /include,$(DEP_PREFIXES))) $(EXTRA_INCS)
# Libraries and library paths
LIBPATH = $(addprefix -L,$(addsuffix /lib,$(DEP_PREFIXES))) $(EXTRA_LIBPATH) $(addprefix -Wl$(comma)-rpath$(comma),$(RPATH_DIRS))
LIBS = -lgsl -lboost_iostreams -lboost_program_options -lgslcblas -lm

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
