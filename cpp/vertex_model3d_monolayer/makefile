EXECUTABLE = organoid
OUT_DIR = .
MKDIR_P = mkdir -p

# Compiler settings
COMP = g++
# gnu++20 rather than c++20: in strict mode, conda's clang does not get INFINITY and NAN from
# the math.h of the macOS 27 SDK
FLAGS = -Wall -Wextra -Wpedantic -std=gnu++20 -O3
comma := ,

# GSL and Boost from a conda environment (pixi.toml) or from Homebrew
ifneq ($(wildcard $(CONDA_PREFIX)/include/gsl),)
  COMP = $(CXX)
  DEP_PREFIXES = $(CONDA_PREFIX)
  # conda's libraries are found through the rpath, on macOS too
  RPATH_DIRS = $(CONDA_PREFIX)/lib
else
  # Any Boost works; prefer the one OpenFPM installs (boost, or boost@1.85 for older OpenFPM).
  BREW_PREFIX ?= $(shell brew --prefix)
  BOOST_PREFIX ?= $(firstword $(wildcard $(BREW_PREFIX)/opt/boost $(BREW_PREFIX)/opt/boost@1.85))
  DEP_PREFIXES = $(BOOST_PREFIX) $(BREW_PREFIX)/opt/libomp $(BREW_PREFIX)/opt/gsl
  EXTRA_INCS = -I/usr/local/include
  EXTRA_LIBPATH = -L/usr/local/lib
  # Homebrew's library directories are not on the Linux loader path (and boost@1.85 is keg-only)
  ifeq ($(shell uname -s),Linux)
    RPATH_DIRS = $(addsuffix /lib,$(DEP_PREFIXES))
  endif
endif
INCS = -I./src/ -I./examples/ -I$(OUT_DIR) -I./src/eigen-3.4.0 $(addprefix -I,$(addsuffix /include,$(DEP_PREFIXES))) $(EXTRA_INCS)
LIBPATH = $(addprefix -L,$(addsuffix /lib,$(DEP_PREFIXES))) $(EXTRA_LIBPATH) $(addprefix -Wl$(comma)-rpath$(comma),$(RPATH_DIRS))
LIBS = -lgsl -lboost_iostreams -lboost_program_options -lgslcblas -lm

HEADERS = src/tools.hpp
BASE = \
		 $(OUT_DIR)/options.o \
		 $(OUT_DIR)/Vertex.o \
		 $(OUT_DIR)/Edge.o \
		 $(OUT_DIR)/Triangle.o \
		 $(OUT_DIR)/Cell.o \
		 $(OUT_DIR)/Tissue.o

.PHONY: all clean cleanall directories debug

all: directories $(EXECUTABLE)

debug: FLAGS += -g -DDEBUG -ggdb3
debug: all

directories: $(OUT_DIR)

$(OUT_DIR):
	$(MKDIR_P) $@

$(OUT_DIR)/%.o: src/%.cpp src/%.hpp src/header.hpp.gch
	$(COMP) -c $< -o $@ $(INCS) $(FLAGS) $(LIBPATH) $(LIBS)

$(EXECUTABLE): $(BASE) examples/$(EXECUTABLE).cpp
	$(COMP) examples/$(EXECUTABLE).cpp $(BASE) -o $@ $(INCS) $(FLAGS) $(LIBPATH) $(LIBS)

src/header.hpp.gch: src/header.hpp
	$(COMP) -x c++-header $< -o $@ $(INCS) $(FLAGS)

clean:
	rm -f $(OUT_DIR)/*.o

cleanall: clean
	rm -rf $(OUT_DIR)
	rm -rf src/header.hpp.gch
	rm -f $(EXECUTABLE)
