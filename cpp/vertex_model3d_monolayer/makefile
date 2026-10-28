EXECUTABLE = organoid
OUT_DIR = .
MKDIR_P = mkdir -p

# Homebrew dependencies. Any Boost works; prefer boost@1.85, which OpenFPM installs.
BREW_PREFIX ?= $(shell brew --prefix)
BOOST_PREFIX ?= $(firstword $(wildcard $(BREW_PREFIX)/opt/boost@1.85 $(BREW_PREFIX)/opt/boost))
DEP_PREFIXES = $(BOOST_PREFIX) $(BREW_PREFIX)/opt/libomp $(BREW_PREFIX)/opt/gsl

# Compiler settings
COMP = g++
FLAGS = -Wall -Wextra -Wpedantic -std=c++20 -O3
INCS = -I./src/ -I./examples/ -I$(OUT_DIR) -I./src/eigen-3.4.0 $(addprefix -I,$(addsuffix /include,$(DEP_PREFIXES))) -I/usr/local/include
LIBPATH = $(addprefix -L,$(addsuffix /lib,$(DEP_PREFIXES))) -L/usr/local/lib
LIBS = -lgsl -lboost_iostreams -lboost_program_options -lgslcblas -lm

# Homebrew's library directories are not on the Linux loader path (and boost@1.85 is keg-only)
ifeq ($(shell uname -s),Linux)
  comma := ,
  LIBPATH += $(addprefix -Wl$(comma)-rpath$(comma),$(addsuffix /lib,$(DEP_PREFIXES)))
endif

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
