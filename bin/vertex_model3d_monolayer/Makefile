include monolayer.mk

EXECUTABLE = organoid
OUT_DIR = .
MKDIR_P = mkdir -p

# Compiler settings
COMP = mpic++  # Ensure this is set to clang++ if using Clang
FLAGS = -Wall -Wextra -Wpedantic -std=c++20
INCS = -I src/ -I examples/ -I$(OUT_DIR) -I./src/eigen-3.4.0 -I$(shell brew --prefix boost)/include -I$(shell brew --prefix libomp)/include -I$(shell brew --prefix gsl)/include -I/usr/local/include
LIBPATH = -L$(shell brew --prefix boost)/lib -L$(shell brew --prefix libomp)/lib -L$(shell brew --prefix gsl)/lib -L/usr/local/lib
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
	$(COMP) -x c++-header -o $@ $< $(INCS) $(FLAGS)

clean:
	rm -f $(OUT_DIR)/*.o

cleanall: clean
	rm -rf $(OUT_DIR)
	rm -rf src/header.hpp.gch
	rm -f $(EXECUTABLE)