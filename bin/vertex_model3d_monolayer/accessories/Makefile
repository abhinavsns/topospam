EXECUTABLE = analyze

INCS = -I src/ -I./../src/ -I examples/ -I build/ -I./../src/eigen-3.4.0

FLAGS = -Wall -Wextra -Wpedantic -std=c++20 -O3 -fopenmp -g

LIBS = # -lboost_program_options

LIBPATH =

COMP = g++

HEADERS = ./../src/tools.hpp

BASE = \
			 build/analyze_options.o \
			 build/Writer.o \
			 build/Vertex.o \
			 build/Triangle.o \
       build/Edge.o \
       build/Cell.o \
			 build/Triangulation.o \
       build/Tissue.o

#

OUT_DIR = build/

MKDIR_P = mkdir -p

.PHONY: all clean cleanall directories debug splat

all: LIBS += -lboost_program_options -lgsl -lgslcblas -lm

all: directories $(EXECUTABLE)

debug: FLAGS += -g -DDEBUG -ggdb3

debug: all

splat: SPLATLIBS = /usr/local/boost/1.75/lib/libboost_program_options.a

splat: LIBPATH += -L/usr/local/boost_1_75_0

splat: directories $(EXECUTABLE)

directories: ${OUT_DIR}

${OUT_DIR}:
	${MKDIR_P} ${OUT_DIR}

src/header.hpp.gch: ./../src/header.hpp
	$(COMP) -c $< -o $@ $(INCS) $(FLAGS) $(LIBPATH) $(LIBS)

build/%.o: src/%.cpp src/%.hpp src/header.hpp.gch
	$(COMP) -c $< -o $@ $(INCS) $(FLAGS) $(LIBPATH) $(LIBS)

build/%.o: ./../src/%.cpp ./../src/%.hpp src/header.hpp.gch
	$(COMP) -c $< -o $@ $(INCS) $(FLAGS) $(LIBPATH) $(LIBS)

$(EXECUTABLE): examples/$(EXECUTABLE).cpp $(BASE) $(HEADERS)
	$(COMP) $^ -o $@ $(INCS) $(FLAGS) $(LIBPATH) $(LIBS) $(SPLATLIBS)

clean:
	rm -f build/*.o

cleanall: clean
	rm -rf build/
	rm -rf src/header.hpp.gch
	rm $(EXECUTABLE)
