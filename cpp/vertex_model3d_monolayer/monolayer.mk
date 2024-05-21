INCS=-I./src/ -I./examples/ -I./build/ -I./src/eigen-3.4.0  -I$(shell brew --prefix boost)/include -I$(shell brew --prefix libomp)/include -I$(shell brew --prefix gsl)/include -I/usr/local/include
LIBPATH= -L$(shell brew --prefix boost)/lib -L$(shell brew --prefix libomp)/lib -L$(shell brew --prefix gsl)/lib -L/usr/local/lib
FLAGS = -Wall -Wextra -Wpedantic -std=c++20 -O3 -g -Wl,-ld_classic
LIBS = -lgsl -lboost_iostreams -lboost_program_options # -lboost_program_options
COMP = mpic++
