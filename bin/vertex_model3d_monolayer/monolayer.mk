INCS= -I/opt/homebrew/opt/boost/include -I/opt/homebrew/opt/libomp/include -I/opt/homebrew/opt/gsl/include -I/opt/homebrew/include
LIBPATH= -L/opt/homebrew/opt/boost/lib -L/opt/homebrew/opt/libomp/lib -L/opt/homebrew/opt/gsl/lib -L/opt/homebrew/lib
FLAGS = -Wall -Wextra -Wpedantic -std=c++20 -O3 -fopenmp -g -Wl,-ld_classic
LIBS = -lgsl -lboost_iostreams -lboost_program_options # -lboost_program_options
COMP = mpic++
