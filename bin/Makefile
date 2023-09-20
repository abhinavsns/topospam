include /usr/local/openfpm/source/.openfpm.mk

CC=mpic++

LDIR =

OBJ = Active2d.o
OBJ2= SpringLattice.o

%.o: %.cpp
	$(CC) -O3 -march=native -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

Active2d: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)
SpringLattice: $(OBJ2)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: Active2d SpringLattice

run: all
	mpirun -np 4 /Active2d
OBJ = main.o

.PHONY: clean all run

clean:
	rm -f *.o *~ core Active2d SpringLattice