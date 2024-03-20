include ~/example.mk

CC=mpic++

LDIR =

OBJ = Active2d.o
OBJ2= SpringLattice.o

%.o: %.cpp
	$(CC) -O0 -g -march=native -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

Active2d: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)
SpringLattice: $(OBJ2)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: Active2d SpringLattice

.PHONY: clean all run

clean:
	rm -f *.o *~ core Active2d SpringLattice