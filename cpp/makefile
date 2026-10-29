include ./example.mk

.DEFAULT_GOAL := all

CC=mpic++
CFLAGS=
LDIR =
# Tune for the build machine by default. Binaries that are copied to other machines
# (e.g. the Docker image) must use a portable value such as ARCH_FLAGS=-march=x86-64-v2.
ARCH_FLAGS ?= -march=native

OBJ = Active2d.o
OBJ2= SpringLattice.o

# -fsigned-char: OpenFPM's DCPSE operators keep a sign in a plain char (Dcpse.hpp), which is
# unsigned on Linux/ARM64 and breaks all even-order derivatives there. x86-64 and macOS
# already use a signed char.
%.o: %.cpp
	$(CC) -O3 -g $(ARCH_FLAGS) -fsigned-char -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

Active2d: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)
SpringLattice: $(OBJ2)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: Active2d SpringLattice

.PHONY: clean all run

clean:
	rm -f *.o *~ core Active2d SpringLattice
