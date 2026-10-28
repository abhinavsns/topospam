# Builds the TopoSPAM C++ simulators. Their dependencies (OpenFPM and friends from Homebrew)
# have to be installed first: ./install.sh does both, see README.md.
#
#   make backend                                build all simulators
#   make backend ARCH_FLAGS=-march=x86-64-v2    portable binaries (e.g. for a Docker image)
#   make clean                                  remove the build products

.PHONY: backend clean

# On Linux, Homebrew installs its own glibc when the system one is older than its bottles need.
# Programs linked against Homebrew's libraries must then be built with Homebrew's GCC as well.
ifeq ($(shell uname -s),Linux)
  BREW_PREFIX ?= $(shell brew --prefix 2>/dev/null)
  ifneq ($(wildcard $(BREW_PREFIX)/opt/glibc),)
    BREW_GXX := $(shell ls $(BREW_PREFIX)/bin/g++-[0-9]* 2>/dev/null | sort -V | tail -n 1)
    ifneq ($(BREW_GXX),)
      OMPI_CXX ?= $(BREW_GXX)
      export OMPI_CXX
      VM_COMPILER := COMP=$(BREW_GXX)
    endif
  endif
endif

backend:
	$(MAKE) -C cpp all
	$(MAKE) -C cpp/vertex_model3d_monolayer $(VM_COMPILER)
	$(MAKE) -C cpp/vertex_model3d_monolayer/accessories $(VM_COMPILER)

clean:
	rm -f cpp/*.o cpp/Active2d cpp/SpringLattice
	rm -f cpp/vertex_model3d_monolayer/*.o cpp/vertex_model3d_monolayer/organoid
	rm -f cpp/vertex_model3d_monolayer/accessories/analyze
	rm -f cpp/vertex_model3d_monolayer/src/header.hpp.gch cpp/vertex_model3d_monolayer/accessories/src/header.hpp.gch
