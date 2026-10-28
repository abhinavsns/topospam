# example.mk
#
# Compiler and linker flags for the OpenFPM-based simulators (Active2d, SpringLattice).
# The dependencies come from a conda environment (pixi.toml), in which OpenFPM is installed
# in $CONDA_PREFIX/openfpm, or from Homebrew (install.sh). Homebrew paths go through its
# stable "opt" links, so nothing has to be regenerated after `brew upgrade`.
SHELL := /bin/bash
comma := ,
OPENFPM_INCLUDE_DIRS := openfpm_numerics/include openfpm_pdata/include/config openfpm_pdata/include \
                        openfpm_data/include openfpm_vcluster/include openfpm_io/include openfpm_devices/include

ifneq ($(wildcard $(CONDA_PREFIX)/openfpm/openfpm_pdata/include),)
  # conda environment: all dependencies are installed in one prefix
  OPENFPM_PREFIX ?= $(CONDA_PREFIX)/openfpm
  DEP_INCLUDE_DIRS := $(CONDA_PREFIX)/include $(CONDA_PREFIX)/include/eigen3 $(CONDA_PREFIX)/include/suitesparse
  DEP_LIB_DIRS := $(CONDA_PREFIX)/lib
  # conda's libraries are found through the rpath, on macOS too
  RPATH_DIRS := $(DEP_LIB_DIRS)
  # OpenFPM's headers first: the environment's include directory is shared by all packages
  INCLUDE_PATH := $(addprefix -I$(OPENFPM_PREFIX)/,$(OPENFPM_INCLUDE_DIRS)) $(addprefix -I,$(DEP_INCLUDE_DIRS))
else
  BREW_PREFIX ?= $(shell brew --prefix 2>/dev/null)
  ifeq ($(BREW_PREFIX),)
    $(error Neither a conda environment with OpenFPM (pixi.toml) nor Homebrew was found. See the installation section of README.md)
  endif
  OPT := $(BREW_PREFIX)/opt
  OPENFPM_PREFIX ?= $(OPT)/openfpm
  # Formulae that provide headers and libraries
  DEPS := petsc metis parmetis libhilbert boost@1.85 vc blitz algoim suite-sparse libomp hdf5-mpi
  DEP_INCLUDE_DIRS := $(wildcard $(addsuffix /include,$(addprefix $(OPT)/,$(DEPS))))
  DEP_LIB_DIRS := $(wildcard $(addsuffix /lib,$(addprefix $(OPT)/,$(DEPS))))
  # OpenFPM 5.0.1 does not compile with Eigen 5 (Homebrew's "eigen" since late 2025): use eigen@3 if present
  EIGEN_INCLUDE := $(firstword $(wildcard $(OPT)/eigen@3/include/eigen3 $(OPT)/eigen/include/eigen3))
  INCLUDE_PATH := $(addprefix -I,$(DEP_INCLUDE_DIRS)) \
                  -I$(EIGEN_INCLUDE) -I$(OPT)/suite-sparse/include/suitesparse \
                  $(addprefix -I$(OPENFPM_PREFIX)/,$(OPENFPM_INCLUDE_DIRS))
  # Homebrew's library directories are not on the Linux loader path (and boost@1.85 is
  # keg-only), so record them in the binaries. macOS uses absolute install names instead.
  ifeq ($(shell uname -s),Linux)
    RPATH_DIRS := $(DEP_LIB_DIRS) $(BREW_PREFIX)/lib
  endif
endif
ifeq ($(wildcard $(OPENFPM_PREFIX)/openfpm_pdata/include),)
  $(error OpenFPM was not found in $(OPENFPM_PREFIX). Install it with ./install.sh or see README.md)
endif
OPENFPM_LIB_DIRS := $(addprefix $(OPENFPM_PREFIX)/,openfpm_devices/lib openfpm_pdata/lib openfpm_vcluster/lib)

# The OpenFPM 5.0.1 install does not ship its generated config.h, which its headers include.
# openfpm_config/ provides a copy and comes last, so a config.h from OpenFPM itself wins.
INCLUDE_PATH += -I$(CURDIR)/openfpm_config

LIBS_PATH := $(addprefix -L,$(DEP_LIB_DIRS) $(OPENFPM_LIB_DIRS)) \
             $(addprefix -Wl$(comma)-rpath$(comma),$(RPATH_DIRS))

# Predefined libraries list
LIBS := -lvcluster -lofpm_pdata -lofpmmemory -lparmetis -lmetis -lboost_iostreams \
        -lboost_program_options -lhdf5 -llibhilbert -lVc -lpetsc -lumfpack -lamd \
        -lbtf -lcamd -lccolamd -lcholmod -lcolamd -lcxsparse -lklu -ldl -lrbio \
        -lspqr -lsuitesparseconfig -lboost_filesystem

# CUDA flag: defaults to NO; override via 'cuda_on_cpu' variable if needed
CUDA_ON_CPU ?= NO
ifeq ($(cuda_on_cpu),YES)
  CUDA_ON_CPU := YES
endif
