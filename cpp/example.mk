# example.mk
#
# Compiler and linker flags for the OpenFPM-based simulators (Active2d, SpringLattice).
# All dependencies are Homebrew formulae (macOS and Linux). Paths go through Homebrew's
# stable "opt" links, so nothing has to be regenerated after `brew upgrade`.
SHELL := /bin/bash

BREW_PREFIX ?= $(shell brew --prefix 2>/dev/null)
ifeq ($(BREW_PREFIX),)
  $(error Homebrew (brew) was not found on PATH. See the installation section of README.md)
endif
OPT := $(BREW_PREFIX)/opt
OPENFPM_PREFIX ?= $(OPT)/openfpm
ifeq ($(wildcard $(OPENFPM_PREFIX)/openfpm_pdata/include),)
  $(error OpenFPM was not found in $(OPENFPM_PREFIX). Install it with ./install.sh or see README.md)
endif

# Formulae that provide headers and libraries
DEPS := petsc metis parmetis libhilbert boost@1.85 vc blitz algoim suite-sparse libomp hdf5-mpi
DEP_INCLUDE_DIRS := $(wildcard $(addsuffix /include,$(addprefix $(OPT)/,$(DEPS))))
DEP_LIB_DIRS := $(wildcard $(addsuffix /lib,$(addprefix $(OPT)/,$(DEPS))))
OPENFPM_LIB_DIRS := $(addprefix $(OPENFPM_PREFIX)/,openfpm_devices/lib openfpm_pdata/lib openfpm_vcluster/lib)
# OpenFPM 5.0.1 does not compile with Eigen 5 (Homebrew's "eigen" since late 2025): use eigen@3 if present
EIGEN_INCLUDE := $(firstword $(wildcard $(OPT)/eigen@3/include/eigen3 $(OPT)/eigen/include/eigen3))

# The OpenFPM 5.0.1 install does not ship its generated config.h, which its headers include.
# openfpm_config/ provides a copy and comes last, so a config.h from OpenFPM itself wins.
INCLUDE_PATH := $(addprefix -I,$(DEP_INCLUDE_DIRS)) \
                -I$(EIGEN_INCLUDE) -I$(OPT)/suite-sparse/include/suitesparse \
                $(addprefix -I$(OPENFPM_PREFIX)/,openfpm_numerics/include openfpm_pdata/include/config \
                  openfpm_pdata/include openfpm_data/include openfpm_vcluster/include \
                  openfpm_io/include openfpm_devices/include) \
                -I$(CURDIR)/openfpm_config

LIBS_PATH := $(addprefix -L,$(DEP_LIB_DIRS) $(OPENFPM_LIB_DIRS))

# Homebrew's library directories are not on the Linux loader path (and boost@1.85 is
# keg-only), so record them in the binaries. macOS uses absolute install names instead.
ifeq ($(shell uname -s),Linux)
  comma := ,
  LIBS_PATH += $(addprefix -Wl$(comma)-rpath$(comma),$(DEP_LIB_DIRS) $(BREW_PREFIX)/lib)
endif

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
