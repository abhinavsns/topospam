# example.mk
SHELL := /bin/bash

# Compute include paths for dependencies
INCLUDE_PATH := $(shell \
  deps="petsc metis parmetis libhilbert boost@1.85 vc blitz algoim suitesparse libomp hdf5-mpi"; \
  inc=""; \
  for dep in $$deps; do \
    prefix=`brew --prefix $$dep`; \
    if [ -d "$$prefix/include" ]; then \
      inc="$$inc -I$$prefix/include"; \
    fi; \
  done; \
  inc="$$inc -I`brew --prefix eigen`/include/eigen3 -I`brew --prefix suitesparse`/include/suitesparse"; \
  openfpm_prefix=`brew --prefix openfpm`; \
  inc="$$inc -I$$openfpm_prefix/openfpm_numerics/include -I$$openfpm_prefix/openfpm_pdata/include/config -I$$openfpm_prefix/openfpm_pdata/include -I$$openfpm_prefix/openfpm_data/include -I$$openfpm_prefix/openfpm_vcluster/include -I$$openfpm_prefix/openfpm_io/include -I$$openfpm_prefix/openfpm_devices/include"; \
  echo $$inc \
)

# Compute library paths for dependencies
LIBS_PATH := $(shell \
  deps="petsc metis parmetis libhilbert boost@1.85 vc blitz algoim suitesparse libomp hdf5-mpi"; \
  libs=""; \
  for dep in $$deps; do \
    prefix=`brew --prefix $$dep`; \
    if [ -d "$$prefix/lib" ]; then \
      libs="$$libs -L$$prefix/lib"; \
    fi; \
  done; \
  openfpm_prefix=`brew --prefix openfpm`; \
  libs="$$libs -L$$openfpm_prefix/openfpm_devices/lib -L$$openfpm_prefix/openfpm_pdata/lib -L$$openfpm_prefix/openfpm_vcluster/lib"; \
  echo $$libs \
)

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
