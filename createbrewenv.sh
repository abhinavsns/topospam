#! /bin/bash

# On Linux, force the MPI C++ compiler to be g++-11 from Linuxbrew.
if [ "$(uname)" = "Linux" ]; then
  export MPICXX=/home/linuxbrew/bin/g++-11
fi

# If MPICXX is not already set (e.g. on macOS or if you unset it), default to the one in PATH.
if [ -z "$MPICXX" ]; then
  MPICXX=$(which mpic++)
fi

# Function to add dependency include and lib paths.
add_dependency_paths() {
  local dep_prefix
  dep_prefix=$(brew --prefix "$1")
  
  # Add include paths if available.
  if [ -d "$dep_prefix/include" ]; then
    INCLUDE_PATH="$INCLUDE_PATH -I$dep_prefix/include"
  fi
  
  # Add library paths if available.
  if [ -d "$dep_prefix/lib" ]; then
    LIBS_PATH="$LIBS_PATH -L$dep_prefix/lib"
  fi
}

# Initialize variables.
INCLUDE_PATH=""
LIBS_PATH=""
LIBS=""

# Array of dependencies installed via Homebrew.
depends=("petsc" "metis" "parmetis" "libhilbert" "boost@1.85" "vc" "blitz" "algoim" "suitesparse" "libomp" "hdf5-mpi")

# Add paths for each dependency.
for dep in "${depends[@]}"; do
  add_dependency_paths "$dep"
done

# Extra include directories.
INCLUDE_PATH="$INCLUDE_PATH -I$(brew --prefix eigen)/include/eigen3/"
INCLUDE_PATH="$INCLUDE_PATH -I$(brew --prefix suitesparse)/include/suitesparse/"

# Additional libraries to link.
LIBS="$LIBS -lvcluster -lofpm_pdata -lofpmmemory -lparmetis -lmetis -lboost_iostreams -lboost_program_options -lhdf5 -llibhilbert -lVc -lpetsc -lumfpack -lamd -lbtf -lcamd -lccolamd -lcholmod -lcolamd -lcxsparse -lklu -ldl -lrbio -lspqr -lsuitesparseconfig -lboost_filesystem"

# Add OpenFPM paths.
dep_prefix=$(brew --prefix openfpm)
INCLUDE_PATH="$INCLUDE_PATH -I$dep_prefix/openfpm_numerics/include -I$dep_prefix/openfpm_pdata/include/config -I$dep_prefix/openfpm_pdata/include -I$dep_prefix/openfpm_data/include -I$dep_prefix/openfpm_vcluster/include -I$dep_prefix/openfpm_io/include -I$dep_prefix/openfpm_devices/include"
LIBS_PATH="$LIBS_PATH -L$dep_prefix/openfpm_devices/lib -L$dep_prefix/openfpm_pdata/lib -L$dep_prefix/openfpm_vcluster/lib"

# Generate example.mk with the paths.
{
  echo "INCLUDE_PATH=$INCLUDE_PATH"
  echo "LIBS_PATH=$LIBS_PATH"
  echo "LIBS=$LIBS"
  # Optionally add CUDA flag.
  CUDA_ON_CPU=NO
  if [ x"$cuda_on_cpu" == x"YES" ]; then
     CUDA_ON_CPU=YES
  fi
  echo "CUDA_ON_CPU=$CUDA_ON_CPU"
} > ./cpp/example.mk

# Setup monolayer-specific paths.
INCS="-I./src/ -I./examples/ -I./build/ -I./src/eigen-3.4.0 "

# Function for monolayer dependencies.
add_dependency_paths_monolayer() {
  local dep_prefix_monolayer
  dep_prefix_monolayer=$(brew --prefix "$1")
  INCS="$INCS -I$dep_prefix_monolayer/include"
  LIBPATH="$LIBPATH -L$dep_prefix_monolayer/lib"
}

dependsMonolayer=("boost@1.85" "libomp" "gsl")

for dep in "${dependsMonolayer[@]}"; do
  add_dependency_paths_monolayer "$dep"
done

# Generate the makefile for the vertex model 3d monolayer.
{
  echo "INCS=$INCS"
  echo "LIBPATH=$LIBPATH"
  echo 'FLAGS = -Wall -Wextra -Wpedantic -std=c++20 -O3 -g'
  echo 'LIBS = -lgsl -lboost_iostreams -lboost_program_options'
  # Include the MPI C++ compiler variable.
  echo "MPICXX = $MPICXX"
  echo "COMP = $MPICXX"
} > ./cpp/vertex_model3d_monolayer/monolayer.mk
