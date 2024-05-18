#! /bin/bash

# Function to find and add a dependency's include and lib paths
add_dependency_paths() {
  local dep_prefix=$(brew --prefix $1)

  # Include paths
  if [ -d "$dep_prefix/include" ]; then
    INCLUDE_PATH="$INCLUDE_PATH -I$dep_prefix/include"
  fi
  
  # Library paths
  if [ -d "$dep_prefix/lib" ]; then
    LIBS_PATH="$LIBS_PATH -L$dep_prefix/lib"
    # Add library link flags as necessary. You might need a specific approach to determine the library names.
  fi
}

# Initialize variables
INCLUDE_PATH=""
LIBS_PATH=""
LIBS=""

# Dependencies array. Ensure these are the correct formula names in Homebrew.
depends=("petsc" "metis" "parmetis" "libhilbert" "boost" "vc" "blitz" "algoim" "suitesparse" "libomp" "hdf5-mpi")

# Add paths for each dependency
for dep in "${depends[@]}"; do
  add_dependency_paths $dep
done

INCLUDE_PATH="$INCLUDE_PATH -I$(brew --prefix eigen)/include/eigen3/"
INCLUDE_PATH="$INCLUDE_PATH -I$(brew --prefix suitesparse)/include/suitesparse/"


# Additional dependencies or paths can be added directly. Adjust as necessary.
INCLUDE_PATH="$INCLUDE_PATH"
LIBS_PATH="$LIBS_PATH"
LIBS="$LIBS -lvcluster -lofpm_pdata -lofpmmemory -lparmetis -lmetis -lboost_iostreams -lboost_program_options -lhdf5 -llibhilbert -lVc -lpetsc -lumfpack -lamd -lbtf -lcamd -lccolamd -lcholmod -lcolamd -lcxsparse -lklu -ldl -lrbio -lspqr -lsuitesparseconfig -lboost_filesystem"

# Assuming `openfpm` is a prefix for related paths, not a dependency to be iterated over.
dep_prefix=$(brew --prefix openfpm)
INCLUDE_PATH="$INCLUDE_PATH -I$dep_prefix/openfpm_numerics/include -I$dep_prefix/openfpm_pdata/include/config -I$dep_prefix/openfpm_pdata/include -I$dep_prefix/openfpm_data/include -I$dep_prefix/openfpm_vcluster/include -I$dep_prefix/openfpm_io/include -I$dep_prefix/openfpm_devices/include"
LIBS_PATH="$LIBS_PATH -L$dep_prefix/openfpm_devices/lib -L$dep_prefix/openfpm_pdata/lib -L$dep_prefix/openfpm_vcluster/lib"

# Generate example.mk
echo "INCLUDE_PATH=$INCLUDE_PATH" > ./bin/example.mk
echo "LIBS_PATH=$LIBS_PATH" >> ./bin/example.mk
echo "LIBS=$LIBS" >> ./bin/example.mk
# Handle CUDA_ON_CPU flag
CUDA_ON_CPU=NO  # Default setting
if [ x"$cuda_on_cpu" == x"YES" ]; then
   CUDA_ON_CPU=YES
fi
echo "CUDA_ON_CPU=$CUDA_ON_CPU" >> ./bin/example.mk

INCS="-I./src/ -I./examples/ -I./build/ -I./src/eigen-3.4.0 " 

add_dependency_paths_monolayer() {
  local dep_prefix_monolayer=$(brew --prefix $1)

  # Include paths
  #if [ -d "$dep_prefix_monolayer/include" ]; then
  
    INCS="$INCS -I$dep_prefix_monolayer/include"
  #fi
  
  # Library paths
  #if [ -d "$dep_prefix_monolayer/lib" ]; then
    LIBPATH="$LIBPATH -L$dep_prefix_monolayer/lib"
  #fi
}

dependsMonolayer=("boost" "libomp" "gsl" " ")

# Add paths for each dependency
for dep in "${dependsMonolayer[@]}"; do
  add_dependency_paths_monolayer $dep
done

echo "INCS=$INCS" > ./bin/vertex_model3d_monolayer/monolayer.mk
echo "LIBPATH=$LIBPATH" >> ./bin/vertex_model3d_monolayer/monolayer.mk

echo "FLAGS = -Wall -Wextra -Wpedantic -std=c++20 -O3 -g >> ./bin/vertex_model3d_monolayer/monolayer.mk

echo "LIBS = -lgsl -lboost_iostreams -lboost_program_options # -lboost_program_options" >> ./bin/vertex_model3d_monolayer/monolayer.mk

echo "COMP = mpic++" >> ./bin/vertex_model3d_monolayer/monolayer.mk