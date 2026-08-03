#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "$0")" && pwd)"
jobs="${JOBS:-8}"

cmake_args=(
  -S "$script_dir"
  -B "$script_dir/build"
  -DCMAKE_BUILD_TYPE=Release
  -DCMAKE_DISABLE_FIND_PACKAGE_MPI=TRUE
  -DCMAKE_INSTALL_PREFIX="$script_dir/install"
)

if [[ "$(uname -s)" == "Darwin" ]]; then
  libomp_prefix="${LIBOMP_PREFIX:-$(brew --prefix libomp)}"
  cmake_args+=(
    "-DOpenMP_CXX_FLAGS=-Xpreprocessor -fopenmp"
    -DOpenMP_CXX_LIB_NAMES=omp
    "-DOpenMP_omp_LIBRARY=$libomp_prefix/lib/libomp.dylib"
    "-DCMAKE_CXX_FLAGS=-I$libomp_prefix/include"
  )
fi

cmake "${cmake_args[@]}"
cmake --build "$script_dir/build" --parallel "$jobs"
cmake --install "$script_dir/build"
