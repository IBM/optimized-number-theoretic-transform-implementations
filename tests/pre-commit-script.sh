#!/bin/bash -ex
# Copyright IBM Inc. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Avoid removing the "build" directory if the script does not run from the 
# package root directory 
basedir=`pwd`
if [[ ! -f "$basedir/tests/pre-commit-script.sh" ]]; then
  >&2 echo "Script does not run from the root directory"
  exit 0
fi

# Clean previous build content
rm -rf build;

mkdir build;
cd build;

# Test clang-format
cmake ..; make format; 
rm -rf *

# Test clang-tidy
CC=clang-12 cmake -DCMAKE_C_CLANG_TIDY="clang-tidy;--format-style=file" ..
make -j20
rm -rf *

for flag in "" "-DASAN=1" "-DUBSAN=1" ; do
  CC=clang-12 cmake $flag ..; 
  make -j20 
  ./ntt-variants
  rm -rf *
done
