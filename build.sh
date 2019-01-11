#!/bin/sh

mkdir build
cd build
cmake ../src/
make -j # if you get Internal Compiler Error(ICE), try make without -j

cd ..
mv run.sh ./build/
