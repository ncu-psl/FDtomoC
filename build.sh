#!/bin/bash

mkdir build
cd build
cmake ..
make -j # if you get Internal Compiler Error(ICE), try make without -j

cd ..
mv run.sh ./build/
