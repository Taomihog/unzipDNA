#!/bin/bash

read -p "Enter look up table size (Default=128): " LUT_SIZE
if [[ -z "$LUT_SIZE" ]]; then
    LUT_SIZE="128"
fi

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BUILD_DIR="$SCRIPT_DIR/build"

if [[ ! -d "$BUILD_DIR" ]]; then
    mkdir "$BUILD_DIR"
fi
cd "$BUILD_DIR"

# basepair energy related stuff
g++ -c ../BPEnergy/BPEnergy.cpp -o BPEnergy.o
ar rcs BPEnergy.a BPEnergy.o

# create LUT
g++ -c ../UnzipLUT/UnzipLUT.cpp -std=c++20 -fconstexpr-ops-limit=100000000000 -DJ_SIZE="$LUT_SIZE" -DEXT_SIZE="$LUT_SIZE" -o UnzipLUT.o
ar rcs UnzipLUT.a UnzipLUT.o

# executable
g++ ../main.cpp -L. -lUnzipLUT -lBPEnergy -o UnzipDNA

./UnzipDNA ../exampleData.txt

cd ..