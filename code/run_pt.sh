#!/usr/bin/env bash
SAMPLES=100
METHOD=pt
cmake -B build
cmake --build build

# Run all testcases. 
# You can comment some lines to disable the run of specific examples.
mkdir -p output

build/PA1 testcases/scene08_smallpt.txt output/scene08.bmp $METHOD $SAMPLES