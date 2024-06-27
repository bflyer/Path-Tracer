#!/usr/bin/env bash
SAMPLES=15
METHOD=pt
cmake -B build
cmake --build build

# Run all testcases. 
# You can comment some lines to disable the run of specific examples.
mkdir -p output

# build/PA1 testcases/scene08_smallpt.txt output/scene08.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene08_smallpt_square.txt output/scene08_square.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene10_smallpt_square.txt output/scene10_square.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene11_box.txt output/scene11_box.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene13_star.txt output/scene13_star.bmp $METHOD $SAMPLES
build/PA1 testcases/scene14_star_yellow.txt output/scene14_star_yellow.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene12.txt output/scene12.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene09_sky.txt output/scene09.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene08_smallpt_allball.txt output/scene08_allball.bmp $METHOD $SAMPLES