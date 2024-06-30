#!/usr/bin/env bash
SAMPLES=1
SAMPLES1000=1000
SAMPLES500=500
SAMPLES200=50
METHOD=pt
cmake -B build
cmake --build build

# Run all testcases. 
# You can comment some lines to disable the run of specific examples.
mkdir -p output

# build/PA1 testcases/scene08_smallpt_square.txt output/scene08_square.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene10_smallpt_star.txt output/scene10_star.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene14_star_yellow.txt output/scene14_star_yellow.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene22_cow_texture.txt output/scene22_cow_texture.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene24_dof.txt output/scene24_dof.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene25_glossy.txt output/scene25_glossy.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene28_bump.txt output/scene28_bump.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene33_kitten.txt output/scene33_kitten.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene34_move.txt output/scene34_move.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene36_hall_of_fame.txt output/scene36_hall_of_fame.bmp $METHOD $SAMPLES
build/PA1 testcases/scene40_vase.txt output/scene40_vase.bmp $METHOD $SAMPLES