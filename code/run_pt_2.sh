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

# build/PA1 testcases/scene08_smallpt.txt output/scene08.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene08_smallpt_square.txt output/scene08_square.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene10_smallpt_square.txt output/scene10_square.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene11_box.txt output/scene11_box.bmp $METHOD $
# build/PA1 testcases/scene13_star.txt output/scene13_star.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene15_cube_and_ball.txt output/scene15_cube_and_ball.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene14_star_yellow.txt output/scene14_star_yellow.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene12.txt output/scene12.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene09_sky.txt output/scene09.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene08_smallpt_allball.txt output/scene08_allball.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene16.txt output/scene16_banny_and_ball.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene17_texture_ball.txt output/scene17_texture_ball.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene18_global.txt output/scene18_global.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene19_bunny_texture.txt output/scene19_bunny_texture.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene20_bunny_bump.txt output/20_bunny_bump.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene21_bunny_star.txt output/scene21_bunny_star.bmp $METHOD $SAMPLES
build/PA1 testcases/scene22_cow_texture.txt output/scene22_cow_texture.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene23_dof.txt output/scene23_dof.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene24_dof.txt output/scene24_dof.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene24_dof_0.txt output/scene24_dof_0.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene25_glossy.txt output/scene25_glossy.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene26_wineglass.txt output/scene26_wineglass.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene27_bump.txt output/scene27_bump.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene28_bump.txt output/scene28_bump.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene29_floor_fire.txt output/scene29_great.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene09_norm.txt output/scene09_norm.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene31_norm.txt output/scene31_norm_2.bmp $METHOD $SAMPLES 
# build/PA1 testcases/scene32_norm.txt output/scene32_norm.bmp $METHOD $SAMPLES 
# build/PA1 testcases/scene33_vase.txt output/scene33_vase.bmp $METHOD $SAMPLES 
# build/PA1 testcases/scene33_kitten.txt output/scene33_kitten.bmp $METHOD $SAMPLES
# build/PA1 testcases/scene24_dof.txt output/scene24_dof.bmp $METHOD $SAMPLES500
# build/PA1 testcases/scene37_player.txt output/scene37_player.bmp $METHOD $SAMPLES