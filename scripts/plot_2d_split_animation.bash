#!/bin/bash

OUTPUT="2d_split_animation.png"
TMP_FILE="tmp_sph"

cd $(dirname $0)/../result
N=$(ls | sed 's/dambreaking2d\.txt\.\([0-9]*\)$/\1/g' | sort -n | tail -n 1)
echo "
  set term pngcairo size 1500,550
  set size ratio -1
  set noborder
  set notics
  unset colorbox
  set xrange [] writeback
  set yrange [] writeback
  # set palette rgbformula 22,13,-31
  # set palette defined ( 1 '#2d5ae0', 2 '#333333' )
  load '../scripts/viridis.pal'
  do for [i=0:$N] {
    set output sprintf('${TMP_FILE}_%d.png', i)
    plot sprintf('dambreaking2d.txt.%d', i) u 1:2 lc 'gray50' pt 7 notitle, \
         sprintf('particle_tree2d.txt.%d', i) u ((\$4+\$2)/2):((\$5+\$3)/2):((\$4-\$2)/2):((\$5-\$3)/2):1 w boxxyerrorbars fs transparent solid 0.3 lw 1.3 palette notitle
    # plot sprintf('dambreaking2d.txt.%d', i) u 1:2 lc 'gray50' pt 7 notitle, \
    #      sprintf('particle_tree2d.txt.%d', i) u ((\$4+\$2)/2):((\$5+\$3)/2):((\$4-\$2)/2):((\$5-\$3)/2) w boxxyerrorbars fs transparent solid 0.3 lw 1.3 notitle
    # plot sprintf('dambreaking2d.txt.%d', i) u 1:2:3 lc palette pt 7 ps 0.3 notitle
    set xrange restore
    set yrange restore
  }
  " | gnuplot -

apngasm $OUTPUT $(ls ${TMP_FILE}_*.png) 1 10
# rm ${TMP_FILE}_*.png
