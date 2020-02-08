#!/bin/bash
set -euo pipefail
export LC_ALL=C
export LANG=C

OUTPUT="2d_animation.gif"

cd $(dirname $0)/../result
N=$(ls | sed 's/dambreaking2d\.txt\.\([0-9]*\)$/\1/g' | sort -n | tail -n 1)

echo "
  set term gif animate optimize delay 10 size 1000,350
  set size ratio -1
  set noborder
  set notics
  unset colorbox
  set xrange [] writeback
  set yrange [] writeback
  set palette defined ( 1 '#2d5ae0', 2 '#333333' )
  set output '${OUTPUT}'
  do for [i=0:$N] {
    plot sprintf('dambreaking2d.txt.%d', i) u 1:2:3 lc palette w dots notitle
    set xrange restore
    set yrange restore
  }
  " | gnuplot -

echo "Generated animation file: ${PWD}/${OUTPUT}"
