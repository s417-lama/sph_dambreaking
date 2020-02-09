#!/bin/bash
set -euo pipefail
export LC_ALL=C
export LANG=C

SPH_2D=1
DOUBLE=1
REUSE_TREE=1
SCALE=8
OUTPUT_INTERVAL=0
MAX_STEP=200

export SPH_CUDA_PARALLEL=1

cd $(dirname $0)/..

echo "## Generating data..."
mkdir -p data
./scripts/gen_data_2d.py $SCALE > data/data2d.txt

echo "## Building..."
make clean
CFLAGS="-DSPH_2D=$SPH_2D -DSPH_DOUBLE=$DOUBLE -DSPH_REUSE_TREE=$REUSE_TREE -DSPH_OUTPUT_INTERVAL=$OUTPUT_INTERVAL -DSPH_DATA_SCALE=$SCALE -DSPH_MAX_STEP=$MAX_STEP" make -j

echo "## Starting..."
mkdir -p result
./sph.out
# nvprof -o profile.nvp ./sph.out

echo "## Generating an animation..."
./scripts/plot_2d_animation.bash
