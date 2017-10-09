#!/bin/bash
#run prisim simulations for vivaldi feed
source activate prisim
zs=( 06 07 08 09 10 11 12 13 14 15 16 )
export TMPDIR=/tmp/
for z in "${zs[@]}"; do 
    mpirun -n 2 python /Users/ewallwic/miniconda2/envs/prisim/bin/run_prisim.py -i /Users/ewallwic/Dropbox_MIT/Science/simulations-hera-eox/data/simulations/hera_19_vivaldi_z${z}.yaml
done
