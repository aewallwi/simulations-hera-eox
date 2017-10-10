#!/bin/bash
#run prisim simulations for vivaldi feed
source activate prisim
export TMPDIR=/tmp/

mpirun -n 2 python /Users/ewallwic/miniconda2/envs/prisim/bin/run_prisim.py -i /Users/ewallwic/Dropbox_MIT/Science/simulations-hera-eox/data/simulations/hera_19_vivaldi_z09_lst0.yaml
mpirun -n 2 python /Users/ewallwic/miniconda2/envs/prisim/bin/run_prisim.py -i /Users/ewallwic/Dropbox_MIT/Science/simulations-hera-eox/data/simulations/hera_19_sinuous_gr_0.5_bh_50_z08_lst0.yaml
mpirun -n 2 python /Users/ewallwic/miniconda2/envs/prisim/bin/run_prisim.py -i /Users/ewallwic/Dropbox_MIT/Science/simulations-hera-eox/data/simulations/hera_19_highband_z09_lst0.yaml
