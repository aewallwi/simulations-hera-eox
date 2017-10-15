#!/bin/bash
#source ~/Python/analysis/bin/activate
source activate prisim
export TMPDIR=/tmp
models=( tall short )
heights=( 4 9 )
zs=( 13 14 15 16 )
#models=( tall )
#heights=( 4 )
#zs=( 13 )
export TMPDIR=/tmp/
for model in "${models[@]}"; do
    for h in "${heights[@]}"; do
	for z in "${zs[@]}"; do 
	    #mpirun -n 2 python ~/Python/analysis/bin/run_prisim.py -i /Users/aaronew/Dropbox_MIT/Science/simulations-hera-eox/data/simulations/hera_19_ctp_${model}_z${z}_h4p${h}m.yaml
	    mpirun -n 2 python /Users/ewallwic/miniconda2/envs/prisim/bin/run_prisim.py -i /Users/ewallwic/Dropbox_MIT/Science/simulations-hera-eox/data/simulations/hera_19_ctp_${model}_z${z}_h4p${h}m.yaml
	done
    done
done


