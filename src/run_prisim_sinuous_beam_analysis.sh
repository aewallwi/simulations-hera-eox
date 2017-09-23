#!/bin/bash
#source ~/Python/analysis/bin/activate
source activate prisim
growthrates=( 0.5 0.8 )
bheigts=( 50 90 )
zs=( 13 14 15 16 )
#models=( tall )
#heights=( 4 )
#zs=( 13 )
export TMPDIR=/tmp/
for gr in "${growthrates[@]}"; do
    for bh in "${bheights[@]}"; do
	for z in "${zs[@]}"; do 
	    mpirun -n 2 python ~/Python/analysis/bin/run_prisim.py -i /Users/aaronew/Dropbox_MIT/Science/simulations-hera-eox/data/simulations/hera_19_ctp_${model}_z${z}_h4p${h}m.yaml
	done
    done
done


