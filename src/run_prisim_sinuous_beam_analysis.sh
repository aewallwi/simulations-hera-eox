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
	if [ $bh != "50" ] || [ ! $gr != "0.8" ]
	then
	    for z in "${zs[@]}"; do 
		mpirun -n 2 python ~/Python/analysis/bin/run_prisim.py -i /Users/aaronew/Dropbox_MIT/Science/simulations-hera-eox/data/simulations/hera_19_sinuous_gr_${gr}_bh_${bh}_z${z}.yaml
	    done
	fi
    done
done


