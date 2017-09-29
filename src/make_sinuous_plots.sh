#!/bin/bash
source activate hera_antenna
zs=( 06 07 08 09 10 11 12 13 14 15 16 )
for z in "${zs[@]}"
do
    python plot_power_spec_sim_config.py --input compare_sinuous_z$z.conf --output compare_sinuous_z$z.png
done
