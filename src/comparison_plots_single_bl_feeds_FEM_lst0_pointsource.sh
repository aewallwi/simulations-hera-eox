#!/bin/bash
source activate hera_antenna

python plot_power_spec_sim_config.py --input compare_feeds_z08_single_bl_FEM_lst0_pointsources.conf --output feeds_fem_z08_bl2.png --b 2
python plot_power_spec_sim_config.py --input compare_feeds_z08_single_bl_FEM_lst0_pointsources.conf --output feeds_fem_z08_bl12.png --b 12
python plot_power_spec_sim_config.py --input compare_feeds_z08_single_bl_FEM_lst0_pointsources.conf --output feeds_fem_z08_bl9.png --b 9
