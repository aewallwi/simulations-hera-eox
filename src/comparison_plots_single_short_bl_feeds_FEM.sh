#!/bin/bash
#Compare feeds with time domain that includes the HERA FEM. 
#
source activate hera_antenna
python plot_power_spec_sim_config.py --input compare_feeds_z06_single_bl_FEM.conf --output feeds_fem_shortbl_z06.png --baseline 2
python plot_power_spec_sim_config.py --input compare_feeds_z07_single_bl_FEM.conf --output feeds_fem_shortbl_z07.png --baseline 2
python plot_power_spec_sim_config.py --input compare_feeds_z08_single_bl_FEM.conf --output feeds_fem_shortbl_z08.png --baseline 2
python plot_power_spec_sim_config.py --input compare_feeds_z09_single_bl_FEM.conf --output feeds_fem_shortbl_z09.png --baseline 2
python plot_power_spec_sim_config.py --input compare_feeds_z10_single_bl_FEM.conf --output feeds_fem_shortbl_z10.png --baseline 2
python plot_power_spec_sim_config.py --input compare_feeds_z11_single_bl_FEM.conf --output feeds_fem_shortbl_z11.png --baseline 2
python plot_power_spec_sim_config.py --input compare_feeds_z12_single_bl_FEM.conf --output feeds_fem_shortbl_z12.png --baseline 2
