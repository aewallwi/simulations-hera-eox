#!/bin/bash
#Compare feeds with time domain that includes the HERA FEM. 
#include convertible feedk
source activate hera_antenna
python plot_power_spec_sim_config.py --input compare_feeds_z06_single_bl_FEM_convertible.conf --output feeds_fem_convertible_z06.png
python plot_power_spec_sim_config.py --input compare_feeds_z07_single_bl_FEM_convertible.conf --output feeds_fem_convertible_z07.png
python plot_power_spec_sim_config.py --input compare_feeds_z08_single_bl_FEM_convertible.conf --output feeds_fem_convertible_z08.png
python plot_power_spec_sim_config.py --input compare_feeds_z09_single_bl_FEM_convertible.conf --output feeds_fem_convertible_z09.png
python plot_power_spec_sim_config.py --input compare_feeds_z10_single_bl_FEM_convertible.conf --output feeds_fem_convertible_z10.png
python plot_power_spec_sim_config.py --input compare_feeds_z11_single_bl_FEM_convertible.conf --output feeds_fem_convertible_z11.png
python plot_power_spec_sim_config.py --input compare_feeds_z12_single_bl_FEM_convertible.conf --output feeds_fem_convertible_z12.png
python plot_power_spec_sim_config.py --input compare_feeds_z13_single_bl_FEM_convertible.conf --output feeds_fem_convertible_z13.png
python plot_power_spec_sim_config.py --input compare_feeds_z14_single_bl_FEM_convertible.conf --output feeds_fem_convertible_z14.png
python plot_power_spec_sim_config.py --input compare_feeds_z15_single_bl_FEM_convertible.conf --output feeds_fem_convertible_z15.png
python plot_power_spec_sim_config.py --input compare_feeds_z16_single_bl_FEM_convertible.conf --output feeds_fem_convertible_z16.png

