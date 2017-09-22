#!/bin/bash
source activate hera_antenna
python plot_power_spec_sim_config.py --input compare_ctp_z13.conf --output single_bl_z13.png
python plot_power_spec_sim_config.py --input compare_ctp_z14.conf --output single_bl_z14.png
python plot_power_spec_sim_config.py --input compare_ctp_z15.conf --output single_bl_z15.png
python plot_power_spec_sim_config.py --input compare_ctp_z16.conf --output single_bl_z16.png

