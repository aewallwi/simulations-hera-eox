#!/bin/bash
#Compare feeds with time domain that includes the HERA FEM. 
#include convertible feedk
source activate hera_antenna
ouputdir=/Users/ewallwic/Dropbox_MIT/Science/simulatins-hera-eox/analysis/
python plot_power_spec_sim_config.py --input compare_feeds_z06_single_bl_FEM_convertible.conf --output ${outputdir}feeds_fem_convertible_z06_FEM.png
python plot_power_spec_sim_config.py --input compare_feeds_z07_single_bl_FEM_convertible.conf --output ${outputdir}feeds_fem_convertible_z07_FEM.png
python plot_power_spec_sim_config.py --input compare_feeds_z08_single_bl_FEM_convertible.conf --output ${outputdir}feeds_fem_convertible_z08_FEM.png
python plot_power_spec_sim_config.py --input compare_feeds_z09_single_bl_FEM_convertible.conf --output ${outputdir}feeds_fem_convertible_z09_FEM.png
python plot_power_spec_sim_config.py --input compare_feeds_z10_single_bl_FEM_convertible.conf --output ${outputdir}feeds_fem_convertible_z10_FEM.png
python plot_power_spec_sim_config.py --input compare_feeds_z11_single_bl_FEM_convertible.conf --output ${outputdir}feeds_fem_convertible_z11_FEM.png
python plot_power_spec_sim_config.py --input compare_feeds_z12_single_bl_FEM_convertible.conf --output ${outputdir}feeds_fem_convertible_z12_FEM.png
python plot_power_spec_sim_config.py --input compare_feeds_z13_single_bl_FEM_convertible.conf --output ${outputdir}feeds_fem_convertible_z13_FEM.png
python plot_power_spec_sim_config.py --input compare_feeds_z14_single_bl_FEM_convertible.conf --output ${outputdir}feeds_fem_convertible_z14_FEM.png
python plot_power_spec_sim_config.py --input compare_feeds_z15_single_bl_FEM_convertible.conf --output ${outputdir}feeds_fem_convertible_z15_FEM.png
python plot_power_spec_sim_config.py --input compare_feeds_z16_single_bl_FEM_convertible.conf --output ${outputdir}feeds_fem_convertible_z16_FEM.png
#now do 100 ohm
python plot_power_spec_sim_config.py --input compare_feeds_z06_single_bl_100Ohm_convertible.conf --output ${outputdir}feeds_fem_convertible_z06_100_Ohm.png
python plot_power_spec_sim_config.py --input compare_feeds_z07_single_bl_100Ohm_convertible.conf --output ${outputdir}feeds_fem_convertible_z07_100_Ohm.png
python plot_power_spec_sim_config.py --input compare_feeds_z08_single_bl_100Ohm_convertible.conf --output ${outputdir}feeds_fem_convertible_z08_100_Ohm.png
python plot_power_spec_sim_config.py --input compare_feeds_z09_single_bl_100Ohm_convertible.conf --output ${outputdir}feeds_fem_convertible_z09_100_Ohm.png
python plot_power_spec_sim_config.py --input compare_feeds_z10_single_bl_100Ohm_convertible.conf --output ${outputdir}feeds_fem_convertible_z10_100_Ohm.png
python plot_power_spec_sim_config.py --input compare_feeds_z11_single_bl_100Ohm_convertible.conf --output ${outputdir}feeds_fem_convertible_z11_100_Ohm.png
python plot_power_spec_sim_config.py --input compare_feeds_z12_single_bl_100Ohm_convertible.conf --output ${outputdir}feeds_fem_convertible_z12_100_Ohm.png
python plot_power_spec_sim_config.py --input compare_feeds_z13_single_bl_100Ohm_convertible.conf --output ${outputdir}feeds_fem_convertible_z13_100_Ohm.png
python plot_power_spec_sim_config.py --input compare_feeds_z14_single_bl_100Ohm_convertible.conf --output ${outputdir}feeds_fem_convertible_z14_100_Ohm.png
python plot_power_spec_sim_config.py --input compare_feeds_z15_single_bl_100Ohm_convertible.conf --output ${outputdir}feeds_fem_convertible_z15_100_Ohm.png
python plot_power_spec_sim_config.py --input compare_feeds_z16_single_bl_100Ohm_convertible.conf --output ${outputdir}feeds_fem_convertible_z16_100_Ohm.png
