#!/bin/bash
source activate hera_antenna
ttdir=/Users/ewallwic/Dropbox_MIT/Science/simulations-hera-eox/data/timetraces/
ttdir1=/Users/ewallwic/Dropbox_MIT/Science/AntennaDataHERA/Cambridge_UK_results/AdditionalData/

chdir=${ttdir}HERA_Convertible_high_resend/
cldir=${ttdir}HERA_Convertible_low_resend/


vdir=${ttdir}Vivaldi/

python standardize_time_series.py -i ${chdir}Plane_wave_major.txt -o ${chdir}o1.txt -n ${chdir}convertible_high_interp.txt
python standardize_time_series.py -i ${cldir}Plane_wave_major.txt -o ${cldir}o1.txt -n ${cldir}convertible_low_interp.txt
python standardize_time_series.py -i ${vdir}Excitation_time_signal-50-250MHz-500ns.txt -o ${vdir}Vivaldi+dish-model3-output-voltage-50-250MHz-newFEM.txt -n ${vdir}vivaldi_interp.txt
python standardize_time_series.py -i ${vdir}Excitation_time_signal-50-250MHz-500ns.txt -o ${vdir}Sinuous+dish-output_port_X-voltage-50-250MHz-newFEM.txt -n ${vdir}sinuous_interp.txt
python standardize_time_series.py -i ${vdir}Excitation_time_signal-50-250MHz-500ns.txt -o ${vdir}HERA4.9m-output-voltage-50-250MHz-newFEM.txt -n ${vdir}hera_interp.txt




