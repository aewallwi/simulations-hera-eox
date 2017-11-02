#!/bin/bash
source activate hera_antenna

ttdir=/Users/ewallwic/Dropbox_MIT/Science/HIRAX/HIRAX_adaptive_mesh_results/

chdir=${ttdir}HERA_Convertible_high_resend/
cldir=${ttdir}HERA_Convertible_low_resend/


vdir=${ttdir}Vivaldi/
ifile=hirax_adaptive_mesh_i1.txt
fbase=hirax_adaptive_mesh
labels=( pass1 pass2 pass3 pass4 pass5 pass6 )
for label in "${labels[@]}";
do
    python standardize_time_series.py -i ${ttdir}/${ifile} -o ${ttdir}${fbase}_${label}_o1.txt -n ${ttdir}${fbase}_${label}_interp.txt
done

python standardize_time_series.py -i ${ttdir}/terminal_excitation_80dB.txt -o ${ttdir}terminal_excitation_80dB.txt -n ${ttdir}terminal_excitation_80dB_interp.txt



