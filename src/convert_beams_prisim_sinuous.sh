#!/bin/bash
#source ~/Python/analysis/bin/activate
#source activate hera_antenna
#python cst2prisim_fits.py -b beamlist_short.txt -f freq_list.txt -o /Users/ewallwic/Dropbox_MIT/Science/simulations-hera-eox/data/simulations/CTP_beams_tall/ 

#python cst2prisim_fits.py -b beamlist_tall.txt -f freq_list.txt -o /Users/ewallwic/Dropbox_MIT/Science/simulations-hera-eox/data/simulations/CTP_beams_short/

python cst2prisim_fits.py -b beamlist_sinuous.txt -f freq_list_sinuous.txt -i False -r False -q True -o /Users/ewallwic/Dropbox_MIT/Science/simulations-hera-eox/data/simulations/Sinuous_Beams
