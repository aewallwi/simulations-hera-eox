#!/bin/bash
source ~/Python/analysis/bin/activate

#python cst2prisim_fits.py -b beamlist_short.txt -f freq_list.txt -o /Users/aaronew/Dropbox_MIT/Science/simulations-hera-eox/data/simulations/CTP_beams_tall/ 

#python cst2prisim_fits.py -b beamlist_tall.txt -f freq_list.txt -o /Users/aaronew/Dropbox_MIT/Science/simulations-hera-eox/data/simulations/CTP_beams_short/

python cst2prisim_fits.py -b beamlist_sinuous.txt -f freq_list_sinuous.txt -o /Users/aaronew/Dropbox_MIT/Science/simulations-hera-eox/data/simulations/Sinuous_Beams/
