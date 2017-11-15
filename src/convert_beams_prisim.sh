#!/bin/bash
#source ~/Python/analysis/bin/activate
source activate hera_antenna
python cst2prisim_fits.py -b beamlist_short.txt -f freq_list.txt -o /Users/ewallwic/Dropbox_MIT/Science/simulations-hera-eox/data/simulations/CTP_beams_tall/ -s True

python cst2prisim_fits.py -b beamlist_tall.txt -f freq_list.txt -o /Users/ewallwic/Dropbox_MIT/Science/simulations-hera-eox/data/simulations/CTP_beams_short/ -s True
