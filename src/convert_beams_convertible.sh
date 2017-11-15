#!/bin/bash
#source ~/Python/analysis/bin/activate
source activate hera_antenna
#convert beams for convertible HERA antenna.
python cst2prisim_fits.py -b beamlist_convertible_low.txt -f freq_list_convertible_low.txt -i True -r False -q False -o /Users/ewallwic/Dropbox_MIT/Science/simulations-hera-eox/data/simulations/Convertible_Beams/ -s True
python cst2prisim_fits.py -b beamlist_convertible_high.txt -f freq_list_convertible_high.txt -i True -r False True -q False -o /Users/ewallwic/Dropbox_MIT/Science/simulations-hera-eox/data/simulations/Convertible_Beams/ -s True
