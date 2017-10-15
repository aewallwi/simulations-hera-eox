#!/bin/bash
#source ~/Python/analysis/bin/activate
source activate hera_antenna
#convert beams for convertible HERA antenna.
python cst2prisim_fits.py -b beamlist_vivaldi.txt -f freq_list_vivaldi.txt -i False -r False True -q False -o /Users/ewallwic/Dropbox_MIT/Science/simulations-hera-eox/data/simulations/Vivaldi_Beams/
