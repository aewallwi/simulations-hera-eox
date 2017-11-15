#!/bin/bash
source activate hera_antenna
python cst2prisim_fits.py -b beamlist_highband.txt -f freq_list_highband.txt -i True -r False -q False -o /Users/ewallwic/Dropbox_MIT/Science/simulations-hera-eox/data/simulations/Beams_HighBand/ -s True
