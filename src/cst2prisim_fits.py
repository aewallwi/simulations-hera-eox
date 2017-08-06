import numpy as np
from farfieldData import Beam
from optparse import OptionParser
import re
usage=('Converts a list of CST beam files into a list of healpix .fits files'
       'that are compatible with the prisim.')
parser=OptionParser(usage=usage)
parser.add_option('-b','--beams',dest='beamdirlist',help="text list of CST beam files")
parser.add_option('-f','--freqs',dest='freqlist',help="text list of frequencies in MHz")
parser.add_option('-o','--output',dest='outputdir',help="output directory")
(options,args)=parser.parse_args()
beamdirs=open(options.beamdirlist).readlines()
fstrlist=open(options.freqlist).readlines()
for beamdir in beamdirs:
    beamname=beamdir.split('/')[-1]
    beam=Beam(pols=['X'],rotateY=True,rotatexz=True,invert=True)
    flist=[re.sub('\n','',(beamdir+'/farfield (f=%s) [1].txt')%(s)) for s in fstrlist]
    beam.read_files(fstrlist,flist)
    beam.export_fits_prisim(options.outputdir+beamname+'.fits')




