import numpy as np
from farfieldData import Beam
from optparse import OptionParser
import re
import glob
usage=('Converts a list of CST beam files into a list of healpix .fits files'
       'that are compatible with the prisim.')
parser=OptionParser(usage=usage)
parser.add_option('-b','--beams',dest='beamdirlist',help="text list of CST beam files")
parser.add_option('-f','--freqs',dest='freqlist',help="text list of frequencies in MHz")
parser.add_option('-o','--output',dest='outputdir',help="output directory")
parser.add_option('-r','--rotatexz',dest='rotatexz',help="rotate xz, default=True",default="True")
parser.add_option('-i','--invert',dest='invert',help="invert beam,default=True",default="True")
parser.add_option('-q','--rotatexy',dest='rotatexy',help="rotate xy, default=True",default="False")
parser.add_option('-s','--seperate',dest='seperate',helkp="True if you want seperate beams for each channel",default=False)
(options,args)=parser.parse_args()
beamdirs=open(options.beamdirlist).readlines()
fstrlist=open(options.freqlist).readlines()

assert options.rotatexz in ["True","False"]
assert options.invert in ["True","False"]
assert options.seperate in ["True","False"]

if options.rotatexz=="True":
    rotatexz=True
else:
    rotatexz=False
if options.invert=="True":
    invert=True
else:
    invert=False
if options.rotatexy=="True":
    rotatexy=True
else:
    rotatexy=False
if options.seperate=="True":
    seperate_freqs=True
else:
    seperate_freqs=False
print('invert=%s'%invert)
print('rotatexz=%s'%rotatexz)
    
for beamdir in beamdirs:
    beamdir=re.sub('\n','',beamdir)
    beamdirstr=beamdir
    flist_unorg=glob.glob(beamdir+'/*')
    #print beamdir
    #print flist_unorg
    beamdir=beamdir.split('/')
    if len(beamdir[-1])<2:
        beamname=beamdir[-2]
    else:
        beamname=beamdir[-1]
    #print('beamname=%s'%beamname)
    beam=Beam(pols=['X'],rotateY=True,rotatexz=rotatexz,invert=invert,rotatexy=rotatexy)
    #flist=[re.sub('\n','',(beamdir+'/farfield (f=%s) [1].txt')%(s)) for s in fstrlist]
    flist=[]
    #reorder flist according to frequencies
    for f in fstrlist:
        for filename in flist_unorg:
            filename=filename.split('/')[-1]
            if str(f)[:-1] in filename and not( float(f[:-1])<100 and ('1'+str(f)[:-1] in filename or '2'+str(f)[:-1] in filename)):
                flist.append(beamdirstr+filename)

    #print fstrlist
    #print flist
    
    beam.read_files(fstrlist,flist)
    freqlist=np.array([float(f) for f in fstrlist])*1e6
    if seperate_freqs:
        for fnum,freq in enumerate(freqlist):
            beam.export_fits_prisim(options.outputdir+beamname+'_'+fstrlist[fnum]+'.fits',['I'],
                                    [freqlist[fnum]])
            beam.export_beam_integrals(options.outputdir+beamname+'_integrals_'+fstrlist[fnum])

    else:
        beam.export_fits_prisim(options.outputdir+beamname+'.fits',['I'],
                                freqlist)
        beam.export_beam_integrals(options.outputdir+beamname+'_integrals')



