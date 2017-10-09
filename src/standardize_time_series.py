#************************************************************
#Aaron Ewall-Wice
#aaron.m.ewall-wice@jpl.nasa.gov
#This script is meant to perform interpolation (patch up) the
#time-domain files 
#************************************************************
import scipy.interpolate as interp
import numpy as np
import argparse
DEBUG=True

#Read CST time trace file
def readTrace(fileName,comment=''):
    dataFile=open(fileName)
    dataLines=dataFile.readlines()
    header=dataLines[:2]
    if('ns' in header[0]):
        tFactor=1.
    if('ms' in header[0]):
        tFactor=1e6
    if('micro' in header[0]):
        tFactor=1e3
    if('sec' in header[0]):
        tFactor=1e9
    inputTrace=[]
    outputTrace1=[]
    outputTrace2=[]
    lNum=0
    hout1=''
    hout2=''
    hin=''
    while lNum <len(dataLines):
        if('o1' in dataLines[lNum] or 'Port1' in dataLines[lNum]):
            thisTrace=outputTrace1
            hout1=dataLines[lNum]+dataLines[lNum+1]
            lNum+=2
        elif('o2' in dataLines[lNum] or 'Port2' in dataLines[lNum]):
            thisTrace=outputTrace2
            hout2=dataLines[lNum]+dataLines[lNum+1]
            lNum+=2
        elif('i1' in dataLines[lNum]):
            thisTrace=inputTrace
            hin=dataLines[lNum]+dataLines[lNum+1]
            lNum+=2
            dtype='Terminal Excitation'
        elif('Plane wave' in dataLines[lNum]):
            thisTrace=inputTrace
            hin=dataLines[lNum]+dataLines[lNum+1]
            lNum+=2
            dtype='PlaneWave Excitation'
        else:
            entry=dataLines[lNum].split()
            if(len(entry)==2):
                thisTrace.append([float(entry[0]),float(entry[1])])
            lNum+=1
    inputTrace=np.array(inputTrace)
    outputTrace1=np.array(outputTrace1)
    outputTrace2=np.array(outputTrace2)
    if len(inputTrace)>0:
        inputTrace[:,0]*=tFactor
    if len(outputTrace1)>0:
        outputTrace1[:,0]*=tFactor
    if(len(outputTrace2)>0):
        outputTrace2[:,0]*=tFactor
    if np.mod(len(inputTrace),2)==1:
        if len(outputTrace1)>0:
            outputTrace1=outputTrace1[:-1,:]
        if len(outputTrace2)>0:
            outputTrace2=outputTrace2[:-1,:]
        if len(inputTrace)>0:
            inputTrace=inputTrace[:-1,:]
    return [inputTrace,outputTrace1,outputTrace2],[hin,hout1,hout2]

parser=argparse.ArgumentParser(description='Standardize interpolation time-series.')
parser.add_argument('--inputsignal','-i',dest='input',type=str,
                    help=('Name of file for input wave'))
parser.add_argument('--outputsignal','-o',dest='output',type=str,
                    help=('Name of file for output wave'),default=None)
parser.add_argument('--outputname','-n',dest='name',type=str,
                    help=('Name of file to write output too'))
args=parser.parse_args()
inputfile=args.input
outputfile=args.output
outputname=args.name
if outputfile is None:
    traces,headers=readTrace(inputfile)
    tracei=traces[0]
    traceo=traces[1]
    headeri=headers[0]
    headero=headers[1]
else:
    tracesi,headersi=readTrace(inputfile)
    traceso,headerso=readTrace(outputfile)
    tracei=tracesi[0]
    traceo=traceso[1]
    headeri=headersi[0]
    headero=headerso[1]
tmin=np.max([tracei[0,0],traceo[0,0]])
tmax=np.min([tracei[-1,0],traceo[-1,0]])
nsamples=len(tracei)
#make sure interpolation axis is even in length.
if np.mod(nsamples,2)==1:
    nsamples+=1
tinterp=np.linspace(tmin,tmax,nsamples,endpoint=False)
#perofmr interpolation
interp_i=interp.interp1d(tracei[:,0],tracei[:,1])(tinterp)
interp_o=interp.interp1d(traceo[:,0],traceo[:,1])(tinterp)
#write everything out to a single file.
ofile=open(outputname,'wb')
ofile.write(headeri)
#if DEBUG:
#    import matplotlib.pyplot as plt
#    plt.plot(tracei[:,0],tracei[:,1],color='k')
#    plt.plot(traceo[:,0],traceo[:,1],color='r')
#    print('max output=%.1f'%(traceo[:,1].max()))
#    print('max output interp=%.1f'%(interp_o.max()))
#    plt.plot(tinterp[::10],interp_i[::10],marker='o',color='k',ls='None')
#    plt.plot(tinterp[::10],interp_o[::10],marker='o',color='r',ls='None')
#    plt.show()
for lnum in range(nsamples):
    ofile.write('%.10e               %.10e\n'%(tinterp[lnum],interp_i[lnum]))
ofile.write(headero)
for lnum in range(nsamples):
    ofile.write('%.10e               %.10e\n'%(tinterp[lnum],interp_o[lnum]))
ofile.close()






    

                    
