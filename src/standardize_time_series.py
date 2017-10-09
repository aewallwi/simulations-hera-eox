#************************************************************
#Aaron Ewall-Wice
#aaron.m.ewall-wice@jpl.nasa.gov
#This script is meant to perform interpolation (patch up) the
#time-domain files 
#************************************************************
import scipy.interpolate as interp
import numpy as np
import argparse


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
        if('o1' in dataLines[lNum]):
            thisTrace=outputTrace1
            hout1=dataLines[lNum]
            lNum+=2
        elif('o2' in dataLines[lNum]):
            thisTrace=outputTrace2
            hout2=dataLines[lNum]
            lNum+=2
        elif('i1' in dataLines[lNum]):
            thisTrace=inputTrace
            hin=dataLines[lNum]
            lNum+=2
            dtype='Terminal Excitation'
        elif('Plane wave' in dataLines[lNum]):
            thisTrace=inputTrace
            hin=dataLines[lNum]
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
    print('len(inputtrace=%d'%(len(inputTrace)))
    print('len(outputtrace=%d'%(len(outputTrace1)))
    inputTrace[:,0]*=tFactor
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
    return [inputTrace,outputTrace1,outputTrace2],dtype

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
    traceso,headerso=readTrace(inputfile)
    tracei=tracesi[0]
    traceo=traceso[1]
    headeri=tracesi[0]
    headero=traceso[1]

tmin=np.max([traci[:,0],traceo[:,0]])
tmax=np.min([traci[:,-1],traceo[:,-1]])
nsamples=len(tmin)
#make sure interpolation axis is even in length.
if np.mod(nsamples,2)==1:
    nsamples+=1
tinterp=np.linspace(tmin,tmax,nsampes,endpoint=False)
#perofmr interpolation
interp_i=interp.interp1d(tinterp,tracei)(tinterp)
interp_o=interp.interp1d(tinterp,traceo)(tinterp)
#write everything out to a single file.
ofile=open(outputname,'w')
ofile.write(headeri)
for lnum in range(nsamples):
    ofile.write('%.10f \t %.10f\n'%(tinterp[lnum],interp_i[lnum]))
ofile.write(headero)
for lnum in range(nsamples):
    ofile.write('%.10f \t %.10f\n'%(tinterp[lnum],interp_o[lnum]))
ofile.close()






    

                    
