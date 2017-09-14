#************************************************************
#Aaron Ewall-Wice
#root.zarantan@gmail.com
#September 11th 2017
#Script to make comparitive power spectrum plots. 
#************************************************************
import numpy as np
import matplotlib.pyplot as plt
import numpy.fft as fft
import scipy.signal as sig
import copy
import scipy
import healpy as hp
import astropy.fits as fits
import getPower
import scipy.optimize as op
import scipy.signal as signal
import scipy.integrate as integrate
import scipy.interpolate as interp
import gainData as gd
import argparse
import delayGridding.py as dg
import glob

PI=np.pi
C=299792458.
LITTLEH=0.68

def ps_line_model(z0,lw,color,label):
    
    

def ps_line(simfile,bandpassfile,beamPfile,z0,deltaz,blindex,lst,ax,
            lw=2,color='k',ls='-',wind='Blackman-Harris',label=None,model=False,horizon=False):
    '''
    function to plot power spectrum line
    Args:
    simfile, string with complete path to PRISIM output.
    bandfile, string with complete path to CST simulation of antenna bandpass
    beamPfile, string with complete path to text file containing beam integrals.
    z0, center redshift of observed band.
    deltaz, redshift width of observation.
    blindex, index of baseline to plot. 
    lst, lst (in hours) to plot. 
    ax, axis to plot line on
    lw, width of plot line
    color, color of plot line.
    ls, style of plot line. 
    label, label associated with plot
    Returns:
    line object. 
    '''
    #load simulation file
    data=np.load(simfile)
    freqs=data['freq']
    df=np.abs(freqs[1]-freqs[0])
    beamIntegrals=np.load(beamPfile)
    lstbin=np.where(np.abs(data['lst']-lst)\
                    ==np.abs(data['lst']-lst).min())[0][0]
    f0=cosmology.z2f(z0)
    uvw_l=data['bl_length'][blindex]*c/f0
    #determine the correct bandwidth based on window function
    window_t=signal.blackmanharris(100)
    if wind='Nithya':
        window_t=signal.fftconvolve(window_t,window_t,mode='full')
        window_t.append(window_t,window_t[-1])
    wnorm=1./np.sqrt(np.mean(np.abs(window_t)**2.))
    bwidth=-cosmology.dz2df(wnorm*delta_z,z0)
    flow=f0-bwidth/2.
    fhigh=f0+bwidth/2.
    select=n.logical_and(freqs>=flow,freqs<fhigh)
    freqs_select=freqs[select]
    nf=len(freqs_select)
    if np.mod(nf,2)==1:
        maxind=np.where(select)[0].max()
        select[maxind+1]=True
        nfs+1
        freqs_select=freqs[select]
    band_select=freqs_select.max()-freqs_select.min()
    delays=np.arange(-nf/2,nf/2)/(nf*df)
    kparas=cosmology.eta2kpara(delays,z0)/LITTLEH
    band_pass=gd.GainData(bandpassfile,fileType='CST_TimeTrace',fMin=0.05,fMax=0.150)
    _,kernel=band_pass.interpolate_subband(nf,df*1e-9,f0*1e-9)
    kernel=np.abs(kernel)**2.
    dtv=dg.delayTransformVisibilities(data['skyvis_freq'][[blindex],select,lstbin],df,kernel=kernel)
    integralchan=np.where(np.abs(beamIntegrals[:,0]*1e6-f0)==np.min(np.abs(beamIntegrals[:,0]*1e6-f0)))[0][0]
    ps=dg.delaySq2Ps(dtv,f0,beamIntegrals[integralchan,2])
    ps_line=ax.plot(kparas,ps,ls=ls,lw=lw,color=color,label=label)
    if horizon:
        horzn=cosmology.eta2kpara(data['bl_length'][blindex]/C,z0)/LITTLEH
        ax.axvline(horzn,color='k',ls='--')
        ax.axvline(-horzn,color='k',ls='--')
    if model:
        kvals=np.sqrt(kparas**2.+(cosmology.u2kperp(data['bl_length'][blindex]/C*f0,z0)/LITTLEH)**2)
        zvals=np.ones_like(kvals)*z0
        model_ps=getPower.ps21(kvals,zvals)
        model_line=ax.plot(kparas,model_ps,color='k',lw=4,label='21cmFAST')
    return 

parser=argpares.ArgumentParser(description='Plot delay power spectrum for PRISIM simulations with and without bandpass.')
parser.add_argument('--input','-i',dest='input',type=str,
                    help=('Name of config file'))
parser.add_argument('--output','-o',dest='output',type=str,
                    help=('Name of output file. Default None will cause no output to be written.'), default=None)

parser.add_argument('--window','-w',dest='window',type=str,
                    help=('Type of window-function to use. (Blackman-Harris,Nithya). Default=Blackman-Harris'),default='Blackman-Harris')
parser.add_argument('--title','-t',dest='title',type=str,default='',help='Title for the plot.')

args=parser.parse_args()

inputfile=args.input
output=args.output
window=args.window
title=args.title

simfiles=[]
bandpassfiles=[]
beamPfiles=[]
zs=[]
deltazs=[]
blindices=[]
lsts=[]
lws=[]
colors=[]
lss=[]
windows=[]
labels=[]
modellist=[]

firstline=True
wdir=''
for line in open(inputfile).readlines():
    if not '#' in line:
        if firstline:
            wdir=line[:-1]
        line_items=line.split(',')
        simfiles.append(wdir+line_items[0])
        bandpassfiles.append(wdir+line_items[1])
        beamPfiles.append(wdir+line_items[2])
        zs.append(float(line_items[3]))
        deltazs.append(float(line_items[4]))
        blindices.append(int(line_items[5]))
        lsts.append(float(line_items[6]))
        lws.append(float(line_items[7]))
        colors.append(line_items[8])
        lss.append(line_items[9])
        windows.append(line_items[10])
        labels.append(line_items[-1])

fig=plt.figure()
ax=fig.add_axes([.1,.1,.8,.8])


for (simfile,bandfile,beamfile,z,deltaz
     ,blindex,lst,lw,color,ls,window,model) in zip(simfiles,bandpassfiles,beamPfiles,zs,deltazs,blindices,lsts,
                                                   lws,colors,lss,windows,modellist):
    simfile=glob.glob(simfile+'/*')[0]+'/simvis/simvis.npz'
    ps_line(simfile,bandpassfile,beamfile,z,deltaz,blindex,lst,ax,lw,color,ls,window,label,model)

ax.grid()
ax.set_xlabel('$k_\\parallel$ ($h$Mpc$^{-1}$)')
ax.set_ylabel('P(k) $h^{-3}Mpc$^3$mK$^2$')
ax.legend(loc='best')

if not title is None:
    ax.set_title(title)

if not output is None:
    plt.savefig(output,bbox_inches='tight')

plt.show()

    



        
    

