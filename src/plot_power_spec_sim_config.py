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
from astropy.io import fits
import getPower
import scipy.optimize as op
import scipy.signal as signal
import scipy.integrate as integrate
import scipy.interpolate as interp
import gainData as gd
import argparse
import delayGridding as dg
import glob
import cosmology

PI=np.pi
C=299792458.
LITTLEH=0.68

def ps_line(simfile,bandpassfile,beamPfile,z0,bwidth,blindex,lst,ax,
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
    #print('lst=%s'%(data['lst']))
    lstbin=np.where(np.abs(data['lst']-lst)\
                    ==np.abs(data['lst']-lst).min())[0][0]
    f0=cosmology.z2f(z0)
    uvw_l=data['bl_length'][blindex]*C/f0
    #determine the correct bandwidth based on window function
    window_t=signal.blackmanharris(100)
    if wind=='Nithya':
        window_t=signal.fftconvolve(window_t,window_t,mode='full')
        window_t=np.append(window_t,window_t[-1])
    #window_t/=window_t.max()
    #wnorm=np.sqrt(np.mean(np.abs(window_t)**2.))
    #print('deltaz=%s'%(deltaz))
    #bwidth=-cosmology.dz2df(wnorm*deltaz,z0)
    flow=f0-bwidth/2.
    fhigh=f0+bwidth/2.
    select=np.logical_and(freqs>=flow,freqs<fhigh)
    freqs_select=freqs[select]
    nf=len(freqs_select)
    if np.mod(nf,2)==1:
        maxind=np.where(select)[0].max()
        select[maxind+1]=True
        nf=nf+1
        freqs_select=freqs[select]
    window_t=signal.blackmanharris(nf)
    if window=='Nithya':
        window_t=signal.blackmanharris(nf/2)
        window_t=signal.fftconvolve(window_t,window_t,mode='full')
        window_t=np.append(window_t,window_t[-1])
        window_t/=window_t.max()
        
    band_select=freqs_select.max()-freqs_select.min()
    delays=np.arange(-nf/2,nf/2)/(nf*df)
    kparas=cosmology.eta2kpara(delays,z0)/LITTLEH
    if not bandpassfile=='':
        band_pass=gd.GainData(bandpassfile,fileType='CST_TimeTrace',fMin=0.05,fMax=0.150)
        _,kernel=band_pass.interpolate_subband(nf,df*1e-9,f0*1e-9)
        kernel=np.abs(kernel)**2.
        figk=plt.figure()
        axk=figk.add_axes([.1,.1,.8,.8])
        axk.plot(freqs_select,kernel)
    else:
        kernel=np.ones(nf)
    print('blindex=%d'%blindex)
    print('bl_len=%.2f'%data['bl_length'][blindex])
    vis=np.array([data['skyvis_freq'][[blindex],select,lstbin]]).T
    #print(vis.shape)
    dtv=np.abs(dg.delayTransformVisibilities(vis,df,kernel=kernel,wind=window))**2.
    integralchan=np.where(np.abs(beamIntegrals[:,0]*1e6-f0)==np.min(np.abs(beamIntegrals[:,0]*1e6-f0)))[0][0]
    ps=dg.delaySq2Ps(dtv,f0,beamIntegrals[integralchan,2],bwidth).squeeze()
    #print('ps=%s'%(ps))
    ps_line=ax.plot(kparas,ps,ls=ls,lw=lw,color=color,label=label)
    if horizon:
        horzn=cosmology.eta2kpara(data['bl_length'][blindex]/C,z0)/LITTLEH
        print('hrzn=%.2f'%horzn)
        ax.axvline(horzn,color='k',ls='--')
        ax.axvline(-horzn,color='k',ls='--')
    if model:
        kvals=np.sqrt(kparas**2.+(cosmology.u2kperp(data['bl_length'][blindex]/C*f0,z0)/LITTLEH)**2)
        zvals=np.ones_like(kvals)*z0
        model_ps=getPower.ps21(kvals,zvals)
        model_line=ax.plot(kparas,model_ps,color='k',lw=4,label='21cmFAST')
    return kernel*window_t

parser=argparse.ArgumentParser(description='Plot delay power spectrum for PRISIM simulations with and without bandpass.')
parser.add_argument('--input','-i',dest='input',type=str,
                    help=('Name of config file'))
parser.add_argument('--output','-o',dest='output',type=str,
                    help=('Name of output file. Default None will cause no output to be written.'), default=None)
parser.add_argument('--title','-t',dest='title',type=str,default='',help='Title for the plot.')

args=parser.parse_args()

inputfile=args.input
output=args.output
title=args.title

simfiles=[]
bandpassfiles=[]
beamPfiles=[]
zs=[]
bwidths=[]
blindices=[]
lsts=[]
lws=[]
colors=[]
lss=[]
windows=[]
labels=[]
modellist=[]
horizons=[]

firstline=True
wdir=''
for line in open(inputfile).readlines():
    if not '#' in line:
        if firstline:
            wdir=line[:-1]
            #print('wdir=%s'%(wdir))
            firstline=False
        else:
            line_items=line.split(',')
            simfiles.append(wdir+line_items[0])
            if not line_items[1]=='':
                bandpassfiles.append(wdir+line_items[1])
            else:
                bandpassfiles.append(line_items[1])
            beamPfiles.append(wdir+line_items[2])
            zs.append(float(line_items[3]))
            bwidths.append(float(line_items[4]))
            blindices.append(int(line_items[5]))
            lsts.append(float(line_items[6]))
            lws.append(float(line_items[7]))
            colors.append(line_items[8])
            lss.append(line_items[9])
            windows.append(line_items[10])
            labels.append(line_items[11])
            horizons.append(line_items[12])
            modellist.append(bool(line_items[13]))

fig=plt.figure()
myax=fig.add_axes([.1,.1,.8,.8])
#print simfiles
#print bandpassfiles
#print beamPfiles
#print zs
#print deltazs
#print blindices
#print lsts
#print lws
#print colors
#print lss
#print windows
#print modellist


for (simfile,bandfile,beamfile,z,bwidth
     ,blindex,lst,lw,color,ls,window,
     model,label,horizon) in zip(simfiles,bandpassfiles,beamPfiles,
                                 zs,bwidths,blindices,lsts,
                                 lws,colors,lss,windows,
                                 modellist,labels,horizons):
    simfile=glob.glob(simfile+'/*')[0]+'/simdata/simvis.npz'
    windowkernel=ps_line(simfile,bandfile,beamfile,z,bwidth,blindex,lst,myax,lw,color,ls,window,label,model,horizon)
    windowkernel/=np.abs(windowkernel).max()

myax.grid(b=True,which='major')
myax.grid(b=True,which='minor')

myax.set_xlim(-.6,.6)
#if all zs in plot line list are the same, create a delay axis.
samez=np.all(np.array(zs)==zs[0])
if samez:
    myax_t=myax.twiny()
    ticks=myax.get_xticks()
    labels=[]
    for tick in ticks:
        labels.append('%d'%(1e9*cosmology.kpara2eta(LITTLEH*tick,zs[0])))
    myax_t.set_xticks(ticks)
    myax_t.set_xticklabels(labels)
    myax_t.set_xlim(myax.get_xlim())
    myax_t.set_ylim(1e3,1e17)

sameb=np.all(np.array(blindices)==blindices[0])
samedz=np.all(np.array(bwidths)==bwidths[0])
windowkernel/=windowkernel.max()
print('mean=%.2f'%np.sqrt(np.mean(windowkernel**2.)))
delta_z_eff=cosmology.df2dz(bwidths[0],zs[0])*np.sqrt(np.mean((windowkernel)**2.))



if sameb and samedz and samez:
    data=np.load(simfile)
    myax.text(-0.6,1e8,'b=%.2f\n$\\Delta z=%.1f$\nLST=%.1f hours\nz=%.2f'%(data['bl_length'][blindex],delta_z_eff,lsts[0]/15.,zs[0]),fontsize=20)
        


    
    
myax.set_yscale('log')
myax.set_xlabel('$k_\\parallel$ ($h$Mpc$^{-1}$)')
myax.set_ylabel('P(k) $h^{-3}$Mpc$^3$mK$^2$')
myax.legend(loc='best')
#add delay axis



if not title is None:
    myax.set_title(title,y=1.05)

if not output is None:
    plt.savefig(output,bbox_inches='tight')

plt.show()

    



        
    

