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

PI=np.pi
C=299792458.
LITTLEH=0.68


def ps_line(simfile,bandfile,z0,deltaz,blindex,lst,lw=2,color='k',ls='-'):
    '''
    function to plot power spectrum line
    '''
    



parser=argpares.ArgumentParser(description='Plot delay power spectrum for PRISIM simulations with and without bandpass.')
parser.add_argument('--input','-i',dest='input',type=str,
                    help=('Name of config file'))
parser.add_argument('--output','-o',dest='output',type=str,
                    help=('Name of output file. Default None will cause no output to be written.'), default=None)

parser.add_argument('--window','-w',dest='window',type=str,
                    help=('Type of window-function to use. (Blackman-Harris,Nithya). Default=Blackman-Harris'),default='Blackman-Harris')
parser.add_argument('--model','-m',dest='model',type=bool,help=('Determine whether you wish to plot power spectrum model. Default=True'),default=True)
args=parser.parse_args()


