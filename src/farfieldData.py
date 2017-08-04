import numpy as n
import healpy as hp
from astropy.io import fits
pi=n.pi
import copy
import scipy.optimize as op
c=299792458


#take a cut through a beam
def hpCut(phi,nPix,data):
    '''
    '''
    nSide=hp.npix2nside(len(data))
    output=n.zeros(nPix)
    thetaVals=n.arange(nPix/2)/(nPix/2.)*pi/2.
    thetaVals=n.hstack([n.flipud(thetaVals),thetaVals,]).T
    phiVals=n.ones(len(thetaVals))
    phi1=phi+pi
    phiVals[:nPix/2]=phi1
    phiVals[nPix/2:]=phi
    output=hp.get_interp_val(data,thetaVals,phiVals)
    return output


#rotate
def rotateBeam(inputMap,rot=[90,0,0]):
    rotator=hp.Rotator(rot=rot)
    npix=len(inputMap)
    nside=hp.npix2nside(npix)
    theta,phi=hp.pix2ang(nside,range(npix))
    newtheta,newphi=rotator(theta,phi)
    output=hp.get_interp_val(inputMap,newtheta,newphi)
    return output
    


class Beam:
    def __init__(self,nside=64,pols=['XX','Y'],rotateY=False,invert=False,rotatexz=False):
        self.nside=nside
        self.npolsOriginal=len(pols)
        self.npols=max(len(pols),2)
        self.nPix=hp.nside2npix(nside)
        self.nSide=nside
        self.invert=invert
        self.pixArea=hp.nside2pixarea(self.nSide)
        self.pols=pols
        self.rotatexz=rotatexz
        if(rotateY):
            pols.append('Y')
            
    def _interp_beam(self,beam_file,pol,chan):
        data=n.loadtxt(beam_file,skiprows=2);
        theta,phi=hp.pix2ang(self.nSide,range(self.nPix))
        theta=n.round(n.degrees(theta)).astype(int)
        phi=n.round(n.degrees(phi)).astype(int)
        self.data[pol,chan,:]=10**((data[:,2].squeeze().reshape(360,181))[phi,theta]/10.)
        self.data[pol,chan,:]/=self.data[pol,chan,:].flatten().max();    
        if self.rotatexz:
            self.data[pol,chan,:]=rotateBeam(self.data[pol,chan,:].flatten(),rot=[0,-90,0])
        if(self.invert):
            self.data[pol,chan,:]=rotateBeam(self.data[pol,chan,:].flatten(),rot=[0,180,0])
        self.data[pol,chan,theta>90.]=0.
        self.solidAngles[pol,chan]=self.pixArea*n.sum(self.data[pol,chan,:])
        self.effArea[pol,chan]=(c/(self.fAxis[chan]))**2./self.solidAngles[pol,chan]
        if(self.npolsOriginal==1):
            self.data[1,chan,:]=rotateBeam(self.data[0,chan,:].flatten())
            self.solidAngles[1,chan]=self.pixArea*n.sum(self.data[1,chan,:])
            self.effArea[1,chan]=(c/(self.fAxis[chan]))**2./(self.solidAngles[1,chan])
        if(len(self.pols)>1 and self.pols[0]=='X' and self.pols[1]=='Y'):
            self.ellipticity[chan]=n.sum((self.data[0,chan]-self.data[1,chan])**2.)/n.sum((self.data[0,chan]+self.data[1,chan])**2.)
    def read_files(self,flist,filelist):
        self.nf=len(flist)
        self.fAxis=n.zeros(self.nf)
        self.data=n.zeros((self.npols,self.nf,self.nPix))   
        self.solidAngles=n.zeros((self.npols,self.nf))
        self.effArea=n.zeros_like(self.solidAngles)
        self.ellipticity=n.zeros(self.nf)

        for m in range(self.nf):
            tempf=flist[m].split('p')
            self.fAxis[m]=float(tempf[0])*1e6
            for p in range(self.npolsOriginal):
                self._interp_beam(filelist[p*self.nf+m],p,m)
    def export_fits_prisim(self,fitsfile,pol_list,freq_list,scheme='RING',nside_out=None):
        '''
        export fits-file at channel chan and polarization pol
        Args:
        fitsfile, str, name of file to save .fits to
        pol_list, list of labels of polarizations to write
        chan_list, list of frequencies to write
        '''
        if nside_out is None:
            nside_out=self.nside
        pol_list=n.array(pol_list)
        freq_list=n.array(freq_list)
        pol_inds=[]
        freq_inds=[]
        for pol,freq in zip(pol_list,freq_list):
            assert pol in self.pols
            pol_inds.append(n.where(n.array(self.pols)==pol)[0][0])
            assert freq in self.fAxis
            freq_inds.append(n.where(n.array(self.fAxis)==freq)[0][0])
        data=self.data[:,freq_inds,:].reshape(-1,1)
        theta_out,phi_out=hp.pix2ang(self.nside,n.arange(hp.nside2npix(nside_out)))
        #freq_col=[fits.Column(name='Frequency [MHz]',format='D',array=n.array(freq_list))]
        #freq_columns=fits.ColDefs(freq_col,ascii=False)
        #freq_tbhdu = fits.BinTableHDU.from_columns(freq_col)
        #freq_tbhdu = fits.BinTableHDU.from_columns(n.array(freq_list))
        
        hduprimary=fits.PrimaryHDU()
        hduprimary.header.set('EXTNAME','PRIMARY')
        hduprimary.header.set('NEXTEN',3)
        hduprimary.header.set('FITSTYPE','IMAGE')
        hduprimary.header['NSIDE']=(nside_out,'NSIDE')
        hduprimary.header['PIXAREA']=(hp.nside2pixarea(nside_out),'pixel solid angle (steradians)')
        hduprimary.header['NEXTEN']=(3,'Number of extensions')
        hduprimary.header['NPOL'] = (len(pol_inds), 'Number of polarizations')
        hduprimary.header['SOURCE'] = ('HERA-CST', 'Source of data')
        hdulist=[hduprimary]
        fits.HDUList(hdulist).writeto(fitsfile,clobber=True)
        for pol in pol_list:
            #freq_tbhdu.header.set('EXTNAME','FREQS_{0}'.format(pol))
            freq_tbhdu=fits.ImageHDU(freq_list,name='FREQS_{0}'.format(pol))
            fits.append(fitsfile,freq_tbhdu.data,freq_tbhdu.header,verify=False)
        data_interp=n.zeros((hp.nside2npix(nside_out),len(freq_inds)))
        
        for polind,pol in zip(pol_inds,pol_list):
            for fi,freqind in enumerate(freq_inds):
                data=self.data[polind,freqind,:].flatten()
                data_interp[:,fi]=hp.get_interp_val(data,theta_out,phi_out)
            imghdu = fits.ImageHDU(data_interp, name='BEAM_{0}'.format(pol))
            imghdu.header['PIXTYPE'] = ('HEALPIX', 'Type of pixelization')
            imghdu.header['ORDERING'] = (scheme, 'Pixel ordering scheme, either RING or NESTED')
            imghdu.header['NSIDE'] = (nside_out, 'NSIDE parameter of HEALPIX')
            imghdu.header['NPIX'] = (hp.nside2npix(nside_out), 'Number of HEALPIX pixels')
            imghdu.header['FIRSTPIX'] = (0, 'First pixel # (0 based)')
            imghdu.header['LASTPIX'] = (len(data_interp)-1, 'Last pixel # (0 based)')
            fits.append(fitsfile,imghdu.data,imghdu.header,verify=False)
            #hdulist += [hdu]
            #hdulist += [freq_tbhdu.data,freq_tbhdu.header]
            #hdulist += [fits.ImageHDU([100e6], name='FREQS_{0}'.format(pol))]
                                
#outhdu = fits.HDUList(hdulist)
#outhdu.writeto(outfile, clobber =True)

                
class FarFieldData:
    def __init__(self,fListFeedOnly,fListFeedOverDish,
                 freqList,nside,
                 pols=['X','Y'],rotateY=False,
                 invertFeedOnly=True,
                 invertDish=True,
                 dDish=14.,dFocus=4.5,rotatexz=False):
        self.nf=len(freqList)
        self.nSide=nside
        self.nPix=hp.nside2npix(nside)
        self.daveFeedEfficiency=n.zeros(self.nf)
        self.daveTaper=n.zeros(self.nf)
        self.davePol=n.zeros(self.nf)
        self.daveBeamEfficiency=n.zeros(self.nf)
        self.dDish=dDish
        self.dFocus=dFocus
        self.thetaEdge=2.*n.arctan(self.dDish/(4.*self.dFocus))
        self.beamFeedAndDish=Beam(nside,copy.deepcopy(pols),rotateY,invertDish,rotatexz=rotatexz)
        self.beamFeedAndDish.read_files(freqList,fListFeedOverDish)
        self.beamFeedOnly=Beam(nside,pols,rotateY,invertFeedOnly)
        self.beamFeedOnly.read_files(freqList,fListFeedOnly)
        self.fAxis=self.beamFeedOnly.fAxis
        self.beamFWHM=n.zeros(len(self.fAxis))
        self.gain=4.*pi/self.beamFeedAndDish.solidAngles[0,:]
        theta,phi=hp.pix2ang(self.nSide,range(self.nPix))


        self.nSide=nside
        self.nPix=hp.nside2npix(self.nSide)
        self.pols=self.beamFeedAndDish.pols
        for mm in range(self.nf):
            dataFeed=self.beamFeedOnly.data[0,mm,:]
            dataDish=self.beamFeedAndDish.data[0,mm,:]
            #compute feed efficiency by integrating feed beam intercepted by dish
            selectionEdge=theta<self.thetaEdge
            self.daveFeedEfficiency[mm]=n.sum(dataFeed[selectionEdge])/n.sum(dataFeed)
            #compute taper by averaging gain on ring at dish edge
            phiRange=n.radians(n.arange(0,360))
            thetaRange=n.ones(len(phiRange))*self.thetaEdge
            self.daveTaper[mm]=1./n.mean(hp.get_interp_val(dataFeed,thetaRange,phiRange))
            #compute polarization mismatch by integrating phi=0 and phi=90
            thetaRange=n.radians(n.arange(0,180))
            phiRange=n.ones(len(thetaRange))*0.
            xArc=hp.get_interp_val(dataDish,thetaRange,phiRange)
            yArc=hp.get_interp_val(dataDish,thetaRange,phiRange+pi/2.)
            self.davePol[mm]=n.sum((n.sqrt(xArc)-n.sqrt(yArc))**2.)/n.sum((n.sqrt(xArc)+n.sqrt(yArc))**2.)
            #compute beam efficiency. First determine FWHM by taking average in 360 degree ring. 
            stddevs=[]
            nth=1000
            tha=n.arange(-nth/2,nth/2)*pi/nth
            for nn in range(0,180,10):
                cut=hpCut(n.radians(nn),nth,dataDish)
                try:
                    stddevs.append(op.curve_fit(lambda x,p: n.exp(-x**2./(2.*p**2.)),tha,cut,p0=[n.radians(15.)])[0])
                except Exception as e:
                    print e
                    print 'error'
                    p.plot(tha,cut)
                    p.yscale('log')
                    p.show()
            stddeves=n.array(stddevs)
            stddev=n.max(stddevs)
            fwhm=2.*n.sqrt(2.*n.log(2.))*stddev
            self.beamFWHM[mm]=fwhm
            selectionFWHM=theta<fwhm/2.
            self.daveBeamEfficiency[mm]=n.sum(dataDish[selectionFWHM])/n.sum(dataDish)
