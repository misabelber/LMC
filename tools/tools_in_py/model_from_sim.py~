from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
from astropy.wcs import WCS
import numpy as np
from matplotlib.colors import LogNorm
import scipy
import argparse

parser = argparse.ArgumentParser(description = "Name of component")
parser.add_argument('--component', '-c', type=str,
                    dest='component',
                    help='name of model component')

parser.add_argument('--nobs', '-n', type=int,
                    dest='nobs',
                    help='Number of observations')

args = parser.parse_args()

"""
def rebin(a, *args):
    shape = a.shape
    lenShape = len(shape)
    factor = asarray(shape)/asarray(args)
    evList = ['a.reshape('] + \
             ['args[%d],factor[%d],'%(i,i) for i in range(lenShape)] + \
             [')'] + ['.mean(%d)'%(i+1) for i in range(lenShape)]
    print ''.join(evList)
    return eval(''.join(evList))
"""
# ======================== #                                                                           
# Main routine entry point #
# ======================== #

if __name__ == '__main__':

    component=args.component
    nobs=args.nobs
    
    suf = ""
    PATH = "/scratch/bernardos/LMC/Obs_"+component+suf+"/" #Para extended srcs, CR(Irf), Diffuse
    #PATH = "/scratch/bernardos/LMC/Obs_PS/" #Para Point Sources
    #PATH = "/scratch/bernardos/LMC/Obs_DM/" #Para DM
    filecntcube=PATH+"cntcube_LMC_"+component+"_functions_KSP_100GeV-100TeV"
    #filecntcube=PATH+"cntcube_dm_LMC_"+component+"_KSP_100GeV-100TeV" #Para DM
    hdul = fits.open(filecntcube+"000.fits")
    hdu = hdul[0]
    mean = np.array(hdu.data)
    mean2 = np.array(mean**2)

    for i in range(1,nobs):
        obstring = "%03d" % i
        newcntcube = np.array(fits.open(filecntcube+obstring+'.fits')[0].data)
        mean = np.array(mean+newcntcube)
        mean2 = np.array(mean2+newcntcube**2)

   
    mean = mean/nobs
    mean2 = mean2/nobs
    sigma2 = mean2-mean*mean
    realsigma2 = sigma2*nobs/(nobs-1)
    
    sigma = np.sqrt(sigma2)
    realsigma=np.sqrt(realsigma2)
    
    hdumean = fits.PrimaryHDU(mean)
    hdusigma = fits.PrimaryHDU(realsigma)
    header = hdu.header
    hdumean.header = header
    ebins = hdul[2]
    newhdul = fits.HDUList([hdumean,ebins])

    newhdul.writeto(PATH+"Model_"+component+"_functions_KSP_100GeV-100TeV.fits",overwrite=True)

    plt.subplot(131)
    cmap = plt.get_cmap('inferno')
    cmap.set_bad((0,0,0))
    plt.imshow(mean.sum(0), cmap=cmap, origin='lower')
    plt.grid(color='white',ls='solid')
    plt.xlabel("Right Ascension J2000")
    plt.ylabel("Declination J2000")
    plt.colorbar(label="# of Counts (Log scale)")
    
    plt.subplot(132)
    cmap = plt.get_cmap('inferno')
    cmap.set_bad((0,0,0))
    plt.imshow(np.sqrt(sigma2.sum(0)), cmap=cmap, origin='lower')
    plt.grid(color='white',ls='solid')
    plt.xlabel("Right Ascension J2000")
    plt.ylabel("Declination J2000")
    plt.colorbar(label="# of Counts (Log scale)")
    
    plt.subplot(133)
    cmap = plt.get_cmap('inferno')
    cmap.set_bad((0,0,0))
    plt.imshow(np.sqrt(realsigma2.sum(0)), cmap=cmap, origin='lower')
    plt.grid(color='white',ls='solid')
    plt.xlabel("Right Ascension J2000")
    plt.ylabel("Declination J2000")
    plt.colorbar(label="# of Counts (Log scale)")

    plt.show()

    """
    for bin in range(0,20):
        
        plt.subplot(131)
        cmap = plt.get_cmap('inferno')
        cmap.set_bad((0,0,0))
        plt.imshow(mean[bin], cmap=cmap, origin='lower')
        plt.grid(color='white',ls='solid')
        plt.xlabel("Right Ascension J2000")
        plt.ylabel("Declination J2000")
        plt.colorbar(label="# of Counts (Log scale)")
        
        plt.subplot(132)
        cmap = plt.get_cmap('inferno')
        cmap.set_bad((0,0,0))
        plt.imshow(np.sqrt(sigma2[bin]), cmap=cmap, origin='lower')
        plt.grid(color='white',ls='solid')
        plt.xlabel("Right Ascension J2000")
        plt.ylabel("Declination J2000")
        plt.colorbar(label="# of Counts (Log scale)")
        
        plt.subplot(133)
        cmap = plt.get_cmap('inferno')
        cmap.set_bad((0,0,0))
        plt.imshow(np.sqrt(realsigma2[bin]), cmap=cmap, origin='lower')
        plt.grid(color='white',ls='solid')
        plt.xlabel("Right Ascension J2000")
        plt.ylabel("Declination J2000")
        plt.colorbar(label="# of Counts (Log scale)")
        
        plt.show()
    """

