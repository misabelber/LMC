from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
from astropy.wcs import WCS
import numpy as np
from matplotlib.colors import LogNorm
import scipy



def rebin(a, *args):
    shape = a.shape
    lenShape = len(shape)
    factor = np.asarray(shape)/np.asarray(args)
    evList = ['a.reshape('] + \
             ['args[%d],factor[%d],'%(i,i) for i in range(lenShape)] + \
             [')'] + ['.mean(%d)'%(i+1) for i in range(lenShape)]
#    print ''.join(evList)
#    return eval(''.join(evList))


# ======================== #                                                                           
# Main routine entry point #
# ======================== #
                                                                             
if __name__ == '__main__':

    #Modcube
    PATH = "/scratch/bernardos/LMC/Obs_Irf+CR+DiffuseSources+PS/"
    PATH_func = "/scratch/bernardos/LMC/Obs_Irf+CR+DiffuseSources+PS_functions/"
    #filemodcube = PATH+"modcube_LMC_Irf+CR+DiffuseSources+PS_KSP_100GeV-100TeV.fits"
    filemodcube = PATH_func+"modcube_LMC_Irf+CR+DiffuseSources+PS_functions_KSP_100GeV-100TeV000.fits"
    #filemodcube = "/scratch/bernardos/Tests/modcube.fits"
    #filecntcube = PATH+"cntcube_LMC_Irf+CR+DiffuseSources+PS_KSP_100GeV-100TeV"
    filecntcube = PATH_func+"cntcube_LMC_Irf+CR+DiffuseSources+PS_functions_KSP_100GeV-100TeV"
    #filecntcube = "/scratch/bernardos/Tests/cntcube"
    hdulmodcube = fits.open(filemodcube)[0]

    modcube = np.array(hdulmodcube.data)
    cntcube = np.array(modcube)
    cntcube.fill(0)
    
    for i in range(0,20):
        obstring = "%03d" % i
        newcntcube = np.array(fits.open(filecntcube+obstring+'.fits')[0].data)
        cntcube = np.array(cntcube+newcntcube)
 
    cntcube = cntcube/20
            
    #residuals = np.sqrt(((cntcube-modcube)/np.sqrt(modcube+0.5))**2)
    residuals = (cntcube-modcube)/np.sqrt(modcube/20+0.5)
    
    
    
#    residuals = rebin(residuals,20,25,25)

    for bin in range(0,20):
        
        plt.subplot(131)
        cmap = plt.get_cmap('inferno')
        cmap.set_bad((0,0,0))
        plt.imshow(modcube[bin], cmap=cmap, origin='lower')
        plt.grid(color='white',ls='solid')
        plt.xlabel("Right Ascension J2000")
        plt.ylabel("Declination J2000")
        plt.colorbar(label="# of Counts (Log scale)")
        
        plt.subplot(132)
        cmap = plt.get_cmap('inferno')
        cmap.set_bad((0,0,0))
        plt.imshow(cntcube[bin], cmap=cmap, origin='lower')
        plt.grid(color='white',ls='solid')
        plt.xlabel("Right Ascension J2000")
        plt.ylabel("Declination J2000")
        plt.colorbar(label="# of Counts (Log scale)")
        
        plt.subplot(133)
        cmap = plt.get_cmap('inferno')
        cmap.set_bad((0,0,0))
        plt.imshow(residuals[bin], cmap=cmap, origin='lower')
        plt.grid(color='white',ls='solid')
        plt.xlabel("Right Ascension J2000")
        plt.ylabel("Declination J2000")
        plt.colorbar(label="# of Counts (Log scale)")
        
        plt.show()

