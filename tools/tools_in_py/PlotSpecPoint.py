import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from matplotlib import rc
import os
from matplotlib import ticker

formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((-2,2))

def plot_Spectra(filename,case,marker):

    data = np.genfromtxt(filename,unpack=True)
    energy = data[0]
    flux = data[1]
    energy = energy[flux>1e-60]
    flux = flux[flux>1e-60]
    plt.plot(energy,flux,marker,label=case,linewidth=2)

    
# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    # Filepath to results
    path = "/afs/ciemat.es/user/b/bernardos/GitHub/LMC/spectra/sources/"

    sourcename = ["J0537-691",
                  "J0524.5-6937",          
                  "J0534.1-6732",                  
                  "J0525.2-6614",                                                              
                  "J0535.3-6559",
                  "J0454.6-6825",  
                  "J0537.0-7113",  
                  "J0535-691",  
                  "J0525-696"]
    
    for source in sourcename:
        filename = "spec_"+source+".dat"
        plot_Spectra(path+filename,source,"-")
        #plt.title("Spectrum of "+source)
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Energy (MeV)')
        plt.ylabel('Flux $phcm^{-2}s^{-1}MeV^{-1}$')
        #plt.title("$\mu\overline{\mu}$ annihilation channel")
        plt.legend()
    plt.show()
        
