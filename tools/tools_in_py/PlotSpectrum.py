import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from matplotlib import rc
import os
from matplotlib import ticker

formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((-2,2))

def plot_Spectra(filename,case,marker,color):

    data = np.genfromtxt(filename,unpack=True)
    energy = data[0]
    flux = data[1]
    plt.plot(energy,flux,marker,label=case,linewidth=2,color=color)

    
# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    # Filepath to results
    path = "/afs/ciemat.es/user/b/bernardos/GitHub/LMC/spectra/sources/"
    sourcename       = ["J0509.9-6418","J0454.6-6825","J0509.9-6418","J0537-691"] 
    suf              = ["","_3FGL","_3FHL","_HESS"]
    filename         =["spec_"+sourcename[0]+suf[0]+".dat",
                       "spec_"+sourcename[1]+suf[1]+".dat",
                       "spec_"+sourcename[2]+suf[2]+".dat",
                       "spec_"+sourcename[3]+suf[3]+".dat"]
    case             =["CTA","3FGL","3FHL","H.E.S.S."]
    marker           =["--","-"]
    color            =["red","white","green"]

    num = 2
    # Call plot function
    plot_Spectra(path+filename[0],case[0],marker[0],color[1])
    plot_Spectra(path+filename[num],case[num],marker[1],color[2])
    #plot_Spectra(path+filename,case[2])
    #plot_Spectra(path+filename,case[3])
    plt.title("Spectrum of "+case[num]+" "+sourcename[num])
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Energy (MeV)')
    plt.ylabel('Flux $phcm^{-2}s^{-1}MeV^{-1}$')
    #plt.title("$\mu\overline{\mu}$ annihilation channel")
    plt.legend()
    plt.show()
    
    # Call plot function
    plot_Spectra(path+filename[0],case[0],marker[0],color[0])
    plot_Spectra(path+filename[num],case[num],marker[1],color[2])
    
    #plot_Spectra(path+filename,case[2])
    #plot_Spectra(path+filename,case[3])
    plt.title("Spectrum of "+case[num]+" "+sourcename[num])
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Energy (MeV)')
    plt.ylabel('Flux $phcm^{-2}s^{-1}MeV^{-1}$')
    #plt.title("$\mu\overline{\mu}$ annihilation channel")
    plt.legend()
    plt.show()

