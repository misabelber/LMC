import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from matplotlib import rc
import os
from matplotlib import ticker

formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((-2,2))

def plotCTAsens(filename):

    data = np.genfromtxt(filename,unpack=True)
    E1 = data[0]
    E2 = data[1]
    Sens = data[2]
    E = (data[0]+data[1])/2.
    
    plt.plot(E,Sens,"o--",label="CTA South Sensitivity 50h",linewidth=0.4)
        
def plotHESSlimit(filename):
    data = np.genfromtxt(filename,unpack=True)
    E = data[0]
    Flux = 1.60218*data[0]*data[1]/10000
    plt.plot(E,Flux,"x",label="H.E.S.S. upper limits on SN 1987A",linewidth=0.4)
    
def plotSpec(filename,srcname):
    data = np.genfromtxt(filename,unpack=True)
    E = data[0]
    Flux = data[1]
    plt.plot(E,Flux,"x",label="",linewidth=0.4)

if __name__ == '__main__':
    
    plotCTAsens("CTAsens50h.dat")
    plotHESSlimit("HESSUpperLimitsonSN1987A.dat")
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Energy (TeV)')
    plt.ylabel('$E^{2}$ $Flux(erg$ $cm^{-2} s^{-1})$')
    plt.legend()
    plt.show()
    
    
