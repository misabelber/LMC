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
    Sens = np.log10(data[2]*6.2415E11)
    E = np.log10((data[0]+data[1])/2.*1E12)
    
    plt.plot(E,Sens,"--",color="red",label="CTA South Sensitivity 50h",linewidth=2)
        
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
    #plotHESSlimit("HESSUpperLimitsonSN1987A.dat")
    #plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel('log(Energy) (eV)')
    plt.ylabel('$log(E)$ $Flux(eV$ $cm^{-2} s^{-1})$')
    yticks = [-1.5,-1,-0.5,0]
    plt.yticks(yticks)
    plt.xlim(xmax=14.5,xmin=10)
    xticks=[10,11,12,13,14]
    plt.xticks(xticks)
    plt.ylim(ymax=0,ymin=-1.5)
    #plt.legend()
    plt.show()
    
    
