import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from matplotlib import rc
import os
from matplotlib import ticker

formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((-2,2))


def plot_Limits(filename,case,thermal):

    ThermalX = 3e-26
        
    data = np.genfromtxt(filename,unpack=True)
    dm_mass = data[0]
    limit = data[1]*ThermalX
    print(filename)
    print(data[0],data[1])
    tx = np.full((data[0].size),ThermalX)
    if thermal==True:
        plt.plot(dm_mass,tx,"--",color='red',label="Thermal$<\sigma v>$" )
    plt.plot(dm_mass,limit,"o-",label=case)
    
    
    
# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    # Filepath to results
    particle="W"
    path = "/afs/ciemat.es/user/b/bernardos/GitHub/LMC/results/"
    filename         =["Limits"+particle+"_jfactorNFW",
                       "Limits"+particle+"_jfactorgamma0.5",
                       "Limits"+particle+"_jfactorgamma1.5"]
    
    case             =["NFW","$\Gamma=0.5$","$\Gamma=1.5$"]
    # Call plot function
    plot_Limits(path+filename[0],case[0],True)
    plot_Limits(path+filename[1],case[1],False)
    plot_Limits(path+filename[2],case[2],False)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Mass (TeV)')
    plt.ylabel('$<\sigma v> cm^{3}s^{-1}$')
    plt.title("$\mu\overline{\mu}$ annihilation channel")
    
    plt.legend()


    plt.show()
