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
    plt.plot(dm_mass,limit,"x--",label=case,linewidth=0.4)

    
# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    # Filepath to results
    particle="W"
    path = "/afs/ciemat.es/user/b/bernardos/GitHub/LMC/results/"
    suf              ="_rebin_0.1x100_Pointin5deg"    
    filename         =["Limits"+particle+suf+".dat",
                       "Limits"+particle+suf+"_gamma0.5.dat",
                       "Limits"+particle+suf+"_gamma1.5.dat"]
    

    '''
    filename         =["Limits"+particle+suf+".dat",
                       "Limits"+particle+suf+"_DiffTRUNC500GeV.dat",
                       "Limits"+particle+suf+"_DiffTRUNC1TeV.dat"]
    '''
    case             =["NFW","$\gamma=0.5$","$\gamma=1.5$"]
    #case             =["Full spectra","Cut 500GeV","Cut 1TeV"]
    # Call plot function
    plot_Limits(path+filename[0],case[0],True)
    plot_Limits(path+filename[1],case[1],False)
    plot_Limits(path+filename[2],case[2],False)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Mass (TeV)')
    plt.ylabel('$<\sigma v> cm^{3}s^{-1}$')
    #plt.title("$\mu\overline{\mu}$ annihilation channel")
    plt.title("$W^{+}W^{-}$ annihilation channel")
    plt.legend()


    plt.show()
