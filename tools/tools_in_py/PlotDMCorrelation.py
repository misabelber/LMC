import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from matplotlib import rc
import os
from matplotlib import ticker

formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((-2,2))


def plot_Limits(filename,components,component):
        
    data = np.genfromtxt(filename,unpack=True)
    dm_mass = data[0]
    cfactor = data[component+1]
    print(filename)
    print(data[0],data[component])
    plt.plot(dm_mass,cfactor,".--",label=components[component],linewidth=0.5)
    
    
    
# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    # Filepath to results
    particle="W"
    path = "/afs/ciemat.es/user/b/bernardos/GitHub/LMC/results/"
    filename         ="Cfactors_b_modelHM.dat"
    
    components = ["CRBkg",
                  "Leptonic",                                                             
                  "Hadronic",                                                             
                  "3FHL_J0500.9-6945e",                                                   
                  "3FHL_J0530.0-6900e",                                                   
                  "3FHL_J0531.8-6639e",
                  "J0537-691",                                                                
                  "J0524.5-6937",                                                             
                  "J0534.1-6732",                                                             
                  "J0525.2-6614",                                                             
                  "J0535.3-6559",                                                             
                  "J0454.6-6825",                                                             
                  "J0537.0-7113",                                                             
                  "J0535-691",                                                                
                  "J0525-696"]
    # Call plot function
    plot_Limits(path+filename,components,0)
    plot_Limits(path+filename,components,1)
    plot_Limits(path+filename,components,2)
    plot_Limits(path+filename,components,3)
    plot_Limits(path+filename,components,4)
    plot_Limits(path+filename,components,5)
    plt.xscale('log')
    plt.xlabel('Mass (TeV)')
    plt.ylabel('Correlation Factor')
    #plt.title("$W^{+}W^{-}$ annihilation channel")
    plt.title("$b\overline{b}$ annihilation channel")
    plt.legend()
    plt.show()
    plot_Limits(path+filename,components,6)
    plot_Limits(path+filename,components,7)
    plot_Limits(path+filename,components,8)
    plot_Limits(path+filename,components,9)
    plot_Limits(path+filename,components,10)
    plot_Limits(path+filename,components,11)
    plot_Limits(path+filename,components,12)
    plot_Limits(path+filename,components,13)
    plot_Limits(path+filename,components,14)
    plt.xscale('log')
    plt.xlabel('Mass (TeV)')
    plt.ylabel('Correlation Factor')
    #plt.title("$W^{+}W^{-}$ annihilation channel")
    plt.title("$b\overline{b}$ annihilation channel")
    plt.legend()
    plt.show()
