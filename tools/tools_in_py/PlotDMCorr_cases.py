import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from matplotlib import rc
import os
from matplotlib import ticker

formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((-2,2))


def plot_Limits(filename,components,component,color,marker):
        
    data = np.genfromtxt(filename,unpack=True)
    dm_mass = data[0]
    cfactor = data[component+1]
    print(filename)
    print(data[0],data[component])
    plt.plot(dm_mass,cfactor,marker,label=components[component],linewidth=0.5,color=color)
    
    
    
# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    # Filepath to results
    particle="W"
    path = "/afs/ciemat.es/user/b/bernardos/GitHub/LMC/results/"
    filenames         =["DMCorrelations_"+particle+".dat",
                        "DMCorrelations_"+particle+"_gamma0.5.dat",
                        "DMCorrelations_"+particle+"_gamma1.5.dat"]
    
    components        =[["Irf",
                         "Leptonic",                                                             
                         "Hadronic",
                         "3FHL_J0500.9-6945e",
                         "3FHL_J0530.0-6900e",
                         "3FHL_J0531.8-6639e"],
                        ["Irf $\gamma$=0.5",
                         "Leptonic $\gamma$=0.5",
                         "Hadronic $\gamma$=0.5",
                         "3FHL_J0500.9-6945e $\gamma$=0.5",
                         "3FHL_J0530.0-6900e $\gamma$=0.5",
                         "3FHL_J0531.8-6639e $\gamma$=0.5"],
                        ["Irf $\gamma$=1.5",
                         "Leptonic $\gamma$=1.5",
                         "Hadronic $\gamma$=1.5",
                         "3FHL_J0500.9-6945e $\gamma$=1.5",
                         "3FHL_J0530.0-6900e $\gamma$=1.5",
                         "3FHL_J0531.8-6639e $\gamma$=1.5"]]
    colors            =["red",
                        "blue",
                        "green",
                        "orange",
                        "salmon",
                        "purple"]
    markers           =[".-",
                        ".--",
                        ".-."]
    for i,filename in enumerate(filenames):
        # Call plot function
        #plot_Limits(path+filename,components[i],0,colors[0],markers[i])
        #plot_Limits(path+filename,components[i],1,colors[1],markers[i])
        #plot_Limits(path+filename,components[i],2,colors[2],markers[i])
        #plot_Limits(path+filename,components[i],3,colors[3],markers[i])
        #plot_Limits(path+filename,components[i],4,colors[4],markers[i])
        plot_Limits(path+filename,components[i],5,colors[5],markers[i])
    plt.xscale('log')
    plt.xlabel('Mass (TeV)')
    plt.ylabel('Correlation Factor')
    plt.title("$W^{+}W^{-}$ annihilation channel")
    handles,labels = plt.gca().get_legend_handles_labels()
    #handles = [handles[0],handles[3],handles[6],handles[1],handles[4],handles[7],handles[2],handles[5],handles[8]]
    #labels = [components[0][0],components[1][0],components[2][0],components[0][1],components[1][1],components[2][1],components[0][2],components[1][2],components[2][2]]
    #plt.legend(handles,labels)
    plt.legend()
    plt.show()
    
