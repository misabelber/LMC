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
    #particle=["W","b","Mu","Tau","e","g","Z","t"]
    particle = "W"
    path = "/afs/ciemat.es/user/b/bernardos/GitHub/LMC/results/"
    #suf              ="_KSP_100GeV-100TeVisomean"    
    suf = ["_KSP_100GeV-100TeV","_KSP_100GeV-100TeVgamma1.5","_KSP_100GeV-100TeVisomean","_KSP_100GeV-100TeVisomin"]
    '''
    filename         =["Limits"+particle[0]+suf+".dat",
                       "Limits"+particle[1]+suf+".dat",
                       "Limits"+particle[2]+suf+".dat",
                       "Limits"+particle[3]+suf+".dat",
                       "Limits"+particle[4]+suf+".dat",
                       "Limits"+particle[5]+suf+".dat",
                       "Limits"+particle[6]+suf+".dat",
                       "Limits"+particle[7]+suf+".dat"]
    '''    

    filename         =["Limits"+particle+suf[0]+".dat",
                       "Limits"+particle+suf[1]+".dat",
                       "Limits"+particle+suf[2]+".dat",
                       "Limits"+particle+suf[3]+".dat",];

    
#    case             =["$W^{+}W^{-}$","$b\overline{b}$",r"$\mu^{+}\mu^{-}$",
#                       r"$\tau^{+}\tau^{-}$",r"$e^{+}e^{-}$","gg","$Z^{+}Z^{-}$","$t\overline{t}$"]
    case             =["NFW","gamma=1.5","Isothermal(mean)","Isothermal(min)"]
    #case             =["Full spectra","Cut 500GeV","Cut 1TeV"]
    # Call plot function
    plot_Limits(path+filename[0],case[0],True)
    plot_Limits(path+filename[1],case[1],False)
    plot_Limits(path+filename[2],case[2],False)
    plot_Limits(path+filename[3],case[3],False)
    #plot_Limits(path+filename[4],case[4],False)
    #plot_Limits(path+filename[5],case[5],False)
    #plot_Limits(path+filename[6],case[6],False)
    #plot_Limits(path+filename[7],case[7],False)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Mass (TeV)')
    plt.ylabel('$<\sigma v> cm^{3}s^{-1}$')
    #plt.title("$\mu\overline{\mu}$ annihilation channel")
    #plt.title("Mean sensitivity curves for Isothermal(mean) profile")
    plt.title("Mean sensitivity curves for $W^{+}W^{-}$ annihilation channel")
    plt.legend()


    plt.show()
