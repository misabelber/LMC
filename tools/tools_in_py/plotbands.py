# =========================================================================== 
# Plot DM constraints ("upper limits")
# 
# =========================================================================== 
import matplotlib.pyplot as plt
import numpy as np
import os


def plotbands(filepath, masses, particle,suf, bestCase=True, tol=0.01):
    """
    Plot DM constraints

    Input
    -----
    filepath : string
        path to results (constraints)
    masses : list
        DM masses
    partile : string
        final annihilation particle
    suf  : string
       suffix on the filename
    bestCase : bool
        True in case limits computed imposing the normalization not to change 
        much
    tol : double
        In case bestCase=True, how much you let the normalized factor to vary

    Output
    ------
        Plot with Mean Expected constraint and 84% and 95% containment
    """
    # Number of masses
    nmasses    = len(masses)
    # Mean
    meanlimit  = np.ones(nmasses)
    # 95% containment bands
    sstdevs    = np.ones(nmasses)
    # 84% containment bands
    sstdevs2   = np.ones(nmasses)
    # Thermal averaged annihilation cross-section
    thermalX   = np.ones(nmasses)
    
    for i in range(nmasses):
        dm_mass = masses[i];
        # Build wisely the name of the file
        if (dm_mass < 1.0):
            masstr = str(int(dm_mass*1000))+"GeV";      
        else:
            masstr = str(int(dm_mass))+"TeV"; 

        if bestCase==True:
            filename = filepath+"Limits_"+particle+masstr+suf+"_bestCase_tol0.01.dat"
        else:
            filename = filepath+"Bands_"+particle+masstr+suf+".dat"

        data   = np.genfromtxt(filename, unpack=True) 
        limits = data*3e-26
        mean   = np.mean(limits)
        std    = np.std(limits)
        # Plot mean and standard deviation
        print(mean/3e-26, "  ", std)

        meanlimit[i] = mean;
        # 95% containment bands
        sstdevs[i]   = 1.96*std
        # 84% containment bands
        sstdevs2[i]  = 1.4*std
        # Thermal averaged cross section
        thermalX[i]  = 3e-26
        
    #fig, ax = plt.subplots()
    
    """
    # Plot mean constraint line
    ax.plot(masses, meanlimit, lw=2., color="k", label="Mean Expected")
    # Plot 95% containment band
    ax.fill_between(masses, meanlimit-sstdevs, meanlimit+sstdevs, 
                    color="#FFFF33", alpha=0.7, label="95$\%$ containment")
    # Plot 84% containment band
    ax.fill_between(masses, meanlimit-sstdevs2, meanlimit+sstdevs2, 
                    color="#00CC33", alpha=0.7, label="84$\%$ containment")
    # Plot thermal average annihilation cross section
    ax.axhline(3e-26, color="red", ls="--", 
               label="Thermal $\\langle\\sigma v \\rangle$")
    
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_ylabel("$\\langle\\sigma v \\rangle$ $\\rm [cm^3s^{-1}]$")
    ax.set_xlabel("DM mass [TeV]")
    ax.grid(True, which="both", ls="--", alpha=0.5)
    ax.legend(frameon=False)
    """
    # Plot mean constraint line
    plt.plot(masses, meanlimit, lw=2., color="k", label="Mean Expected")
    # Plot 95% containment band
    plt.fill_between(masses, meanlimit-sstdevs, meanlimit+sstdevs, 
                    color="#FFFF33", alpha=0.7, label="95$\%$ containment")
    # Plot 84% containment band
    plt.fill_between(masses, meanlimit-sstdevs2, meanlimit+sstdevs2, 
                    color="#00CC33", alpha=0.7, label="84$\%$ containment")
    # Plot thermal average annihilation cross section
    plt.axhline(3e-26, color="red", ls="--", 
               label="Thermal $\\langle\\sigma v \\rangle$")
        
    #plt.show()

# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    # DM masses
    #masses           = [0.200,0.300,0.400,0.500,0.600,0.800,1,4,5,8,10,40,50,80,100]
    #masses           = [0.2,0.5,1,2,5,10,20,50,100]
    masses           = [0.2,0.5,1,5,10,50,100]
    # Annihilated final particle
    particle         = "Tau";
    # Filepath to results
    filepath_to_here = os.getcwd()
    filepath         = filepath_to_here + "/../../results/"
    suf = "_modelHM"
    # Call plotbands function
    plotbands(filepath, masses, particle, suf,False)
    #plt.title("$W^{+}W^{-}$ annihilation channel")
    #plt.title("$b\overline{b}$ annihilation channel")
    plt.title(r"$\tau\overline{\tau}$ annihilation channel")
    
    """
    masses           = [0.2,0.5,1,5,10,50,100]
    # Annihilated final particle
    particle         = "b";
    # Filepath to results
    filepath_to_here = os.getcwd()
    filepath         = filepath_to_here + "/../../results/"
    suf = "_modelHM"
    # Call plotbands function
    plt.legend(frameon=False)
    plotbands(filepath, masses, particle, suf,False)
    """ 
    plt.yscale("log")
    plt.xscale("log")
    plt.ylabel("$\\langle\\sigma v \\rangle$ $\\rm [cm^3s^{-1}]$")
    plt.xlabel("DM mass [TeV]")
    plt.legend(frameon=False)
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.show()


