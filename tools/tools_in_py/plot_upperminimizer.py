# =========================================================================== 
# Plot outputs from Upper_Minimizer()
# 
# =========================================================================== 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from matplotlib import rc
import os
rc('font', family='times new roman', size=16)
rc('text', usetex=True)


def plot_results_upperminimizer(filepath, bar_comp):
    """
    Plot results from Upper_Minimizer

    Input
    -----
    filepath : string
        path to results 
    bar_comp : string
        Baryonic component
    Output
    ------
        Plot with 2D 2*(maxLogL - LogL) with random paths obtained in
        trial and iter for loops. Highlighting the paths accepted.
    """
    # Load difference = 2*(maxLogL - LogL)
    data  = np.genfromtxt(filepath + "difference.dat", unpack=True)
    x_dif = data[0]
    y_dif = data[1]
    dif   = data[2]
    # Grid difference
    xgrid        = np.linspace(x_dif.min(), x_dif.max(), 100)
    ygrid        = np.linspace(y_dif.min(), y_dif.max(), 100)
    xgrid, ygrid = np.meshgrid(xgrid, ygrid)
    zgrid        = griddata(x_dif, y_dif, dif, xgrid, ygrid, interp="linear")
    data   = np.genfromtxt(filepath + "normalization.dat", unpack=True)
    x_norm = data[0]
    y_norm = data[1]
    trial  = data[2]
    del data
    trial_acc = np.genfromtxt(filepath + "trial_acc.dat", unpack=True)
    # Save each random path within dictionary
    paths ={}
    for i in range(len(trial)):
        paths[i] = [[], []]
        for j in range(len(trial)):
            if abs(i-trial[j]) < 1e-1:
                paths[i][0].append(x_norm[j])
                paths[i][1].append(y_norm[j])
    # Plot
    k       = 0
    fig, ax = plt.subplots()
    levels  = np.arange(0, 3., 0.4)
    cax     = ax.contourf(xgrid, ygrid, zgrid, levels=levels)
    fig.colorbar(cax)
    for i in range(100):
        if abs(i-trial_acc[k])<1e-1 and k<len(trial_acc)-1:
            ax.plot(paths[i][0], paths[i][1], lw=3., color="k")
            ax.plot(paths[i][0][-1], paths[i][1][-1], "ko")
            k = k+1
        elif abs(i-trial_acc[k])<1e-1 and k==len(trial_acc)-1:
            ax.plot(paths[i][0], paths[i][1], lw=3., color="red")
            ax.plot(paths[i][0][-1], paths[i][1][-1], "ro")
        else:
            ax.plot(paths[i][0], paths[i][1], color="grey", alpha=0.7)
    ax.set_xlabel("normalization DM")
    ax.set_ylabel("normalization" + bar_comp)
    plt.show()

# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    # Filepath to results
    filepath_to_here = os.getcwd()
    filepath         = filepath_to_here + "/../../results/tests/"
    bar_comp         = "Leptonic"
    # Call plotbands function
    plot_results_upperminimizer(filepath, bar_comp)
