import numpy as np
import matplotlib.pyplot as plt

energies  = np.logspace(4.477121254719663,8.0,num=1000)
values = 5.7e-16*pow((energies/0.3e6),-2.48)
bins = np.linspace(0,950,20)
integral=[]
ebins =[]
error=[]
limits = np.loadtxt("/afs/ciemat.es/user/b/bernardos/GitHub/LMC/pipelines/pipes_in_C/fluxlimits_base.dat")

for fbin in bins:
    bin =  int(fbin)
    result = np.trapz(values[bin:bin+49],energies[bin:bin+49])
    ebin = (energies[bin]+energies[bin+49])/2
    error.append(energies[bin+49]-ebin)
    ebins.append(ebin)
    integral.append(result/ebin)
    
print integral

plt.errorbar(ebins,integral*limits*ebins*ebins,xerr=error,fmt="o",markersize=4,color='black')
plt.xscale('log')                                                                              
plt.yscale('log')
plt.show()
