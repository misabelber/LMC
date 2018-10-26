import numpy as np
import matplotlib.pyplot as plt
import os

energies = np.logspace(4,8,100)

#====================#
#J0537.9-6909/J0537-691
#====================#

Fermi = np.array([])
Hess = np.array([])
Broken = np.array([])

f = 5.637e-9*(1e6)**-1.502
h = 8.2e-2*(1e6)**-2.8
factor = f/h/2

for e in energies:
    f = 5.637e-9*e**-1.502
    Fermi = np.append(Fermi,f) 
    h = 8.2e-2*e**-2.8
    Hess = np.append(Hess,h)

    b = (5.637e-9*e**-1.502*(1+(e/1e6)**((-1.502+2.8)/0.2))**-0.2)/factor
    Broken = np.append(Broken,b)
plt.plot(energies[energies<=1e6],Fermi[energies<=1e6])
plt.plot(energies[energies>=1e6],Hess[energies>=1e6])
plt.plot(energies,Broken)
plt.xscale('log')
plt.yscale('log')
plt.show()

data = np.genfromtxt("/afs/ciemat.es/user/b/bernardos/GitHub/LMC/spectra/sources/spec_J0530.0-6900e.dat",unpack=True)
energy = data[0]
flux = data[1]

#====================#
#J0534.1-6732/J0536-675
#====================#

Fermi = np.array([])
Hess = np.array([])
Broken = np.array([])

#Difference at 1 TeV

f = 0.13*(1e6)**-2.8
h = 2e-4*(1e6)**-2.5

factor = f/h/2

for e in energies:
    f = 0.13*e**-2.8
    Fermi = np.append(Fermi,f) 
    h = 2e-4*e**-2.5
    Hess = np.append(Hess,h)

    b = (0.13/factor*1.5*e**-2.8*(1+(e/1e6)**((-2.8+2.5)/-0.1))**0.1)
    Broken = np.append(Broken,b)
plt.plot(energies[energies<=1e6],Fermi[energies<=1e6])
plt.plot(energies[energies>=1e6],Hess[energies>=1e6])
plt.plot(energies,Broken)
plt.xscale('log')
plt.yscale('log')
plt.show()

#====================#
#J0524.5-6937
#====================#

Fermi = np.array([])
CutOff = np.array([])

for e in energies:
    f = 0.626*1e-14*(e/7467.335)**-1.633
    Fermi = np.append(Fermi,f) 
    
    c = 0.626*1e-14*(e/7467.335)**-1.633*np.exp(-(e/1e6)**3)
    CutOff = np.append(CutOff,c)

plt.plot(energies,Fermi)
plt.plot(energies[CutOff>1e-30],CutOff[CutOff>1e-30])
plt.xscale('log')
plt.yscale('log')
plt.show()

#====================#
#J0524.5-6937
#====================#

Fermi = np.array([])
CutOff = np.array([])

for e in energies:
    f = 0.626*1e-14*(e/7467.335)**-1.633
    Fermi = np.append(Fermi,f) 
    
    c = 0.626*1e-14*(e/7467.335)**-1.633*np.exp(-(e/1e6)**-3)
    CutOff = np.append(CutOff,c)

plt.plot(energies,Fermi)
plt.plot(energies[CutOff>1e-30],CutOff[CutOff>1e-30])
plt.xscale('log')
plt.yscale('log')
plt.show()
