from scipy.interpolate import interp1d
import numpy as np
#from io importpython StringIO
#from io import StringIO
from matplotlib import pyplot as plt
import io
import math
import sys

filename = 'AtProductionNoEW_gammas.dat'
finalstate = sys.argv[1]  # choose particle W, b, etc
suf = sys.argv[2]

with open(filename) as f:
    lines = (line for line in f if not line.startswith('#'))
    data = np.genfromtxt (lines, names = True ,dtype = None)

massvals = data["mDM"]
TeVtoGeV=1000
GeVtoMeV=1000
TeVtoMeV=1e6

fileresults = "/afs/ciemat.es/user/b/bernardos/GitHub/LMC/pipelines/pipes_in_py/Ulimitsresults_empty_"+finalstate+suf+".txt"
with open(fileresults) as fr:
    results = np.genfromtxt(fr,unpack=True)

masses = np.array([])
limits = np.array([])

for i in range(results.shape[1]):

    mass = results[0][i]*TeVtoGeV
    eref = results[1][i]*TeVtoGeV
    dFdE = results[2][i]
    integral_flux = results[3][i]
    
    skip = np.array([])
    
    if mass in skip:
        continue

    emin = results[5][i]*TeVtoGeV                                                                 
    emax = results[6][i]*TeVtoGeV
    
    index = np.where(np.abs( (massvals - mass) / mass) < 1.e-3)
    xvals = 10**(data["Log10x"][index])

    flux = data[finalstate][index]/(np.log(10)*xvals)
    flux = flux[:180]

    try:
        loadspec = interp1d(xvals,flux)
    
    except:
        continue

    def dNdx(x):
        fluxval = loadspec(x)
        if (x >1 or fluxval<0):
            return 0
        else:
            return fluxval
    
    jfactor=1e26  # MeV ^ 2 /cm ^ 5
    sv=3.0*1e-26 #cm3s-1
    
    evals = np.linspace(emin,emax,1000)
    dNdEvals = np.array([])
    '''
    for e in evals:
        dNdE = dNdx(e/mass)/(e*GeVtoMeV)
        dNdEvals = np.append(dNdEvals,dNdE)
    '''
    xrange = xvals[xvals*mass >= emin]
    xmin = xrange[0]
    xmax = xrange[xrange.shape[0]-1]
    new_x = np.linspace(xmin,xmax,1000)

    for x in new_x:
        dNdE = dNdx(x)/mass/GeVtoMeV
        dNdEvals = np.append(dNdEvals,dNdE)                                                        
        
       
    dNdE_integral = np.trapz(dNdEvals,x=new_x*mass*GeVtoMeV)
    dNdE = dNdx(eref/mass)/(eref*GeVtoMeV)        
    #limit = (dFdE*8*math.pi*(mass*GeVtoMeV)**2)/(jfactor*sv*dNdE)
    limit = (integral_flux*8*math.pi*(mass*GeVtoMeV)**2)/(jfactor*sv*dNdE_integral)

    masses = np.append(masses,mass)
    limits = np.append(limits,limit*sv)
    print(mass/TeVtoGeV,limit)

plt.plot(masses,limits,
         label = 'Limit on $\\langle\\sigma v \\rangle$' )
plt.axhline(y=sv, color='r', linestyle='--',
            label="Thermal $\\langle\\sigma v \\rangle$")
plt.yscale('log')
plt.xscale('log')
plt.xlabel("DM mass [TeV]")
plt.ylabel("$\\langle\\sigma v \\rangle$ $\\rm [cm^3s^{-1}]$")
plt.legend()
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.title("$b\overline{b}$ annihilation channel")
plt.show()
