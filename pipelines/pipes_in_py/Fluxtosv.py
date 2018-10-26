from scipy.interpolate import interp1d
import numpy as np
#from io import StringIO
from StringIO import StringIO
from matplotlib import pyplot as plt
import io
import math

filename = '../../tools/tools_in_py/AtProduction_gammas.dat'

with open(filename) as f:
    lines = (line for line in f if not line.startswith('#'))
    data = np.genfromtxt (lines, names = True ,dtype = None)

GeVtoMeV = 1000
TeVtoGeV = 1000

massvals = data["mDM"]

#Define DM model:

jfactor = 1e20 #GeV^2/cm^5
finalstate = "b"
sv=3.0*1e-26 #cm^3/s
jfactor = jfactor*GeVtoMeV*GeVtoMeV #MeV^2/cm^5
#Get flux from file

fileflux = "Ulimitsresults.txt"
datalimits = np.genfromtxt(fileflux,unpack=True)
masses = datalimits[0].tolist()
fluxes = datalimits[2]#MeV/cm^2/s
masses = [0.1.0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100]
#masses = [0.1]

for mass,fluxulim in zip(masses,fluxes):

    mass = mass*TeVtoGeV
    #Get dN/dE from Cirelli et. al
    index = np.where(np.abs( (massvals - mass) / mass) < 1.e-3)
    xvals = 10**(data["Log10x"][index])
    
    flux = data[finalstate][index]/(np.log(10)*xvals)
    loadspec = interp1d(xvals,flux)

    def dNdx(x):
        fluxval = loadspec(x)
        if (x>1 or fluxval<0):
            return 0
        else:
            return fluxval
            
    #all in MeV
    dNdE = dNdx(0.9)/mass/GeVtoMeV
    
    mass = mass*GeVtoMeV
    print mass,fluxulim,jfactor,dNdE
    svulim = 8*np.pi*mass*mass*fluxulim/jfactor/dNdE
    print svulim
