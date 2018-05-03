
#Mab: This program substracts the Dark Matter spectrum for certain final state (particle) and DM masses from the Cirelli et al. file. 

#!/usr/bin/python

# Import required modules
#from astropy.io import ascii
from scipy.interpolate import interp1d
import numpy as np
#from io import StringIO
from StringIO import StringIO
from matplotlib import pyplot as plt
import io
import math

#Define output path to store DM spectra

OUT_PATH = "../../spectra/DM/"

filename = 'AtProduction_gammas.dat'

finalstate = "W"  # choose particle W, b, etc

with open(filename) as f:
    lines = (line for line in f if not line.startswith('#'))
    #data = np.genfromtxt (StringIO(lines), names = True ,dtype = None)
    #data = np.genfromtxt (BytesIO(lines), names = True ,dtype = None)
    data = np.genfromtxt (lines, names = True ,dtype = None)
#masses = [0.300,0.400,0.600,0.700,0.800,0.900,2,3,4,5,6,7,8,9,10]
masses = [0.100,0.200,0.300,0.400,0.500,0.600,0.800,1,4,5,8,10,40,50,80,100]     
#masses = [20,30,40,50,60,70,80,90,100] #Mab: Masses of dark matter particle in TeV

massvals = data["mDM"]
#units in GeV
#mass = 100000  # mass of dark matter eq 1 TeV
#units in MeV for ctools
GeVtoMeV=1000

for mass in masses:
    mass = mass*1000 #Mab: Convert mass from TeV to GeV
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
    #Mab: I don't want to plot anything here, so I commented it out.    
    #plt.plot(mass*GeVtoMeV*xvals,[np.log10(dNdx(x)/mass/GeVtoMeV) for x in xvals],color='dimgray')
    #plt.title('Ex. 3: Gamma Direct Spectrum into '+finalstate+ ' for $m_{\chi}=$'+massname,fontsize=14)
    #plt.xscale('log')
    #plt.xlabel('$E (MeV)$', fontsize=18)
    #plt.ylabel('$dN / dE (MeV^{-1})$', fontsize=18)
    #plt.ylim([-8,-1])
    #plt.xlim([10**2,10**8])
    #plt.legend(fontsize=12,loc=2)
    #plt.show()

    #Mab: Build wisely the name of the file so it can be easily used by other programs
    #____________________________________________________________________________
    if mass<1000:
        massname=str(int(mass))+"GeV"
    else:
        massname=str(int(mass)/1000)+"TeV"
    filename = OUT_PATH+"flux"+finalstate+massname+".txt"
    foutfinal = open(filename,"w")
    #____________________________________________________________________________

    jfactor=1e20  # GeV ^ 2 /cm ^ 5
    sv=3.0*1e-26

    #do not know how nditer put elements into an array
    vector=[]
    for x in np.nditer(xvals):
        vector.append(x)
    
        # dN/dE  = dNdx(x)/mass
        # prints two lines E , dark matter flux (multiplied by J_tot)
                
    zeros=np.logspace(math.log10(mass*GeVtoMeV),5.5+math.log10(GeVtoMeV),num=10) #Mab: Bins that will be filled with zero flux.

    #Mab: HERE TWO OPTIONS:
    # 1. Write the full spectrum as it is and add 10 more bins with zeros up to more than 100TeV so ctools doesn't have to extrapolate. This will keep the "bump" at the end of the "W" spectra, which causes weird behaviour for several masses in the upper limits plots. 
    #2.  Erase the last 6 bins of "W" spectra to eliminate the bump, letting ctools extrapolate this bins, and then add again 10 more bins with zeros up to more than 100 TeV.

    #OPTION 1: ADD ZEROS
    """
    with foutfinal as fff:
        for x in vector[118:179]:
            #for x in vector[119:180]:
            print >> fff, x*mass*GeVtoMeV, dNdx(x)/mass/mass/mass/GeVtoMeV/8./3.14*sv*jfactor 
            #print >> fff, x*mass*GeVtoMeV, dNdx(x)/mass/GeVtoMeV   # E  and dN/dE                
        
        for i in range(len(zeros)-1):
            print >> fff,zeros[i+1],1e-300 #Mab: Add zeros
    """
    # OPTION 2: ERASE THE BUMP FOR W+ADD ZEROS
    nbin=1
    with foutfinal as fff:
        for x in vector[118:179]:
            #for x in vector[119:180]:
            if nbin < 56 or (nbin >= 56 and (finalstate!= 'W' and finalstate!='Z')): #Mab: Don't print the last bins of 'W' spectra (the bump)
                print >> fff, x*mass*GeVtoMeV, dNdx(x)/mass/mass/mass/GeVtoMeV/8./3.14*sv*jfactor
                #print >> fff, x*mass*GeVtoMeV, dNdx(x)/mass/GeVtoMeV   # E  and dN/dE                
            nbin=nbin+1
        for i in range(len(zeros)-1):
            print >> fff,zeros[i+1],1e-300
                        
    foutfinal.close()

