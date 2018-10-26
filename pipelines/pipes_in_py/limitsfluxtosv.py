from scipy.interpolate import interp1d
import numpy as np
#from io import StringIO
from StringIO import StringIO
from matplotlib import pyplot as plt
import io
import math

filename = '../../tools/tools_in_py/AtProduction_gammas.dat'
finalstate = "W"  # choose particle W, b, etc

with open(filename) as f:
    lines = (line for line in f if not line.startswith('#'))
    data = np.genfromtxt (lines, names = True ,dtype = None)

masses = np.linspace(0.1,0.9,9)
masses = np.append(masses,np.linspace(1,9,9))
masses = np.append(masses,np.linspace(10,90,9))
masses = np.append(masses,np.linspace(100,300,3))

massvals = data["mDM"]

GeVtoMeV = 1000

for mass in masses:
    


