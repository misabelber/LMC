import gammalib
import ctools
import config
import numpy as np
import sys
from scipy import optimize
def powerlaw(x,k,index,pivot):
    return k*(x/(pivot))**(-index)

def fitspectra(filename):
    data = np.genfromtxt(filename,unpack=True)
    energy = data[0]
    flux = data[1]
    energy = energy[flux>1e-60]
    flux = flux[flux>1e-60]
    params, params_covariance = optimize.curve_fit(powerlaw, 
                                                   energy[energy>=1e5], 
                                                   flux[energy>=1e5],
                                                   p0= [1.,2.,1.e3])
    return params
#mass = 1 #Mass in TeV
#eref = 0.9*mass

particle=sys.argv[1]
suf = sys.argv[2]
f = open("Ulimitsresults_empty_"+particle+suf+".txt","a")

#masses = [0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100]
masses = [0.2,0.3,0.5,0.6,0.8,1.,2.,3.,5.,6.,8.,10.,20.,30.,50.,100.]
#masses = [0.2,0.5,1.,5.,10.,100.]
#masses = [0.4,0.6,0.8,0.9,2,3,4,6,8,9,20,30,40,60,80,90,100]

for mass in masses:
    #Build wisely the strings for the filenames
    #______________________________________________
    if mass < 1:                                                        
        masstr = str(int(mass*1000))+'GeV'         
    else:
        masstr = str(int(mass))+'TeV'
    eref = 0.8*mass
    
    if eref > mass:
        eref = mass
 
    emin = 0.1
    emax = mass

    caldb_ = 'prod3b-v1' #Calibration files for ctobssim
    irf_ = 'South_z40_average_50h'

    Component = 'test_ctools' #Components that we are simulating: IRF, DM, Leptonic, Leptonic+Irf, etc. This is needed to read and write consistent filenames to be used by other programs.

    PATH_HERE = "../pipes_in_py" #Path where we are running
    PATH_MODEL = "../../models/test_ctools/" #Path where the model to simulate is stored
    PATH_OBS = "/home/bernardos/LMC/test_ctools/" 
#    suf = "_KSP_100GeV-100TeV"
    suf = "_empty"
    expcube = PATH_OBS+"expcube_"+'LMC_'+Component+'_'+suf+'.fits'
    psfcube = PATH_OBS+"psfcube_"+'LMC_'+Component+'_'+suf+'.fits'
    bkgcube = PATH_OBS+"bkgcube_"+'LMC_'+Component+'_'+suf+'.fits'
    cntcube = PATH_OBS+"cntcube_"+'LMC_'+Component+'_'+suf+'.fits'
    
    bkgmodel = gammalib.GModels(PATH_MODEL+"LMC_test_ctools_empty_bkgmodel.xml")
    
    jfactor ='/pnfs/ciemat.es/data/cta/mabel/LMC/DM/jfactor/annihil_LMC2D_FOVdiameter10.0deg_alphaint0.10deg_nside1024NFW-JFACTOR-Jsmooth-image.fits'
    
    dm_flux = '../../spectra/DM/flux'+particle+masstr+'.txt'
    
    spatial = gammalib.model.GModelSpatialDiffuseMap(jfactor)
    spectral = gammalib.model.GModelSpectralFunc(gammalib.GFilename(dm_flux),1)

    #Create a model from a power law fit
    '''
    p = fitspectra(dm_flux)
    print(p)
    spectral = gammalib.model.GModelSpectralPlaw()
    model["Prefactor"].value(p[0])
    model["Prefactor"].scale(1)
    model["Prefactor"].fix()
    model["Index"].value(p[1])
    model["Index"].scale(-1)
    model["Index"].fix()
    model["PivotEnergy"].value(p[2])
    model["PivotEnergy"].scale(1)
    model["PivotEnergy"].fix()
    '''
    model = gammalib.GModelSky(spatial,spectral)
    model.name('DM')
    bkgmodel.append(model)
        
    bkgmodel.save("fitmodel.xml")
    
    ulimit = ctools.ctulimit()
    ulimit["inobs"] = cntcube
    ulimit["expcube"] = expcube
    ulimit["psfcube"] = psfcube
    ulimit["bkgcube"] = bkgcube
    ulimit["inmodel"] = "fitmodel.xml"
    ulimit["caldb"]=caldb_
    ulimit["irf"]=irf_
    ulimit["srcname"] = 'DM'
    ulimit["eref"] = eref
    ulimit["emin"] = emin
    ulimit["emax"] = emax
    ulimit["debug"] = True
    ulimit.execute()
    
    print >> f,mass,eref,ulimit.diff_ulimit(),ulimit.flux_ulimit(),ulimit.eflux_ulimit(),emin,emax
f.close()



