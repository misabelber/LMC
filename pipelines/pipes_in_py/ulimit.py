import gammalib
import ctools
import config
import numpy as np

#mass = 1 #Mass in TeV
#eref = 0.9*mass

masses = np.linspace(0.1,0.9,9)
masses = np.append(masses,np.linspace(1,9,9))
masses = np.append(masses,np.linspace(10,90,9))
masses = np.append(masses,np.linspace(100,300,3))

f = open("Ulimitsresults.txt","a")

for mass in masses:

    eref = 0.9*mass

    if eref > mass:
        eref = mass

    caldb_ = 'prod3b-v1' #Calibration files for ctobssim
    irf_ = 'South_z40_average_50h'

    Component = 'Irf+CR+DiffuseSources+PS' #Components that we are simulating: IRF, DM, Leptonic, Leptonic+Irf, etc. This is needed to read and write consistent filenames to be used by other programs.

    PATH_HERE = "../pipes_in_py" #Path where we are running
    PATH_MODEL = "../../models/" #Path where the model to simulate is stored
    PATH_OBS = config.DATA_PATH+"/Obs_"+Component+"/" 
    suf = "_KSP_100GeV-100TeV"

    expcube = PATH_OBS+"expcube_"+'LMC_'+Component+suf+'.fits'
    psfcube = PATH_OBS+"psfcube_"+'LMC_'+Component+suf+'.fits'
    bkgcube = PATH_OBS+"bkgcube_"+'LMC_'+Component+suf+'.fits'
    cntcube = PATH_OBS+"cntcube_"+'LMC_'+Component+suf+'.fits'
    
    bkgmodel = gammalib.GModels(PATH_OBS+"bkgmodel_"+'LMC_'+Component+suf+'000.xml')

    jfactor = config.DATA_PATH+'/DM/jfactor/annihil_LMC2D_FOVdiameter10.0deg_alphaint0.10deg_nside1024NFW-JFACTOR-Jsmooth-image.fits'

    spatial = gammalib.model.GModelSpatialDiffuseMap(jfactor)
    spectral = gammalib.model.GModelSpectralPlaw()
    prefactor = spectral[0]
    index = spectral[1]
    pivot = spectral[2]

    prefactor.scale(1e-16)
    pivot.scale(1e6)
    prefactor.free()
    pivot.fix()
    index.fix()
    
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
    ulimit["debug"] = True
    ulimit.execute()
    
    print >> f,mass,eref,ulimit.diff_ulimit(),ulimit.flux_ulimit(),ulimit.eflux_ulimit()
f.close()



