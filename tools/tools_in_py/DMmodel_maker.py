
#This program is to produce '.xml' models for different dark matter masses and final states (particles).

import gammalib
import math
import numpy as np

#PATHS

LMC_PATH = "/afs/ciemat.es/user/b/bernardos/GitHub/LMC/"
DATA_PATH = "/pnfs/ciemat.es/data/cta/mabel/LMC/"
PATH_SPEC=LMC_PATH+'/spectra/DM/' #Path where the spectra (filefunctions) are stored.
PATH_MODELS=LMC_PATH+'/models/'  #Path to store the models created (.xml)

#Open model file
models = gammalib.GModels(PATH_MODELS+"example_model.xml") #Open an example model file to modify. Its a diffuse source where the spatial part is a counts cube (J Factor fits file) and the spectral part is a FileFunction. 

particle = 'b' #Final state particle (W, b, Mu...)

#Masses of dark matter particle in TeV
masses = [0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100]
#masses = [0.300,0.400,0.600,0.700,0.800,0.900,2,3,4,5,6,7,8,9,10]

#jfactorname = "DM/jfactor/jfactorLMC.fits"

jfactorname = "DM/jfactor/annihil_LMC2D_FOVdiameter10.0deg_alphaint0.10deg_nside1024NFW-JFACTOR-Jsmooth-image.fits"

for mass in masses:
    #Standarized name for the produced model so the other programs can use it easily:
    #_______________________________________________________________________________
    if mass < 1:
        masstr = str(int(mass*1000))+'GeV'
    else:
        masstr = str(mass)+'TeV'
    specname = 'flux'+particle+masstr+'.txt'
    modelname = 'dm_LMC_'+particle+masstr
    suf = '_jfactorNFW'
    #_______________________________________________________________________________

    for model in models:
        modelcontainer = gammalib.GModels()
        model.spectral().filename(PATH_SPEC+specname) #Put the spectrum for the specific mass in the model file
        spectral = model.spectral()
        spatial = gammalib.GModelSpatialDiffuseMap(DATA_PATH+jfactorname)
        model = gammalib.GModelSky(spatial,spectral)
        modelcontainer.append(model)
    modelcontainer.save(PATH_MODELS+modelname+suf+'.xml') #Save the model
