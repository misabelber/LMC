import gammalib
import ctools
import config

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
bkgmodel = PATH_OBS+"bkgmodel_"+'LMC_'+Component+suf+'.xml'

NObs = 100

#f = open("TSresults.txt","w")
f = open("TSresults.txt","a")
f.write("\n")

for i in range(6,NObs):

    obstring = "%03d" % i
    cntcube = PATH_OBS+"cntcube_"+'LMC_'+Component+suf+obstring+'.fits'
    result = PATH_OBS+"bkgmodel_"+'LMC_'+Component+suf+obstring+'.xml'

    like = ctools.ctlike()
    like["inobs"] = cntcube
    like["expcube"] = expcube
    like["psfcube"] = psfcube
    like["bkgcube"] = bkgcube
    like["inmodel"] = bkgmodel
    like["caldb"]=caldb_
    like["irf"]=irf_
    like["outmodel"] = result
    like["debug"] = True
    like["outcovmat"] = "covmatrix.fits"
    like.execute()

    model = gammalib.GModels(result)
    
    
    if (i==0):
        for component in model:
            Source = component.name()
            print >> f,Source,
    f.write("\n")
    for component in model:    
        TS = component.ts()
        print >> f,TS,

f.close()
