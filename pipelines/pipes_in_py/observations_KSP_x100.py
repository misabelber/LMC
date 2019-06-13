
#This program simulate observations with a Pointing Pattern similar to the one suggested by CTA KSP for LMC.
#6 pointings around LMC center separated 2-3 degrees.

import gammalib
import ctools
import cscripts
import numpy as np
import multiprocessing
import math
import shutil
import os
import random
import config


#Set Observation Variables

centerx = 80.0
centery = -69.5

#Define pointings
r = 2.0
ra_list=np.zeros(6)
dec_list=np.zeros(6)
for i in range(0,6):
    angle=(i-1)*2*math.pi/6
    ra_list[i] = r*math.cos(angle)+centerx
    dec_list[i] = r*math.sin(angle)+centery
    
rad = 3 #Radius of ROI
emin = 0.1 #Minimum energy in TeV
emax = 100.0 #Maximum energy in TeV
tstart = 0.0 #Starting time
duration = 204000 #Duration of each observation
deadc = 0.95 #Dead time

#Binning
binsz = 0.1 #Spatial binning
nxpix = 100
nypix = 100
enumbins = 20 #Energy bins

#I need to define the calibration files in two ways because "ctobssim" can read simply this, but gammalib GObservation class needs the full path
caldb = gammalib.GCaldb(config.CTOOLS_PATH+"/share/caldb/data/cta/prod3b-v1/bcf/South_z40_average_50h") #Calibration Files for gammalib class
irf = "irf_file.fits"

caldb_ = 'prod3b-v1' #Calibration files for ctobssim
irf_ = 'South_z40_average_50h'

debug = True

Component = 'test_ctools' #Components that we are simulating: IRF, DM, Leptonic, Leptonic+Irf, etc. This is needed to read and write consistent filenames to be used by other programs.

#Define Paths
PATH_HERE = "../pipes_in_py" #Path where we are running
PATH_MODEL = "../../models/" #Path where the model to simulate is stored
PATH_OBS = "/home/bernardos/LMC/test_ctools/" #Path to store the observation files. I create a different directory to store  each component (or set of components) data. 
if not os.path.exists(PATH_OBS): #If the observation path doesn't exist, create it.
    os.makedirs(PATH_OBS)

#Set Filenames wisely:

suf = "_hesssources"
input_model = PATH_MODEL+'test_ctools/LMC_'+Component+suf+'.xml' 

outfile = PATH_OBS+'observations_'+'LMC_'+Component+suf+'.xml' #List of Observations file that will be produced ('.xml')
file = open(outfile,'w') 


rndseed = random.randint(1,300000)


NObs = 10

for i in range(0,NObs):

    Obs_list = gammalib.GObservations()
    xml = gammalib.GXml(outfile)

    obstring = "%03d" % i

    cntcube = PATH_OBS+"cntcube_"+'LMC_'+Component+suf+obstring+'.fits'
    modcube = PATH_OBS+"modcube_"+'LMC_'+Component+suf+obstring+'.fits'
    
    expcube = PATH_OBS+"expcube_"+'LMC_'+Component+suf+obstring+'.fits'
    psfcube = PATH_OBS+"psfcube_"+'LMC_'+Component+suf+obstring+'.fits'
    bkgcube = PATH_OBS+"bkgcube_"+'LMC_'+Component+suf+obstring+'.fits'
    
    bkgmodel = PATH_OBS+"bkgmodel_"+'LMC_'+Component+suf+'.xml'

    number = 1

    for ra in ra_list:
    
        dec = dec_list[number-1]
        # SIMULATE EACH OBSERVATION
        filename = 'events_'+'LMC_'+Component+'_'+str(duration/3600)+'h_'+'0'+str(number)+'.fits'
        eventfile = gammalib.GFilename(filename)

        sim = ctools.ctobssim()
        sim["seed"]=i+number+rndseed
        sim["inmodel"]=input_model
        sim["outevents"]=filename
        sim["caldb"]=caldb_
        sim["irf"]=irf_
        sim["ra"]=ra
        sim["dec"]=dec
        sim["rad"]=rad
        sim["tmin"]=tstart
        sim["tmax"]=tstart+duration
        sim["emin"]=emin
        sim["emax"]=emax
        sim["debug"]=True
        sim.execute()
        
        #STORE THE OBSERVATION IN GObservations Class and store it in the .xml output file.
        #Allocate CTA observation
        obs = gammalib.GCTAObservation()
        #Set pointing direction
        pntdir = gammalib.GSkyDir()
        pntdir.radec_deg(ra,dec)
        
        pnt = gammalib.GCTAPointing()
        pnt.dir(pntdir)
        obs.pointing(pnt)
        
        #Set ROI
        
        roi = gammalib.GCTARoi()
        instdir = gammalib.GCTAInstDir()
        instdir.dir(pntdir)
        roi.centre(instdir)
        roi.radius(rad)
        
        #Set GTI
        
        gti = gammalib.GGti()
        start = gammalib.GTime(tstart)
        stop = gammalib.GTime(tstart+duration)
        gti.append(start,stop)
        
        #Set Energy Boundaries
        ebounds = gammalib.GEbounds()
        e_min = gammalib.GEnergy(emin,'TeV')
        e_max = gammalib.GEnergy(emax,'TeV')
        ebounds.append(e_min,e_max)
        
        #Allocate event list
        events = gammalib.GCTAEventList(eventfile)
        obs.eventfile(eventfile)
        events.roi(roi)
        events.gti(gti)
        events.ebounds(ebounds)
        obs.events(events)
        
        #Set instrument response
        obs.response(irf,caldb)

        #Set ontime, livetime, and deadtime correction factor
        obs.ontime(duration)
        obs.livetime(duration*deadc)
        obs.deadc(deadc)
        
        obs.id(str(number))
        obs.name('events_'+'LMC_'+Component+'_'+str(duration/3600)+'h'+'_'+str(number))
        
        
        Obs_list.append(obs)
        
        number=number+1
        
    Obs_list.models(input_model)
    Obs_list.save(outfile)

    #Move observations to destination folder.Observations are simulated and stored in the present folder, but I prefer to store them somewhere else because of storaging space, so I move them there.

    for src_dir,dirs,files in os.walk(PATH_HERE):
        dst_dir = PATH_OBS
        if not os.path.exists(dst_dir):
            os.makedirs(dst_dir)
        for file in files:
            src_file = os.path.join(PATH_HERE,file)
            dst_file = os.path.join(PATH_OBS,file)
            if file.startswith('events_'+'LMC_'+Component):
                if os.path.exists(dst_file):
                    os.remove(dst_file)
                shutil.move(src_file,PATH_OBS)

    # BIN THE DATA
    binn = ctools.ctbin()
    binn["inobs"]=outfile
    binn["outcube"]=cntcube
    binn["coordsys"] = "CEL"
    binn["proj"] = "CAR"
    binn["ebinalg"] = "LOG"
    binn["xref"]=centerx
    binn["yref"]=centery
    binn["nxpix"] = nxpix
    binn["nypix"] = nypix
    binn["binsz"] = binsz
    binn["enumbins"] = enumbins
    binn["emin"]=emin
    binn["emax"]=emax
    binn["debug"]=True
    binn.execute()
    
    if i==0:
        
        #Pre-comute the binned response
    
        exp = ctools.ctexpcube()
        exp["inobs"] = outfile
        exp["incube"] = cntcube
        exp["caldb"]=caldb_
        exp["irf"]=irf_
        exp["outcube"] = expcube 
        exp["debug"]=True 
        #exp.execute()
    
        psf = ctools.ctpsfcube()
        psf["inobs"] = outfile
        psf["incube"] = cntcube
        psf["caldb"]=caldb_
        psf["irf"]=irf_
        psf["outcube"] = psfcube 
        psf["debug"]=True 
        #psf.execute()
    
        bkg = ctools.ctbkgcube()
        bkg["inobs"] = outfile
        bkg["incube"] = cntcube
        bkg["caldb"]=caldb_
        bkg["irf"]=irf_
        bkg["inmodel"]= input_model
        bkg["outcube"] = bkgcube
        bkg["outmodel"] = bkgmodel
        bkg["debug"]=True 
        #bkg.execute()
    
        #PRODUCE MODELCUBE
        mod = ctools.ctmodel()
        mod["inobs"]= outfile
        mod["inmodel"]= input_model
        mod["outcube"] = modcube
        mod["incube"] = cntcube
        mod["expcube"] = expcube 
        mod["psfcube"] = psfcube
        mod["bkgcube"] = "NONE"
        mod["caldb"]=caldb_
        mod["irf"]=irf_
        mod["debug"]=True
        #mod.execute()
        
    
    
