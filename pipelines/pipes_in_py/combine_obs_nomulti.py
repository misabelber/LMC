
#This program allows the simulation of long exposure observations (> 50h) and combines observations in an observation .xml file, so it requires less memory

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

ra = 80.0
dec = -69.5 

rad = 5 #Radius of ROI
emin = 0.03 #Minimum energy in TeV
emax = 100.0 #Maximum energy in TeV
duration = 180000 #Duration of each observation
deadc = 0.95 #Dead time

Total_time =1080000 #1440000 # 300h, Total amount of time required
nobs = int(Total_time/duration) #Number of observations of 50h that will be simulated and combined
 
#I need to define the calibration files in two ways because "ctobssim" can read simply this, but gammalib GObservation class needs the full path
caldb = gammalib.GCaldb(config.CTOOLS_PATH+"/gamma/share/caldb/data/cta/1dc/bcf/South_z20_50h") #Calibration Files for gammalib class
irf = "irf_file.fits"

caldb_ = '1dc' #Calibration files for ctobssim
irf_ = 'South_z20_50h'

debug = True
Component = 'Irf' #Components that we are simulating: IRF, DM, Leptonic, Leptonic+Irf, etc. This is needed to read and write consistent filenames to be used by other programs.

#Define Paths
PATH_HERE = "../pipelines/pipes_in_py" #Path where we are running
PATH_MODEL = "../models/" #Path where the model to simulate is stored
PATH_OBS = config.DATA_PATH+"/Obs_"+Component+"/" #Path to store the observation files. I create a different directory to store  each component (or set of components) data. 
if not os.path.exists(PATH_OBS): #If the observation path doesn't exist, create it.
    os.makedirs(PATH_OBS)

#Set Filenames wisely:
#This program is very useful to simulate Irf or Irf+Other things (for the observation file that we need to calculate limits), because Irf is what causes more memory problems.
#Other diffuse sources(alone) can be simulated directly with one ctobssim, but I use this program anyways.
#DM(alone) simulation is not implemented in this program, and is not needed.
if Component != 'DM': #DM Not implemented. Anycase, to simulate dark matter is not necessary to split observations. Use pipeline "producedata.py" instead.
    input_model = PATH_MODEL+'LMC_'+Component+'.xml'
    
outfile = PATH_OBS+'observations_'+'LMC_'+Component+'_'+str(nobs*duration/3600)+'h.xml' #List of Observations file that will be produced ('.xml')
file = open(outfile,'w') 

Obs_list = gammalib.GObservations()
xml = gammalib.GXml(outfile)

rndseed = random.randint(1,300000)

for i in range(nobs):

    # SIMULATE EACH 50h OBSERVATION
    filename = 'events_'+'LMC_'+Component+'_'+str(duration/3600)+'h_'+'0'+str(i+1)+'.fits'
    eventfile = gammalib.GFilename(filename)
    tstart = i*duration

    sim = ctools.ctobssim()
    sim["seed"]=i+rndseed
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

    obs.id(str(i+1))
    obs.name('events_'+'LMC_'+Component+'_'+str(duration/3600)+'h'+'_'+str(i+1))


    Obs_list.append(obs)

Obs_list.models(input_model)
Obs_list.save(outfile)

#Move observations to destination folder.Observations are simulated and stored in the present folder, but I prefer to store them somewhere else, so I move them there.

for src_dir,dirs,files in os.walk(PATH_HERE):
    dst_dir = PATH_OBS
    if not os.path.exists(dst_dir):
        os.makedirs(dst_dir)
    for file in files:
        src_file = os.path.join(PATH_HERE,file)
        dst_file = os.path.join(PATH_OBS,file)
        if file.endswith('.fits'):
            if os.path.exists(dst_file):
                os.remove(dst_file)
            shutil.move(src_file,PATH_OBS)


