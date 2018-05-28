
# This program simulates, bins and creates model cubes for several Dark Matter models

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

PATH_MODEL = "../../models/" #Path where models are stored
PATH_OBS = config.DATA_PATH+"/Obs_PS/" #Path to store resulting .fits files.
PATH_HERE = "../pipes_in_py/" #Path where we are running


if not os.path.exists(PATH_OBS): #If the observation path doesn't exist, create it.
    os.makedirs(PATH_OBS)
#Set Observation Variables

centerx = 80.0
centery = -69.5

r = 3.0
ra_list=np.zeros(6)
dec_list=np.zeros(6)
for i in range(0,6):
    angle=(i-1)*2*math.pi/6
    ra_list[i] = r*math.cos(angle)+centerx
    dec_list[i] = r*math.sin(angle)+centery


#Observation variables

rad = 5 #Radius of ROI
emin = 0.03 #Minimum energy in TeV
emax = 100.0 #Maximum enery in TeV
tstart = 0.0 #Starting time
duration = 180000 #Ending time
deadc = 0.95 #Dead time

binsz = 0.5 #Spatial binning
nxpix = 20
nypix = 20
enumbins = 20

caldb = gammalib.GCaldb(config.CTOOLS_PATH+"/share/caldb/data/cta/1dc/bcf/South_z20_50h") #Calibration Files for gammalib class
irf = "irf_file.fits"

caldb_= "1dc"
irf_="South_z20_50h"

#Get model name
models = gammalib.GModels(PATH_MODEL+"LMC_closer_files.xml")
for model in models:
    modelcontainer = gammalib.GModels()
    modelname = model.name()
    srctype = model.type()
    par = model.spectral().at(0)
    if srctype=='PointSource' and par.is_free():
        modelcontainer.append(model)
        filename = 'LMC_'+modelname+'.xml'
        
    
        model = PATH_MODEL+'LMC_'+modelname+".xml"
        
        cntcube = PATH_OBS+"cntcube_"+"LMC_"+modelname+'_KSPpointing_v2_.fits'
        modcube = PATH_OBS+"modcube_"+"LMC_"+modelname+'_KSPpointing_v2_.fits'
        outfile = PATH_OBS+'observations_'+'LMC_'+modelname+'_KSPpointing_v2_.xml' #List of Observations file that will be produced ('.xml')
        file = open(outfile,'w')

        #_______________________________________________
    
        Obs_list = gammalib.GObservations()
        xml = gammalib.GXml(outfile)
        
        rndseed = random.randint(1,300000)
        number = 1

        for ra in ra_list:

            dec = dec_list[number-1]
            # SIMULATE EACH 50h OBSERVATION
            filename = 'events_'+'LMC_'+modelname+'_'+str(duration/3600)+'h_'+'0'+str(number)+'.fits'
            eventfile = gammalib.GFilename(filename)
            
            sim = ctools.ctobssim()
            sim["inmodel"]=model
            sim["seed"] = number+rndseed
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
            obs.name('events_'+'LMC_'+modelname+'_'+str(duration/3600)+'h'+'_'+str(number))

            
            Obs_list.append(obs)

            number=number+1

        Obs_list.models(model)
        Obs_list.save(outfile)

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

        #PRODUCE MODELCUBE
        mod = ctools.ctmodel()
        mod["inobs"]= outfile
        mod["inmodel"]= model
        mod["outcube"] = modcube
        mod["incube"] = cntcube
        mod["expcube"] = "NONE" 
        mod["psfcube"] = "NONE"
        mod["bkgcube"] = "NONE"
        mod["caldb"]=caldb_
        mod["irf"]=irf_
        mod["debug"]=True
        mod.execute()
