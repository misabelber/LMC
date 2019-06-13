# This program simulates, bins and creates model cubes for several Dark Matter models using a Pointing pattern similar to the one propopsed in the CTA KSP for LMC.

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
import sys

PATH_MODEL = "../../models/" #Path where models are stored
PATH_OBS = "/home/bernardos/LMC/test_ctools/test_ctools_call/DMobs/" #Path to store resulting .fits files.
PATH_HERE = "../pipes_in_py/" #Path where we are running
# If the observation path doesn't exist, create it.
if not os.path.exists(PATH_OBS): 
    os.makedirs(PATH_OBS)

particle = sys.argv[1] #Final state particle

#Pointing
centerx = 80.0
centery = -69.5

r = 2.0
#r = 5.0
ra_list=np.zeros(6)
dec_list=np.zeros(6)
for i in range(0,6):
    angle=(i-1)*2*math.pi/6
    ra_list[i] = r*math.cos(angle)+centerx
    dec_list[i] = r*math.sin(angle)+centery

#Observation variables

rad = 3 #Radius of ROI
emin = 0.1 #Minimum energy in TeV
emax = 100.0 #Maximum enery in TeV
tstart = 0.0 #Startind time
duration = 204000 #Ending time
deadc = 0.95 #Dead time

binsz = 0.1 #Spatial binning
nxpix = 100
nypix = 100

enumbins = 20

caldb = gammalib.GCaldb(config.CTOOLS_PATH+"/share/caldb/data/cta/prod3b-v1/bcf/South_z40_average_50h") #Calibration Files for gammalib class
irf = "irf_file.fits"

caldb_ = 'prod3b-v1' #Calibration files for ctobssim
irf_ = 'South_z40_average_50h'

#masses = [0]
#masses = [0.100,0.200,0.500,1,5,10,50,100]
#masses = [0.200,0.300,0.400,0.500,0.600,0.800,1,4,5,8,10,40,50,80,100]
masses = [0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100]
#masses = [0]

NObs = 10

for mass in masses:
    rndseed = random.randint(1,300000)
    #Build wisely the strings for the filenames
    #______________________________________________
    if mass < 1:
        masstr = str(int(mass*1000))+'GeV'
    else:
        masstr = str(mass)+'TeV'
    specname = 'flux'+particle+masstr+'.txt'
    modelname = 'dm_LMC_'+particle+masstr
    suf = "_jfactorNFW"
    #modelname = 'dm_LMC_Crab.xml'
    
    model = PATH_MODEL+modelname+suf+'.xml'
    time = str(int(duration/3600))
    outfile = PATH_OBS+'observations_'+'LMC_'+modelname+suf+'.xml' #List of Observations file that will be produced ('.xml')

    rndseed = random.randint(1,300000)

    for i in range(0,NObs):

        obstring = "%03d" % i

        cntcube = PATH_OBS+"cntcube_"+modelname+suf+obstring+'.fits'
        modcube = PATH_OBS+"modcube_"+modelname+suf+obstring+'.fits'
        expcube = PATH_OBS+"expcube_"+modelname+suf+obstring+'.fits'
        psfcube = PATH_OBS+"psfcube_"+modelname+suf+obstring+'.fits'
        #_______________________________________________
    
        file = open(outfile,'w')
        Obs_list = gammalib.GObservations()
        xml = gammalib.GXml(outfile)
    
        # RUN SIMULATION
        number=1
        for ra in ra_list:
            dec = dec_list[number-1]
            filename = 'events_'+'LMC_'+modelname+'_KSP'+'0'+str(number)+'.fits'
            eventfile = gammalib.GFilename(filename);

            sim = ctools.ctobssim()
            sim["inmodel"]=model
            sim["seed"] = number+rndseed
            sim["outevents"]=filename
            sim["caldb"]=caldb_
            sim["irf"]=irf_
            sim["ra"]=ra
            sim["dec"]=dec
            sim["rad"]=rad
            sim["tmin"]=gammalib.GTime(tstart)
            sim["tmax"]=gammalib.GTime(tstart+duration)
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
            obs.name('events_'+'LMC_'+modelname+'_KSP'+'0'+str(number))
            
            Obs_list.append(obs)    
            number=number+1

        Obs_list.models(model)
        Obs_list.save(outfile)

        #Move observations to destination folder.Observations are simulated and stored in the present folder, but I prefer to store them somewhere else, so I move them there.

        for src_dir,dirs,files in os.walk(PATH_HERE):
            dst_dir = PATH_OBS
            if not os.path.exists(dst_dir):
                os.makedirs(dst_dir)
            for file in files:
                src_file = os.path.join(PATH_HERE,file)
                dst_file = os.path.join(PATH_OBS,file)
                if file.startswith('events_'+'LMC_'+modelname):
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
            exp = ctools.ctexpcube()
            exp["inobs"] = outfile
            exp["incube"] = cntcube
            exp["caldb"]=caldb_
            exp["irf"]=irf_
            exp["outcube"] = expcube 
            exp["debug"]=True 
            exp.execute()
            
            psf = ctools.ctpsfcube()
            psf["inobs"] = outfile
            psf["incube"] = cntcube
            psf["caldb"]=caldb_
            psf["irf"]=irf_
            psf["outcube"] = psfcube 
            psf["debug"]=True 
            psf.execute()
            
            #PRODUCE MODELCUBE
            mod = ctools.ctmodel()
            mod["inobs"]= outfile
            mod["inmodel"]= model
            mod["outcube"] = modcube
            mod["incube"] = cntcube
            mod["expcube"] = expcube 
            mod["psfcube"] = psfcube
            mod["bkgcube"] = "NONE"
            mod["caldb"]=caldb_
            mod["irf"]=irf_
            mod["debug"]=True
            mod.execute()
   
