import gammalib
import math
import numpy as np

#Open model file

models = gammalib.GModels("../../models/LMC_model.xml")

centerRA = 80.0
centerDEC = -69.5
centerRA = math.radians(centerRA)
centerDEC = math.radians(centerDEC)

for model in models:
    name = model.name()
    tipo = model.type()
    if tipo=="PointSource":
        RA = model.spatial()["RA"].value()
        dec = model.spatial()["DEC"].value()
        RA = math.radians(RA)
        dec = math.radians(dec)
        DA = math.sin(dec)*math.sin(centerDEC)+math.cos(dec)*math.cos(centerDEC)*math.cos(RA-centerRA)
        DA = math.degrees(math.acos(DA))
        print DA
        if DA > 6.:
            for par in model:
                par.fix()
                
models.save("../models/LMC_closer_files.xml")
