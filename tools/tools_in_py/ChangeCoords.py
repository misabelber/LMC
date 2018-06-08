import numpy as np
from astropy.io import fits
from astropy import wcs
import sys

DATA_PATH = "/home/queenmab/DATA/LMC/DM/jfactor/"

filename = DATA_PATH+"annihil_LMC2D_FOVdiameter10.0deg_alphaint0.10deg_nside1024NFW-JFACTOR-Jsmooth-image.fits"

hdul = fits.open(filename)

#Create a new WCS object, number of axes set from start

w = wcs.WCS(naxis=2)

# Set up an "Airy's zenithal" projection
# Vector properties may be set with Python lists, or Numpy arrays
w.wcs.crpix = [210.0, 210.0]
w.wcs.cdelt = np.array([0.02863961813842482, 0.02863961813842482])
w.wcs.crval = [80.0, -69.5]
w.wcs.ctype = ["RA---CAR", "DEC--CAR"]

w.wcs.set_pv([(2, 1, 45.0)])

# Now, write out the WCS object as a FITS header
header = w.to_header()

# header is an astropy.io.fits.Header object.  We can use it to create a new
# PrimaryHDU and write it to a file.
hdu = fits.PrimaryHDU(header=header,data=hdul[0].data)
# Save to FITS file
hdu.writeto(filename,overwrite=True)
