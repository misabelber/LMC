import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import astropy_mpl_style
from astropy.wcs import WCS
from matplotlib.patches import Circle
import math
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm

plt.style.use(astropy_mpl_style)


image_file = fits.open('/scratch/bernardos/LMC/Obs_CR+DiffuseSources+PS/skymap_LMC_Irf+CR+DiffuseSources+PS_KSP_100GeV-100TeV.fits')

hdu = image_file[0]
image_data = hdu.data

wcs = WCS(hdu.header)

ax = plt.subplot(projection=wcs)

#Plot the SkyMap

cmap = plt.get_cmap('inferno')
cmap.set_bad((0,0,0))
plt.imshow(image_data, cmap=cmap, origin='lower',norm=LogNorm())
plt.grid(color='white',ls='solid')
plt.xlabel("Right Ascension J2000")
plt.ylabel("Declination J2000")
plt.colorbar(label="# of Counts (Log scale)")

#Plot the Pointings
centerx = 80.0
centery = -69.5
rad = 3
r = 2.0                                                                 
for i in range(0,6):                                                                                  
    angle=(i-1)*2*math.pi/6                                                                           
    ra = r*math.cos(angle)+centerx                                                            
    dec = r*math.sin(angle)+centery
    c = Circle((ra,dec),rad,transform=ax.get_transform('fk5'),edgecolor='green',linestyle='--', facecolor='none')
    if i==0:
        c.set_label("Pointing")
    ax.add_patch(c)

#Plot Sources position
filename = 'pointsourcespositions.txt'                                                                                   
dat = pd.read_table(filename,delim_whitespace=True)

ax.scatter(dat['RA'],dat['dec'], transform=ax.get_transform('fk5'), s=300,
           edgecolor='white', facecolor='none',label='Point sources')

for label, x, y in zip(dat['Sourcename'], dat['RA'], dat['dec']):

    if label!='J0525-696' and label!='J0537-691':
        ax.text(x-0.3,y+0.2,label,transform=ax.get_transform('fk5'),color='white')
    else:
        ax.text(x+4,y,label,transform=ax.get_transform('fk5'),color='white')

#Plot Extended sources

filename = 'extsourcespositions.txt'
dat = pd.read_table(filename,delim_whitespace=True)
s = (8000,4000,5000)

ax.scatter(dat['RA'],dat['dec'], transform=ax.get_transform('fk5'),
           edgecolor='cyan', facecolor='none',linestyle="-.",label='Extended sources', s=s)

for label, x, y in zip(dat['Sourcename'], dat['RA'], dat['dec']):
    
        ax.text(x-0.3,y+0.2,label,transform=ax.get_transform('fk5'),color='cyan')
    


plt.legend(markerscale=0.2)
plt.show()
