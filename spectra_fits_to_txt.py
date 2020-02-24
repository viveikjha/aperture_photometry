from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import csv,sys


filename='spec-4401-55510-0440.fits' # You can use glob to load multiple files here.
data=fits.open(filename)
lam=10**data[1].data['loglam']        # observed wavelength
flux=data[1].data['flux']             # observed fluc
err=1./np.sqrt(data[1].data['ivar'])  # 1 sigma error
z=data[2].data['z'][0]  # redshift
ra=data[2].data['PLUG_RA'][0] # RA 
dec=data[2].data['PLUG_DEC'][0] # DEC
plateid = data[0].header['plateid']   # SDSS plate ID
mjd = data[0].header['mjd']           # SDSS MJD
fiberid = data[0].header['fiberid'] 

mod_file=filename.replace('fits','')
with open (mod_file+'txt', "w") as f:
    writer= csv.writer(f, delimiter=" ")
    writer.writerows(zip(lam,flux,err))
    sys.exit()
