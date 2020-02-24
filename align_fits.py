import astroalign as aa
from glob import glob
from astropy.io import fits
path='/home/vivek/phot_trial/i_band/'

nfiles=sorted(glob(path+'j0947*.fits'))
data=fits.open(path+'j0947_i_1.fits')
ref=data[0].data
for i in range(len(nfiles)):
    data=fits.open(nfiles[i])
    source=data[0].data
    im_al,footprint=aa.register(source,ref)
    fits.writeto(path+'j0947_corrected_%i.fits'%i,im_al,overwrite=True)
    print('No. %i done'%i)
