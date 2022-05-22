from astropy.wcs import WCS
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astroquery.astrometry_net import AstrometryNet
from astropy import units as u
from astropy.coordinates import SkyCoord
ast=AstrometryNet()
ast.api_key= 'iqmqwvazpvolmjmn'
filename='sa92_standard/sa92_b4.fits'
wcs_header=ast.solve_from_image(filename)
wcs=WCS(wcs_header)
# print(wcs.to_header)
# hdu=fits.open(filename)[0]
# image=hdu.data
# header=hdu.header
# hdu.header.update(wcs.to_header)
# hdu.writeto('new.fits')
ra=wcs.wcs.crval[0]
dec=wcs.wcs.crval[1]
# print(ra,dec)
# px,py=wcs.wcs_pix2world(250,250,1)
# c=SkyCoord(ra=px*u.degree,dec=py*u.degree)
# print(c.to_string('hmsdms'))


### Now querying the surveys for the source (for calibration)

from astroquery.vizier import Vizier
import astropy.coordinates as coord
result = Vizier.query_region(coord.SkyCoord(ra=ra, dec=dec,unit=(u.deg, u.deg),frame='icrs'),width="18m",catalog=["NOMAD", "UCAC"])

bmag=result[2]['Bmag']
print(np.max(bmag))

