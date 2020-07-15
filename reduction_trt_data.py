#!/usr/bin/env python
# coding: utf-8

# # Reduction for the TRT data

# The script to work with the TRT data. This code will perform the reduction and susequently aperture photometry on the TRT data.

# Following packages are necessary to import:

import numpy as np
import matplotlib as mp
import ccdproc,os,sys,time,random
from astropy import units as u
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.time import Time
from glob import glob
#from astroquery.astrometry_net import AstrometryNet
from astropy.wcs import WCS
import astroalign as aa
from astropy.nddata import CCDData
from astropy.stats import sigma_clipped_stats, SigmaClip
from astropy.visualization import ImageNormalize, LogStretch
from matplotlib.ticker import LogLocator
from astropy.stats import SigmaClip, mad_std
import os
import matplotlib.animation as animation
from do_photometry import view_image, do_aperture_photometry

image_path='/mnt/eac99553-b108-449e-bf41-0716c975df8b/TRT-data/NGC 5947/band_wise/B_band/aligned/'
calibration_path='/mnt/eac99553-b108-449e-bf41-0716c975df8b/TRT-data/calibration/'
filename='flat_files/MF-1-B*fits'
flats=sorted(glob(os.path.join(calibration_path,filename)))
print(len(flats))
'''
for i in range(2,3):
    flat_file=fits.open(flats[i])
    flat_image=flat_file[0].data
    flat_image=flat_image/np.mean(flat_image)
    fits.writeto('normalised_flat.fits',flat_image,overwrite=True)
    view_image('normalised_flat.fits',10)
'''
ims=[]
sfilename='*fits'


s=sorted(glob(os.path.join(image_path,sfilename)))
print(len(s))
for k in range(70,len(s)):
    image=ccdproc.CCDData.read(s[k],unit='adu')
    head=image.header
    jd=head['JD']
    #view_image(s[k],1)
    do_aperture_photometry(s[k],k,6,jd)
    '''
    #view_image(s[k],1)
    dark=ccdproc.CCDData.read('mdark.fits',unit='adu')
    bias=ccdproc.CCDData.read('mbias.fits',unit='adu')
    flat=ccdproc.CCDData.read('normalised_flat.fits',unit='adu')

    bias_corrected = ccdproc.subtract_bias(image,bias)
    #dark_subtracted = ccdproc.subtract_dark(bias_corrected,dark,exposure_time=300,exposure_unit=u.second)#,scale=True)
    flat_corrected=ccdproc.flat_correct(bias_corrected,flat)
    #cr_cleaned = ccdproc.cosmicray_lacosmic(flat_corrected,readnoise=1, sigclip=5,satlevel=35500,niter=5,cleantype='meanmask',gain_apply=True)
    mean,std=np.mean(flat_corrected),np.std(flat_corrected)
    clean_file=s[k].replace('.fits','')
    fits.writeto('{}_corrected.fits'.format(clean_file),flat_corrected,header=head,overwrite=True)
    print('no. {} file written'.format(k+1))
    #print(mean,std)
    #view_image(s[k])
    #do_aperture_photometry(flat_corrected,k,6,jd)
    #print(jd,mag_back)
#plt.show()

    #plt.scatter(323,344, s=1000,edgecolor='red',linewidth=1.5, facecolor='none')
    #plt.imshow(nbf,cmap='gray_r')
    '''
