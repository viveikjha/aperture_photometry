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


image_path='/mnt/eac99553-b108-449e-bf41-0716c975df8b/TRT-data/NGC 5947/band_wise/I_band/aligned/'
calibration_path='/mnt/eac99553-b108-449e-bf41-0716c975df8b/TRT-data/calibration/'
filename='flat_files/MF-1-I*fits'
flats=sorted(glob(os.path.join(calibration_path,filename)))
print(len(flats))
for i in range(5,6):
    flat_file=fits.open(flats[i])
    flat_image=flat_file[0].data
    flat_image=flat_image/np.median(flat_image)
    fits.writeto('normalised_flat.fits',flat_image,overwrite=True)

sfilename='*fits'
s=sorted(glob(os.path.join(image_path,sfilename)))
print(len(s))
for k in range(0,80):
    image=ccdproc.CCDData.read(s[k],unit='adu')
    dark=ccdproc.CCDData.read('mdark.fits',unit='adu')
    bias=ccdproc.CCDData.read('mbias.fits',unit='adu')
    flat=ccdproc.CCDData.read('normalised_flat.fits',unit='adu')

    bias_corrected = ccdproc.subtract_bias(image,bias)
    #dark_subtracted = ccdproc.subtract_dark(bias_corrected,dark,exposure_time=300,exposure_unit=u.second)#,scale=True)
    flat_corrected=ccdproc.flat_correct(bias_corrected,flat)
    #cr_cleaned = ccdproc.cosmicray_lacosmic(flat_corrected,readnoise=1, sigclip=5,satlevel=35500,niter=5,cleantype='meanmask',gain_apply=True)
    mean,std=np.mean(flat_corrected),np.std(flat_corrected)
    print(mean,std)
    plt.imshow(flat_corrected[700:1350,700:1350],cmap='gray_r',vmin=mean-std,vmax=mean+std)
    plt.scatter(323,344, s=1000,edgecolor='red',linewidth=1.5, facecolor='none')
    #plt.imshow(nbf,cmap='gray_r')
    plt.show()