#!/usr/bin/env python
# coding: utf-8



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
#from pyraf import iraf
#from iraf import noao,imred,specred
from astropy.nddata import CCDData
from astropy.stats import sigma_clipped_stats, SigmaClip
from astropy.visualization import ImageNormalize, LogStretch
from matplotlib.ticker import LogLocator
from photutils.background import Background2D, MedianBackground, MeanBackground,SExtractorBackground
from photutils.aperture import  CircularAperture, CircularAnnulus, aperture_photometry
from photutils.centroids import centroid_2dg
from astropy.stats import SigmaClip, mad_std
from photutils.utils import calc_total_error
from photutils.detection import find_peaks,DAOStarFinder, IRAFStarFinder
from astropy.coordinates import SkyCoord
from astropy import coordinates as coord
from astropy import units as u
from astroquery.simbad import Simbad
import pandas as pd
import os, csv

import warnings
warnings.filterwarnings("ignore")

# Not really required usually

#import math
#import re
#import astroalign as aa   # use this for alignment function.
#from astropy.visualization import ImageNormalize, LogStretch
#from matplotlib.ticker import LogLocator

### Cleaning images

path=os.getcwd()
filename='crr_ALHe090282_astro_r.fits'


def view_image(name,time):
    '''
    This module is meant for viewinga large number of images. The images
    can be seen as an animation. Will be helpful for large datasets as it doesn't require a lot of pointing and clicking.
    INPUT:
    name: The name of the image file.
    time: No. of seconds you'd like to see the image on screen

    OUTPUT:

    the images seen as animation, with a gap of 1 second between successive images.
    '''

    ims=[]
    data=fits.open(name)
    image=data[0].data
    head=data[0].header

    mean,std=np.mean(image),np.std(image)
    #image=image/np.median(image)
    im=plt.imshow(image,cmap='gray_r',origin='lower',vmin=mean-2*std,vmax=mean+2*std)
    #plt.scatter(124, 146, s=400,edgecolor='red',linewidth=1, facecolor='none')
    ims.append([im])
    plt.pause(time)
    plt.clf()

def clean_the_images(path):
    '''
    This module is meant for cleaning the images. The tasks to be included are: bias correction,
    flat correction, trimming, overscan as well as the cosmic ray removal from the science cases.
    (For the time we are skipping the overscan and trimming part and only performing the bias, flat and cosmic ray correction.)

    INPUT:
    path: The directory where the images are kept (string)
    filename: The first few characters and the extension of the images (string). Example:
    j0946*fits, HD1104*fits etc.

    OUTPUT:

    cleaned images in the new directory: path/cleaned
    '''
    dir = path
    gain = 2 * u.electron / u.adu  # gain and readout noise are properties of the CCD and will change for different CCDs.
    readnoise = 7.5 * u.electron

    bias_files = sorted(glob(os.path.join(dir,'B-*fits')))
    biaslist = []
    for i in range (0,len(bias_files)):
        data= ccdproc.CCDData.read(bias_files[i],unit='adu')
        #data = ccdproc.create_deviation(data, gain=gain, readnoise=readnoise)
        #data= data-(data.uncertainty.array)
        biaslist.append(data)
    masterbias = ccdproc.combine(biaslist,method='average',sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
                             sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std)
    print('Step1:Master Bias generated.')
    print(" Mean and median of the masterbias: ",np.mean(masterbias), np.median(masterbias))


    flat_files=sorted(glob(os.path.join(dir,'F-*.fits')))
    flatlist = []
    for j in range(0,len(flat_files)):
        flat=ccdproc.CCDData.read(flat_files[j],unit='adu')
        flat_bias_removed=ccdproc.subtract_bias(flat,masterbias)
        flatlist.append(flat_bias_removed)

        def inv_median(a):
            return 1 / np.median(a)

    masterflat = ccdproc.combine(flatlist,method='median', scale=inv_median,
                                 sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
                                 sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std)
    print('Step 2: Master flat generated')
    print(" Mean and median of the masterflat: ",np.mean(masterflat), np.median(masterflat))

    files=sorted(glob(os.path.join(dir,'S-*.fits')))
    for i in range(0,len(files)): 
        image=ccdproc.CCDData.read(files[i],unit='adu')
        header=fits.getheader(files[i],0)
        bias_subtracted = ccdproc.subtract_bias(image, masterbias)
        flat_corrected = ccdproc.flat_correct(bias_subtracted, masterflat)#, 
        cr_cleaned = ccdproc.cosmicray_lacosmic(flat_corrected,readnoise=7.5, sigclip=5, verbose=True)
        print('Step 3: Cosmic rays removed.')
        clean_file=files[i].replace('.fits','')
        fits.writeto(clean_file+'_cleaned.fits',cr_cleaned, header=header,overwrite=True) 

def align_the_images(path,filename,ref_image):



    '''
    This function is meant for alignment of the images with respect to a reference image. To do this task we are using the astro-align package.

    INPUT:

    path: The directory where the images are kept (string)
    filename: The first few characters and the extension of the images (string). Example:
    j0946*fits, HD1104*fits etc.
    ref_image: Reference image which will be used to align the images.

    OUTPUT:

    aligned images.

    '''



    nfiles=sorted(glob(path+filename))
    image_data=fits.open(path+ref_image)
    reference_image=image_data[0].data
    for i in range(len(nfiles)):
        image_data=fits.open(nfiles[i])
        source_image=image_data[0].data
        header=image_data[0].header
        image_aligned,footprint=aa.register(source_image,reference_image)

        aligned_file=nfiles[i].replace('.fits','')
        fits.writeto(aligned_file+'_aligned'+'.fits',image_aligned,header,overwrite=True)

        print('Alignment of image no. %i  done'%i)

def do_photometry(path,filename,fwhm=4,sigma=3):
    '''
    To identify the sources based on either IRAFStarFinder, or Star/DAO finder.
    Inputs:
    path,filenam, guess fwhm and sigma
    Output:

    Returns a file with the instrumental magnitudes
    '''
    data= fits.open(os.path.join(path,filename))
    image=data[0].data
    header=data[0].header
    jd=header['JD']
    wcs = WCS(data[0].header)
    sigma_clip = SigmaClip(sigma=2, maxiters=10)
    bkg_estimator = SExtractorBackground()
    bkg = Background2D(image, (10,10), filter_size=(3, 3),sigma_clip=sigma_clip, bkg_estimator=bkg_estimator) 
    back=bkg.background
    mean, median, std = np.mean(image), np.median(image), np.std(image)
    mean, median, std = sigma_clipped_stats(image, sigma=3.0)  
    iraf_find = IRAFStarFinder(fwhm=fwhm, threshold=sigma*std, sigma_radius=5, minsep_fwhm=3, sharplo=0.1, sharphi=100.0, roundlo=0.0, roundhi=100.0,)
    sources=iraf_find(image-back) # customizable
    print(len(sources),' sources have been found.')
    fwhm_0=sources['fwhm']
    fwhm=np.median(fwhm_0)
    positions = [(sources['xcentroid'][i], sources['ycentroid'][i]) for i in range(len(sources))]
    ra , dec = wcs.all_pix2world(sources['xcentroid'],sources['ycentroid'],0)
    radii=[ fwhm,2*fwhm, 3*fwhm,4*fwhm,5*fwhm]
    apertures = [CircularAperture(positions, r=r) for r in radii]
    an_ap = CircularAnnulus(positions, r_in=6*fwhm, r_out=6.2*fwhm)
    print('Aperture selected...')
    exposure=header['EXPTIME']
    effective_gain=exposure # Need to verify this one. Why effective gain= exposure?
    error=calc_total_error(image,back,effective_gain)
    print('Doing photometry..')
    phot_table = aperture_photometry(image-back, apertures,error=error)
    phot_table2=aperture_photometry(image-back,an_ap)

        
    bkg_mean = phot_table2['aperture_sum'] / an_ap.area
    bkg_sum = bkg_mean * an_ap.area


    final_sum0=phot_table['aperture_sum_0']-bkg_sum
    final_sum1=phot_table['aperture_sum_1']-bkg_sum
    final_sum2=phot_table['aperture_sum_2']-bkg_sum
    final_sum3=phot_table['aperture_sum_3']-bkg_sum
    final_sum4=phot_table['aperture_sum_4']-bkg_sum
            
    mag_back=-2.5*np.log10(bkg_mean/exposure)
    mag_0=-2.5*np.log10(final_sum0/exposure)
    mag_1=-2.5*np.log10(final_sum1/exposure)
    mag_2=-2.5*np.log10(final_sum2/exposure)
    mag_3=-2.5*np.log10(final_sum3/exposure)
    mag_4=-2.5*np.log10(final_sum4/exposure)

    try:
        wcs = WCS(header)

        fig = plt.figure()
        ax = plt.subplot(projection=wcs)
        ax.imshow(image, cmap='gray', origin='lower', vmin=mean-1.5*std, vmax=mean+1.5*std)
        ax.set_xlabel('Right Ascension (J2000)')
        ax.set_ylabel('Declination (J2000)')
        #colors=['red','green','yellow','blue','black']
        #for i in range(len(apertures)):
        #    apertures[i].plot(color=colors[i], alpha=0.7) 
        for aperture in apertures:
            positions = aperture.positions
            plt.scatter(positions[:, 0], positions[:, 1], marker='x', color='red', s=10) 
        #an_ap.plot(color='red', alpha=0.7) 

    except (IOError, KeyError):
        fig = plt.figure()
        plt.imshow(image, cmap='gray', origin='lower', vmin=mean-2*std, vmax=mean+2*std)


    plt.show()


    flux_err_0=phot_table['aperture_sum_err_0']
    mag_err_0=1.09*flux_err_0/final_sum0

    flux_err_1=phot_table['aperture_sum_err_1']
    mag_err_1=1.09*flux_err_1/final_sum1

    flux_err_2=phot_table['aperture_sum_err_2']
    mag_err_2=1.09*flux_err_2/final_sum2

    flux_err_3=phot_table['aperture_sum_err_3']
    mag_err_3=1.09*flux_err_3/final_sum3

    flux_err_4=phot_table['aperture_sum_err_4']
    mag_err_4=1.09*flux_err_4/final_sum4 

    jd_list = [jd] * len(sources['xcentroid'])
    new_name=filename.replace('_cleaned.fits','_instrumental')
    with open(path+"{}.csv".format(new_name), 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['JD' ,'RA', 'DEC', 'Magnitude0', 'Magnitude_error0', 'Magnitude1', 'Magnitude_error1', 'Magnitude2', 'Magnitude_error2', 'Magnitude3', 'Magnitude_error3', 'Magnitude4', 'Magnitude_error4'])
        writer.writerows(zip(jd_list, ra,dec, mag_0, mag_err_0, mag_1, mag_err_1, mag_2, mag_err_2, mag_3, mag_err_3, mag_4, mag_err_4))
    
    print('Magnitude file generated')
    return path+"{}.csv".format(new_name)
    #eturn jd,ra,dec,mag_0, mag_err_0,
    #mag_1,mag_err_1,mag_2,mag_err_2 ,mag_3,mag_err_3,mag_4,mag_err_4


#def optimal_aperture(filename)

def calibrate_magnitudes(instrumental_magnitude_file):
    start = time.time()
    # Read RA and DEC from the text file
    star_coordinates = []
    with open(instrumental_magnitude_file) as file: #
        reader = csv.reader(file)
        next(reader)  # Skip the header row
        for row in reader:
            ra, dec = row[1], row[2]
            star_coordinates.append((float(ra), float(dec)))


    # Query Simbad for magnitudes of stars
    Simbad.reset_votable_fields()
    Simbad.add_votable_fields('flux(R)', 'flux_error(R)')

    ra_list = []
    dec_list = []
    magnitude_list = []
    magnitude_error_list = []
    print('Searching SIMBAD for co-ordinates')
    for ra, dec in star_coordinates:
        coords = SkyCoord(ra=ra, dec=dec, unit=u.deg, frame='icrs')
        result_table = Simbad.query_region(coords, radius=5 * u.arcsec)
        if result_table is not None:
            magnitude = result_table['FLUX_R'][0]
            magnitude_error = result_table['FLUX_ERROR_R'][0]
            ra_list.append(ra)
            dec_list.append(dec)
            magnitude_list.append(magnitude)
            magnitude_error_list.append(magnitude_error)
        else:
            ra_list.append(ra)
            dec_list.append(dec)
            magnitude_list.append('Notfound')
            magnitude_error_list.append('Notfound')
    # Save the retrieved magnitudes to a .csv file
    output_file = 'simbad_magnitudes3.csv'
    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['RA', 'DEC', 'Magnitude', 'Magnitude_error'])
        writer.writerows(zip(ra_list, dec_list, magnitude_list, magnitude_error_list))

    print(f"Star magnitudes saved to '{output_file}'.")
    file1 = pd.read_csv(instrumental_magnitude_file) #reference frame file
    file2 = pd.read_csv(output_file) #simbad code generated file

    jd=file1['JD']
    RA=file1['RA']
    Dec=file1['DEC']


    # Extract the desired columns
    x_values0 = file1['Magnitude0']  # Replace 'desired_column_name' with the actual column name from file1
    x_error0 = file1['Magnitude_error0']
    y_values = file2['Magnitude']  # Assuming the third column in file2 is named 'Column3'
    y_error = file2['Magnitude_error']  # Assuming the fourth column in file2 is named 'Column4'

    x_values1 = file1['Magnitude1']  # Replace 'desired_column_name' with the actual column name from file1
    x_error1 = file1['Magnitude_error1']
    x_values2 = file1['Magnitude2']  # Replace 'desired_column_name' with the actual column name from file1
    x_error2 = file1['Magnitude_error2']
    x_values3 = file1['Magnitude3']  # Replace 'desired_column_name' with the actual column name from file1
    x_error3 = file1['Magnitude_error3']
    x_values4 = file1['Magnitude4']  # Replace 'desired_column_name' with the actual column name from file1
    x_error4 = file1['Magnitude_error4']
    # Convert y_values and y_error to numeric, handling invalid or missing values
    y_values = pd.to_numeric(y_values, errors='coerce')
    y_error = pd.to_numeric(y_error, errors='coerce')

    # Filter out invalid or missing values from x_values, y_values, and y_error
    valid_indices = (~np.isnan(x_values0)) & (~np.isnan(y_values)) & (~np.isnan(y_error))
    ra_calib= RA[valid_indices]
    dec_calib= Dec[valid_indices]
    x_values0_new = x_values0[valid_indices]
    x_values1_new = x_values1[valid_indices]
    x_values2_new = x_values2[valid_indices]
    x_values3_new = x_values3[valid_indices]
    x_values4_new = x_values4[valid_indices]
    y_values = y_values[valid_indices]
    y_error = y_error[valid_indices]
    x_error0_new = x_error0[valid_indices]
    x_error1_new = x_error1[valid_indices]
    x_error2_new = x_error2[valid_indices]
    x_error3_new = x_error3[valid_indices]
    x_error4_new = x_error4[valid_indices]
   
    plt.plot(x_values2_new, y_values,'ko',markersize=4)
    plt.errorbar(x_values2_new, y_values,xerr=x_error2_new, yerr=y_error, color='black',fmt='o',capsize=3, label='Data')

    slope2, intercept2 = np.polyfit(x_values2_new, y_values, 1)
    best_fit_line2 = slope2 * x_values2_new + intercept2

    # Plot the best-fit line
    plt.plot(x_values2_new, best_fit_line2, color='darkgreen', lw=1,label='Best Fit Line')
    equation2 = f'y = {slope2:.2f}x + {intercept2:.2f}'
    plt.text(0.5, 0.85, equation2,color = 'black', ha='center', va='center', transform=plt.gca().transAxes)
    plt.xlabel('Instrumental Magnitude')  # Replace with the desired X axis label
    plt.ylabel('Simbad Magnitude')  # Replace with the desired Y axis label
    plt.title('Magnitude Calibration')  # Replace with the desired plot title
    plt.legend()
    plt.show()
    with open(path+"simbad_matched_magnitudes.csv", 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow([ 'RA','Dec','Instrumental_Mag', 'Instrumental_Mag_err', 'Simbad_Mag', 'Simbad_Mag_error'])
        writer.writerows(zip(ra_calib,dec_calib,x_values2_new, x_error2_new,y_values, y_error))



    calib_mag0=slope2 * x_values0 + intercept2
    calib_mag1=slope2 * x_values1 + intercept2
    calib_mag2=slope2 * x_values2 + intercept2
    calib_mag3=slope2 * x_values3 + intercept2
    calib_mag4=slope2 * x_values4 + intercept2
    calib_name=instrumental_magnitude_file.replace('_instrumental.csv','_calibrated')
    with open("{}.csv".format(calib_name), 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['JD' ,'RA', 'DEC', 'Magnitude0', 'Magnitude_error0', 'Magnitude1', 'Magnitude_error1', 'Magnitude2', 'Magnitude_error2', 'Magnitude3', 'Magnitude_error3', 'Magnitude4', 'Magnitude_error4'])
        writer.writerows(zip(jd, RA,Dec, calib_mag0, x_error0, calib_mag1,x_error1, calib_mag2, x_error2, calib_mag3, x_error3, calib_mag4, x_error4))
    
    print('Calibrated magnitude file generated')
    return "{}.csv".format(calib_name)


def plot_the_histogram(calibrated_magnitude_file):
    data = pd.read_csv(calibrated_magnitude_file)
    magnitude = data['Magnitude2']  # Replace 'Column3' with the actual column name from your CSV file

    # Plotting the histogram
    plt.hist(magnitude, bins=10,color='blue',edgecolor='white' )  # You can adjust the number of bins as per your preference
    #plt.axvline(x=17, color='r', linestyle='--', label='x=17')
    # Add labels and title to the plot
    plt.xlabel('Calibrated Magnitude')  # Replace with the desired label for the x-axis
    plt.ylabel('No. of Sources')  # Replace with the desired label for the y-axis
    plt.title('Histogram of sources obtained with aperture photometry')  # Replace with the desired title for the plot
    #plt.legend()
    # Display the histogram
    plt.show()

def main():
    start = time.time()
    #view_image(filename,4)
    #clean_the_images(path)
    photometry_list=do_photometry(path,filename,4,3)
    calib_magnitudes=calibrate_magnitudes(photometry_list)
    #plot_the_histogram(calib_magnitudes)
    end = time.time()
    print('Execution took: ',end - start,' seconds')




if __name__ == '__main__':
    main()
