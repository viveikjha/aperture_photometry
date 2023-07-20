#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib as mp
import math
import re
import matplotlib.pyplot as plt
import ccdproc,os,sys,time,random, csv
import astroalign as aa
from glob import glob
from astropy import units as u
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from astropy.nddata import CCDData
from astropy.stats import sigma_clipped_stats, SigmaClip
from astropy.visualization import ImageNormalize, LogStretch
from matplotlib.ticker import LogLocator
from astropy.stats import SigmaClip, mad_std
from photutils.background import Background2D, MeanBackground,SExtractorBackground
from photutils import find_peaks, CircularAperture, CircularAnnulus, aperture_photometry
from photutils.centroids import centroid_2dg
from photutils import Background2D, MedianBackground, DAOStarFinder
from photutils.utils import calc_total_error
from astropy.coordinates import SkyCoord
from astropy import coordinates as coord
from astropy import units as u
#from pyraf import iraf




# CREATING MAGNITUDE APERTURE FILE FOR single FRAME
files=sorted(glob(os.path.join('/home/aries/NGC2506_3/20211220_DFOT/cleaned_ngc/NGC*cleaned.fits')))
f = 0
filename=files[f]
start = time.time()
data_0,header_0=fits.getdata(files[f],header=True)
print(files[f])
source_0 = source(data_0 , header_0)
print('No. of sources: ',len(source_0))
print(source_0)
end = time.time()
print('Execution took:',end - start,' seconds')



radii = [5, 6, 7, 8, 9]
positions = [(source_0['xcentroid'][i], source_0['ycentroid'][i]) for i in range(len(source_0))]
apertures = [CircularAperture(positions, r=r) for r in radii] # Can we be flexible with variable apertures?
an_ap = CircularAnnulus(positions, r_in=14, r_out=15) # Is it needed if we get the 2D backgound  map?

start = time.time()
mag_0, mag_err_0, mag_1, mag_err_1, mag_2, mag_err_2, mag_3, mag_err_3, mag_4, mag_err_4 = magnitude(data_0, apertures, an_ap)
jd = header_0['JD']
jd_list = [jd] * len(source_0['xcentroid'])
h = fits.open(files[f])
wcs = WCS(h[0].header)
ra , dec = wcs.all_pix2world(source_0['xcentroid'],source_0['ycentroid'],0)
print(ra, dec)
#data = np.array([jd_list, ra,dec, mag_0, mag_err_0, mag_1, mag_err_1, mag_2, mag_err_2, mag_3, mag_err_3, mag_4, mag_err_4])
#transposed_data = np.transpose(data)

with open(f"/home/aries/NGC2506_3/20211220_DFOT/zipped_data_{f}.csv", 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['JD', 'RA', 'DEC', 'Magnitude0', 'Magnitude_error0', 'Magnitude1', 'Magnitude_error1', 'Magnitude2', 'Magnitude_error2', 'Magnitude3', 'Magnitude_error3', 'Magnitude4', 'Magnitude_error4'])
    writer.writerows(zip(jd_list, ra,dec, mag_0, mag_err_0, mag_1, mag_err_1, mag_2, mag_err_2, mag_3, mag_err_3, mag_4, mag_err_4))

end = time.time()
print('Execution took:', end - start, 'seconds')


# In[20]:


# PLOTING TO CHOSE OPTIMUM APERTURE SIZE
x=radii
start = time.time()
plt.figure(figsize = (8,8))
for i in range(0,350):
    y=[]
    y.append(mag_0[i])
    y.append(mag_1[i])
    y.append(mag_2[i])
    y.append(mag_3[i])
    y.append(mag_4[i])  
    plt.plot(x,y, 'r--',lw=0.5)
    plt.xlabel('Aperture diameter (pixels)', fontsize=15)
    plt.ylabel('Largest aperture mag - aperture mag', fontsize = 15)
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
plt.show()
    
end = time.time()
print('Execution took: ',end - start,' seconds')


# In[22]:


# CREATING MAGNITUDE APERTURE FILE FOR ALL FRAMES

files=sorted(glob(os.path.join('/home/aries/NGC2506_3/20211220_DFOT/cleaned_ngc/NGC*cleaned.fits')))
radii = [5,6,7,8,9]
start = time.time()
for f in range(len(files)):
    filename=files[f]
    data_0,header_0=fits.getdata(files[f],header=True)
    print(files[f])
    source_0 = source(data_0 , header_0)
    print('No. of sources:',len(source_0))

    positions = [(source_0['xcentroid'][i], source_0['ycentroid'][i]) for i in range(len(source_0))]
    apertures = [CircularAperture(positions, r=r) for r in radii]
    an_ap = CircularAnnulus(positions, r_in=14, r_out=15)

    start = time.time()
    mag_0, mag_err_0, mag_1, mag_err_1, mag_2, mag_err_2, mag_3, mag_err_3, mag_4, mag_err_4 = magnitude(data_0, apertures, an_ap)
    jd = header_0['JD']
    jd_list = [jd] * len(source_0['xcentroid'])
    h = fits.open(files[f])
    wcs = WCS(h[0].header)
    ra , dec = wcs.all_pix2world(source_0['xcentroid'],source_0['ycentroid'],0)
    #print(ra, dec)
    #data = np.array([jd_list, ra,dec, mag_0, mag_err_0, mag_1, mag_err_1, mag_2, mag_err_2, mag_3, mag_err_3, mag_4, mag_err_4])
    #transposed_data = np.transpose(data)

    with open(f"/home/aries/NGC2506_3/20211220_DFOT/zipped_five_mag_{f}.csv", 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['JD' ,'RA', 'DEC', 'Magnitude0', 'Magnitude_error0', 'Magnitude1', 'Magnitude_error1', 'Magnitude2', 'Magnitude_error2', 'Magnitude3', 'Magnitude_error3', 'Magnitude4', 'Magnitude_error4'])
        writer.writerows(zip(jd_list, ra,dec, mag_0, mag_err_0, mag_1, mag_err_1, mag_2, mag_err_2, mag_3, mag_err_3, mag_4, mag_err_4))
    
    print('file',f,'generated')
end = time.time()
print('Execution took: ',end - start,' seconds')


# In[21]:


#SIMBAD CODE

import astropy.units as u
import csv
from astropy.coordinates import SkyCoord
from astroquery.simbad import Simbad
start = time.time()
# Read RA and DEC from the text file
star_coordinates = []
with open('zipped_data_0.csv', 'r') as file: #
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
end = time.time()
print('Execution took: ',end - start,' seconds')


# NOTE: RUNNING SIMBAD CODE IS LUCK, SOMETIMES IT DOESN'T WORK AND MANYTIMES IT WORKS BUT BELOW CODE DON'T GIVE LINEAR PLOT! BE PATIENT...!

# In[ ]:





# In[20]:


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read the CSV files
file1 = pd.read_csv('zipped_data_0.csv') #reference frame file
file2 = pd.read_csv('simbad_magnitudes.csv') #simbad code generated file

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
x_values0 = x_values0[valid_indices]
x_values1 = x_values1[valid_indices]
x_values2 = x_values2[valid_indices]
x_values3 = x_values3[valid_indices]
x_values4 = x_values4[valid_indices]
y_values = y_values[valid_indices]
y_error = y_error[valid_indices]
x_error0 = x_error0[valid_indices]
x_error1 = x_error1[valid_indices]
x_error2 = x_error2[valid_indices]
x_error3 = x_error3[valid_indices]
x_error4 = x_error4[valid_indices]
# Plotting
# plt.errorbar(x_values0, y_values,xerr=x_error0, yerr=y_error, fmt='.', )
# plt.errorbar(x_values1, y_values,xerr=x_error1, yerr=y_error, fmt='.', )
plt.errorbar(x_values2, y_values,xerr=x_error2, yerr=y_error, color='blue',fmt='.', label='Data')
# plt.errorbar(x_values3, y_values,xerr=x_error3, yerr=y_error, fmt='.', )
# plt.errorbar(x_values4, y_values,xerr=x_error4, yerr=y_error, fmt='.',)

# Calculate the best-fit line parameters
# slope0, intercept0 = np.polyfit(x_values0, y_values, 1)
# best_fit_line0 = slope0 * x_values0 + intercept0

# Plot the best-fit line
# plt.plot(x_values0, best_fit_line0, color='r', label='Best Fit Line 0')
# equation0 = f'y = {slope0:.2f}x + {intercept0:.2f}'

# slope1, intercept1 = np.polyfit(x_values1, y_values, 1)
# best_fit_line1 = slope1 * x_values1 + intercept1

# # Plot the best-fit line
# plt.plot(x_values1, best_fit_line1, color='g', label='Best Fit Line 1')
# equation1 = f'y = {slope1:.2f}x + {intercept1:.2f}'

slope2, intercept2 = np.polyfit(x_values2, y_values, 1)
best_fit_line2 = slope2 * x_values2 + intercept2

# Plot the best-fit line
plt.plot(x_values2, best_fit_line2, color='blue', label='Best Fit Line 2')
equation2 = f'y = {slope2:.2f}x + {intercept2:.2f}'

# slope3, intercept3 = np.polyfit(x_values3, y_values, 1)
# best_fit_line3 = slope3 * x_values3 + intercept3

# # Plot the best-fit line
# plt.plot(x_values3, best_fit_line3, color='black', label='Best Fit Line 3')
# equation3 = f'y = {slope3:.2f}x + {intercept3:.2f}'

# slope4, intercept4 = np.polyfit(x_values4, y_values, 1)
# best_fit_line4 = slope4 * x_values4 + intercept4

# # Plot the best-fit line
# plt.plot(x_values4, best_fit_line4, color='purple', label='Best Fit Line 4')
# equation4 = f'y = {slope4:.2f}x + {intercept4:.2f}'


# Add equation annotation to the plot
# plt.text(0.5, 0.95, equation0, color = 'r',ha='center', va='center', transform=plt.gca().transAxes)
# plt.text(0.5, 0.9, equation1,color = 'g', ha='center', va='center', transform=plt.gca().transAxes)
plt.text(0.5, 0.85, equation2,color = 'blue', ha='center', va='center', transform=plt.gca().transAxes)
# plt.text(0.5, 0.8, equation3, color = 'black',ha='center', va='center', transform=plt.gca().transAxes)
# plt.text(0.5, 0.75, equation4,color = 'purple', ha='center', va='center', transform=plt.gca().transAxes)

plt.xlabel('Magnitude from code')  # Replace with the desired X axis label
plt.ylabel('Magnitude from simbad')  # Replace with the desired Y axis label
plt.title('Calibration of magnitude')  # Replace with the desired plot title
plt.legend()
plt.show()

mag=y_values - x_values2
plt.hist(mag)
plt.xlabel('Calibration Magnitude')
plt.title('Histogram to correctly pick calibration magnitude')
plt.show()


# In[39]:


#Calibrating reference frame (single)

files=sorted(glob(os.path.join('/home/aries/NGC2506_2/20210215_DFOT/try/NGC*cleaned.fits')))
f = 0
filename=files[f]
start = time.time()
data_0,header_0=fits.getdata(files[f],header=True)
print(files[f])
source_0 = source(data_0 , header_0)
print('No. of sources: ',len(source_0))
end = time.time()
print('Execution took: ',end - start,' seconds')

intercept= 22.10

radii = [7]
positions = [(source_0['xcentroid'][i], source_0['ycentroid'][i]) for i in range(len(source_0))]
apertures = [CircularAperture(positions, r=r) for r in radii]
an_ap = CircularAnnulus(positions, r_in=14, r_out=15)

start = time.time()
mag_0, mag_err_0 = magnitude(data_0, apertures, an_ap)
calibrated_mag = mag_0 + intercept
jd = header_0['JD']
jd_list = [jd] * len(source_0['xcentroid'])
h = fits.open(files[f])
wcs = WCS(h[0].header)
ra , dec = wcs.all_pix2world(source_0['xcentroid'],source_0['ycentroid'],0)
print(ra, dec)
data = np.array([jd_list, ra,dec, mag_0, mag_err_0])
transposed_data = np.transpose(data)

with open(f"Calibrated_zipped_data_{f}.csv", 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['JD', 'RA', 'DEC', 'Calibrated_Magnitude', 'Magnitude_error0'])
    writer.writerows(zip(jd_list, ra,dec, calibrated_mag, mag_err_0))

end = time.time()
print('Execution took:', end - start, 'seconds')


# In[45]:


#PLOTING HISTOGRAM ALONG WITH BRIGHT SOURCE BOUNDARY
import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV file
data = pd.read_csv('Calibrated_zipped_data_0.csv')

# Extract the data from the 3rd column
column_3_data = data['Calibrated_Magnitude']  # Replace 'Column3' with the actual column name from your CSV file

# Plotting the histogram
plt.hist(column_3_data, bins=10)  # You can adjust the number of bins as per your preference
plt.axvline(x=17, color='r', linestyle='--', label='x=17')
# Add labels and title to the plot
plt.xlabel('Values')  # Replace with the desired label for the x-axis
plt.ylabel('Frequency')  # Replace with the desired label for the y-axis
plt.title('Histogram of CALIBRATED magnitude (ref) Data')  # Replace with the desired title for the plot
plt.legend()
# Display the histogram
plt.show()


# In[48]:


# CREATING MAGNITUDE APERTURE FILE FOR ALL FRAMES

files=sorted(glob(os.path.join('/home/aries/NGC2506_2/20210215_DFOT/try/NGC*cleaned.fits')))
radii = [7]
intercept= 22.10
start = time.time()
for f in range(1,len(files)):
    filename=files[f]
    data_0,header_0=fits.getdata(files[f],header=True)
    print(files[f])
    source_0 = source(data_0 , header_0)
    print('No. of sources: ',len(source_0))

    positions = [(source_0['xcentroid'][i], source_0['ycentroid'][i]) for i in range(len(source_0))]
    apertures = [CircularAperture(positions, r=r) for r in radii]
    an_ap = CircularAnnulus(positions, r_in=14, r_out=15)

    mag_0, mag_err_0 = magnitude(data_0, apertures, an_ap)
    calibrated_mag = mag_0 + intercept
    jd = header_0['JD']
    jd_list = [jd]*len(source_0['xcentroid'])
    h = fits.open(files[f])
    wcs = WCS(h[0].header)
    ra , dec = wcs.all_pix2world(source_0['xcentroid'],source_0['ycentroid'],0)
    #data = np.array([jd_list, ra,dec, mag_0, mag_err_0, mag_1, mag_err_1, mag_2, mag_err_2, mag_3, mag_err_3, mag_4, mag_err_4])
    #transposed_data = np.transpose(data)

    with open(f"Calibrated_zipped_data_{f}.csv", 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['JD' ,'RA', 'DEC', 'Calibrated_Magnitude', 'Magnitude_error0'])
        writer.writerows(zip(jd_list, ra,dec, calibrated_mag, mag_err_0))
end = time.time()
print('Execution took: ',end - start,' seconds')


# In[16]:


#Sorting calibrated magnitude files  wrt calibrated magnitude
import pandas as pd
f=0
csv_file = f'/home/aries/NGC2506_2/20210215_DFOT/try/calibrated_zipped_data/Calibrated_zipped_data_{f}.csv'
# Read the CSV file
data = pd.read_csv(csv_file)

# Sort the data based on column 4
sorted_data = data.sort_values(by='Calibrated_Magnitude')

# Save the sorted data to a new CSV file
sorted_data.to_csv( f"Calibrated_zipped_data_{f}.csv", index=False) 

print("CSV file sorted successfully.")


# ### Ploting light curves

# In[3]:


#Ploting light curve for single star

import pandas as pd
import matplotlib.pyplot as plt

#f'/home/aries/NGC2506_2/20210215_DFOT/try/calibrated_zipped_data/Calibrated_zipped_data_{i}.csv'
start = time.time()
# Read the reference file
reference_file = pd.read_csv('/home/aries/NGC2506_2/20210215_DFOT/try/calibrated_zipped_data/Calibrated_zipped_data_0.csv')

# Get the RA and DEC values from the reference file
reference_ra = reference_file['RA'].iloc[0]
reference_dec = reference_file['DEC'].iloc[0]

# Initialize lists to store the matching target files' JD and magnitude values
matching_jd = []
matching_magnitude = []

# Iterate over the target files
for i in range(1, 375):  # Assuming your target files are named as 'target_2.csv', 'target_3.csv', ..., 'target_375.csv'
    target_file = pd.read_csv(f'/home/aries/NGC2506_2/20210215_DFOT/try/calibrated_zipped_data/Calibrated_zipped_data_{i}.csv')
    
    # Check if any RA and DEC values in the target file match with the reference values within the tolerance
    matches = ((target_file['RA'] - reference_ra).abs() <= 0.0001) & ((target_file['DEC'] - reference_dec).abs() <= 0.0001)
    
    # If there is a match, append the JD and magnitude values to the respective lists
    if matches.any():
        matching_jd.append(target_file.loc[matches, 'JD'].iloc[0])
        matching_magnitude.append(target_file.loc[matches, 'Calibrated_Magnitude'].iloc[0])

# Plot the magnitudes of the matching target files against their JD values
plt.plot(matching_jd, matching_magnitude, 'o')
plt.xlabel('JD')
plt.ylabel('Magnitude')
plt.title('Light curve for 1st star')
plt.show()

end = time.time()
print('Execution took: ',end - start,' seconds')


# In[7]:


#Ploting  light cuves for desired number of stars

import pandas as pd
import matplotlib.pyplot as plt

# Read the reference file
reference_file = pd.read_csv('/home/aries/NGC2506_2/20210215_DFOT/try/calibrated_zipped_data/Calibrated_zipped_data_0.csv')
start = time.time()
start_index = 0  # Index of the first row to consider
end_index = 10  # Index of the last row to consider
# Iterate over the RA and DEC values in the reference file
for index, row in reference_file.iloc[start_index:end_index].iterrows():
    reference_ra = row['RA']
    reference_dec = row['DEC']

    # Initialize lists to store the matching target files' JD and magnitude values
    matching_jd = []
    matching_magnitude = []

    # Iterate over the target files
    for i in range(1, 375):  # Assuming your target files are named as 'target_2.csv', 'target_3.csv', ..., 'target_375.csv'
        target_file = pd.read_csv(f'/home/aries/NGC2506_2/20210215_DFOT/try/calibrated_zipped_data/Calibrated_zipped_data_{i}.csv')

        # Check if any RA and DEC values in the target file match with the reference values within the tolerance
        matches = (
            (target_file['RA'] - reference_ra).abs() <= 0.0001
        ) & (
            (target_file['DEC'] - reference_dec).abs() <= 0.0001
        )

        # If there is a match, append the JD and magnitude values to the respective lists
        if matches.any():
            matching_jd.append(target_file.loc[matches, 'JD'].iloc[0])
            matching_magnitude.append(target_file.loc[matches, 'Calibrated_Magnitude'].iloc[0])

    # Plot the magnitudes of the matching target files against their JD values
    plt.plot(matching_jd, matching_magnitude, 'o')
    plt.xlabel('JD')
    plt.ylabel('Magnitude')
    plt.title(f'Magnitudes of star with RA={reference_ra} and DEC={reference_dec}')
    plt.show()
end = time.time()
print('Execution took: ',end - start,' seconds')


# In[5]:


#Ploting light curves

import pandas as pd
import matplotlib.pyplot as plt


start = time.time()
#f'/home/aries/NGC2506_2/20210215_DFOT/try/calibrated_zipped_data/Calibrated_zipped_data_{i}.csv'

# Read the reference file
reference_file = pd.read_csv('/home/aries/NGC2506_2/20210215_DFOT/try/calibrated_zipped_data/Calibrated_zipped_data_0.csv')
for n in range(0,10):
# Get the RA and DEC values from the reference file
    reference_ra = reference_file['RA'].iloc[n]
    reference_dec = reference_file['DEC'].iloc[n]

# Initialize lists to store the matching target files' JD and magnitude values
    matching_jd = []
    matching_magnitude = []

# Iterate over the target files
    for i in range(1, 375):  # Assuming your target files are named as 'target_2.csv', 'target_3.csv', ..., 'target_375.csv'
        target_file = pd.read_csv(f'/home/aries/NGC2506_2/20210215_DFOT/try/calibrated_zipped_data/Calibrated_zipped_data_{i}.csv')
    
    # Check if any RA and DEC values in the target file match with the reference values within the tolerance
        matches = ((target_file['RA'] - reference_ra).abs() <= 0.0001) & ((target_file['DEC'] - reference_dec).abs() <= 0.0001)
    
    # If there is a match, append the JD and magnitude values to the respective lists
        if matches.any():
            matching_jd.append(target_file.loc[matches, 'JD'].iloc[0])
            matching_magnitude.append(target_file.loc[matches, 'Calibrated_Magnitude'].iloc[0])
    
    mag_mean = np.mean(matching_magnitude)
    plt.axhline(y=mag_mean, color='r', linestyle='--', label='mean')
# Plot the magnitudes of the matching target files against their JD values
    plt.plot(matching_jd, matching_magnitude, 'o')
    plt.xlabel('JD')
    plt.ylabel('Magnitude')
    plt.title(f'Light curve for {n} star' )
    plt.show()

end = time.time()
print('Execution took: ',end - start,' seconds')


# In[ ]:


#PLOTING LIGHT CURVE AND HISTOGRAM OF DEVIATION FROM MEAN - FOR ANY SINGLE STAR

import random

csv_files = [f for f in os.listdir(directory_path) if f.startswith('Mag_file_star_')]
csv_files.sort(key=lambda x: int(x.split('_')[3].split('.')[0]))

random_number = random.randint(0, 375)
df = pd.read_csv(csv_files[80])
print(random_number)
magni_tude = [df.iloc[:,1]]
print(np.mean(magni_tude))
#print(magni_tude)
diff1=[]
for i in range(374):
    diff =(df.iloc[i, 1] - np.mean(magni_tude))
    diff1.append(diff)
    
    
plt.plot(df.iloc[:, 0], df.iloc[:, 1], 'o')
  # Adjust the width and height as needed
plt.xlabel('JD')
plt.ylabel('Magnitude')
plt.title('Light curve')
plt.axhline(y=np.mean(magni_tude), color='r', linestyle='--', label='mean')
plt.legend()
plt.show()


plt.hist(diff1, bins=50)
plt.xlabel('values')
plt.ylabel('frequency')
plt.title('Deviation from mean')
#plt.axvline(x=np.mean(magni_tude), color='r', linestyle='--',label='mean')

plt.show()


# ### Transpose of data: now as many files as stars (not frames)

# In[95]:


import pandas as pd
import numpy as np

# Read the reference file
reference_file = pd.read_csv('/home/aries/NGC2506_2/20210215_DFOT/try/calibrated_zipped_data2/Calibrated_zipped_sorted_0.csv')

# Filter reference file by 'Calibrated_magnitude' < 17
reference_file = reference_file[reference_file['Calibrated_Magnitude'] < 17]
start = time.time()
# Iterate through each row in the reference file
start_index = 0  # Index of the first row to consider
end_index = 2000  # Index of the last row to consider
# Iterate over the RA and DEC values in the reference file
for index, row in reference_file.iloc[start_index:end_index].iterrows():
    reference_ra = row['RA']
    reference_dec = row['DEC']
    
    # Create an empty DataFrame to store the results for each matching file
    matching_results = pd.DataFrame(columns=['JD', 'Calibrated_Magnitude', 'Magnitude_error0'])
    
    # Iterate through each file (excluding the reference file)
    for i in range(1, 375):
        filename = f'/home/aries/NGC2506_2/20210215_DFOT/try/calibrated_zipped_data2/Calibrated_zipped_sorted_{i}.csv'  # Assuming the file names are in the format '1.csv', '2.csv', ...
        file_data = pd.read_csv(filename)
        
        # Filter the file data based on matching 'RA' and 'DEC' values with the reference file
        matching_data = file_data[(np.isclose(file_data['RA'], reference_ra, atol=0.00001)) &
                                  (np.isclose(file_data['DEC'], reference_dec, atol=0.00001))]
        
        # Append the matching data to the results DataFrame
        matching_results = pd.concat([matching_results, matching_data[['JD', 'Calibrated_Magnitude', 'Magnitude_error0']]], ignore_index=True)
    
    # Save the matching results to a CSV file
    matching_results.to_csv(f'Mag_file_star_{index}.csv', index=False)
    
end = time.time()
print('Execution took: ',end - start,' seconds')


# ### Minimum difference using transposed data

# In[ ]:


#TO CHECK IF ALL TRANSPOSED STAR FILES HAVE SAME jd COLUMN: hence to DIRECTLY PLOT MAG/ MAG_DIFF WITH jd
import csv
import os

def check_first_column(csv_files):
    first_column_values = set()
    for file in csv_files:
        with open(file, 'r') as csv_fil:
            reader = csv.reader(csv_fil)
            for row in reader:
                first_column_values.add(row[0])
                break  # Only check the first row

    return len(first_column_values) == 1

# Directory path containing the CSV files
directory_path = '/home/aries/NGC2506_2/20210215_DFOT/try/sortjd_cal_star/'

# Get all CSV files in the directory
#csv_files = [f for f in os.listdir(directory_path) if f.startswith('Mag_') and f.endswith('.csv')]
csv_files =[]
for filename in os.listdir(directory_path):
    if filename.startswith('Mag_file_star') and filename.endswith('.csv'):
        file_path = os.path.join(directory_path, filename)
        csv_files.append(file_path)
# Check if all CSV files have the same first column
if check_first_column(csv_files):
    print("All CSV files have the same first column.")
else:
    print("CSV files have different values in the first column.")


# In[ ]:


# sorting calibrated_zipped_data wrt any column

import pandas as pd
import os

# Directory containing the CSV files
directory = '/home/aries/NGC2506_2/20210215_DFOT/try/'
output_directory =  '/home/aries/NGC2506_2/20210215_DFOT/try/sortjd_cal_star'
os.makedirs(output_directory, exist_ok=True)
# Iterate over each file in the directory
for filename in os.listdir(directory):
    if filename.startswith('Mag_file_star') and filename.endswith('.csv'):
        file_path = os.path.join(directory, filename)
        print('here')
        # Read the CSV file
        df = pd.read_csv(file_path)

        # Sort the DataFrame based on the first column
        sorted_df = df.sort_values(by=df.columns[0])
        
        # Generate the new file name
        new_filename = os.path.splitext(filename)[0] + '_sorted.csv'
        new_file_path = os.path.join(output_directory, new_filename)

        # Write the sorted DataFrame to the new file
        sorted_df.to_csv(new_file_path, index=False)


# In[8]:


#TO FILTER STAR THAT LIE IN MAG RANGE (+-1)

import pandas as pd
import glob

# Step 1: Read all CSV files
file_list =glob.glob('/home/aries/NGC2506_2/20210215_DFOT/try/Mag_file_star_*.csv')
file_list = sorted(file_list, key=lambda x: int(x.split('_')[5].split('.')[0]))
df_list = []
start = time.time()

for file_name in file_list:
    df = pd.read_csv(file_name)
    df_list.append(df)

selected_file = file_list[30] # can change this as per center of bin

# Step 3: Read selected file to get the first value in Calibrated_Magnitude column (x)
selected_df = pd.read_csv(selected_file)
x = selected_df['Calibrated_Magnitude'].iloc[0]
print(x)
# Step 4: Filter files where all values in Calibrated_Magnitude lie within range (x-1, x+1)
filtered_files = []

for file_name in file_list:
    df = pd.read_csv(file_name)
    min_val = df['Calibrated_Magnitude'].min()
    max_val = df['Calibrated_Magnitude'].max()

    if min_val >= (x - 1) and max_val <= (x + 1):
        filtered_files.append(file_name)

print(filtered_files)


# In[5]:


#TO Find the minimum difference between Calibrated_Magnitude of two files
min_diff = float('inf')
min_diff_files = []

for i in range(len(filtered_files)):
    for j in range(i + 1, len(filtered_files)):
        if i!=j:
            file1 = filtered_files[i]
            file2 = filtered_files[j]
            df1 = pd.read_csv(file1)
            df2 = pd.read_csv(file2)
            if len(df1)==len(df2)==374:
                #print(file2)
                diff = abs(df1['Calibrated_Magnitude'].values[:len(df2)] - df2['Calibrated_Magnitude'].values)
                if diff.min() < min_diff:
                    min_diff = diff.min()
                    min_diff_files = [file1, file2]
                    
print(min_diff_files)
# Step 6: Print the filenames and the minimum difference
for file_name in min_diff_files:
    df = pd.read_csv(file_name)
    print("File Name:", file_name)
    print("Calibrated Magnitude:", df['Calibrated_Magnitude'])
    print()

print("Minimum Difference:", min_diff)


# ### PLOTING DIFFERENTIAL LIGHT CURVES (within +-1 range)

# In[25]:


#PLOTING LIGHT CURVES OF STARS WITH MAG RANGE (+-1)
import pandas as pd
import matplotlib.pyplot as plt


start = time.time()

# Read the first CSV file
df1 = pd.read_csv('/home/aries/NGC2506_2/20210215_DFOT/try/Mag_file_star_67.csv')
# Read the second CSV file
df2 = pd.read_csv('/home/aries/NGC2506_2/20210215_DFOT/try/Mag_file_star_70.csv')
for i in range(8,91):
    
    df3 = pd.read_csv(f'/home/aries/NGC2506_2/20210215_DFOT/try/Mag_file_star_{i}.csv')
# Extract the first and second columns from both dataframes
    x = df1.iloc[:, 0]  # 1st column of file1 (assuming index starts from 0) 
    if len(df1)==len(df2)==len(df3):
        diff1 = df3.iloc[:, 1] - df1.iloc[:, 1]  # Difference of the 2nd columns
        diff2 = df3.iloc[:, 1] - df2.iloc[:, 1]
        diff3 = df2.iloc[:, 1] - df1.iloc[:, 1]
# Plot the scatter plot
        
        #plt.axhline(y=np.mean(diff3), color='r', linestyle='--', label='mean of C1-C2')
        plt.plot(x, diff1, '.', label='T-C1')
        plt.ylabel('T-C1 ')

        plt.plot(x, diff2,'.',label='T-C2')
        plt.ylabel('T-C2 ')

        plt.plot(x, diff3 + diff1,'.',label='C1-C2')
        plt.ylabel('C1-C2 ')
# Add labels and title

        plt.xlabel('JD')
        plt.title(f'Difference of Magnitude {i}')
        plt.legend()
# Display the plot
        plt.show()
    
end = time.time()
print('Execution took: ',end - start,' seconds')


# In[ ]:




