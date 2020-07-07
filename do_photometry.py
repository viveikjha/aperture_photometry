''''
The modules required for the reduction of astronomical images and performing aperture photometry on the selected sources. All the modules except one use pure PYTHON packages, mostly derived from the astropy project. One of the modules still uses IRAF which we'll upgrade soon.

This program contains functions that perform:
1. BIAS correction.
2. Flat correction.
3. Cosmic Ray removal.
4. Finding the FWHM for a few sources, to determine the aperture radius.
5. Selecting the source and comparision stars for performing aperture photometry.
6. Conversion of time to Julian Dates (JD) or MJD.
7. Align the images with respect to a reference image (to make sure the co-ordinates are constant throughout the sample)
8. Calculate and subtract the global rms background.
9. Finds the sources in an image based on the DAOFIND algorithm. IRAFStarFinder, which uses IRAF algorithm can also be used.
10. Generates multiple aprtures on the locations selected by the user.
11. Generates annulus to estimate the local backgound and subtract it later.
12. Gets the total flux within the apertures and instrumental magnitudes.
13. Estimates the error on the flux and magnitudes assuming poisson noise.
14. Writes the output to ASCII file.
15. Plot the light curves using matplotlib.

written by:
Vivek, Kiran and Dimple at ARIES Nainital.

This was written keeping the 1.3m DFOT in mind and the customizations are needed for various telescopes/datasets if required.
'''
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


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
from photutils.background import Background2D, MeanBackground,SExtractorBackground
from photutils import find_peaks, CircularAperture, CircularAnnulus, aperture_photometry
from photutils.centroids import centroid_2dg
from astropy.stats import SigmaClip, mad_std
from photutils import Background2D, MedianBackground, DAOStarFinder
from photutils.utils import calc_total_error
from photutils.detection import findstars

def view_image(name):
    '''
    This module is meant for viewinga large number of images. The images
    can be seen as an animation. Will be helpful for large datasets as it doesn't
    require a lot of pointing and clicking.
    INPUT:
    name: The name of the image file.

    OUTPUT:

    the images seen as animation, with a gap of 1 second between successive images.
    '''

    data=fits.open(name)
    image=data[0].data
    head=data[0].header

    mean,std=np.mean(image),np.std(image)
    #image=image/np.median(image)
    im=plt.imshow(image,cmap='gray_r',origin='lower',vmin=mean-2*std,vmax=mean+2*std)
    #plt.scatter(124, 146, s=400,edgecolor='red',linewidth=1, facecolor='none')
    ims.append([im])
    plt.pause(1)
    plt.clf()




def clean_the_images(path,filename):
    '''
    This module is meant for cleaning the images. The tasks to be included are: bias correction,
    flat correction, trimming, overscan as well as the cosmic ray removal from the science cases.
    (For the time we are skipping the overscan and trimming part.

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

    ra=input('Enter the RA of the source:   ')
    dec=input('Enter the DEC of the source: ')


    bias_files = sorted(glob(os.path.join(dir,'bias*.fits')))
    biaslist = []
    for i in range (0,len(bias_files)):
        data= ccdproc.CCDData.read(bias_files[i],unit='adu')
        #data = ccdproc.create_deviation(data, gain=gain, readnoise=readnoise)
        #data= data-(data.uncertainty.array)
        biaslist.append(data)
    masterbias = ccdproc.combine(biaslist,method='average',sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
                             sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std)
    masterbias.write('masterbias.fits', overwrite=True)
    mbias=ccdproc.CCDData.read('masterbias.fits',unit='adu')
    print('Master bias generated')
    print(" Mean and median of the masterbias: ",np.mean(masterbias), np.median(masterbias))



    flat_files=sorted(glob(os.path.join(dir,'flat*.fits')))
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
    masterflat.write('masterflat.fits', overwrite=True)
    mflat=ccdproc.CCDData.read('masterflat.fits',unit='adu')
    print('Master flat generated')
    print(" Mean and median of the masterflat: ",np.mean(masterflat), np.median(masterflat))



    file_names = sorted(glob(os.path.join(dir,filename)))
    for i in range(0,len(file_names)):
        image=ccdproc.CCDData.read(file_names[i],unit='adu')
        header=fits.getheader(file_names[i],0)
        bias_subtracted = ccdproc.subtract_bias(image, masterbias)
        flat_corrected = ccdproc.flat_correct(bias_subtracted, masterflat)
        cr_cleaned = ccdproc.cosmicray_lacosmic(flat_corrected,readnoise=7.5, sigclip=5,satlevel=65535,niter=20,cleantype='meanmask',gain_apply=True)
        #print('Cosmic rays removed')
        clean_file=file_names[i].replace('.fits','')


        fits.writeto(clean_file+'_cleaned.fits',cr_cleaned,header,overwrite=True)
        print('Image no-%i has been cleaned'%i)


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

        print('No. %i done'%i)




def time_to_jd(path,filename):

    '''
    This function is used to take the observation epoch from the header and convert it to Julian dates. A header value 'JD' is written to the fits header. For the kinetic mode images, we take the time from the first frame and slice the images first and then add the intregration time to the frames to get the JD for each frame.

    INPUT:

    path: The directory where the images are kept (string)
    filename: The first few characters and the extension of the images (string). Example:
    j0946*fits, HD1104*fits etc.


    OUTPUT:

    1. The images sliced in the case of kinetic mode images.
    2. The Julian date updated to the header of the fits files.


    '''

    dir=path

    files=sorted(glob(os.path.join(dir,filename)))
    nof=np.zeros(len(files))
    for i in range(0,len(files)):
        data=fits.open(files[i])
        header=data[0].header
        image=data[0].data
        k=np.shape(image)
        nof[i]=k[0]

        check_header=header['ACQMODE']

        if (check_header=='Single Scan'):
            jd_up=image
            time=header['DATE']
            t=Time(time,format='isot',scale='utc')
            time_jd=t.jd
            header.insert(15,('JD',time_jd))
            files[i]
            mod_file_1=files[i].replace('.fits','')
            fits.writeto(mod_file_1+'_sliced'+'.fits',jd_up,header,overwrite=True)

            #print(files[i],t.jd,t.mjd,'single scan image')



        elif (check_header=='Kinetics'):

            print('kinetic mode image with no. of files:',nof[i])

            name_of_file=files[i]
            mod_file=name_of_file.replace('.fits','')
            time=header['DATE']
            t=Time(time,format='isot',scale='utc')
            t_jd=t.jd
            temp=int(nof[i])
            mod_jd=np.zeros(temp)
            exp_time=header['ACT']/86400  # for the 'day' from seconds calculation.
            mod_jd[0]=t_jd   # first frame time from the header and the next ones from the 'ACT'
            for j in range(1,temp):
                mod_jd[j]=mod_jd[j-1]+exp_time


            for k in range(0,len(mod_jd)):
                mod_header=header
                sliced_image=image[k]
                time_jd=mod_jd[k]
                mod_header.insert(15,('JD',time_jd))
                fits.writeto(mod_file+'_sliced_%g'%k+'.fits',sliced_image,mod_header,overwrite=True)
                print(mod_file+'_sliced_%g'%k+'.fits has been written')
                continue


        #os.system("mv *sliced* cleaned/")

def iraf_fwhm(filename):

    '''
    This function calculates the fwhm of the sources we'll select from the list of sources available to us. This uses a combination of IRAF and DS9 software at the moment, and we'll replace it with suitable python packages in time.

    INPUT:

    filename: The name of the file. Be sure to chose the best image among the group of images that you have. This image will also be the reference image for future alignment.


    OUTPUT:

    1.The mean  FWHM calculated from the sources selected by the user.
    2. The filename, which will be used as a reference image.



    '''
    print(' We need to get the FWHM using IRAF tasks.')
    print(' Opening DS9 for image display.')
    os.system("ds9 &")
    #filename=input('Enter the filename:')
    input('Please press enter here!!!')
    iraf.display(filename,1)
    print ('Please do the following:\n1. Press enter here and then click on any source in the DS9 window. \n2. Press (comma) in the middle of source(s) to get FWHM.\n3. Press q to quit. \n ')
    imx=iraf.imexam(Stdout=1)
    sum1=0
    if(imx[1].split()[10]=='INDEF'):
        fwhm=8.0
        print ('Fwhm:', 'INDEF', 'but taken as 8.0')
    else:
        for i in range(1,len(imx)):
            token=i
            sum1+=eval(imx[i].split()[10])
        fwhm=sum1/token
        print("Average FWHM in pixels:",fwhm, ' in arc seconds:',0.53*fwhm)
    return(fwhm,filename)





def source_list(filename):
    '''
    This function selects the sources based on the IRAF imexa task and dispayed through DS9. We point at the centre of the sources and the return is the x,y co-ordinates. Be careful to click in the centre. You can use zoom in function to get the centre. The apertures will be selected keeping the same locations as the centre.


    INPUT:
    filename: The name of the file. Be sure to chose the best image among the group of images that you have. This image will also be the reference image for future alignment.


    OUTPUT:

    The x,y co-ordinates selected by the user.

    '''
    print('*****************************')
    print('now selecting the sources')
    os.system("ds9 &")
    input("Please press enter")
    iraf.display(filename,1)
    print ('Please do the following:\n1. Press enter here and then click on any source in the DS9 window. \n2. Press (comma) in the middle of source(s) to get FWHM.\n3. Press q to quit. \n ')
    imx=iraf.imexam(Stdout=1)
    xval=np.zeros(len(imx))
    yval=np.zeros(len(imx))
    for i in range(1,len(imx)):
        xval[i]=eval(imx[i].split()[0])
        yval[i]=eval(imx[i].split()[1])
    with open('co_ordinates_list.txt', 'w') as f:
        for i in range(1,len(imx)):
            print(xval[i],'\t',yval[i],file=f)

    return(xval,yval)




def py_source_list(filename):
    '''
    This function selects the sources based on the python tasks task and dispayed through matplotlib. We point at the centre of the sources and the return is the x,y co-ordinates. Be careful to click in the centre. You can use zoom in function to get the centre. The apertures will be selected keeping the same locations as the centre.


    INPUT:
    filename: The name of the file. Be sure to chose the best image among the group of images that you have. This image will also be the reference image for future alignment.


    OUTPUT:

    The x,y co-ordinates selected by the user.
    '''

    #data1=fits.open(filename)
    #image1=data1[0].data

    def tellme(s):
        data=fits.open(filename)
        image=data[0].data
        mean,std=np.mean(image),np.std(image)
        print(s)
        plt.imshow(image,origin='lower',cmap='gray_r',vmin=mean-std,vmax=mean+std)
        plt.title(s, fontsize=16)
        #plt.draw()

    plt.clf()
    tellme('In this interactive window, we select sources')

    def input_source(val):

        val=int(val)
        print(val)
        return val

    vals=input('How many sources you want to select?')
    nsources=input_source(vals)

    plt.setp(plt.gca(), autoscale_on=True)


    plt.waitforbuttonpress()

    pts = []


    while len(pts) < nsources:
        tellme('Select 3 corners with mouse')
        pts = np.asarray(plt.ginput(nsources, timeout=-1))
        #print(pts[0],pts[1])

    #raise SystemExit
    plt.clf()
    for i in range(len(pts)):
        k=pts[i]
        x=k[0]
        y=k[1]
        new_x=np.arange(x-40,x+40,1).astype(int)
        new_y=np.arange(y-40,y+40,1).astype(int)
        mod_image=image1[new_x[0]:new_x[78],new_y[0]:new_y[78]]
        plt.plot(mod_image)
        plt.show()

def do_aperture_photometry(filename, fwhm,date):

    #fwhm,files=iraf_fwhm()

    #xpix,ypix=source_list(files)

    #ast=AstrometryNet()
    #ast.api_key= 'iqmqwvazpvolmjmn'
    '''
    choice=input("Enter the mode. Please use the keywords\n 'single' for single image, 'multiple'  for multiple images: \n\n")

    if (choice=='single'):
        print('single mode')
    elif (choice ==' multiple'):
        print(' multiple image mode')
    else:
        print('Enter valid choice!!!')
    '''
    data,header=fits.getdata(filename,header=True)
    #exposure=header['EXPOSURE']
    exposure=300
    #print('Exposure is',exposure)

    sigma_clip = SigmaClip(sigma=3, maxiters=10)
    bkg_estimator = SExtractorBackground()
    bkg = Background2D(data, (10,10), filter_size=(3, 3),sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    back=bkg.background # this is the background we need for the background subtraction.
    back2=np.median(bkg.background)
    print('median background is',back2)

    mask = data == 0
    unit = u.electron / u.s


    xdf_image = CCDData(data, unit=unit, meta=header, mask=mask)
    norm_image = ImageNormalize(vmin=1e-4, vmax=5e-2, stretch=LogStretch(), clip=False)
    xdf_image_clipped = np.clip(xdf_image, 1e-4, None)

    mean, median, std = sigma_clipped_stats(xdf_image.data, sigma=3.0, maxiters=20, mask=xdf_image.mask)

    print('Finding the sources')

    #daofind = DAOStarFinder(fwhm=fwhm, threshold=5*std) # 3 sigma above the background.
    #sources = daofind(data-back)

    #sources_findpeaks = find_peaks(xdf_image.data, mask=xdf_image.mask, threshold=30.*std, box_size=30, centroid_func=centroid_2dg)

    #print('We have found:',len(sources),' sources')
    #print(sources)

    #print(sources['xcentroid'], sources['ycentroid'],sources['fwhm'])
    #positions=sources['xcentroid'], sources['ycentroid']
    positions=np.genfromtxt('co_ordinates_list.txt',unpack=True,usecols=(0,1))
    print(positions)
    radii=[ fwhm,2*fwhm, 3*fwhm,4*fwhm,5*fwhm]

    #positions=(sources['xcentroid'], sources['ycentroid'])
    apertures = [CircularAperture(positions, r=r) for r in radii]

    an_ap = CircularAnnulus(positions, r_in=6*fwhm, r_out=6.2*fwhm)
    #apers = [apertures, annulus_apertures]


    #bkg_sigma=mad_std(data)
    effective_gain=exposure
    error=calc_total_error(data,back,effective_gain)


    #error=0.1*data
    phot_table = aperture_photometry(data-back, apertures,error=error)
    phot_table2=aperture_photometry(data-back,an_ap)


    bkg_mean = phot_table2['aperture_sum'] / an_ap.area
    bkg_sum = bkg_mean * an_ap.area


    final_sum0=phot_table['aperture_sum_0']-bkg_sum
    final_sum1=phot_table['aperture_sum_1']-bkg_sum
    final_sum2=phot_table['aperture_sum_2']-bkg_sum
    final_sum3=phot_table['aperture_sum_3']-bkg_sum
    final_sum4=phot_table['aperture_sum_4']-bkg_sum

    mag_back=-2.5*np.log10(bkg_mean/exposure)+22
    mag_0=-2.5*np.log10(final_sum0/exposure)+22
    mag_1=-2.5*np.log10(final_sum1/exposure)+22
    mag_2=-2.5*np.log10(final_sum2/exposure)+22
    mag_3=-2.5*np.log10(final_sum3/exposure)+22
    mag_4=-2.5*np.log10(final_sum4/exposure)+22

    print(mag_back,mag_0,mag_1,mag_2,mag_3,mag_4)


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

    fig=plt.figure()
    plt.imshow(data,cmap='gray',origin='lower',vmin=mean-4*std,vmax=mean+4*std)
    colors=['red','salmon','yellow','blue','cyan']
    for i in range(len(apertures)):
        apertures[i].plot(color=colors[i], alpha=0.7)

    an_ap.plot(color='green', alpha=0.7)
    plt.show()

    with open ('{}.dat'.format(date),'w') as r:
        for i in range (len(phot_table)):
            #print(final_sum0[i],final_sum1[i],final_sum2[i],final_sum3[i],final_sum4[i],final_sum5[i],file=r)
            print(mag_back[i],mag_0[i],mag_err_0[i],mag_1[i],mag_err_1[i],mag_2[i],mag_err_2[i],mag_3[i],mag_err_3[i],mag_4[i],mag_err_4[i],file=r)



    return mag_back
    '''



    mag=-2.5*np.log10(final_sum/30)+25

    flux=final_sum
    flux_err=phot_table['aperture_sum_err_0']
    mag_err=1.09*flux_err/flux

    x=[phot.value for phot in phot_table['xcenter']]
    y=[phot.value for phot in phot_table['ycenter']]

    #with open('result.dat', 'w') as f:
    #with open('out.txt', 'w') as f:
    for i in range(len(x)):
        print(x[i],y[i],'\t',mag[i],mag_err[i])


    outfile=' '
    for i in range (len(phot_table)):
        outfile+=x[i]+ " "+ y[i]+" "+ mag[i]+" " +mag_err[i]
        outfile+='\n'

    out=open('result.txt','w')
    out.write(outfile,overwrite=True)
    out.close()

    '''

    '''

    phot_table['xcenter']

    wcs_header=ast.solve_from_image('bfc.fits')
    wcs=WCS(wcs_header)
    xpix=[phot.value for phot in phot_table['xcenter']]
    ypix=[phot.value for phot in phot_table['ycenter']]
    ra=np.zeros(len(xpix))
    dec=np.zeros(len(ypix))
    for i in range(len(ra)):
        ra[i],dec[i]=wcs.all_pix2world(xpix[i],ypix[i],0)



    plt.imshow(data,cmap='gray',vmin=30,vmax=300)

    plt.show()

    #print(phot_table)

    '''
def plot_light_curve(filename,list):
#~~~~~~~~~~~~~~~~~~~~~~for plotting purposes~~~~~~~~#
    mp.rcParams['font.family']='serif'
    mp.rcParams['xtick.major.size']=10
    mp.rcParams['xtick.major.width']=2
    mp.rcParams['xtick.minor.size']=7
    mp.rcParams['xtick.minor.width']=2
    mp.rcParams['ytick.major.size']=10
    mp.rcParams['ytick.major.width']=2
    mp.rcParams['ytick.minor.size']=7
    mp.rcParams['ytick.minor.width']=2
    mp.rcParams['axes.linewidth']=1.5
    mp.rcParams['xtick.labelsize']=36
    mp.rcParams['ytick.labelsize']=36
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    filelist=list
    files=sorted(glob(filelist))
    julian_date=np.zeros(len(files))
    for i in range(0,len(files)):
        data=fits.open(files[i])
        header=data[0].header
        julian_date[i]=header['JD']

    # This is a dummy for flux
    julian_date-=2458870
    fluxd=np.zeros(len(files))
    for i in range(len(files)):
        fluxd[i]=random.randint(5,8)

    flux_errd=0.05*fluxd
    #################

    fig,ax=plt.subplots()
    plt.plot(julian_date,fluxd,'ko',markersize=25,label='I band light curve')
    plt.errorbar(julian_date,fluxd,yerr=flux_errd,capsize=5,fmt=' ')
    plt.ylabel("Magnitudes (instrumental)",fontsize=36)
    plt.xlabel("Julian Dates (2458870+)",fontsize=36)
    #plt.title("Deviation in lag values with respect to the variability",fontsize=36)
    ax.tick_params(axis="both",which='minor',direction="in")
    ax.tick_params(axis="both",which='major',direction="in")
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.minorticks_on()
    ax.legend(fontsize=30)
    plt.show()
