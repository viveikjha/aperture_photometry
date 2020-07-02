import numpy as np
import ccdproc,os
from astropy import units as u
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.time import Time
from glob import glob
#from astroquery.astrometry_net import AstrometryNet
from astropy.wcs import WCS
import astroalign as aa


def clean_the_images(path,filename):
    #ast=AstrometryNet()
    #ast.api_key= 'iqmqwvazpvolmjmn'
    dir = path
    gain = 2 * u.electron / u.adu
    readnoise = 7.5 * u.electron
    ra=input('Enter the RA of the source:   ')
    dec=input('Enter the DEC of the source: ')
    '''
    wcs_header=ast.solve_from_image(path+filename)
    wcs=WCS(wcs_header)
    ran,decn=wcs.all_pix2world(1024,1024,0)
    print(ran,decn)
    '''
    file_name = os.path.join(dir,filename)
    image=ccdproc.CCDData.read(file_name,unit='adu')
    header=fits.getheader(file_name,0)

    time=header['DATE']
    t=Time(time,format='isot',scale='utc')
    print(t.jd,t.mjd)
    header.insert(15,('RA',ra))
    header.insert(16,('DEC',dec))


    a = sorted(glob(os.path.join(dir,'bias*.fits')))
    biaslist = []
    for i in range (0,len(a)):
        data= ccdproc.CCDData.read(a[i],unit='adu')
        #data = ccdproc.create_deviation(data, gain=gain, readnoise=readnoise)
        #data= data-(data.uncertainty.array)
        biaslist.append(data)
    combiner = ccdproc.Combiner(biaslist)
    masterbias = combiner.median_combine()
    masterbias.write('masterbias.fit', overwrite=True)
    mbias=ccdproc.CCDData.read('masterbias.fit',unit='adu')
    #masterbias.meta=image.meta
    print('master bias generated')
    print(np.mean(masterbias), np.median(masterbias))



    c=sorted(glob(os.path.join(dir,'flat*.fits')))
    flatlist = []
    for j in range(0,len(c)):
        flat=ccdproc.CCDData.read(c[j],unit='adu')
        #flat= ccdproc.create_deviation(flat, gain=gain, readnoise=readnoise)
        flat=ccdproc.subtract_bias(flat,masterbias)
        flatlist.append(flat)
    combiner = ccdproc.Combiner(flatlist)
    masterflat = combiner.median_combine()
    masterflat.write('masterflat.fits', overwrite=True)
    mflat=ccdproc.CCDData.read('masterflat.fits',unit='adu')
    print('master flat generated')
    print(np.mean(masterflat), np.median(masterflat))

    #masterflat.meta=image.meta


    bias_subtracted = ccdproc.subtract_bias(image, masterbias)
    flat_corrected = ccdproc.flat_correct(bias_subtracted, masterflat)
    cr_cleaned = ccdproc.cosmicray_lacosmic(flat_corrected,readnoise=7.5, sigclip=5)
    print('cosmic ray removed')



    fits.writeto(dir+'j_0947_i_1_clean.fits',cr_cleaned,header,overwrite=True)


    print('image cleaned')

# To align multiple images with respect to one image we use the astroalign package.

def align_the_images(path,filename,ref_image):
    nfiles=sorted(glob(path+filename))
    image_data=fits.open(path+ref_image)
    reference_image=image_data[0].data
    for i in range(len(nfiles)):
        image_data=fits.open(nfiles[i])
        source_image=image_data[0].data
        image_aligned,footprint=aa.register(source_image,reference_image)
        fits.writeto(path+'j0947_corrected_%i.fits'%i,image_aligned,overwrite=True)
        print('No. %i done'%i)


def time_to_jd(path,filename):
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
            fits.writeto(mod_file_1+'_sliced_'+'.fits',jd_up,header,overwrite=True)
            #print(files[i],t.jd,t.mjd,'single scan image')



        elif (check_header=='Kinetics'):
            exposure=header['EXPOSURE']
            print('kinetic mode image with no. of files:',files[i])

            name_of_file=files[i]
            mod_file=name_of_file.replace('.fits','')
            time=header['DATE']
            #print(time)
            t=Time(time,format='isot',scale='utc')
            tim=t.jd
            temp=int(nof[i])
            mod_jd=np.zeros(temp)
            exp_time=header['EXPOSURE']
            exp_time=exp_time/86400  # for the 'day' from seconds calculation.
            mod_jd[0]=tim
            for j in range(1,temp):
                mod_jd[j]=mod_jd[j-1]+exp_time

            for k in range(0,len(mod_jd)):
                sliced_image=image[k]
                time_jd=mod_jd[k]
                header.insert(15,('JD',time_jd))
                fits.writeto(mod_file+'_sliced_%g'%k+'.fits',sliced_image,header,overwrite=True)
                print(mod_file+'_sliced_%g'%k+'.fits has been written')
