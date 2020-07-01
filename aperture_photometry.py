
from astropy.nddata import CCDData
from astropy.stats import sigma_clipped_stats, SigmaClip
from astropy.visualization import ImageNormalize, LogStretch
from matplotlib.ticker import LogLocator
from photutils.background import Background2D, MeanBackground
from photutils import find_peaks, CircularAperture, CircularAnnulus, aperture_photometry
from photutils.centroids import centroid_2dg
from astropy.stats import SigmaClip, mad_std
from photutils import Background2D, MedianBackground, DAOStarFinder
from photutils.utils import calc_total_error
from photutils.detection import findstar s


def do_aperture_photometry(filename):

    fwhm,filename=iraf_fwhm()

    xpix,ypix=source_list(filename)

    #ast=AstrometryNet()
    #ast.api_key= 'iqmqwvazpvolmjmn'


    data,header=fits.getdata(filename,header=True)
    exposure_time=header['EXPOSURE']

    sigma_clip = SigmaClip(sigma=3., maxiters=10) 
    bkg_estimator = MedianBackground() 
    bkg = Background2D(data, (10,10), filter_size=(3, 3),sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    back=bkg.background # this is the background we need for the background subtraction.
    back2=np.median(bkg.background)


    mask = data == 0
    unit = u.electron / u.s


    xdf_image = CCDData(data, unit=unit, meta=header, mask=mask)
    norm_image = ImageNormalize(vmin=1e-4, vmax=5e-2, stretch=LogStretch(), clip=False)
    xdf_image_clipped = np.clip(xdf_image, 1e-4, None)

    mean, median, std = sigma_clipped_stats(xdf_image.data, sigma=3.0, maxiters=20, mask=xdf_image.mask)

    print('Finding the sources')

    #daofind = DAOStarFinder(fwhm=fwhm, threshold=5*std) # 3 sigma above the background.
    #sources = daofind(data - median)

    #sources_findpeaks = find_peaks(xdf_image.data, mask=xdf_image.mask, threshold=30.*std, box_size=30, centroid_func=centroid_2dg)

    #print('We have found:',len(sources),' sources')
    #print(sources)

    #print(sources['xcentroid'], sources['ycentroid'],sources['fwhm'])
    #positions=sources['xcentroid'], sources['ycentroid']
    positions=np.genfromtxt('co_ordinates_list_1.txt',unpack=True,usecols=(0,1))
    #print(positions)
    radii=[ fwhm,2*fwhm, 4*fwhm, 6*fwhm]

    #positions=(sources['xcentroid'], sources['ycentroid'])
    apertures = [CircularAperture(positions, r=r) for r in radii] 

    an_ap = CircularAnnulus(positions, r_in=8*fwhm, r_out=10*fwhm)
    #apers = [apertures, annulus_apertures]


    #bkg_sigma=mad_std(data)
    effective_gain=exposure_time
    error=calc_total_error(data,back,effective_gain)


    #error=0.1*data
    phot_table = aperture_photometry(data, apertures,error=error)
    phot_table2=aperture_photometry(data,an_ap)


    bkg_mean = phot_table2['aperture_sum'] / an_ap.area()
    bkg_sum = bkg_mean * an_ap.area()


    final_sum0=phot_table['aperture_sum_0']-bkg_sum
    final_sum1=phot_table['aperture_sum_1']-bkg_sum
    final_sum2=phot_table['aperture_sum_2']-bkg_sum
    final_sum3=phot_table['aperture_sum_3']-bkg_sum



    mag_0=-2.5*np.log10(final_sum0/exposure_time)+25
    mag_1=-2.5*np.log10(final_sum1/exposure_time)+25
    mag_2=-2.5*np.log10(final_sum2/exposure_time)+25
    mag_3=-2.5*np.log10(final_sum3/exposure_time)+25

    fig=plt.figure()
    plt.imshow(data,cmap='gray',origin='lower',vmin=50,vmax=400)
    colors=['red','green','yellow','blue']
    for i in range(len(apertures)):
        apertures[i].plot(color=colors[i], alpha=0.7) 
    
    an_ap.plot(color='green', alpha=0.7) 
    plt.show()

    for i in range (len(phot_table)):
        print(mag_0[i],mag_1[i],mag_2[i],mag_3[i])



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








