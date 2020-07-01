from astropy.io import fits
from pylab import *
import glob
#import astroalign as aa
#from do_photometry import do_aperture_photometry

files=sorted(glob.glob('/MCG-10-16-052/band_wise/U_band/*U.fits'))
image_data=fits.open(files[0])
reference_image=image_data[0].data
'''
for i in range(3,len(files)):
    image_data=fits.open(files[i])
    source_image=image_data[0].data
    header=image_data[0].header
    image_aligned,footprint=aa.register(source_image,reference_image)

    aligned_file=files[i].replace('.fits','')
    fits.writeto(aligned_file+'_aligned'+'.fits',image_aligned,header,overwrite=True)

    print('No. %i done'%i)
'''
def view_image(name):
    data=fits.open(name)
    image=data[0].data
    head=data[0].header
    mean,std=np.mean(image),np.std(image)
    imshow(image[800:1100,800:1100],cmap='gray_r',origin='lower',vmin=mean-0.3*std,vmax=mean+0.3*std)
    plt.scatter(128, 185, s=400,edgecolor='red',linewidth=1, facecolor='none')
    #plt.text(600, 600,'B band image', fontsize=25)
    show()

for i in range(0,len(files)):
    name_of_file=files[i]
    view_image(name_of_file)
    #do_aperture_photometry(name_of_file,5)
