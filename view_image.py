from astropy.io import fits
from pylab import *
import glob
import matplotlib.animation as animation
#import astroalign as aa
#from do_photometry import do_aperture_photometry

files=sorted(glob.glob('/mnt/eac99553-b108-449e-bf41-0716c975df8b/TRT-data/calibration/flat_files/MF-1-I*fits'))
print(len(files))
ims = []
image_data=fits.open(files[0])

def view_image(name):
    data=fits.open(name)
    image=data[0].data
    head=data[0].header
    mean,std=np.mean(image),np.std(image)
    image=image/np.median(image)
    im=imshow(image,cmap='gray_r',origin='lower')#,vmin=mean-0.2*std,vmax=mean+0.2*std)
    #plt.scatter(128, 185, s=400,edgecolor='red',linewidth=1, facecolor='none')
    ims.append([im])
    plt.pause(1)
    plt.clf()
    #plt.text(600, 600,'B band image', fontsize=25)

for i in range(0,len(files)):#len(files)):
    name_of_file=files[i]
    print(name_of_file)
    im=view_image(name_of_file)

plt.show()
    #do_aperture_photometry(name_of_file,5)
