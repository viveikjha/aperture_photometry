import time
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

def find_the_sources():
    filename='data/j0947_i_7_sliced.fits'
    data1=fits.open(filename)
    image1=data1[0].data

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
        #
        #plt.plot(new_x,new_y)
        #
    '''
    tellme('All Done!')
    plt.show()
    plt.pause(1)

    '''
