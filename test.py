import os,sys
import numpy as np
from astropy.io import fits,ascii
from pyraf import iraf
from iraf import noao,imred,specred
from stsci import tools   			    
# The cosmisc.py is a python module used for the cosmic ray correction. I obtained it from: http://obswww.unige.ch/~tewes/cosmics_dot_py/

'''
print("***********************************************************")
print("Starting spectroscopic reduction of your data!!!!.")
print("***********************************************************")
dpath= input('Please enter the path to the directory where your raw data files are=====>')
os.chdir(dpath)


#For Bias Correction
os.system("ls *Bias* > bias.in")
iraf.imcombine('@bias.in', 'mbias', combine="median",rdnoise=7.5,gain=2.0)
os.system("ls FeAr* J105829* halogen* GD248* > subbias.in")
os.system("sed s/.fits/b.fits/g subbias.in > subbias.out")
iraf.imarith('@subbias.in','-','mbias.fits','@subbias.out')
print("***********************************************************")
print("Bias correction done!")


#For Flat correction.
os.system("ls halogen*b.fits > flat.in")
iraf.imcombine('@flat.in', 'mflat.fits', combine="median",rdnoise=7.5,gain=2.0)

#making master flat.
iraf.response("mflat","mflat", "nmflat", interactive="yes", threshold="INDEF", sample="*", naverage=1,
      function="spline3", order=20, low_reject=0., high_reject=0., niterate=1,
      grow=0., graphics="stdgraph", cursor="")

iraf.imstat('nmflat.fits')    
os.system("ls J105829*b.fits GD248*b.fits > flatf.in")
os.system("sed s/b.fits/bf.fits/g flatf.in > flatf.out")
iraf.imarith('@flatf.in','/','nmflat.fits','@flatf.out')
print("***********************************************************")
print("Flat correction done!")



# Removal of cosmic rays from the images
array,header=cosmics.fromfits("J105829_250_3500_4_Grism_7_2015_06_17_yf170018bf.fits")
c= cosmics.cosmicsimage(array, gain=1.22, readnoise=4.8, sigclip=4.5, sigfrac=0.3, objlim=4.5)
c.run(maxiter=4)
cosmics.tofits("J105829_250_3500_4_Grism_7_2015_06_17_yf170018bfc.fits", c.cleanarray, header)

array, header=cosmics.fromfits("GD248_standard_250_3500_4_Grism_7_yf170031bf.fits")
c= cosmics.cosmicsimage(array, gain=1.22, readnoise=4.8, sigclip=4.5, sigfrac=0.3, objlim=4.5)
c.run(maxiter=4)
cosmics.tofits("GD248_standard_250_3500_4_Grism_7_yf170031bfc.fits", c.cleanarray, header)

print("***********************************************************")
print("Cosmic rays have been removed! Want to see the change (yes/no) ???")

c=raw_input()

if  c=="yes":
	os.system("ds9  J105829_250_3500_4_Grism_7_2015_06_17_yf170018bfc.fits") 
	print("The raw image:")
	os.system("ds9 J105829_250_3500_4_Grism_7_2015_06_17_yf170018.fits") 




#---------------------------------------------------------------------------------------------------------------------------------------------------
# At this point the image is bias subtracted, flat corrected and the cosmic rays have been removed.
#----------------------------------------------------------------------------------------------------------------------------------------------------

'''


def iraf_fwhm():
    print(' We need to get the FWHM using IRAF tasks.')
    print(' Opening DS9 for image display.')
    os.system("ds9 &")
    filename=input('Enter the filename:')
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
    print('*****************************')
    print('now selecting the sources')
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


