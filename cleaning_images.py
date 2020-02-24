'''
This python code is written to perform the spectroscopic data reduction. Initially, the bias subtraction, flat correction and cosmic ray removal from the raw image are done.
Further, the spectroscopy related IRAF tasks are executed step by step. The arguments for these tasks can be changed depending on the requirements.
The final result obtained is the spectrum of the source with calibrated wavelength. 
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
Note:It is a work in progress and the code segments can be replaced with better/precise codes whenever they are found.
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
Written by:
Vivek K Jha,
JRF @ ARIES Nainital
27/03/2018
----------------------------------------------------------------------------------------------------------------------------------------------------------------------
'''
from astropy.io import fits,ascii
from pyraf import iraf
from iraf import noao,imred,specred
from stsci import tools
import numpy as np
import os
import sys
sys.path.append("~/Documents/cosmics.py_0.4/") #the path where cosmisc.py will be located.
import cosmics    			    
# The cosmisc.py is a python module used for the cosmic ray correction. I obtained it from: http://obswww.unige.ch/~tewes/cosmics_dot_py/


print("***********************************************************")
print("Starting spectroscopic reduction of your data!!!!.")
print("***********************************************************")
print('Please enter the path to the directory where your raw data files are=====>')
dpath = raw_input()
os.chdir(dpath)


#For Bias Correction
os.system("ls *Bias* > bias.in")
iraf.imcombine('@bias.in', 'mbias', combine="median",rdnoise=4.8,gain=1.22)
os.system("ls FeAr* J105829* halogen* GD248* > subbias.in")
os.system("sed s/.fits/b.fits/g subbias.in > subbias.out")
iraf.imarith('@subbias.in','-','mbias.fits','@subbias.out')
print("***********************************************************")
print("Bias correction done!")


#For Flat correction.
os.system("ls halogen*b.fits > flat.in")
iraf.imcombine('@flat.in', 'mflat.fits', combine="median",rdnoise=4.8,gain=1.22)

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

print("Starting apall task on the source image")

os.system("dispaxis=2")
iraf.apall(input="J105829_250_3500_4_Grism_7_2015_06_17_yf170018bfc.fits", output="J105829_250_3500_4_Grism_7_2015_06_17_yf170018bfcap.fits", 
         apertures="", format="multispec", references="", profiles="",
        interactive="yes", find="yes", recenter="yes", resize="yes", edit="yes", trace="yes",
        fittrace="yes", extract="yes", extras="yes", review="yes", line="INDEF", nsum=20,
        lower=-5., upper=5., apidtable="", b_function="chebyshev", b_order=1,
        b_sample="-10:-6,6:10", b_naverage=-3, b_niterate=0, b_low_reject=3.,
        b_high_rejec=3., b_grow=0., width=5., radius=10., threshold=0., minsep=5.,
        maxsep=1000., order="increasing", aprecenter="", npeaks="INDEF", shift="yes",
        llimit="INDEF", ulimit="INDEF", ylevel=0.1, peak="yes", bkg="yes", r_grow=0.,
        avglimits="no", t_nsum=10, t_step=10, t_nlost=3, t_function="legendre",
        t_order=2, t_sample="*", t_naverage=1, t_niterate=0, t_low_reject=3.,
        t_high_rejec=3., t_grow=0., background="fit", skybox=1, weights="variance",
        pfit="fit1d", clean="yes", saturation="INDEF", readnoise=4.8, gain=1.22,
        lsigma=4., usigma=4., nsubaps=1)
print("apall done on source")


print("Starting the apall on LAMP spectra")

iraf.apall (input="FeAr_250_3500_4_Grism_7_yf170019b.fits",output="FeAr_250_3500_4_Grism_7_2015_06_17_yf170019bap.fits", 
	apertures="", format="multispec", references="J105829_250_3500_4_Grism_7_2015_06_17_yf170018bfc.fits", profiles="",
        interactive="no", find="no", recenter="no", resize="no", edit="no", trace="no",
        fittrace="no", extract="yes", extras="no", review="yes", line="INDEF", nsum=20,
        lower=-5., upper=5., apidtable="", b_function="chebyshev", b_order=1,
        b_sample="-10:-6,6:10", b_naverage=-3, b_niterate=0, b_low_reject=3.,
        b_high_rejec=3., b_grow=0., width=5., radius=10., threshold=0., minsep=5.,
        maxsep=1000., order="increasing", aprecenter="", npeaks="INDEF", shift="yes",
        llimit="INDEF", ulimit="INDEF", ylevel=0.1, peak="yes", bkg="yes", r_grow=0.,
        avglimits="no", t_nsum=10, t_step=10, t_nlost=3, t_function="legendre",
        t_order=2, t_sample="*", t_naverage=1, t_niterate=0, t_low_reject=3.,
        t_high_rejec=3., t_grow=0., background="none", skybox=1, weights="variance",
        pfit="fit1d", clean="no", saturation='INDEF', readnoise=4.8, gain=1.22,
        lsigma=4., usigma=4., nsubaps=1)

print("Starting the identify task")

iraf.identify(images="FeAr_250_3500_4_Grism_7_2015_06_17_yf170019bap.fits", section="middle line", database="database",
	coordlist="/home/aries/Music/Atsoa2018/fear_new.dat", units="", nsum="10", match=10.,
	maxfeatures=50, zwidth=100., ftype="emission", fwidth=4., cradius=5.,
	threshold=0., minsep=2., function="spline3", order=1, sample="*", niterate=0,
	low_reject=3., high_reject=3., grow=0., autowrite="no", graphics="stdgraph",
	cursor="", crval="", cdelt="", aidpars="")

print('header editing')

iraf.hedit (images="J105829_250_3500_4_Grism_7_2015_06_17_yf170018bfcap.fits", fields="refspec1", value="FeAr_250_3500_4_Grism_7_2015_06_17_yf170019bap.fits", add="yes", addonly="no",
	delete="no", verify="no", show="yes", update="yes")


print('Starting dispersion correction')

iraf.dispcor(input="J105829_250_3500_4_Grism_7_2015_06_17_yf170018bfcap.fits",output="J105829_250_3500_4_Grism_7_2015_06_17_yf170018bfcapw.fits",     
        linearize="yes", database="database", table="", w1="3000",w2="8000", dw="INDEF", nw="INDEF", log="no", flux="yes", blank=0., samedisp="no",
	 ignoreaps="no", confirm="yes", listonly="no", verbose="no", logfile="yes")

#os.system("ls *apbfcw.fits> air.in")
#iraf.asthedit(' @air.in /home/aries/Music/Atsoa2018/airmass.dat')


iraf.standard(input= "J105829_250_3500_4_Grism_7_2015_06_17_yf170018bfcapw.fits", output="std", samesta= "yes", beam_sw= "no", apertur="",
          extinct="/home/aries/Music/Atsoa2018/iao_extinction.dat", bandwid= "INDEF", bandsep="INDEF", fnuzero=  3.6800000000000E-20,
          caldir="/scisoft/share/iraf/iraf/noao/lib/onedstds/oke1990/",  interac= "yes", graphic= "stdgraph",
          star_name="g191b2b", airmass =1.99, exptime =900, answer = "yes")


iraf.sensfunc(standard="std",  sensitiv="sens", apertures="", ignoreaps="no", logfile="logfile",
	    extinction="/home/aries/Music/Atsoa2018/iao_extinction.dat", newextinctio="extinct.dat", observatory=" ",
	    function="spline3", order=6, interactive="yes", graphs="sr",
	     marks="plus cross box", colors="2 1 3 4", cursor="", device="stdgraph", answer="yes")


iraf.calibrate (input="J105829_250_3500_4_Grism_7_2015_06_17_yf170018apbfcw.fits", output="J105829_250_3500_4_Grism_7_2015_06_17_yf170018apbfcwf.fits", 
	 extinct="yes", flux="yes", extinction="/home/aries/Music/Atsoa2018/iao_extinction.dat", 
           observatory=" ", ignoreaps="yes", sensitivity="sens", fnu="no", airmass=2.3, exptime=2700 ) 






 









