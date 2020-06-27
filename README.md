# aperture_photometry
A code to perform aperture photometry on data obtained through the 1.3m Devasthal  fast optical Telescope. We are in the process of developing a photometry pipeline and this code is an integral part of the same. We aim to make this package work on data from any telescope with minor tweaks here and there. The salient features are:

1. This is being developed as a pure PYTHON software.
2. None of the operations depend on IRAF tasks. We reiterate here that the aim is not to discourage IRAF and/or related pacakges, but with time since the packages have become unsupported and PYTHON packages have devloped enough to perform all related tasks.
For image display and related purposes, we use the GINGA package instead of DS9.
3. To reduce the images we use CCDPROC alongwith our custom written tools.
4. The PYTHON package PHOTUTILS is being used for the aperture photometry.


We built all our codes on Python 3.6. The required packages: 

(It is recommended to install packages though pip and keep them to the latest versions.  )

1. Numpy
2. Astropy and its affiliated packages: (astropy 4.0 and above)
3. Photutils  (photutils 0.6 and above)
4. CCDProc
5. Astroalign
6. astroquery
7. Matplotlib
8. Ginga

We aim it to develop it as a fast and easy to handle alternative to the users dependent on IRAF/DAOPHOT kind of software.

Currently, we are a team of 3 people, all racing towards our PhD thesis and working on this as a side project:

1. Vivek kumar Jha
2. Dimple
3. Kiran Wani
