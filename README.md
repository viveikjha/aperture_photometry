# aperture_photometry
A code to perform aperture photometry on data obtained through the 1.3m Devasthal  fast optical Telescope. We are in the process of developing a photometry pipeline and this code is an integral part of the same. This is a work in progress as of now and we expect it to develop fully in a couple of months time. The salient features are:

1. This is being developed as a pure PYTHON software.
2. None of the operations depend on IRAF tasks. For image display and related purposes, we use the GINGA package instead of DS9.
3. To reduce the images we use CCDPROC alongwith our custom written tools.
4. The PYTHON package PHOTUTILS is being used for the photometry.

We aim it to develop it as a fast and easy to handle alternative to the users dependent on IRAF/DAOPHOT kind of software.

Currently, we are a team of 3 people, all racing towards our PhD thesis and working on this as a side project:

1. Vivek kumar Jha
2. Dimple
3. Kiran Wani
