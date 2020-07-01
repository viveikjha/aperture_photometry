import numpy as np
import math
from scipy import optimize
from astropy.io import fits
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
def gaussian(height, center_x, center_y, width_x, width_y):
	"""Returns a gaussian function with the given parameters"""
	width_x = float(width_x)
	width_y = float(width_y)
	return lambda x,y: height*np.exp(-(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)



def moments(data):
	"""Returns (height, x, y, width_x, width_y)
	the gaussian parameters of a 2D distribution by calculating its
	moments """
	total = data.sum()
	X, Y = np.indices(data.shape)
	x = (X*data).sum()/total
	y = (Y*data).sum()/total
	col = data[:, int(y)]
	width_x = np.sqrt(np.abs((np.arange(col.size)-y)**2*col).sum()/col.sum())
	row = data[int(x), :]
	width_y = np.sqrt(np.abs((np.arange(row.size)-x)**2*row).sum()/row.sum())
	height = data.max()
	print(width_x,width_y,height)
	return height, x, y, width_x, width_y


def fitgaussian(data):
	"""Returns (height, x, y, width_x, width_y)
	the gaussian parameters of a 2D distribution found by a fit"""
	params = moments(data)
	errorfunction = lambda p: np.ravel(gaussian(*p)(*np.indices(data.shape)) -data)
	p, success = optimize.leastsq(errorfunction, params)
	return p

data = fits.open('bf.fits')
image=data[0].data
psf = image[1200:1550,1500:1850]
print(np.shape(psf))
result = moments(psf)

print (result[3],result[4])

print ("seeing = ", (result[3]+result[4])/2.0*0.0439)

fig=plt.figure()
ax=plt.plot(psf)
plt.show()
#ax.plot3D(psf)
# 1.35"
