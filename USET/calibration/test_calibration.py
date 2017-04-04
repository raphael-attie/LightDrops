import os
import glob
from astropy.io import fits
import matplotlib
import matplotlib.pyplot as plt
import calibration.uset_calibration as uset
import sunpy.cm as cm
from skimage.transform import warp
from skimage.transform import SimilarityTransform
from skimage.transform import rotate
import numpy as np


# Set the directory where the FITS are. Here USETDATA is an environment variable. data_dir can be set to a string
# with a more explicit path if it is not possible to set an environment variable.
data_dir    =  '/Users/rattie/Data/USET/campaign/calibration/HALPHA'
# Get the list of files, change it according to where your data files are and how are they are named.
file_list   = glob.glob(os.path.join(data_dir, '*.FTS'))
# Take 1st one
i = 0
file = file_list[i]
hdu = fits.open(file, ignore_missing_end=True)
# Load header and image data from the 1st data unit: hdu[0]
header = hdu[0].header
image = hdu[0].data
# Transform the image
tform = SimilarityTransform(translation=(-400, -400))
centered_image = warp(image, tform, order=3)

rimage = rotate(image, 30, order=3)
rimage2 = rotate(image, -30, order=3)

# Display
# plt.figure(0)
# plt.subplot(131)
# plt.imshow(image, cmap='gray', origin='lower')
# plt.subplot(132)
# plt.imshow(rimage, cmap='gray', origin='lower')
# plt.subplot(133)
# plt.imshow(rimage2, cmap='gray', origin='lower')





print('done')