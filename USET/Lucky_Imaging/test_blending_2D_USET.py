import glob
from math import sqrt
from astropy.io import fits
import numpy as np
# import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.ndimage.morphology import distance_transform_edt
import ra_graphics as ra
import Lucky_Imaging.ra_lucky_imaging as li
from Calibration.uset_calibration import *

plt.interactive(False)
mpl.rc('image', origin='lower', cmap='gray')

outdir = '/Users/raphaela/Dropbox/Python/PycharmProjects/USET/plots'
files = glob.glob("/Users/raphaela/Data/USET/HAlpha/continuous/2016_Aug_17/registered/*.fits")
hdu = fits.open(files[0])
# Image header
header = hdu[0].header
# Get image dimensions
naxis1 = header['NAXIS1']
naxis2 = header['NAXIS2']
nimages = len(files)
images = np.zeros([naxis2, naxis1, 4])

block_size = 256

for i in range(0, nimages):
    hdu = fits.open(files[i])
    images[:, :, i] = hdu[0].data

# Lucky imaging

# Step 1 of block processing: Get the binned matrix of the quality metrics
nbest = 5
binning = 2
binnedBlkSize = block_size / binning
qBinnedArrays = li.block_processing_setup(images, binning)

# Crop 4 overlapping subfields of size 128 x 128 px

# Reference Coordinates at which to extract the first one.
xr1, yr1 = 1450, 1050
xr2, yr2 = xr1 + block_size, yr1 + block_size
# Make the coordinates of the 4 subfields just dependent on these reference coordinates
# This way just change the coordinate above to change the area on the sun where to benchmark the overlap test.
x1, y1, x2, y2 = np.zeros(4, int), np.zeros(4, int), np.zeros(4, int), np.zeros(4, int)
# 1st subfield (bottom left side)
x1[0], y1[0] = xr1, yr1
# 2nd subfield (bottom right side)
x1[1], y1[1] = xr1 + block_size / 2, yr1
# 3rd subfield (top left side)
x1[2], y1[2] = xr1, yr1 + block_size / 2
# 4th subfield (top right side)
x1[3], y1[3] = xr1 + block_size / 2, yr1 + block_size / 2
# The end coordinates (x2, y2) are just the start coordinates (x1, y1) + block size
x2, y2 = x1 + block_size, y1 + block_size

# Extract the subfields from the data:
# Use an index s -> BL = Bottom Left = 0, BR = Bottom Right = 1, TL = Top Left = 2, TR = Top Right = 3
s = 0
stack_BL, shifts_BL = li.make_aligned_stack(images, qBinnedArrays, nbest, block_size, binnedBlkSize, binning, x1[s],
                                            y1[s])
subfield_BL = np.median(stack_BL, 2)
# subfield_BL = images[y1[s]:y2[s], x1[s]:x2[s], s]
s = 1
stack_BR, shifts_BR = li.make_aligned_stack(images, qBinnedArrays, nbest, block_size, binnedBlkSize, binning, x1[s],
                                            y1[s])
subfield_BR = np.median(stack_BR, 2)
# subfield_BR = images[y1[s]:y2[s], x1[s]:x2[s], s]
s = 2
stack_TL, shifts_TL = li.make_aligned_stack(images, qBinnedArrays, nbest, block_size, binnedBlkSize, binning, x1[s],
                                            y1[s])
subfield_TL = np.median(stack_TL, 2)
# subfield_TL = images[y1[s]:y2[s], x1[s]:x2[s], s]
s = 3
stack_TR, shifts_TR = li.make_aligned_stack(images, qBinnedArrays, nbest, block_size, binnedBlkSize, binning, x1[s],
                                            y1[s])
subfield_TR = np.median(stack_TR, 2)
# subfield_TR = images[y1[s]:y2[s], x1[s]:x2[s], s]

dmin, dmax = 2000, 3000

fig = plt.figure(1)
plt.subplot(221)
plt.imshow(subfield_TL, vmin=dmin, vmax=dmax)
plt.subplot(222)
plt.imshow(subfield_TR, vmin=dmin, vmax=dmax)
plt.subplot(223)
plt.imshow(subfield_BL, vmin=dmin, vmax=dmax)
plt.subplot(224)
plt.imshow(subfield_BR, vmin=dmin, vmax=dmax)

shape = [block_size, block_size]
# Initialize the canvas, with a shape that is enough to contain the 4 blocks that are overlaping by 50%
canvas_shape = np.array([1.5 * block_size, 1.5 * block_size], int)
canvas = np.zeros(canvas_shape)
weight_plane = np.zeros(canvas_shape)

# Pre-compute the blending weights
# The weights are the distance of each pixel to the circle that encompasses the block corners
# radius = block diagonal = shape[0] *  sqrt(2) (where shape are the dimenions of square block => shape[0] = shape[1]

# # Linear
# Use distance transform to calculate distance from the 0-level edges
# weights_1D      = np.ones(shape[0])
# weights_1D[0]   = 0
# weights_1D[-1]  = 0
# weights_1D      = distance_transform_edt(weights_1D)
# weights_2D      = np.tile(weights_1D, [shape[0], 1])

# Radial weights
# x = np.arange(0, block_size) - block_size/2 + 0.5
# xx, yy              = np.meshgrid(x, x)
# dist_to_center      = np.sqrt(xx**2 + yy**2)
# radius              = block_size/2 * np.sqrt(2)
# weights_2D          = radius - dist_to_center
# weights_2D[weights_2D < 0] = 0
#

# Arithmetic average
weights_2D = np.ones(shape)

# position of image A on the canvas [here x1, x2, ... are two-compoment vectors]
# "xc" stands for "x canvas", as the coordinates in the canvas reference frame
xc1, yc1 = x1 - xr1, y1 - yr1
xc2, yc2 = xc1 + block_size, yc1 + block_size
# Bottom Left on canvas
s = 0
canvas[yc1[s]:yc2[s], xc1[s]:xc2[s]] += subfield_BL * weights_2D
weight_plane[yc1[s]:yc2[s], xc1[s]:xc2[s]] += weights_2D
# Bottom Right on canvas
s = 1
canvas[yc1[s]:yc2[s], xc1[s]:xc2[s]] += subfield_BR * weights_2D
weight_plane[yc1[s]:yc2[s], xc1[s]:xc2[s]] += weights_2D
# Top Left on canvas
s = 2
canvas[yc1[s]:yc2[s], xc1[s]:xc2[s]] += subfield_TL * weights_2D
weight_plane[yc1[s]:yc2[s], xc1[s]:xc2[s]] += weights_2D
# Top Right on canvas
s = 3
canvas[yc1[s]:yc2[s], xc1[s]:xc2[s]] += subfield_TR * weights_2D
weight_plane[yc1[s]:yc2[s], xc1[s]:xc2[s]] += weights_2D

# Divide by the weight_plane (check for zeros)
weight_plane[weight_plane == 0] = 1
canvas /= weight_plane
#
#
# # # Display figures

# Plot the blended images
fig = plt.figure(2)
ax1 = plt.subplot(121)
plt.title('2D distance weights')
plt.imshow(weights_2D, cmap='gray', vmin=0, vmax=1)

ax2 = plt.subplot(122)
plt.title('Blended images')
plt.imshow(canvas, cmap='gray', vmin=dmin, vmax=dmax)
plt.xlabel('X')
plt.ylabel('Y')
#
# #plt.savefig(outdir + '/blend_2D_arithmetic_mean.png')
# # plt.savefig(outdir + '/blend_2D_linear.png')
# plt.savefig(outdir + '/blend_2D_radial_sideway.png')
# # plt.savefig(outdir + '/blend_2D_radial2_sideway.png')
# # plt.savefig(outdir + '/blend_2D_bidrectionnal_sideway.png')
# # plt.savefig(outdir + '/blend_2D_linear_two-component_sideway.png')

plt.show()
