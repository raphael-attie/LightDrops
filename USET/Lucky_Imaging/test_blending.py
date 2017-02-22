import glob
from math import sqrt
from astropy.io import fits
import numpy as np
#import matplotlib
#matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
from scipy.ndimage.morphology import distance_transform_edt

plt.interactive(False)

import ra_graphics as ra
import ra_lucky_imaging as li
from uset_calibration import *

# files   = glob.glob('/Users/raphaela/Data/USET/HAlpha/continuous/2016_Aug_17/registered/*.fits')
# f1      = files[0]
# f2      = files[100]
# hdu1    = fits.open(f1)
# hdu2    = fits.open(f2)
# image1  = hdu1[0].data
# image2  = hdu2[0].data
shape = [512, 512]
canvas_shape = [512, 1024]
canvas1  = np.zeros(canvas_shape) + 50
canvas2  = np.zeros(canvas_shape) + 50
count_plane1  = np.zeros(canvas_shape)
count_plane2  = np.zeros(canvas_shape)
weight_plane = np.zeros(canvas_shape)
image1  = np.zeros(shape)
image2  = np.ones(shape)*255

# Pre-compute the blending curve(s)
# Use distance transform to calculate distance from the 0-level edges
weights_1D      = np.ones(shape[0])
weights_1D[0]   = 0
weights_1D[-1]  = 0
weights_1D      = distance_transform_edt(weights_1D)
weights_2D      = np.tile(weights_1D, [shape[0], 1])

# position of image A on the canvas
x1a = 128
x2a = x1a + shape[1]
# Populate canvas with image A (Same for version with and without blending)
canvas1[:, x1a:x2a]         = image1
canvas2[:, x1a:x2a]         = image1
# Increment the count plane
count_plane1[:, x1a:x2a]    = 1
count_plane2[:, x1a:x2a]    = 1
# Populate the weight_plane with the weights for image1
weight_plane[:, x1a:x2a]    = weights_2D

# Position of the image B on the canvas
x1b = x2a - shape[1]/2
x2b = x1b + shape[1]
# Fill the canvas1 (version without blending) with the 2nd image
canvas1[:, x1b:x2b] += image2
count_plane1[:, x1b:x2b]  += 1
count_plane1[count_plane1 == 0] = 1
canvas1 /= count_plane1

# Increment the count plane where image2 is going (used for canvas2, version with blending)
count_plane2[:, x1b:x2b]  += 1
# Overlaping region in the subfield reference frame: use the coordinate x1b, x2b
count_plane_subfield    = count_plane2[:, x1b:x2b]
overlap_mask_subfield   = count_plane_subfield == 2
# Overlap region in the canvas reference frame
overlap_mask            = count_plane2 == 2
# "Bottom" data from the canvas that are in the overlap region
canvas_overlap          = canvas2[overlap_mask]
canvas_overlap_weights  = weight_plane[overlap_mask]
# Fill canvas2 with image2: use the coordinate x1b, x2b
canvas2[:, x1b:x2b] = image2

# Subfield values and weights in the overlap region
subfield_overlap = image2[overlap_mask_subfield]
subfield_weights = weights_2D[overlap_mask_subfield]
# Now calculate the blended values
blend = (canvas_overlap * canvas_overlap_weights + subfield_overlap * subfield_weights) / (canvas_overlap_weights + subfield_weights)
# Put that back into the canvas
canvas2[overlap_mask] = blend


# Display figures
x_a = np.arange(x1a, x2a)
x_b = np.arange(x1b, x2b)
weight1 = np.zeros(canvas_shape[1])
weight2 = np.zeros(canvas_shape[1])
weight1[x_a] = weights_1D
weight2[x_b] = weights_1D
weight_summed = weight1 + weight2
weight_summed[weight_summed == 0] = np.nan



fig = plt.figure(1)
# Plot the non-blended and overlay the distance-based weights
canvasA             = np.zeros(canvas_shape) + 50
canvasB             = np.zeros(canvas_shape) + 50
canvasA[:, x1a:x2a] = image1
canvasB[:, x1b:x2b] = image2

ax1 = plt.subplot(211)
plt.title('1D distance weight')
plt.imshow(canvasA, cmap='gray', vmin=0, vmax=255)
ax1.set_xlim(0, canvas_shape[1])
ax1.invert_yaxis()
# overlay weight individually
ax1_2 = ax1.twinx()
ax1_2.plot(weight1, 'b-')
ax1_2.set_ylabel('Weights')

ax2 = plt.subplot(212, sharex=ax1)
plt.imshow(canvasB, cmap='gray', vmin=0, vmax=255)
ax2.set_xlim(0, canvas_shape[1])
ax2.invert_yaxis()
# overlay weight individually
ax2_2 = ax2.twinx()
ax2_2.plot(weight2, 'g-')
ax2_2.set_ylabel('Weights')

fig = plt.figure(2)

# Plot the blends
ax3 = plt.subplot(211, sharex=ax1)
plt.title('Blended images with arithmetic mean')
plt.imshow(canvas1, cmap='gray')
plt.gca().invert_yaxis()
ax3.set_xlim(0, canvas_shape[1])

ax4 = plt.subplot(212, sharex=ax1)
plt.title('Blended images with 1D distance weights')
plt.imshow(canvas2, cmap='gray')
plt.gca().invert_yaxis()

# Overlay the weights
ax4.set_xlim(0, canvas_shape[1])
plt.xlabel('X')
plt.ylabel('Y')

ax4_2 = ax4.twinx()
ax4_2.plot(weight1/weight_summed)
ax4_2.plot(weight2/weight_summed)
ax4_2.set_xlim(0, canvas_shape[1])
ax4_2.set_ylabel('Weights')
ax4_2.set_ylim(0, 1.2)


plt.show()
