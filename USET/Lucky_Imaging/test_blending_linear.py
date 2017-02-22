import glob
from math import sqrt
from astropy.io import fits
import numpy as np
# import matplotlib
import matplotlib.pyplot as plt
from scipy.ndimage.morphology import distance_transform_edt

plt.interactive(False)

import ra_graphics as ra
import Lucky_Imaging.ra_lucky_imaging as li
from Calibration.uset_calibration import *

outdir = '/Users/raphaela/Dropbox/Python/PycharmProjects/USET/plots'

shape = [512, 512]
canvas_shape = [512, 1024]
canvas = np.zeros(canvas_shape)
weight_plane = np.zeros(canvas_shape)
image1 = np.ones(shape) * 80
image2 = np.ones(shape) * 255

# Pre-compute the blending curve(s)
# Use distance transform to calculate distance from the 0-level edges
weights_1D = np.ones(shape[0])
weights_1D[0] = 0
weights_1D[-1] = 0
weights_1D = distance_transform_edt(weights_1D)
weights_2D = np.tile(weights_1D, [shape[0], 1])

# position of image A on the canvas
x1a = 128
x2a = x1a + shape[1]
# Add the weights to the weight_plane at coordinates of image A
weight_plane[:, x1a:x2a] += weights_2D
# Add weighted image A
canvas[:, x1a:x2a] += image1 * weights_2D

# Position of the image B on the canvas
x1b = x2a - shape[1] / 2
x2b = x1b + shape[1]
# Add the weights to the weight_plane at coordinates of image B
weight_plane[:, x1b:x2b] += weights_2D
# Add weighted image B
canvas[:, x1b:x2b] += image2 * weights_2D

# Divide by the weight_plane (check for zeros)
weight_plane[weight_plane == 0] = 1
canvas /= weight_plane

# Display figures
x_a = np.arange(x1a, x2a)
x_b = np.arange(x1b, x2b)
weight1 = np.zeros(canvas_shape[1])
weight2 = np.zeros(canvas_shape[1])
weight1[x_a] = weights_1D
weight2[x_b] = weights_1D
weight_summed = weight1 + weight2
weight_summed[weight_summed == 0] = np.nan

# Plot the non-blended and overlay the distance-based weights
fig = plt.figure(1)
canvasA = np.zeros(canvas_shape)
canvasB = np.zeros(canvas_shape)
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

# Plot the blended images
fig = plt.figure(2)
plt.title('Blended images with 1D distance weights')
plt.imshow(canvas, cmap='gray')
ax = plt.gca()
ax.invert_yaxis()
ax.set_xlim(0, canvas_shape[1])
plt.xlabel('X')
plt.ylabel('Y')
ax.axis('tight')

# Overlay the weights
ax2 = ax.twinx()
ax2.plot(weight1 / weight_summed)
ax2.plot(weight2 / weight_summed)
ax2.set_xlim(0, canvas_shape[1])
ax2.set_ylabel('Weights')
ax2.set_ylim(0, 1.0)

plt.savefig(outdir + '/blend_1D_linear.png')

plt.show()
