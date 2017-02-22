"""
Test the blending of 2 synthetic constant images using various blending methods.
Does not require USET images.
"""

import glob
from math import sqrt
from astropy.io import fits
import numpy as np
# import matplotlib
import matplotlib.pyplot as plt
from scipy.ndimage.morphology import distance_transform_edt
import ra_graphics as ra

plt.interactive(False)

outdir = '/Users/rattie/PycharmProjects/USET/plots'


shape = np.array([512, 512])
canvas_shape = [1024, 1024]
canvas = np.zeros(canvas_shape)
count_plane = np.zeros(canvas_shape)
weight_plane = np.zeros(canvas_shape)
imageA = np.ones(shape) * 80
imageB = np.ones(shape) * 255

# Pre-compute the blending weights
# The weights are the distance of each pixel to the circle that encompasses the block corners
# radius = block diagonal = shape[0] *  sqrt(2) (where shape are the dimensions of square block => shape[0] = shape[1]

# Radial weights
x = np.arange(0, shape[0]) - shape[0] / 2 + 0.5
xx, yy = np.meshgrid(x, x)
dist_to_center = np.sqrt(xx ** 2 + yy ** 2)
radius = shape[0] / 2 #* np.sqrt(2)
weights_2D = radius - dist_to_center
# weights_2D[weights_2D < 0] = 0

# Bidirectionnal
# x = np.arange(0, shape[0]) - shape[0]/2 + 0.5
# x[0:256]    = np.flipud(x[0:256])
# x[256:512]  = np.flipud(x[256:512])
# xx, yy      = np.abs(np.meshgrid(x, x))
# weights_2D  = np.sqrt(xx**2 + yy**2)

# Linear
# Use distance transform to calculate distance from the 0-level edges
# weights_1D      = np.ones(shape[0])
# weights_1D[0]   = 0
# weights_1D[-1]  = 0
# weights_1D      = distance_transform_edt(weights_1D)
# weights_2D      = np.tile(weights_1D, [shape[0], 1])

# Arithmetic average
# weights_2D  = np.ones(shape)


# position of image A on the canvas [here x1, x2, ... are two-compoment vectors]
x1a = np.array([128, 256])
x2a = x1a + shape
# Add the weights to the weight_plane at coordinates of image A
weight_plane[x1a[1]:x2a[1], x1a[0]:x2a[0]] = weights_2D
# Add image A on canvas
canvas[x1a[1]:x2a[1], x1a[0]:x2a[0]] = imageA * weights_2D

# Position of image B on the canvas
x1b = np.array([x1a[0] + shape[0] / 2, x1a[1]])
x2b = x1b + shape
# Add the weights to the weight_plane at coordinates of image B
weight_plane[x1b[1]:x2b[1], x1b[0]:x2b[0]] += weights_2D
# Add weighted image B
canvas[x1b[1]:x2b[1], x1b[0]:x2b[0]] += imageB * weights_2D

# Divide by the weight_plane (check for zeros)
weight_plane[weight_plane == 0] = 1
canvas /= weight_plane

# # Display figures

# Plot the non-blended
canvasA = np.zeros(canvas_shape) + 50
canvasB = np.zeros(canvas_shape) + 50
canvasA[x1a[1]:x2a[1], x1a[0]:x2a[0]] = imageA
canvasB[x1b[1]:x2b[1], x1b[0]:x2b[0]] = imageB

fig = plt.figure(1)
# Block A
ax1 = plt.subplot(121)
plt.title('Block A')
plt.imshow(canvasA, cmap='gray', vmin=0, vmax=255)
ax1.invert_yaxis()
ax1.set_xlim(0, canvas_shape[1])
ax1.set_ylim(0, canvas_shape[0])
# Block B
ax2 = plt.subplot(122)
plt.title('Block B')
plt.imshow(canvasB, cmap='gray', vmin=0, vmax=255)
ax2.invert_yaxis()
ax1.set_xlim(0, canvas_shape[1])
ax1.set_ylim(0, canvas_shape[0])

plt.savefig(outdir + '/non-blended.png')

# Plot the blended images
fig = plt.figure(2)
fig.clear()
ax1 = plt.subplot(121)
plt.title('2D distance weights')
plt.imshow(weights_2D, cmap='gray', vmin=0, vmax=255)

ax2 = plt.subplot(122)
plt.title('Blended images')
plt.imshow(canvas, cmap='gray', vmin=0, vmax=255)
ax2.invert_yaxis()
plt.xlabel('X')
plt.ylabel('Y')

# plt.savefig(outdir + '/blend_2D_arithmetic_mean.png')
# plt.savefig(outdir + '/blend_2D_linear.png')
plt.savefig(outdir + '/blend_2D_radial_sideway.png')
# plt.savefig(outdir + '/blend_2D_radial2_sideway.png')
# plt.savefig(outdir + '/blend_2D_bidrectionnal_sideway.png')
# plt.savefig(outdir + '/blend_2D_linear_two-component_sideway.png')

plt.show()

print('done')