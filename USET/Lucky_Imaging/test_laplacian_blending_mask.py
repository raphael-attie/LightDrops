import sys
import os
import numpy as np
import cv2
import scipy
from scipy.stats import norm
from scipy.signal import convolve2d
import math


# generate a 5x5 kernel
def generating_kernel(a):
    w_1d = np.array([0.25 - a / 2.0, 0.25, a, 0.25, 0.25 - a / 2.0])
    return np.outer(w_1d, w_1d)


# reduce image by 1/2
def ireduce(image):
    out = None
    kernel = generating_kernel(0.4)
    outimage = scipy.signal.convolve2d(image, kernel, 'same')
    out = outimage[::2, ::2]
    return out


# expand image by factor of 2
def iexpand(image):
    out = None
    kernel = generating_kernel(0.4)
    outimage = np.zeros((image.shape[0] * 2, image.shape[1] * 2), dtype=np.float64)
    outimage[::2, ::2] = image[:, :]
    out = 4 * scipy.signal.convolve2d(outimage, kernel, 'same')
    return out


# create a gaussain pyramid of a given image
def gauss_pyramid(image, levels):
    output = []
    output.append(image)
    tmp = image
    for i in range(0, levels):
        tmp = ireduce(tmp)
        output.append(tmp)
    return output


# build a laplacian pyramid
def lapl_pyramid(gauss_pyr):
    output = []
    k = len(gauss_pyr)
    for i in range(0, k - 1):
        gu = gauss_pyr[i]
        egu = iexpand(gauss_pyr[i + 1])
        if egu.shape[0] > gu.shape[0]:
            egu = np.delete(egu, (-1), axis=0)
        if egu.shape[1] > gu.shape[1]:
            egu = np.delete(egu, (-1), axis=1)
        output.append(gu - egu)
    output.append(gauss_pyr.pop())
    return output


# Blend the two laplacian pyramids by weighting them according to the mask
def blend(lapl_pyr_white, lapl_pyr_black, gauss_pyr_mask):
    blended_pyr = []
    k = len(gauss_pyr_mask)
    for i in range(0, k):
        p1 = gauss_pyr_mask[i] * lapl_pyr_white[i]
        p2 = (1 - gauss_pyr_mask[i]) * lapl_pyr_black[i]
        blended_pyr.append(p1 + p2)
    return blended_pyr


# Reconstruct the image based on its laplacian pyramid
def collapse(lapl_pyr):
    output = None
    output = np.zeros((lapl_pyr[0].shape[0], lapl_pyr[0].shape[1]), dtype=np.float64)
    for i in range(len(lapl_pyr) - 1, 0, -1):
        lap = iexpand(lapl_pyr[i])
        lapb = lapl_pyr[i - 1]
        if lap.shape[0] > lapb.shape[0]:
            lap = np.delete(lap, (-1), axis=0)
        if lap.shape[1] > lapb.shape[1]:
            lap = np.delete(lap, (-1), axis=1)
        tmp = lap + lapb
        lapl_pyr.pop()
        lapl_pyr.pop()
        lapl_pyr.append(tmp)
        output = tmp
    return output


# Main

shape = np.array([512, 512], dtype=int)
A = np.ones(shape) * 80
B = np.ones(shape) * 255

# Colorize inside each side
A[200:300, 0:100] = 255
B[200:300, 400:] = 0

# Define large image frames acting as co-spatial canvas, hosting image A and B down below
CA = np.zeros(2 * shape)
CB = np.zeros(2 * shape)
mask = np.zeros(2 * shape)

AX1, AX2, AY1, AY2 = 0, shape[1], 0, shape[0]
BX1, BX2, BY1, BY2 = shape[1]/2, shape[1] / 2 + shape[1], shape[0]/2, shape[0] / 2 + shape[0]

CA[AY1:AY2, AX1:AX2] = A
CB[BY1:BY2, BX1:BX2] = B
mask[BY1:BY2, BX1:BX2] = 1

A = CA
B = CB

cv2.imwrite('CA.jpg', CA)
cv2.imwrite('CB.jpg', CB)

# Automatically figure out the size
min_size = min(A.shape)
depth = int(math.floor(math.log(min_size, 2))) - 4  # at least 16x16 at the highest level.

gauss_pyr_mask = gauss_pyramid(mask, depth)
gauss_pyr_image1 = gauss_pyramid(A, depth)
gauss_pyr_image2 = gauss_pyramid(B, depth)

lapl_pyr_image1 = lapl_pyramid(gauss_pyr_image1)
lapl_pyr_image2 = lapl_pyramid(gauss_pyr_image2)

outpyr = blend(lapl_pyr_image2, lapl_pyr_image1, gauss_pyr_mask)

outimg = collapse(blend(lapl_pyr_image2, lapl_pyr_image1, gauss_pyr_mask))

# blending sometimes results in slightly out of bound numbers.
outimg[outimg < 0] = 0
outimg[outimg > 255] = 255
outimg = outimg.astype(np.uint8)

cv2.imwrite('blended_with_mask.jpg', outimg)

# Direct blending

mask2 = mask.copy()
mask2[mask == 0] = 1
mask2[mask == 1] = 0

direct_blend = CA*mask2 + CB*mask
cv2.imwrite('direct_blend.jpg', direct_blend)

mask *= 255
mask2 *= 255
cv2.imwrite('mask1.jpg', mask)
cv2.imwrite('mask2.jpg', mask2)
