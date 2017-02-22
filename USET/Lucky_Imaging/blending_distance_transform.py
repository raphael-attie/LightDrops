import cv2
import numpy as np,sys
import scipy.misc
from scipy.ndimage.morphology import distance_transform_edt, distance_transform_cdt, distance_transform_bf
import scipy.ndimage.filters as fi
import matplotlib.pyplot as plt


def gauss_kern(ksize, sigma):
    """Returns a 2D Gaussian kernel array."""
    # create nxn zeros
    inp = np.zeros((ksize, ksize))
    # set element at the middle to one, a dirac delta
    inp[ksize//2, ksize//2] = 1
    # gaussian-smooth the dirac, resulting in a gaussian filter mask
    return fi.gaussian_filter(inp, sigma)


def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * sig**2))

shape = np.array([512, 512], dtype=int)
A = np.ones(shape) * 80
B = np.ones(shape) * 255


# Colorize inside each side
A[200:300, 0:100] = 255
B[200:300, 400:] = 0

# Define large image frames acting as co-spatial canvas, hosting image A and B down below
CA = np.zeros(2 * shape)
CB = np.zeros(2 * shape)
mA = np.zeros(2 * shape)
mB = np.zeros(2 * shape)

offset = 10
AX1, AX2, AY1, AY2 = offset, shape[1]+offset, offset, shape[0]+offset
BX1, BX2, BY1, BY2 = int(shape[1]/2), int(shape[1] / 2 + shape[1]), int(shape[0]/2), int(shape[0] / 2 + shape[0])

CA[AY1:AY2, AX1:AX2] = A
CB[BY1:BY2, BX1:BX2] = B

mA[AY1:AY2, AX1:AX2] = 1
mB[BY1:BY2, BX1:BX2] = 1

A = CA
B = CB
# scipy.misc.imsave('CA.jpg', CA)
# scipy.misc.imsave('CB.jpg', CB)

dmA = distance_transform_edt(mA)
dmB = distance_transform_edt(mB)
blend = (CA * dmA + CB * dmB) / (dmA + dmB)
blend[np.invert(np.isfinite(blend))] = 0

# scipy.misc.imsave('blend_distance_edt.jpg', blend)
# scipy.misc.imsave('d_transformA_edt.jpg', dmA)
# scipy.misc.imsave('d_transformB_edt.jpg', dmB)

# Calculate Gaussian values at the distance transform.
gausswA = gaussian(dmA, shape[0]/2, shape[0]/8)
gausswB = gaussian(dmB, shape[0]/2, shape[0]/8)
# Calculate the weight plane filled with the values above only where the images exist
gauss_wplaneA = np.zeros(2 * shape)
gauss_wplaneB = np.zeros(2 * shape)
gauss_wplaneA[AY1:AY2, AX1:AX2] = gausswA[AY1:AY2, AX1:AX2]
gauss_wplaneB[BY1:BY2, BX1:BX2] = gausswB[BY1:BY2, BX1:BX2]

gauss_blend = (CA * gauss_wplaneA + CB * gauss_wplaneB) / (gauss_wplaneA + gauss_wplaneB)
gauss_blend[np.invert(np.isfinite(gauss_blend))] = 0

scipy.misc.imsave('gauss_blend.jpg', gauss_blend)
scipy.misc.imsave('gausswA.jpg', gausswA)

#
# plt.figure(1)
#
# plt.subplot(321)
# plt.imshow(CA, cmap='gray')
# plt.title('CA')
#
# plt.subplot(322)
# plt.imshow(CB, cmap='gray')
# plt.title('CB')
#
# plt.subplot(323)
# plt.imshow(dmA, cmap='gray')
# plt.title('dmA')
#
# plt.subplot(324)
# plt.imshow(dmB, cmap='gray')
# plt.title('dmB')
#
# plt.subplot(325)
# plt.imshow(gauss_wplaneA, cmap='gray')
# plt.title('gauss_wplaneA')
#
# plt.subplot(326)
# plt.imshow(gauss_wplaneB, cmap='gray')
# plt.title('gauss_wplaneB')
#
# plt.figure(2)
#
# plt.subplot(121)
# plt.imshow(blend, cmap='gray')
# plt.title('distance blend')
#
# plt.subplot(122)
# plt.imshow(gauss_blend, cmap='gray')
# plt.title('gauss_blend')
