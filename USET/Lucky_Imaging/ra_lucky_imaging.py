#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 13:43:46 2016

@author: raphaela
"""

from math import sqrt
from astropy.io import fits
import numpy as np
from skimage.feature import register_translation
from scipy.ndimage import fourier_shift
from scipy.signal import convolve2d
#from astropy.convolution import convolve
#from astropy.convolution.kernels import CustomKernel
#from scipy.ndimage.filters import laplace
from scipy.ndimage.morphology import distance_transform_edt


# Notes on vocabulary: frame ~ image ~ array (all synonyms)
# When written with plural: arrays, frames, images ~ they designate 3D data cube.
# When using indexing, the word "1st" in a comment refers to the 0th index.

def rebin(a, new_shape, operation='sum'):
    """
    Resizes a 2d array by summing elements,
    new dimensions must be integral factors of original dimensions
    Parameters
    ----------
    a : array_like
        Input array.
    new_shape : tuple of int
        Shape of the output array
    Returns
    -------
    rebinned_array : ndarray
        If the new shape is smaller of the input array, the data are summed or averaged,
        if the new shape is bigger array elements are repeated
    Examples
    --------
    a = np.array([[0, 1], [2, 3]])
    b = rebin(a, (4, 6)) #upsize
    b
    array([[0, 0, 0, 1, 1, 1],
           [0, 0, 0, 1, 1, 1],
           [2, 2, 2, 3, 3, 3],
           [2, 2, 2, 3, 3, 3]])
    c = rebin(b, (2, 3)) #downsize
    c
    array([[ 0. ,  0.5,  1. ],
           [ 2. ,  2.5,  3. ]])
    """
    if not operation.lower() in ['sum', 'mean', 'average', 'avg']:
        raise ValueError("Operation {} not supported.".format(operation))

    M, N = a.shape
    m, n = new_shape
    if m < M:
        if operation.lower() == "sum":
            return a.reshape([m, M / m, n, N / n]).sum(3).sum(1)
        elif operation.lower() in ["mean", "average", "avg"]:
            return a.reshape([m, M / m, n, N / n]).mean(3).mean(1)
    else:
        return np.repeat(np.repeat(a, m / M, axis=0), n / N, axis=1)


def rebin2(a):
    # Rebin the input image "a" with 2x2 summing
    b1 = a[0::2, 0:] + a[1::2, 0:]
    b2 = b1[0:, 0::2] + b1[0:, 1::2]
    return b2


def rebin_nd(ndarray, new_shape, operation='sum'):
    """
    From https://gist.github.com/derricw/95eab740e1b08b78c03f
    Bins an ndarray in all axes based on the target shape, by summing or
        averaging.
    Number of output dimensions must match number of input dimensions.
    Example
    -------
    m = np.arange(0,100,1).reshape((10,10))
    n = rebin_nd(m, new_shape=(5,5), operation='sum')
    print(n)
    [[ 22  30  38  46  54]
     [102 110 118 126 134]
     [182 190 198 206 214]
     [262 270 278 286 294]
     [342 350 358 366 374]]
    """
    if not operation.lower() in ['sum', 'mean', 'average', 'avg']:
        raise ValueError("Operation {} not supported.".format(operation))
    if ndarray.ndim != len(new_shape):
        raise ValueError("Shape mismatch: {} -> {}".format(ndarray.shape,
                                                           new_shape))
    compression_pairs = [(d, c // d) for d, c in zip(new_shape,
                                                     ndarray.shape)]
    flattened = [l for p in compression_pairs for l in p]
    ndarray = ndarray.reshape(flattened)
    for i in range(len(new_shape)):
        if operation.lower() == "sum":
            ndarray = ndarray.sum(-1 * (i + 1))
        elif operation.lower() in ["mean", "average", "avg"]:
            ndarray = ndarray.mean(-1 * (i + 1))
    return ndarray


def block_processing_setup(arrays, binning):
    """
    Create a binned version of the full sun image, filled with values of sharpness metrics
    Uses a custom kernel, 2nd-degree, Laplacian-like

    :param arrays: image series (data cube)
    :param binning: binning factor.
    :return: Binned array of size = arrays size / binning.
    """
    # dimensions of the binned array
    naxis1, naxis2, nframes = arrays.shape
    nbaxis1 = int(naxis1 / binning)
    nbaxis2 = int(naxis2 / binning)

    kernel = np.array([[1, 4, 1],
                      [4, -20, 4],
                      [1, 4, 1]])

    # Initialize quality array - "q" for quality
    qBinnedArrays = np.zeros([nbaxis2, nbaxis1, nframes])
    for k in range(0, nframes):
        # Get a binned version of the arrays
        frame           = np.squeeze(arrays[:, :, k])
        # binnedFrame     = rebin(frame, new_shape=(nbaxis2, nbaxis1), operation='sum')
        binnedFrame = rebin2(frame)
        qualityFrame = convolve2d(binnedFrame, kernel, mode='same', boundary='symm')  # laplace(binnedFrame)
        # print 'qualityFrame ='
        # print qualityFrame[0:10, 0:10]
        qBinnedArrays[:, :, k] = qualityFrame

    return qBinnedArrays


def make_aligned_stack(arrays, qbinned_arrays, nbest, blk_size, binned_blk_size, binning, x, y):
    # Position of the block in the qbinned_arrays
    xB = int(x / binning)
    yB = int(y / binning)
    # Extract a series of small block from the quality arrays
    binned_blks = qbinned_arrays[yB: yB + binned_blk_size, xB: xB + binned_blk_size, :]
    binned_blks = binned_blks.reshape(binned_blks.shape[0]*binned_blks.shape[1], binned_blks.shape[2])
    # Quality metric is variance of laplacian, unbiased.
    quality = np.var(binned_blks, 0, ddof=1)

    # Get the sorting indices that sort the quality in descending order (use ::-1 for flipping the vector)
    sort_idx = np.argsort(quality)[::-1]
    bestIndices = sort_idx[0:nbest]

    offset = 4

    # Extract best blocks and reference blocks 
    bestBlks0   = arrays[y - offset: y + blk_size + offset, x - offset: x + blk_size + offset, bestIndices]
    bestBlks    = arrays[y: y + blk_size, x: x + blk_size, bestIndices]
    refBlk      = bestBlks0[:, :, 0]
    # Get the misalignment of each block to the reference block using phase correlation in Fourier space
    shifts = np.zeros([2, nbest-1])
    for i in range(1, nbest):
        blk = bestBlks0[:, :, i]
        shift, error, diffPhase = register_translation(refBlk, blk)
        shifts[:, i-1] = shift
        # align the blocks
        # shifted_blk = fourier_shift(np.fft.fftn(blk), shift)
        # shifted_blk = np.fft.ifftn(shifted_blk)
        # bestBlks[:, :, i] = shifted_blk[offset:offset + blk_size, offset:offset + blk_size]

        # Below, precision of the shift is downgraded to 0.5 px at worse due to rounding
        # We have to subtract the shift instead of summing it, because we are cropping a new block in the original array
        # and not actually shifting the "blk" used above

        xs = int(x - round(shift[1]))
        ys = int(y - round(shift[0]))
        bestBlks[:, :, i] = arrays[ys: ys + blk_size, xs: xs + blk_size, bestIndices[i]]

    return bestBlks, shifts, bestIndices

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * sig**2))

def lucky_imaging(images, globalRefImage, blk_size, nbest, binning):

    binnedBlkSize = int(blk_size / binning)
    # Define some starting coordinates where the block processing will start
    x = 0
    y = 0
    # Initialize Block number
    k = 0
    # offset space for registration
    offset = 4

    # Bound the block processing, regardless of the radius offset. This defines 
    # at which coordinates the while loop ends
    naxis1 = images.shape[1]
    naxis2 = images.shape[0]

    # Offset Coordinate from which the processing actually starts,
    # this avoids going indexing off the image dimensions
    offsetX = blk_size
    offsetY = blk_size
    naxis1End = naxis1 - 1 - 3 * blk_size
    naxis2End = naxis2 - 1 - 3 * blk_size
    # Step 1 of block processing: Get the binned matrix of the quality 2nd order metrics
    qBinnedArrays = block_processing_setup(images, binning)
    # qBinnedArrays is identical as its alternate in Lightdrops / C++

    # Initialize the canvas and weights that will contain the lucky-imaged array and weights
    canvas3D = np.zeros([naxis2, naxis1, nbest])
    weightCanvas = np.zeros([naxis2, naxis1, nbest])
    # Distance-based blending curves
    shape = np.array([blk_size, blk_size], dtype=int)
    weights = np.zeros(shape)
    weights[1:-1, 1:-1] = 1
    dist_weights = distance_transform_edt(weights)
    gauss_dist_weights = gaussian(dist_weights, shape[0] / 2, shape[0] / 8)


    shift_list = []

    ## Start main loop. Account for 50% overlap between consecutive blocks
    while y < naxis2End or x < naxis1End:
        # Coordinates of the current block
        x = int(offsetX + (k * blk_size / 2 % naxis1End))
        y = int(offsetY + ((k * blk_size / 2) / naxis1End) * blk_size / 2)

        stackedBlks, shifts, best_indices = make_aligned_stack(images, qBinnedArrays, nbest, blk_size, binnedBlkSize, binning, x, y)
        blkSlice = stackedBlks[:, :, 0]

        ## Find the best alignment of the stacked Block onto the global reference image
        # Reference block is the one at expected position in the global reference image
        refBlk = globalRefImage[y: y + blk_size, x: x + blk_size]
        #refBlk = globalRefImage[y-offset: y + blk_size + offset, x - offset: x + blk_size + offset]

        # Again, use phase correlation in fourier space. 
        shift, error, diffPhase = register_translation(refBlk, blkSlice)
        # The shifted positions here are in the reference frame of the stacked block,
        # so we need to add the shift instead of subtracting it. We are indeed truly shifting 
        # the position of the block
        # shift_residue = shift - np.round(shift)

        # for i in range(0, nbest):
        #     blk = stackedBlks[:, :, i]
        #     # Shift the block by the shift residue in fourier space.
        #     shifted_blk = fourier_shift(np.fft.fftn(blk), shift_residue)
        #     shifted_blk = np.fft.ifftn(shifted_blk)
        #     stackedBlks[:, :, i] = shifted_blk

        xs = int(x + round(shift[1]))
        ys = int(y + round(shift[0]))

        # Fill the canvas with the stackedBlk, at the shifted position
        # With arithmetic average
        canvas3D[ys: ys + blk_size, xs: xs + blk_size, :] += stackedBlks
        weightCanvas[ys: ys + blk_size, xs: xs + blk_size, :] += 1
        # With distance-based weighted average (blending)
        # canvas3D[ys: ys + blk_size, xs: xs + blk_size, :] += stackedBlks * gauss_dist_weights
        # weightCanvas[ys: ys + blk_size, xs: xs + blk_size, :] += gauss_dist_weights

        # if abs(shift[0]) > 1 or abs(shift[1]) > 1:
        #      print "[k, x, y] = [%d, %d, %d] ; shift= [%.2f, %.2f]" % (k, x, y, shift[1], shift[0])

        # Store the shift vectors
        shift_list.append(shift)
        # Increment the block number
        k += 1
        # end of while loop

    weightMask                  = weightCanvas == 0
    weightCanvas[weightMask]    = 1
    canvas3D                    /= weightCanvas
    canvas                      = np.median(canvas3D, 2)

    # Need to fill the array near the edges where the images were not processed
    fillEnd                = int(round(1.5 * blk_size))
    canvas[0:fillEnd, :]   = globalRefImage[0:fillEnd, :]
    canvas[:, 0:fillEnd]   = globalRefImage[:, 0:fillEnd]
    canvas[naxis2End:, :]  = globalRefImage[naxis2End:, :]
    canvas[:, naxis1End:]  = globalRefImage[:, naxis1End:]

    return canvas, shift_list


def lucky_imaging_single_block(images, globalRefImage, blk_size, nbest, binning, x, y):
    # Size of the binned block
    binnedBlkSize = blk_size / binning
    # Step 1 of block processing: Get the binned matrix of the quality metrics
    # qBinnedArrays is identical as its alternate in Lightdrops / C++
    qBinnedArrays = block_processing_setup(images, binning)

    stackedBlks, shifts, best_indices = make_aligned_stack(images, qBinnedArrays, nbest, blk_size, binnedBlkSize, binning, x, y)
    blkSlice = stackedBlks[:, :, 0]

    ## Find the best alignment of the stacked Block onto the global reference image
    # Reference block is the one at expected position in the global reference image
    refBlk = globalRefImage[y: y + blk_size, x: x + blk_size]
    # refBlk = globalRefImage[y-offset: y + blk_size + offset, x - offset: x + blk_size + offset]

    # Again, use phase correlation in fourier space.
    shift, error, diffPhase = register_translation(refBlk, blkSlice)
    # The shifted positions here are in the reference frame of the stacked block,
    # so we need to add the shift instead of subtracting it. We are indeed truly shifting
    # the position of the block

    stacked_blk = np.median(stackedBlks, 2)
    xs = int(x + round(shift[1]))
    ys = int(y + round(shift[0]))

    print('done')

    return stacked_blk, shift, best_indices, xs, ys
