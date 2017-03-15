#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 13:43:46 2016

@author: Raphael Attie
"""
import sys
import os
import time
import math
from astropy.io import fits
import numpy as np
from skimage.feature import register_translation
from skimage.measure import block_reduce
from scipy.ndimage import fourier_shift
import scipy.ndimage as ndimage
from scipy.signal import convolve2d
from scipy.ndimage.morphology import distance_transform_edt
import cv2
from skimage.filters.rank import entropy
from skimage.morphology import disk
import calibration.uset_calibration as uset

# Notes on vocabulary: frame ~ image ~ array (all synonyms)
# When written with plural: arrays, frames, images ~ they designate 3D data cube.
# When using indexing, the word "1st" in a comment refers to the 0th index.


def rebin(image, binning):
    """
    Rebin image using skimage block_reduce()

    :param image: image to rebin
    :param binning: binning factor
    :param block_func: block processing function
    :return: rebinned image
    """
    r_image = block_reduce(image, block_size =(binning, binning), func=np.mean)
    return r_image

def block_processing_setup(arrays, binning, blk_size, qmetric):
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
    qbinned_arrays = np.zeros([nbaxis2, nbaxis1, nframes])
    for k in range(0, nframes):
        # Get a binned version of the arrays
        frame           = np.squeeze(arrays[:, :, k])
        # binnedFrame     = rebin(frame, new_shape=(nbaxis2, nbaxis1), operation='sum')
        #frame = ndimage.gaussian_filter(frame, sigma=(3, 3), order=0)

        binned_frame = frame.copy()

        if qmetric.lower() == 'laplace':
            if binning != 1:
                binned_frame = rebin(binned_frame, binning)
        # Custom "Laplace-like" quality array
            qframe = convolve2d(binned_frame, kernel, mode='same', boundary='symm')  # laplace(binnedFrame)
        # If using Gradient-entropy
        elif qmetric.lower() == 'entropy':
            # Clip background values so they do not contaminate too much the entropy at the limb
            binned_frame[binned_frame < 700] = 700
            if binning != 1:
                binned_frame = rebin(binned_frame, binning)

            gy, gx = np.gradient(binned_frame)
            binned_frame = np.sqrt(gx ** 2 + gy ** 2)
            binned_frame *= 255.0 / np.max(binned_frame)
            binned_frame = binned_frame.astype(np.uint8) #binned_frame.astype(np.uint16)
            qframe = entropy(binned_frame, disk(blk_size/binning))
        elif qmetric.lower() == 'rentropy':
            # Clip background values so they do not contaminate too much the entropy at the limb
            binned_frame[binned_frame < 700] = 700
            if binning != 1:
                binned_frame = rebin(binned_frame, binning)

            gy, gx = np.gradient(binned_frame)
            binned_frame = np.sqrt(gx ** 2 + gy ** 2)
            binned_frame *= 255.0 / np.max(binned_frame)
            binned_frame = binned_frame.astype(np.uint8)  # binned_frame.astype(np.uint16)
            # Do not calculate entropy here. This is done later for each block
            qframe = binned_frame


        qbinned_arrays[:, :, k] = qframe

    return qbinned_arrays


def make_aligned_stack(arrays, qbinned_arrays, nbest, blk_size, binned_blk_size, binning, qmetric, x, y):
    """
    Create a co-aligned series of quality-sorted blocks (aka subfields).
    Sorting uses the variance of the quality matrix.
    Alignment of the subfield uses phase correlation

    :param arrays: unsorted subfields
    :param qbinned_arrays: quality matrix
    :param nbest: number of best subfields to keep in a series
    :param blk_size: size of the subfield (px)
    :param binned_blk_size: (size of the rebinned subfield)
    :param binning: amount of binning (typically 2 or 4)
    :param qmetric: quality metric
    :param x: x-coordinate of the bottom left corner of the subfield
    :param y: y-coordinate of the bottom left corner of the subfield
    :return: arrays of quality-sorted subfields
    """
    # Position of the block in the qbinned_arrays
    xB = int(x / binning)
    yB = int(y / binning)
    nframes = arrays.shape[2]
    # Extract a series of small block from the quality arrays
    binned_blks = qbinned_arrays[yB: yB + binned_blk_size, xB: xB + binned_blk_size, :]
    binned_blks1D = binned_blks.reshape(binned_blks.shape[0]*binned_blks.shape[1], nframes)


    quality = np.arange(nframes)
    if qmetric.lower() == 'laplace':
    # Quality metric is variance of laplacian, unbiased.
        quality = np.var(binned_blks1D, 0, ddof=1)
    elif qmetric.lower() == 'entropy':
    # If using skimage entropy
        quality = np.sum(binned_blks1D, 0)
    elif qmetric.lower() == 'rentropy':
        # If using custom Shannon entropy and not skimage entropy
        for i in range(0, nframes):
            quality[i] = array_entropy(binned_blks1D[:, i], 256)

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
        shift, error, diffPhase = register_translation(refBlk, blk, 4)
        shifts[:, i-1] = shift
        # align the blocks
        # shifted_blk = fourier_shift(np.fft.fftn(blk), shift)
        # shifted_blk = np.fft.ifftn(shifted_blk)
        # bestBlks[:, :, i] = shifted_blk[offset:offset + blk_size, offset:offset + blk_size]

        # Below, precision of the shift is downgraded to 0.5 px at worse due to rounding
        # We have to subtract the shift instead of adding it, because we are cropping a new block in the original array
        # and not actually shifting the "blk" used above
        xs = int(x - round(shift[1]))
        ys = int(y - round(shift[0]))
        bestBlks[:, :, i] = arrays[ys: ys + blk_size, xs: xs + blk_size, bestIndices[i]]
        # Subpixel accuracy makes no difference and introduces more artifacts than it removes.
        # For enabling it, uncomment what the block below. And comment the line above
        # res = shift - np.round(shift)
        # temp_blk = arrays[ys: ys + blk_size, xs: xs + blk_size, bestIndices[i]]
        # aligned_blk = fourier_shift(np.fft.fftn(temp_blk), res)
        # bestBlks[:, :, i] = np.fft.ifftn(aligned_blk)




    return bestBlks, shifts, bestIndices

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * sig**2))

def lucky_imaging(images, globalRefImage, blk_size, nbest, binning, qmetric, blend_mode='aavg'):

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
    time1 = time.time()
    qBinnedArrays = block_processing_setup(images, binning, blk_size, qmetric)
    # qBinnedArrays is identical as its alternate in Lightdrops / C++
    time2 = time.time()
    print('block_processing_setup() time: %0.1f s' % (time2 - time1))

    # Initialize the canvas and weights that will contain the lucky-imaged array and weights
    canvas3D = np.zeros([naxis2, naxis1, nbest])
    weightCanvas = np.zeros([naxis2, naxis1, nbest])
    # Distance-based blending curves
    shape = np.array([blk_size, blk_size], dtype=int)
    weights = np.zeros(shape)
    weights[1:-1, 1:-1] = 1
    dist_weights = distance_transform_edt(weights)
    gauss_dist_weights = gaussian(dist_weights, shape[0] / 2, shape[0] / 8)
    gauss_dist_weights_3D = np.tile(gauss_dist_weights.reshape(blk_size, blk_size, 1), (1, 1, nbest))


    shift_list = []

    ## Start main loop. Account for 50% overlap between consecutive blocks
    while y < naxis2End or x < naxis1End:
        # Coordinates of the current block
        x = int(offsetX + (k * blk_size / 2 % naxis1End))
        y = int(offsetY + ((k * blk_size / 2) / naxis1End) * blk_size / 2)

        stackedBlks, shifts, best_indices = make_aligned_stack(images, qBinnedArrays, nbest, blk_size, binnedBlkSize, binning, qmetric, x, y)
        blkSlice = stackedBlks[:, :, 0]

        ## Find the best alignment of the stacked Block onto the global reference image
        # Reference block is the one at expected position in the global reference image
        refBlk = globalRefImage[y: y + blk_size, x: x + blk_size]
        #refBlk = globalRefImage[y-offset: y + blk_size + offset, x - offset: x + blk_size + offset]

        # Again, use phase correlation in fourier space. 
        shift, error, diffPhase = register_translation(refBlk, blkSlice, 4)
        # Subpixel accuracy makes no difference and introduces more artifacts than it removes.
        # For enabling it, uncomment what the block below.
        # res = shift - np.round(shift)
        # for i in range(0, nbest):
        #     blk = stackedBlks[:, :, i]
        #     shifted_blk = fourier_shift(np.fft.fftn(blk), res)
        #     shifted_blk = np.fft.ifftn(shifted_blk)
        #     stackedBlks[:, :, i] = shifted_blk

        # Shifted positions in the reference frame of the stacked block.
        # Need to add (+) the shift instead of subtracting (-) it.


        xs = int(x + round(shift[1]))
        ys = int(y + round(shift[0]))

        if blend_mode == 'aavg':
            # Fill the canvas with the stackedBlk, at the shifted position
            # With arithmetic average
            canvas3D[ys: ys + blk_size, xs: xs + blk_size, :] += stackedBlks
            weightCanvas[ys: ys + blk_size, xs: xs + blk_size, :] += 1
        elif blend_mode == 'gblend':
            # With distance-based weighted average (blending)
            canvas3D[ys: ys + blk_size, xs: xs + blk_size, :] += stackedBlks * gauss_dist_weights_3D
            weightCanvas[ys: ys + blk_size, xs: xs + blk_size, :] += gauss_dist_weights_3D

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

def lucky_imaging_wrapper(files, outdir, outdir_jpeg, nImages, interval, nbest, binning, blk_size, qmetric, blend_mode, preview=True):

    print('interval = ' + str(interval))
    print('blend_mode = ' + blend_mode)

    ## Get general information from one sample
    hdu = fits.open(files[0], ignore_missing_end=True)
    # Image header
    header = hdu[0].header
    # Get image dimensions
    naxis1 = header['NAXIS1']
    naxis2 = header['NAXIS2']

    # In case we loop, we'd iterate over this "i", from 0 to the total number of images
    # for i in range(0, nImages / interval):
    i = 0

    # Range of files in the list
    fileRange = np.arange(0, interval) + interval * i
    # Select the file in that range
    sFiles = np.take(files, fileRange)

    ## Import the fits files into an image series
    # First get some info from the 1st fits file header
    hdu = fits.open(sFiles[0], ignore_missing_end=True)
    # Get image header
    header = hdu[0].header
    # Initialize the numpy array
    images = np.zeros([naxis2, naxis1, interval])
    # Load exactly a series of size "interval" of fits files into the 3D array "images"
    for ii in range(0, interval):
        hdu = fits.open(sFiles[ii], ignore_missing_end=True)
        images[:, :, ii] = hdu[0].data

    globalRefImage = np.median(images, 2)
    new_image, shifts = lucky_imaging(images, globalRefImage, blk_size, nbest, binning, qmetric, blend_mode=blend_mode)

    # Export the canvas to fits files. Rount to integers. Float is useless here.
    rImage = np.rint(new_image)
    intImage = np.int16(rImage)
    hdu = fits.PrimaryHDU(intImage)
    # Build file names
    basename = 'lucky_%d_total%d_best%d_blk%d_binning%d_blend_%s_qmetric_%s' % \
               (i, interval, nbest, blk_size, binning, blend_mode, qmetric)
    basename_fits = basename + '.fits'
    fname = os.path.join(outdir, basename_fits)
    uset.write_uset_fits(intImage, header, fname)

    # Rescale the image for preview using the percentile-thresholds.
    new_image = uset.rescale_image_by_histmax(new_image)

    if preview:
        # load a colormap
        #cmap = cm.get_cmap('irissjiFUV')

        for i in range(0, nbest):
            images[:,:,i] = uset.rescale_image_by_histmax(images[:,:,i])

        fov = 512
        fovx1 = 100
        fovx2 = fovx1 + fov
        fovy1 = 100
        fovy2 = fovy1 + fov

        # Export to jpeg2000
        # Normalize
        new_image *= 255.0 / new_image.max()
        basename_jp2 = basename + '.jp2'
        fname = os.path.join(outdir_jpeg, basename_jp2)
        cv2.imwrite(fname, new_image[fovy1:fovy2, fovx1:fovx2])

        # Sample from original

        sample = images[fovy1:fovy2, fovx1:fovx2, 0]
        # Export to jpeg2000
        # Normalize
        sample *= 255.0 / sample.max()
        fname = os.path.join(outdir_jpeg, 'sample.jp2')
        cv2.imwrite(fname, sample)

        # median over first nbest images
        samples = images[fovy1:fovy2, fovx1:fovx2, 0:nbest]
        median_sample = np.median(samples, 2)
        # Export to jpeg2000
        # Normalize
        median_sample *= 255.0 / median_sample.max()
        basename_jp2 = 'median_' + str(nbest) + '.jp2'
        fname = os.path.join(outdir_jpeg, basename_jp2)
        cv2.imwrite(fname, median_sample)

    return shifts


def array_entropy(array, bins):
    """
    Calculate the Shannon entropy of an image

    :param array: input image (2D)
    :param bins: number of bins when doing the histogram
    :return: scalar value of entropy
    """

    pdf, bins = np.histogram(array, bins=bins, density=True)
    pdf[pdf == 0] = 1
    e = -np.sum(pdf * np.log2(pdf))

    return e


def gradient_entropy(array, bins):
    """
    Calculate the Shannon entropy of an image using the neighbouring pixels

    :param array: input image (2D)
    :param bins: number of bins when doing the histogram
    :return: scalar value of entropy
    """
    gy, gx = np.gradient(array)
    g = np.sqrt(gx**2 + gy**2)

    e= array_entropy(g, bins)

    return e