#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 18:52:35 2016

@author: raphaela
"""

import glob 
from math import sqrt
from astropy.io import fits
import numpy as np
#import matplotlib
#matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
from skimage.feature import register_translation

import ra_graphics as ra
import ra_lucky_imaging as li
from uset_calibration import *



# Number of images to consider in the list of files
nImages = 120
# Get the list of files    
files   = glob.glob("/Users/raphaela/Data/USET/HAlpha/Laurent/20150513/*.FTS")[0:nImages]

## Get general information from one sample
hdu     = fits.open(files[0])
# Image header 
header  = hdu[0].header
# Get image dimensions
naxis1  = header['NAXIS1']
naxis2  = header['NAXIS2']

## Parameter for limb fitting
numPts = 64
## Define parameters for the block processing
# Image binning 
# nb of images between each final processed image = nb of images that get sorted in the quality matrix
interval        = 20
# Number of best images to keep in the block processing.
nbest           = 3
# Binning for the quality matrix. Do not change for USET unless you know what you are doing.
binning         = 2
# Buffer space: extra marging around a reference block when aligning a template block
# bufferSpace     = 8
# Block size and binned block size
blk_size         = 32 
# Length (px) in addition to the radius, beyond which the block processing does not occur
# and the values are simply copied from the global averaged image.
# This is a radius assuming the center of the disk is at the center of the array
# radiusOffset    = 100

outdir          = '/Users/raphaela/Data/USET/HAlpha/Laurent/20150513/'

#for i in range(0, 1):
for i in range(0, nImages/interval):
    strIndex = '%02d' % i
    fileRange = np.arange(0, interval) + interval * i
    # Select the file in that range
    sFiles = np.take(files, fileRange)
    
    ## Import the fits files into an image series
    # First get some info from the 1st fits file header
    hdu     = fits.open(sFiles[0])
    # Get image header 
    header  = hdu[0].header
    # Initialize the numpy array
    images  = np.zeros([naxis2, naxis1, interval])
    # Variable for storing the radius from limb fitting
    radius  = np.zeros([interval])
    # Load exactly a series of size "interval" of fits files into the 3D array "images"
    for ii in range(0, interval):
        hdu              = fits.open(sFiles[ii])
        image            = hdu[0].data
        # Rescale the image using the percentile at 99.97% of histogram
        r_image = rescale_image_by_histmax(image)
        # Align the solar disk in the center of the image reference frame
        centeredImage, xc, yc, Rm, limbPts, limbPts2, limbPts3 = uset_limbfit_align(r_image)
        radius[ii]       = Rm
        images[:, :, ii] = centeredImage

        ## Write aligned images as FITS files
        strIndex2 = '%02d' % ii
        if ii == 0:
            rImage = np.rint(centeredImage)
            intImage = np.int16(rImage)
            hdu = fits.PrimaryHDU(intImage)
            fname = outdir + 'python_aligned/' + 'Aligned_' + strIndex + '_' + strIndex2 + '.fits'
            fits.writeto(fname, intImage, header)

        ## Plot the figure with the fitted circle
        # fig = plt.figure(1)
        # plt.cla()
        # plt.autoscale(enable=True, tight=True)
        # plt.imshow(image, cmap='gray')
        # plt.xlabel('X')
        # plt.ylabel('Y')
        # plt.gca().invert_yaxis()
        # ax = plt.gca()
        # for j in range(0, limbPts.shape[0]):
        #     plt.plot(limbPts[j, 0], limbPts[j, 1], 'ro', markersize = 6)
        # for j in range(0, limbPts2.shape[0]):
        #     plt.plot(limbPts2[j, 0], limbPts2[j, 1], 'go', markersize = 4)
        # for j in range(0, limbPts3.shape[0]):
        #     plt.plot(limbPts3[j, 0], limbPts3[j, 1], 'bo', markersize = 2)
        # ## Draw the fitted circle
        # circle1 = plt.Circle((xc, yc), Rm, color='y', fill=False, linewidth=0.5)
        # ax.add_artist(circle1)
        # fname = outdir + 'Python_Aligned/' + 'Plot_limbfit_' + strIndex + '.pdf'
        # fig.savefig(fname)
        globalRefImage  = np.median(images, 2)

    newImage, shift_list = li.lucky_imaging(images, globalRefImage, blk_size, nbest, binning)

    # Export the canvas to fits files. Rount to integers. Float is useless here.
    rImage      = np.rint(newImage)
    intImage    = np.int16(rImage)
    hdu         = fits.PrimaryHDU(intImage)
    fname       = outdir + 'python_processed/' + 'li_' + strIndex + '.fits'
    fits.writeto(fname, intImage, header)
