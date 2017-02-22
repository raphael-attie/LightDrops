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
nImages = 20
# Get the list of files    
files   = glob.glob("/Users/raphaela/Data/USET/HAlpha/ludicrous/s2/aligned_normalized/selection_1fps/*.fits")

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
nbest           = 5
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

outdir          = '/Users/raphaela/Data/USET/HAlpha/ludicrous/s2/aligned_normalized/selection_1fps/python/'

i = 0

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
    images[:, :, ii] = hdu[0].data

globalRefImage  = np.median(images, 2)


# Investigate a chosen block. Give the coordinate of its bottom left corner:
x = 1200
y = 1200

qBinnedArrays = li.block_processing_setup(images, binning)


#stack_block, shift, best_indices, xs, ys = li.lucky_imaging_single_block(images, globalRefImage, blk_size, nbest, binning, x, y)
