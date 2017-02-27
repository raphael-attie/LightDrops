#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Script testing the lucky imaging. Assumes already co-aligned images

@author: raphaela
"""

import os
import sys
sys.path.append('/Users/rattie/Dev/LightDrops/USET/calibration')
import matplotlib
matplotlib.use('Agg')
import uset_calibration as uset
import glob
from astropy.io import fits
import numpy as np
import scipy
import sunpy.cm as cm
import ra_lucky_imaging as li



# Get the list of files
data_dir    = '/Volumes/SDobo-A/Raphael/USET/H_Alpha/UPH20161215_short_exp/calibrated_level1.1/'
files       = glob.glob(os.path.join(data_dir, '*.fits'))
outdir      = os.path.join(data_dir, 'lucky')
outdir_jpeg = os.path.join(outdir , 'jpeg')
# Check if output directory exists. Create if not.
if not os.path.isdir(outdir):
    os.makedirs(outdir)
if not os.path.isdir(outdir_jpeg):
    os.makedirs(outdir_jpeg)

print_preview_fullsun = True

# Number of images to consider in the list of files
nImages = 80
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

new_image, shifts = li.lucky_imaging(images, globalRefImage, blk_size, nbest, binning)

# Export the canvas to fits files. Rount to integers. Float is useless here.
rImage      = np.rint(new_image)
intImage    = np.int16(rImage)
hdu         = fits.PrimaryHDU(intImage)
fname       = outdir + '/python_lucky_NEW_' + strIndex + '.fits'
uset.write_uset_fits(intImage, header, fname)

if print_preview_fullsun:
    # load a colormap
    cmap = cm.get_cmap('irissjiFUV')
    basename_jpg = uset.get_basename(fname) + '_lucky.jpeg'
    fname = os.path.join(outdir_jpeg, basename_jpg)
    # Plot and export preview of limb fitting results
    #uset.export_preview(new_image, fname, cmap)
    scipy.misc.imsave(fname, new_image)
    # Sample from original
    sample = images[:, :, 0]
    basename_jpg = 'sample.jpeg'
    fname = os.path.join(outdir_jpeg, basename_jpg)
    scipy.misc.imsave(fname, sample)


