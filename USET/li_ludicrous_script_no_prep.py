#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Script testing the lucky imaging. It assumes already co-aligned images

@author: Raphael Attie
"""

import os
import sys
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import glob
import Lucky_Imaging.ra_lucky_imaging as li


# Disable interactive mode so figure goes directory to a file and does not show up on screen
plt.ioff()

# Set the directory where the FITS are. Here USETDATA is an environment variable of a parent directory of all USET data
# In my case it's either D:\USET\DATA on my windows PC or /Users/blah/USET/DATA on my Mac.
# data_dir can be set to a string with a more explicit path if it is not possible to set an environment variable.
data_dir    =  os.path.join(os.environ.get('USETDATA'), 'HAlpha', 'UPH20161215_Short_Exp', 'calibrated_level1.1')
# path to master dark to subtract to the data
master_dark_path = '/home/uset/calibration_HAlpha/darks/master_dark.fits'
# Set the output directories (calibrated FITS, jpeg, ...)
outdir      = os.path.join(data_dir, 'lucky')
outdir_jpeg = os.path.join(outdir , 'jpeg')
# Get the list of files
files       = glob.glob(os.path.join(data_dir, '*.fits'))
# Check if output directory exists. Create if not.
if not os.path.isdir(outdir):
    os.makedirs(outdir)
if not os.path.isdir(outdir_jpeg):
    os.makedirs(outdir_jpeg)

print_preview_fullsun = False

# [IGNORED FOR NOW] Number of images to consider in the list of files. Total lucky images = nImages/interval
### IGNORED AT THE MOMENT. THIS SCRIPT ONLY USES "INTERVAL" FOR THE NUMBER OF IMAGES TO CONSIDER
nImages = 20

## Define parameters for the block processing

# Image binning
# Binning for the quality matrix. Do not change for USET unless you know what you are doing.
binning         = 2
# Buffer space: extra marging around a reference block when aligning a template block
# bufferSpace   = 8
# Block size and binned block size
blk_size        = 64
# Nb of best images to keep the stack of subfields
nbest           = 5
# Quality metric
qmetric = 'Laplace'
# blend mode: 'aavg' (arithmetic average), 'gblend' (distance-transform + gaussian).
blend_mode = 'gblend'
# Nb of images out of which the nbest subfields are taken. Typically, nbest << interval
interval = 40

print('Running lucky imaging...')
time1 = time.time()
shifts = li.lucky_imaging_wrapper(files, outdir, outdir_jpeg, nImages,
                                  interval, nbest, binning, blk_size, qmetric, blend_mode,
                                  print_preview_fullsun)
time2 = time.time()
print('Total elapsed time: %0.1f s' % (time2-time1))


## nb of images between each final processed image = nb of images that get sorted in the quality matrix
# intervals        = [20, 60]

# for interval in intervals:
#
#     blend_mode = 'aavg'
#     li.lucky_imaging_wrapper(files, outdir, outdir_jpeg, nImages,
#                             interval, nbest, binning, blk_size, blend_mode,
#                             print_preview_fullsun)
    # blend_mode = 'gblend'
    # li.lucky_imaging_wrapper(files, outdir, outdir_jpeg, nImages,
    #                          interval, nbest, binning, blk_size, blend_mode,
    #                          print_preview_fullsun)

# Change the number of best images

# nbest           = 10
# for interval in intervals:
#
#     blend_mode = 'aavg'
#     li.lucky_imaging_wrapper(files, outdir, outdir_jpeg, nImages,
#                             interval, nbest, binning, blk_size, blend_mode,
#                             print_preview_fullsun)
#     blend_mode = 'gblend'
#     li.lucky_imaging_wrapper(files, outdir, outdir_jpeg, nImages,
#                              interval, nbest, binning, blk_size, blend_mode,
#                              print_preview_fullsun)
print('done')
