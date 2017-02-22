#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This script demonstrates the calibration pipeline for USET.
This script produces also jpeg for a quick check of the limb fitting, can be disabled for faster execution.
The actual functions doing the processing are in the "uset_calibration" module.
Refer to it for more detailed documentation of what each function do.

@author: Raphael Attie
"""

import os
import glob
from astropy.io import fits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import uset_calibration as uset


# Disable interactive mode so figure goes directory to a file and does not show up on screen
plt.ioff()

# Get the list of files, change it according to where your data files are
file_list = glob.glob("/Volumes/SDobo-A/Raphael/USET/USET/H_Alpha/UPH20161215_short_exp/*.FTS")
# Output parent directory
outdir = '/Volumes/SDobo-A/Raphael/USET/USET/H_Alpha/UPH20161215_short_exp'
# Two directories for two configuration levels (1.0 & 1.1), and 1 directory for jpeg image to check limb-fit.
# Calibration level 1.0
outdir_10       = outdir + '/calibrated_level1.0/'
outdir_10_jpeg  = outdir_10 + '/limb_fit_jpeg/'
# Calibration level 1.1
outdir_11       = outdir + '/calibrated_level1.1/'

# Check if output directories exist. Create if not.
if not os.path.isdir(outdir_10):
    os.makedirs(outdir_10)
if not os.path.isdir(outdir_10_jpeg):
    os.makedirs(outdir_10_jpeg)
if not os.path.isdir(outdir_11):
    os.makedirs(outdir_11)

# Enable/Disable agressive clean-up of limb points
limb_cleanup                = True
# Print preview image files of the results of limb fitting
print_preview_fullsun       = True
print_preview_quadrants     = True
# Write level 1.0 fits files
write_fits_level_10         = True
# Write level 1.1 fits files (solar disk aligned to center of FOV)
write_fits_level_11         = True

num_files = len(file_list)
# Loop through all the files
for i in range(0, 2):

    file = file_list[i]
    hdu = fits.open(file, ignore_missing_end=True)
    print(file)

    # Load header and image data from the 1st data unit: hdu[0]
    header = hdu[0].header
    image = hdu[0].data

    centeredImage, xc, yc, Rm, pts, pts2, pts3 = uset.uset_limbfit_align(image, limb_cleanup)
    #taubin = uset.uset_align(image, 'Taubin')

    if write_fits_level_10:
        # Update header with information from limb fitting using WCS.
        # Level 1.0 does not align the image to center of FOV. Level > 1.0 probably would.
        # Here we export level 1.0. With CRPIX1 = xc , CRPIX2 = yc, CRVAL1 = 0, CRVAL2 = 0
        header_10 = uset.set_header_calibrated(header, xc, yc, Rm)
        # Export the original image and new header to a FITS file. Level 1.0
        fname_fits = uset.get_basename(file) + '_level1.0.fits'
        fname = os.path.join(outdir_10, fname_fits)
        uset.write_uset_fits(image, header_10, fname)

    if write_fits_level_11:
        # Update header with information from limb fitting using WCS.
        # Here Level 1.0 does not align the image to center of FOV.
        # Export level 1.0. With CRPIX1 = xc , CRPIX2 = yc, CRVAL1 = 0, CRVAL2 = 0
        header_11 = uset.set_header_calibrated(header, image.shape[1]/2, image.shape[0]/2, Rm)
        # Export the original image and new header to a FITS file. Level 1.0
        fname_fits = uset.get_basename(file) + '_level1.1.fits'
        fname = os.path.join(outdir_11, fname_fits)
        uset.write_uset_fits(centeredImage, header_11, fname)

    if print_preview_quadrants:
        fname_png = uset.get_basename(file) + '_level1.0_quadrants.jpeg'
        fname = os.path.join(outdir_10_jpeg, fname_png)
        # Plot and export preview of limb fitting results
        uset.export_limbfit_preview_quadrants(image, xc, yc, Rm, pts, pts3, fname)

    if print_preview_fullsun:
        fname_png = uset.get_basename(file) + '_level1.0_fullsun.jpeg'
        fname = os.path.join(outdir_10_jpeg, fname_png)
        # Plot and export preview of limb fitting results
        uset.export_limbfit_preview_fullsun(image, xc, yc, Rm, pts, pts3, fname)

print('done')