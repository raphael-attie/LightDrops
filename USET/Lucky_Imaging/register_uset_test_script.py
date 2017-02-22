#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 15:13:12 2016

@author: raphaela
"""

import os
import glob
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import uset_calibration as uset

# Disable interactive mode so figure goes directory to a file and does not show up on screen
plt.ioff()

# Get the list of files
file_list   = glob.glob("/Users/attie/Dropbox/USET/Data/HA/*.FTS")
outdir_10  = os.path.dirname(file_list[0]) + '/calibrated_level1.0/'
outdir_11  = os.path.dirname(file_list[0]) + '/calibrated_level1.1/'

# Check if output directories exist. Create if not.
if not os.path.isdir(outdir_10):
    os.makedirs(outdir_10)
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

# Loop through all the files
for i in range(0, 1):

    file = file_list[i]
    hdu = fits.open(file)

    # Load header and image data from the 1st data unit: hdu[0]
    header = hdu[0].header
    image = hdu[0].data

    # Limb-fit and alignment
    image2 = np.zeros([2048, 2048])
    x1 = 512

    image2[:, x1:] = image[:, 0:(2048-x1)]
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
        # Level 1.0 does not align the image to center of FOV. Level > 1.0 probably would.
        # Here we export level 1.0. With CRPIX1 = xc , CRPIX2 = yc, CRVAL1 = 0, CRVAL2 = 0
        header_11 = uset.set_header_calibrated(header, image.shape[1]/2, image.shape[0]/2, Rm)
        # Export the original image and new header to a FITS file. Level 1.0
        fname_fits = uset.get_basename(file) + '_level1.1.fits'
        fname = os.path.join(outdir_11, fname_fits)
        uset.write_uset_fits(centeredImage, header_11, fname)

    if print_preview_quadrants:
        fname_png = uset.get_basename(file) + '_level1.0_quadrants.jpeg'
        fname = os.path.join(outdir_10, fname_png)
        # Plot and export preview of limb fitting results
        uset.export_limbfit_preview_quadrants(image, xc, yc, Rm, pts, pts3, fname)

    if print_preview_fullsun:
        fname_png = uset.get_basename(file) + '_level1.0_fullsun.jpeg'
        fname = os.path.join(outdir_10, fname_png)
        # Plot and export preview of limb fitting results
        uset.export_limbfit_preview_fullsun(image, xc, yc, Rm, pts, pts3, fname)

print('done')