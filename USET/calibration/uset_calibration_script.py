#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This script demonstrates the calibration pipeline for USET.
It produces two levels of calibration:

- Level 0: only headers filled with the results of limb fitting and WCS-related keywords.
The latter can be used to register (co-align) the image with user's own method.

- Level 1: center of disk into center of FOV, and optionnally rotate by P angle.

This script produces also jpeg for a quick check of the limb fitting, can be disabled for faster execution.
The actual functions doing the processing are in the "uset_calibration" module.
Refer to it for more detailed documentation of what each function does.

Use this script as model/demo, not actual calibration product.


@author: Raphael Attie. raphael.attie@nasa.gov
"""

import os
import glob
from astropy.io import fits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import calibration.uset_calibration as uset
import sunpy.cm as cm


# Disable interactive mode so figure goes directory to a file and does not show up on screen
plt.ioff()

# Set the directory where the FITS are.
data_dir    =  '/Users/rattie/Data/USET/campaign/calibration/HALPHA'
# Get the list of files, change it according to where your data files are and how are they are named.
file_list   = glob.glob(os.path.join(data_dir, '*.FTS'))
# Two output directories for two configuration levels (1.0 & 1.1), and 1 directory for jpeg image to check limb-fit.
# Beware of how to use os.path.join() to create paths of subdirectories. Do not use "/" or "\" at the
# Alternative: use explicit path. (e.g 'C:\USET\blah1\blah2 etc...) instead of os.path.join(...)
# Calibration level 1.0
outdir_0       = os.path.join(data_dir, 'calibrated_lev0')
outdir_0_jpeg  = os.path.join(outdir_0, 'limb_fit_jpeg')
# Calibration level 1.1
outdir_1       = os.path.join(data_dir, 'calibrated_lev1')
outdir_1_jpeg  = os.path.join(outdir_1, 'preview_jpeg')

# Enable / Disable agressive clean-up of limb points
limb_cleanup                    = True
# Enable / disable rotation to align images to solar north
rotate_to_solar_north           = True
# Use "RICE" compression of FITS file
compressed_fits                 = False
# Write level 1.0 fits files
write_fits_lev_0                = True
# Write level 1.1 fits files (solar disk aligned to center of FOV)
write_fits_lev_1                = True
# Toggle preview image of the results of limb fitting. Image does not display, but printed to jpeg file.
print_preview_limbfit_quadrants = False
# Toggle preview of limb fitting over full disk
print_preview_fullsun_lev_0   = False
# Toggle preview of centered image over full disk
print_preview_fullsun_lev_1   = False

# Check if output directories exist. Create if not.
if not os.path.isdir(outdir_0):
    os.makedirs(outdir_0)
if not os.path.isdir(outdir_0_jpeg):
    os.makedirs(outdir_0_jpeg)
if not os.path.isdir(outdir_1):
    os.makedirs(outdir_1)
if not os.path.isdir(outdir_1_jpeg):
    os.makedirs(outdir_1_jpeg)

num_files = len(file_list)
# Loop through all the files
# for i in range(0, num_files):
for i in range(0, 1):

    file = file_list[i]
    hdu = fits.open(file, ignore_missing_end=True)
    print(file)

    # Load header and image data from the 1st data unit: hdu[0]
    header = hdu[0].header
    image = hdu[0].data

    # TODO: subtract master dark
    centered_image, xc, yc, Rm, pts, pts2, pts3 = uset.uset_limbfit_align(image, limb_cleanup)


    if write_fits_lev_0:
        # Update header with information from limb fitting using WCS.
        # Level 0 does not align the image to center of FOV.
        # Here we export level 1.0. With CRPIX1 = xc , CRPIX2 = yc, CRVAL1 = 0, CRVAL2 = 0
        header_0 = uset.set_header_calibrated(header, xc, yc, Rm)
        # Export the original image and new header to a FITS file. Level 1.0
        fname_fits = uset.get_basename(file) + '_lev0.fits'
        fname = os.path.join(outdir_0, fname_fits)
        uset.write_uset_fits(image, header_0, fname, compressed=compressed_fits)

    if write_fits_lev_1:
        # Level 1 aligns the image to center of FOV with rotation to align y-axis to solar north
        # Export level 1.0. With CRPIX1 = xc , CRPIX2 = yc, CRVAL1 = 0, CRVAL2 = 0
        # Coordinate of disk center in centered image
        xc2 = image.shape[0]/2 - 0.5
        yc2 = image.shape[1]/2 - 0.5
        header_1 = uset.set_header_calibrated(header, xc2 , yc2, Rm)
        if rotate_to_solar_north:
            centered_image = uset.align_to_solar_north(centered_image, header_1)

        # Export the aligned image and new header to a FITS file.
        fname_fits = uset.get_basename(file) + '_lev1.fits'
        fname = os.path.join(outdir_1, fname_fits)
        uset.write_uset_fits(centered_image, header_1, fname, compressed=compressed_fits)

    if print_preview_limbfit_quadrants:
        fname_jpeg = uset.get_basename(file) + '_lev1_quadrants.jpeg'
        fname = os.path.join(outdir_0_jpeg, fname_jpeg)
        # Plot and export preview of limb fitting results
        uset.export_limbfit_preview_quadrants(image, xc, yc, Rm, pts, pts3, fname)

    if print_preview_fullsun_lev_0:
        basename_jpg = uset.get_basename(file) + '_fullsun_lev0.jpeg'
        fname = os.path.join(outdir_0_jpeg, basename_jpg)
        # Plot and export preview of limb fitting results
        uset.export_limbfit_preview_fullsun(image, xc, yc, Rm, pts, pts3, fname)

    if print_preview_fullsun_lev_1:
        # load a colormap

        # cmap=list([cm.get_cmap('hinodexrt'),
        #            cm.get_cmap('sdoaia1700'),
        #            cm.get_cmap('irissji1400'),
        #            cm.get_cmap('irissji1600')])

        cmap = cm.get_cmap('irissjiFUV')
        basename_jpg = uset.get_basename(file) + '_fullsun_lev1.jpeg'
        fname = os.path.join(outdir_1_jpeg, basename_jpg)
        # Plot and export preview of limb fitting results
        uset.export_preview(image, fname, cmap)




print('done')