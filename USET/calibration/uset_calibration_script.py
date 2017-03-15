#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This script demonstrates the calibration pipeline for USET.

It produces two levels of calibration:

- Level 1.0 with only headers filled with the results of limb fitting and WCS-related keywords.
The latter can be used to register (co-align) the image with user's own method.

- Level 1.1 are rigidly centered but not rotated by P angle.
Rigidly = centered using integer number of pixels (rounded values of limb fitting results), and no interpolation.

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

# Set the directory where the FITS are. Here USETDATA is an environment variable. data_dir can be set to a string
# with a more explicit path if it is not possible to set an environment variable.
data_dir    =  os.path.join(os.environ.get('USETDATA'), 'UPH20161215111606.060_Short_Exp')
# Get the list of files, change it according to where your data files are and how are they are named.
file_list   = glob.glob(os.path.join(data_dir, '*.FTS'))
# Two output directories for two configuration levels (1.0 & 1.1), and 1 directory for jpeg image to check limb-fit.
# Beware of how to use os.path.join() to create paths of subdirectories. Do not use "/" or "\" at the
# Alternative: use explicit path. (e.g 'C:\USET\blah1\blah2 etc...) instead of os.path.join(...)
# Calibration level 1.0
outdir_10       = os.path.join(data_dir, 'calibrated_level1.0')
outdir_10_jpeg  = os.path.join(outdir_10, 'limb_fit_jpeg')
# Calibration level 1.1
outdir_11       = os.path.join(data_dir, 'calibrated_level1.1')
outdir_11_jpeg  = os.path.join(outdir_11, 'preview_jpeg')

# Check if output directories exist. Create if not.
if not os.path.isdir(outdir_10):
    os.makedirs(outdir_10)
if not os.path.isdir(outdir_10_jpeg):
    os.makedirs(outdir_10_jpeg)
if not os.path.isdir(outdir_11):
    os.makedirs(outdir_11)
if not os.path.isdir(outdir_11_jpeg):
    os.makedirs(outdir_11_jpeg)

# Enable/Disable agressive clean-up of limb points
limb_cleanup                    = True
# Toggle preview image of the results of limb fitting. Image does not display, but printed to jpeg file.
print_preview_limbfit_quadrants = True
# Toggle preview of limb fitting over full disk
print_preview_fullsun_level10   = True
# Toggle preview of centered image over full disk
print_preview_fullsun_level11   = True
# Write level 1.0 fits files
write_fits_level_10             = True
# Write level 1.1 fits files (solar disk aligned to center of FOV)
write_fits_level_11             = True

num_files = len(file_list)
# Loop through all the files
for i in range(0, num_files):

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

    if print_preview_limbfit_quadrants:
        fname_jpeg = uset.get_basename(file) + '_level1.0_quadrants.jpeg'
        fname = os.path.join(outdir_10_jpeg, fname_jpeg)
        # Plot and export preview of limb fitting results
        uset.export_limbfit_preview_quadrants(image, xc, yc, Rm, pts, pts3, fname)

    if print_preview_fullsun_level10:
        basename_jpg = uset.get_basename(file) + '_fullsun_level1.0.jpeg'
        fname = os.path.join(outdir_10_jpeg, basename_jpg)
        # Plot and export preview of limb fitting results
        uset.export_limbfit_preview_fullsun(image, xc, yc, Rm, pts, pts3, fname)

    if print_preview_fullsun_level11:
        # load a colormap

        # cmap=list([cm.get_cmap('hinodexrt'),
        #            cm.get_cmap('sdoaia1700'),
        #            cm.get_cmap('irissji1400'),
        #            cm.get_cmap('irissji1600')])

        cmap = cm.get_cmap('irissjiFUV')
        basename_jpg = uset.get_basename(file) + '_fullsun_level1.1.jpeg'
        fname = os.path.join(outdir_11_jpeg, basename_jpg)
        # Plot and export preview of limb fitting results
        uset.export_preview(image, fname, cmap)


print('done')