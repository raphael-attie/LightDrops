# -*- coding: utf-8 -*-
"""
This module contains all the functions related to the alignment of USET images using limb-fitting.
- uset_limbfit_align() is a convenience function that applies the limb-fitting and align the images.
  It outputs the centered image, sun center coordinates, radius and other output of the limb-fitting algorithm.
  The image is centered such that the limb-fitted disk center is at the middle of the image axes (naxis*/2)
- set_header_calibrated() updates the header with WCS keywords
- See the documentation on "detectLimbPoints()" for more information on how the limb points are detected


@author: Raphael Attie (attie.raphael@gmail.com)
"""
import glob
import os
import shutil
from astropy.io import fits
import numpy as np
from astropy.convolution import convolve, Box1DKernel
from sunpy.sun import solar_north
import sunpy.cm as cm
from scipy import optimize
from skimage.transform import warp
from skimage.transform import SimilarityTransform
from skimage.transform import rotate
import matplotlib.pyplot as plt
import cv2

def calibrate(data_dir, output_dir, fits_extension, dark_path, limb_cleanup, level, preview):
    """
    produces two levels of calibration:

    :param data_dir: directory containing the FITS files to calibrate (FITS, FITS.gz, FTS)
    :param output_dir: directory where the calibrated FITS files are written
    :param extension: file extension of the input fits files (fits, FITS.gz, or FTS) in data_dir
    :param dark_path: path to master dark file to subtract
    :param limb_cleanpup: toggle aggressive limb detection (True / False)
    :param level: calibration level.
    - Level 0: only headers filled with the results of limb fitting and WCS-related keywords.
    The latter can be used to register (co-align) the image with user's own method.
    - Level 1: center of disk into center of FOV, and optionnally rotate by P angle.
    :param preview: print a jpeg preview of the calibrated image
    :return: None. Data are written to disk
    """
    # Get the list of files, change it according to where your data files are and how are they are named.
    file_list = glob.glob(os.path.join(data_dir, '*.' + fits_extension))

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    # Number of frames to process
    num_files = len(file_list)

    for i in range(0, num_files):
        file = file_list[i]
        hdu = fits.open(file, ignore_missing_end=True)
        # Load header and image data from the 1st data unit: hdu[0]
        header = hdu[0].header
        image = hdu[0].data
        naxis1 = image.shape[0]
        naxis2 = image.shape[1]

        if dark_path:
            hdu_dark = fits.open(dark_path)
            image = image - hdu_dark[0].data

        centered_image, xc, yc, Rm, pts, pts2, pts3 = uset_limbfit_align(image, limb_cleanup)
        # Coordinate of disk center in centered image
        xc2 = naxis1 / 2 - 0.5
        yc2 = naxis2 / 2 - 0.5

        if level == 0:
            # Update header with information from limb fitting using WCS.
            # Level 0 does not align the image to center of FOV.
            # Here we export level 1.0. With CRPIX1 = xc , CRPIX2 = yc, CRVAL1 = 0, CRVAL2 = 0
            header = set_header_calibrated(header, xc, yc, Rm)
            # Export the calibrated image and new header into a FITS file.
            fname_fits = get_basename(file) + '_level_' + str(level) + '.fits'
            fname = os.path.join(output_dir, fname_fits)
            write_uset_fits(image, header, fname, compressed=False)

        if level == 1:
            # Level 1 aligns the image to center of FOV with rotation to align y-axis to solar north
            # Export level 1.0. With CRPIX1 = xc2 , CRPIX2 = yc2, CRVAL1 = 0, CRVAL2 = 0

            header = set_header_calibrated(header, xc2, yc2, Rm)
            # Export the calibrated image and new header into a FITS file.
            fname_fits = get_basename(file) + '_level_' + str(level) + '.fits'
            fname = os.path.join(output_dir, fname_fits)
            write_uset_fits(centered_image, header, fname, compressed=False)


        if preview:
            # Apply radial filter to boost off-limb prominences

            # Define grid of coordinates
            xgrid1, ygrid1 = np.meshgrid(np.arange(naxis1), np.arange(naxis2))
            # Distance of each pixel to disk center
            pixel_r = np.sqrt((xgrid1 - xc2) ** 2 + (ygrid1 - yc2) ** 2)

            # apply a radius-based scaling to boost the limb
            radius = header['SOLAR_R']
            # Mask of off-limb pixels
            pixels_off_limb = pixel_r > radius
            # Copy image into new array
            centered_image2 = np.array(centered_image)

            # Apply mask and enhance intensity beyond limb
            centered_image2[pixels_off_limb] = centered_image2[pixels_off_limb] * 2
            # Get the maximum intensity for the rescaled image based on 99.99% percentile
            new_max = compute_intensity_high(centered_image, 99.999)
            # Apply rotation of solar rotation axis to the image y-axis (top-bottom axis)
            centered_image2 = align_to_solar_north(centered_image2, header)
            # load a colormap for the preview function (could be passed as argument to this function...)
            cmap = cm.get_cmap('irissjiFUV')
            # Define file names of png file
            basename_png = get_basename(file) + '_level_1.png'
            fname = os.path.join(output_dir, basename_png)
            # Plot and export preview of limb fitting results
            export_preview(centered_image2, fname, new_max, cmap)

    return



def sort_calibration_files(data_dir, extension):
    """
    Sort calibration files by moving them to directories depending on their nature: light, dark or bias image

    :param data_dir: directory where all files are
    :param extension: file extension. Can be 'FTS', 'fits', 'FTS.gz', 'fits.gz'
    :return: none.
    """
    file_list = glob.glob(os.path.join(data_dir, '*.' + extension))
    # List of tested gain values
    gains = [1000, 2000]

    gain_dirs = [os.path.join(data_dir, 'Gain' + str(gain)) for gain in gains]
    lights_dirs = [os.path.join(gdir, 'lights') for gdir in gain_dirs]
    darks_dirs = [os.path.join(gdir, 'darks') for gdir in gain_dirs]
    bias_dirs = [os.path.join(gdir, 'bias') for gdir in gain_dirs]

    for i in range(0, len(gains)):
        # Check if output directories exist. Create if not.
        if not os.path.isdir(lights_dirs[i]):
            os.makedirs(lights_dirs[i])
        if not os.path.isdir(bias_dirs[i]):
            os.makedirs(bias_dirs[i])
        if not os.path.isdir(darks_dirs[i]):
            os.makedirs(darks_dirs[i])

    for i in range(0, len(file_list)):

        file = file_list[i]
        hdu = fits.open(file, ignore_missing_end=True)
        # Load header and image data from the 1st data unit: hdu[0]
        h = hdu[0].header
        img = hdu[0].data
        # Get gain index of current image
        g = gains.index(h['GAIN'])
        # Move file to directory for according to current gain value, as bias, dark or lights
        if h['EXPTIME'] < 1E-4:
            # Consider image as bias
            shutil.move(file, bias_dirs[g])
        elif img.mean() < 50:
            # Consider image as dark
            shutil.move(file, darks_dirs[g])
        else:
            # if none of the above, consider image as light
            shutil.move(file, lights_dirs[g])
    return


def create_master_frame(calibration_dir, extension, master_path, nframes):
    """
    Take the median of a series of images to produce a master calibration frame (master dark, bias, ...)

    :param calibration_dir: path to calibration images (fits files)
    :param extension: file extension. Can be 'FTS', 'fits', 'FTS.gz', 'fits.gz'
    :return: Master image. Files are printed in dark_dir.
    """
    file_list = glob.glob(os.path.join(calibration_dir, '*.' + extension))
    hdu = fits.open(file_list[0], ignore_missing_end=True)
    header = hdu[0].header
    naxis1 = header['NAXIS1']
    naxis2 = header['NAXIS2']
    # If the number of files in the directory is smaller than the user input number of frames,
    # need to take the latter.
    nframes = np.min([len(file_list), nframes])
    frames = np.zeros([naxis1, naxis2, nframes])

    print('Creating master frame from %d frames...' %nframes)

    for i in range(0, nframes):
        hdu = fits.open(file_list[i], ignore_missing_end=True)
        frames[:,:,i] = hdu[0].data

    master_frame = np.median(frames, 2)

    # Write to disk if path is set for the master frame
    if master_path:
        write_uset_fits(master_frame, header, master_path)
        print('Wrote master into %s ' %master_path)

    return master_frame



def detectLimbPoints(img, numPts, smooth, limb_cleanup):
    """
    Detect points around the limb using the maximum of 1D gradients in slices taken through the disk.
    The limb points are meant to be fitted using Least-square circle fitting.

    The algorithm separates the image in 4 quadrants meant to contain East-, West-, North-, and Souht-limb.

    An optional "limb cleanup" is implemented:
    An image-dependent threshold is calculated for each image. It sets the minimum gradient to consider when detecting
    limb points. If the gradient of the points that are supposed to be at the limb are below that threshold, they are
    rejected. This "limb cleanup" can be enabled/disabled (limb_cleanup = True/False) for less agressive point detection,
    it might be necessary to disable it with cloudy or foggy images whose statistical properties would be too different from
    a good-seeing image.

    Using the above heuristic,  this method tolerates some off-pointing and thus can still fit the limb if the disk
    partially lies outside the FOV: cropping can be up to 25% in each direction. More than that will impact accuracy,
    but it can still be used for rough estimate of pointing corrections of the equatorial table.


    :param img: input solar image
    :param numPts: number of points per quadrant
    :param smooth: width of the boxcar average (1D)
    :param limb_cleanup: True or False. Enable or disable the cleanup of limb points.
    :return: series of x, y coordinates of points lying on the solar limb. [nub of points, coordinates, max gradient]
    """


    naxis1 = img.shape[1]
    naxis2 = img.shape[0]

    # Image is divided in 4 quadrants with a square size of length = naxis_q
    naxis_q = int(naxis1/2)

    # Define array to store the limb points
    limbPts = np.zeros([4 * numPts, 3])

    # The limb points are extracted slices within 1/4 and 3/4 each image axis.
    for ii in range(0, numPts):

        # Coordinates used to extract slices of the image that contain a  part of the limb. 
        x = int(naxis1 / 4 + ii * naxis1 / (2 * numPts))
        y = int(naxis2 / 4 + ii * naxis2 / (2 * numPts))

        ## Extract horizontal and vertical slice 
        sliceH = np.squeeze(img[y: y + 1, :])
        sliceV = np.squeeze(img[:, x: x + 1])

        if smooth > 0:
            # Smooth the slice.
            sliceH = convolve(sliceH, Box1DKernel(smooth), boundary='extend')
            sliceV = convolve(sliceV, Box1DKernel(smooth), boundary='extend')

        # Central derivative of the slices

        # Along X-axis
        gSliceH = np.abs(np.gradient(sliceH))
        # Slice at East side (left) of limb: stops at 1/4 of the axis 
        gSliceLeft = gSliceH[0: naxis_q]
        # Slice at West side (right) of limb: starts at 3/4 of the axis
        gSliceRight = gSliceH[(naxis1 - naxis_q): naxis1]

        # Along Y-axis
        gSliceV = np.abs(np.gradient(sliceV))
        # at North side (Top) of limb: stops at 1/4 of the axis
        gSliceTop = gSliceV[0: naxis_q]
        # at Bottom side (Bottom) of limb: starts at 3/4 of the axis
        gSliceBot = gSliceV[(naxis2 - naxis_q): naxis2]

        # Get the location of the maximum in those slices
        # and convert that location in the reference frame of the image

        # East limb
        limbPts[ii, 0] = np.argmax(gSliceLeft)
        limbPts[ii, 1] = y
        limbPts[ii, 2] = np.max(gSliceLeft)
        # West limb
        limbPts[ii + numPts, 0] = np.argmax(gSliceRight) + naxis1 - naxis_q
        limbPts[ii + numPts, 1] = y
        limbPts[ii + numPts, 2] = np.max(gSliceRight)
        # North limb
        limbPts[ii + 2 * numPts, 0] = x
        limbPts[ii + 2 * numPts, 1] = np.argmax(gSliceTop)
        limbPts[ii + 2 * numPts, 2] = np.max(gSliceTop)
        # South limb
        limbPts[ii + 3 * numPts, 0] = x
        limbPts[ii + 3 * numPts, 1] = np.argmax(gSliceBot) + naxis2 - naxis_q
        limbPts[ii + 3 * numPts, 2] = np.max(gSliceBot)

    # Filter the above points to keep only the ones with a gradient that is at least half the standard deviation
    # of the disk intensity.
    # For this we consider only the intensity of the pixels that are on the disk

    if (limb_cleanup):
        # Reject pixels that would be too much inside the disk.

        # Find background intensity based on the 20%-percentile of the image histogram
        # This is a safe value of off-disk intensity
        bkg_intensity = compute_intensity_percentile(img, 20)
        # Calculate the standard deviation of the background.
        bkg_std = np.std(img[img < bkg_intensity])
        # We consider that pixels on the disk have intensity > at least the background level + 20x the background standard dev.
        disk_threshold = bkg_intensity + 20 * bkg_std
        # Take the mean of the intensity above that
        disk_mean = np.mean(img[img > disk_threshold])
        disk_std = np.std(img[img > disk_mean])
        # The minimum gradient that we consider for a valid limb point is set at half that standard deviation
        # This can be a bit harsh if you need to recover cloudy images. Comment it out if too much valid points are rejected
        grad_min = disk_std / 2.0
        keep_mask = limbPts[:, 2] >= grad_min
        limbPts = limbPts[keep_mask, :]

    return limbPts


def calc_R(x, y, xc, yc):
    """ calculate the distance of each 2D points from the center (xc, yc) """
    return np.sqrt((x - xc) ** 2 + (y - yc) ** 2)


def f_2(c, x, y):
    """ calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc) """
    Ri = calc_R(x, y, *c)
    return Ri - Ri.mean()


def circle_fit(x, y, method):
    """
    Fit a circle with data points (x,y) using a 3-pass least-square fit with sigma-clipping
    Uses method 2 from http://scipy-cookbook.readthedocs.io/items/Least_Squares_Circle.html

    I added some sigma-clipping-based three-passes cleanup for more robustness,
    works well with turbulent or cloudy weather

    :param x: numpy 1D array of the x-coordinates of points to fit
    :param y: numpy 1D array of the y-coordinates of points to fit
    :return: coordinates and radius of fitted circle for each pass
    """

    # coordinates of the barycenter
    xm = np.mean(x)
    ym = np.mean(y)

    xc, yc, xc2, yc2, xc3, yc3, Rm, Rm2, Rm3 = 0, 0, 0, 0, 0, 0, 0, 0, 0
    R_dist, R_dist2, R_dist3 = 0, 0, 0

    if method == 'lsq':
        center_estimate = np.array([xm, ym])
        center = optimize.leastsq(f_2, center_estimate, (x, y))[0]
        xc, yc = center
        R_dist = calc_R(x, y, xc, yc)
        Rm = np.mean(R_dist)
    elif method == 'Taubin':
        xc, yc, Rm = circle_fit_by_Taubin(x, y)[0:3]
        R_dist = calc_R(x, y, xc, yc)

    ## 2nd pass - sigma-clipping of the data points:
    # Eliminate the data points whose distance to fitted-center are > 2*std
    Rmedian = np.median(R_dist)
    keep_mask = np.abs(R_dist - Rmedian) < R_dist.std()
    x2 = x[keep_mask]
    y2 = y[keep_mask]
    limbPts2 = np.array([x2, y2]).transpose()

    # coordinates of the barycenter
    xm = np.mean(x2)
    ym = np.mean(y2)

    if method == 'lsq':
        center_estimate = np.array([xm, ym])
        center = optimize.leastsq(f_2, center_estimate, (x2, y2))[0]
        xc2, yc2 = center
        R_dist2 = calc_R(x2, y2, xc2, yc2)
        Rm2 = np.mean(R_dist2)
    elif method == 'Taubin':
        xc2, yc2, Rm2 = circle_fit_by_Taubin(x2, y2)[0:3]
        R_dist2 = calc_R(x2, y2, xc2, yc2)

    ## 3rd pass
    Rmedian = np.median(R_dist2)
    keep_mask = np.abs(R_dist2 - Rmedian) < R_dist2.std()
    x3 = x2[keep_mask]
    y3 = y2[keep_mask]
    limbPts3 = np.array([x3, y3]).transpose()

    # coordinates of the barycenter
    xm = np.mean(x3)
    ym = np.mean(y3)

    if method == 'lsq':
        center_estimate = np.array([xm, ym])
        center = optimize.leastsq(f_2, center_estimate, (x3, y3))[0]
        xc3, yc3 = center
        R_dist3 = calc_R(x3, y3, xc3, yc3)
        Rm2 = np.mean(R_dist2)
    elif method == 'Taubin':
        xc3, yc3, Rm3 = circle_fit_by_Taubin(x3, y3)[0:3]
        R_dist3 = calc_R(x3, y3, xc3, yc3)

    Rm3 = np.mean(R_dist3)

    return xc3, yc3, Rm3, xc2, yc2, Rm2, xc, yc, Rm, limbPts2, limbPts3



def center_disk(image, xc, yc):
    """
    Shift the image so the solar disk center (xc, yc) is at the image center (axes size / 2).
    E.g. if image size is NAXIS1 x NAXIS2. Center of the disk will be at [NAXIS1 / 2, NAXIS2 / 2].
    Others may prefer to use [ (NAXIS1 + 1)/2 , (NAXIS2  + 1)/2 ].

    This does not use interpolation. Image is shifted "rigidly" using indexing & copying.
    To that end, actual center coordinates from limb fitting are round to the nearest integer.
    You can use the residue of the rounding for accuracy better than 0.5 pixel.

    :param image: input image that needs to be centered
    :param xc: x-coordinate of limb-fitted sun center
    :param yc: y-coordinate of limb-fitted sun center
    :return: Image centered to the middle of the axes.
    """
    naxis1 = image.shape[1]
    naxis2 = image.shape[0]
    x0 = naxis1 / 2 - 0.5
    y0 = naxis2 / 2 - 0.5
    dx = xc - x0
    dy = yc - y0

    # canvas = np.zeros([naxis2, naxis1])
    # # Size of cropped axis
    # newAxis1 = naxis1 - abs(dx)
    # newAxis2 = naxis2 - abs(dy)
    #
    # if (dx >= 0) and (dy >= 0):
    #     canvas[0:newAxis2, 0:newAxis1] = image[dy:, dx:]
    #
    # if (dx >= 0) and (dy < 0):
    #     canvas[abs(dy):, 0:newAxis1] = image[0:newAxis2, dx:]
    #
    # if (dx < 0) and (dy > 0):
    #     canvas[0:newAxis2, abs(dx):] = image[dy:, 0:newAxis1]
    #
    # if (dx < 0) and (dy < 0):
    #     canvas[abs(dy):, abs(dx):] = image[0:newAxis2, 0:newAxis1]

    # define the translation using skimage warp transform
    tform = SimilarityTransform(translation=(dx, dy))
    centered_image = warp(image.astype(np.float), tform, order=3)

    return centered_image


def uset_limbfit_align(image, limb_cleanup):
    """
    Wrapper around the circle-fitting algorithm to fit the limb of USET images.

    :param image: input image to fit. The solar disk may partially lies outside the FOV.
    :param limb_cleanup: Enable/Disable agressive clean-up of limb points
    :return: image centered at the middle of each axis (image.shape / 2), coordinates of disk center,
             radius of the third pad, limb points used at each of the three passes.
             See function circle_fit() for more information.
    """
    # Fix for USET images [not needed with Emil's update, where the faulty Suncap's images are now fixed]
    image[0, :] = image[1, :]
    image[:, 0] = image[:, 1]

    # Default method for limb fitting is least-square ('lsq') from python scipy.optimize package.
    # Alternate method is 'Taubin', as used in LightDrops in its original C++ form,
    # and at Kanzelehoe observatory. I ported here in Python.
    # Discrepancy between both methods is of the order of 1E-3 pixels or even less.
    method = 'lsq'

    # Number of points per limb quadrant. Totol number of points will be typically 4 * numPts if the disk is fully
    # encompossed in the FOV.
    numPts = 128
    smooth = 3

    limbPts = detectLimbPoints(image, numPts, smooth, limb_cleanup)

    limbX = np.squeeze(limbPts[:, 0])
    limbY = np.squeeze(limbPts[:, 1])
    # Fit the circle and output center coordinates and radius
    xc3, yc3, Rm3, xc2, yc2, Rm2, xc, yc, Rm, limbPts2, limbPts3 = circle_fit(limbX, limbY, method)

    # Create a new image with the solar disk centered in the image.
    centered_image = center_disk(image, xc3, yc3)

    return centered_image, xc3, yc3, Rm3, limbPts, limbPts2, limbPts3


def align_to_solar_north(image, header):
    """
    Rotate image by CROTA1 = -solar P angle (from the P B0 L0 solar ephemerids) so the Y-axis (top) aligns to solar north
    CROTA1 must be in the header.
    Rotation is around the center of the image: [naxis1, naxis2] / 2 - 0.5

    :param image: image to be rotated by CROTA1
    :param header: header of the image as produced by set_header_calibrated()
    :return: rotated image using bi-cubic interpolation
    """

    CROTA1 = header['CROTA1']
    rotated_image = rotate(image, CROTA1, order=3)

    return rotated_image


def compute_intensity_percentile(image, percentile):
    """
    Compute the image percentile-based intensity

    :param image: input image
    :param percentile: percentile of the histogram [e.g: 0, 1, ... , 99.5, ..., 100]
    :return: intensity at given percentile
    """
    nbins = 256
    im_hist, bin_edges = np.histogram(image, nbins)
    hist_width = bin_edges[1] - bin_edges[0]
    intensity = calc_threshold(percentile, im_hist, hist_width)

    return intensity

def compute_intensity_high(image, cutoff_high):
    """
    Compute the image percentile-based high intensity threshold up to which we stretch (rescale) the intensity.
    This is specific to USET and the percentile is hardcoded at 99.97%

    :param image: input image
    :param image: percentile. Typical values within [99.97 - 99.99]
    :return: maximum thresholding intensity
    """
    nbins = 256
    im_hist, bin_edges = np.histogram(image, nbins)
    hist_width = bin_edges[1] - bin_edges[0]
    intensity_high = calc_threshold(cutoff_high, im_hist, hist_width)

    return intensity_high


def calc_threshold(cutoff, im_hist, hist_width):
    """
    Calculate the intensity at the input cutoff-percentage of total pixels given an image histogram and bin width
    This can be used to calculate, for example, the min or max threshold for rescaling the image intensity

    :param cutoff: cutoff percentage of total pixels in the histogram
    :param im_hist: image histogram
    :param hist_width: bin width of the image histogram
    :return:
    """
    npixels = im_hist.sum()
    nbins = im_hist.shape[0]

    cdf = 0.0
    i = 0
    count = 0
    while i < nbins:
        cdf += im_hist[i]
        if 100.0 * cdf / npixels > cutoff:
            count = i
            break
        i += 1

    if i == nbins:
        count = nbins - 1

    intensity = hist_width * count

    return intensity


def rescale_image(image, new_min, new_max):
    """
    Rescale linearly an image between new_min and new_max.
    This is also called "contrast stretching".

    :param image:
    :param new_min: minimum intensity value
    :param new_max: maximum intensity value
    :return: rescaled image
    """
    data_max = np.max(image)
    data_min = np.min(image)
    data_range = data_max - data_min
    new_range = new_max - new_min
    alpha = new_range / data_range
    beta = -data_min * new_range / data_range
    rescaled_image = image * alpha + beta

    return rescaled_image


def rescale_image_by_histmax(image):
    """
    Recale linearly an image using the rescale_image() function using new_min = 0 and new_max based on the
    99.97% percentile of the image histogram. This percentile will saturate about the 1200 brightest pixels

    :param image: input image
    :return: rescaled image
    """

    new_min = 0
    new_max = compute_intensity_high(image, 99.97)
    rescaled_image = rescale_image(image, new_min, new_max)

    return rescaled_image


def set_header_calibrated(header, xc, yc, Rm):
    """Return new header with sun center and radius calculated with limb fitting"""

    new_header = header.copy()
    pixel_size = 1.09  # arcsec
    # DATE-OBS	The date of the observation, in the format specified in the FITS Standard.
    # Format: 'yyyy-mm-ddTHH:MM:SS[.sss]'.
    if len(header['DATE-OBS']) < 19: # quick and dirty way to check if the date-obs field only contains a date (before January 2017, USET produced this metadata)
        date = header['DATE-OBS']
        time = header['TIME']
        date_obs = date[6:10] + '-' + date[3:5] + '-' + date[0:2] + 'T' + time
    else:
        date_obs = header['DATE-OBS'] # We should be able to just use the field as it is (already correctly formatted)

    new_header['DATE-OBS'] = date_obs

    # Solar P angle: angle between geocentric north pole and solar rotational north pole
    # CCD vertical axis is assumed to be parallel to the geocentric north pole.
    solar_P = solar_north(date_obs)
    CROTA1 = -solar_P.degree
    new_header.append(('SOLAR_P', solar_P.degree, '[degrees]'))
    # Solar radius from limb fitting
    new_header.append(('SOLAR_R', Rm, '[px] Sun radius from limb fitting'))
    # Standard FITS keywords with WCS
    new_header.append(('WCSNAME', 'Helioprojective-cartesian'))
    new_header.append(('CRPIX1', xc, '[px] Sun center x-coordinate'))
    new_header.append(('CRPIX2', yc, '[px] Sun center y-coordinate'))
    new_header.append(('CRVAL1', 0))
    new_header.append(('CRVAL2', 0))
    new_header.append(('CDELT1', pixel_size, '[arcsec] plate scale'))
    new_header.append(('CDELT2', pixel_size, '[arcsec] plate scale'))
    new_header.append(('CROTA1', CROTA1, '[degrees]'))
    # WCS rotation PC matrix
    PC = np.array([[1, 0, 0],
                   [0, np.cos(CROTA1), -np.sin(CROTA1)],
                   [0, np.sin(CROTA1), np.cos(CROTA1)]])

    new_header.append(('PC1_1', PC[0, 0]))
    new_header.append(('PC2_2', PC[1, 1]))
    new_header.append(('PC2_3', PC[1, 2]))
    new_header.append(('PC3_2', PC[2, 1]))
    new_header.append(('PC3_3', PC[2, 2]))

    return new_header


def write_uset_fits(image, header, fname, compressed=False):
    """Export the USET image to a fits file. Intensity is rounded back to original type."""

    rimage = np.rint(image)
    int_image = np.int16(rimage)

    if not compressed:
        try:
            fits.writeto(fname, int_image, header=header, output_verify='exception', overwrite=True)
        except TypeError:
            fits.writeto(fname, int_image, header=header, output_verify='exception', checksum=True, clobber=True)
    else:
        chdu = fits.CompImageHDU(data=image, header=header, compression_type='RICE_1')
        chdu.writeto(fname, overwrite=True)


def get_basename(filepath):
    """
    return the basename (without file extension) of a file from its full path.

    :param filepath: full file path, with directory and file extension
    :return: file basename without exension

    :Example:
    basename('/Users/blah/blahblah/USET/Data/HA\\UCH20160311105740.FTS')
    'UCH20160311105740'

    """
    path_basename = os.path.splitext(filepath)[0]

    return os.path.basename(path_basename)


def export_preview(image, fname, new_max, colormap='gray'):
    """
    Print preview of the given image with axis labels into an image file

    :param image: image to be written to disk
    :param fname: full path of image file with extension (.png, .jpeg, ...)
    :param colormap: color map for plotting the image. Default is 'gray'.

    :return: figure written on disk
    """

    # Image dimension
    naxis1 = image.shape[1]
    naxis2 = image.shape[0]

    fig = plt.figure(1)
    fig.clear()
    # Plot image
    plt.imshow(image, vmin=0, vmax=new_max, cmap = colormap, origin='lower')
    current_frame = plt.gca()
    current_frame.axes.get_xaxis().set_visible(False)
    current_frame.axes.get_yaxis().set_visible(False)

    plt.savefig(fname, dpi=180, bbox_inches='tight')


def export_limbfit_preview_fullsun(image, xc, yc, Rm, pts, pts3, fname, overlay_pts=False):
    """
    Print preview of the results of limb fitting into an image file. No screen display.

    :param image: limb-fitted image
    :param xc: x-coordinate of disk center [pixels]
    :param yc: y-coordinate of disk center [pixels]
    :param Rm: disk radius (mean)
    :param pts: 1st pass series of points coordinates that were fitted. Dimensions: [# of points, 2]
    :param pts3: 3rd pass series of points coordinates that were fitted. Dimensions: [# of points, 2]
    :param fname: full path of image file with extension (.png, .jpeg, ...)
    :param overlay_pts: toggle overplotting the detected limb points (use supplied pts and pts3).

    :return: figure written on disk
    """

    # Color of symbols in pyplot.plot sym1 for the 1st pass. sym2 for 3rd pass
    sym1 = 'co'
    sym3 = 'ro'
    # Color of limb-fitted circle
    ccolor = 'y'

    # Image dimension
    naxis1 = image.shape[1]
    naxis2 = image.shape[0]

    fig = plt.figure(1)
    fig.clear()
    # Plot image
    plt.imshow(image, cmap = 'gray', origin='lower')
    if overlay_pts:
        # overlay limb points, after 1st pass (sym1) and 3rd pass (sym2)
        for i in range(0, pts.shape[0]):
           plt.plot(pts[i,0], pts[i,1], sym1, ms=2)
        for i in range(0, pts3.shape[0]):
           plt.plot(pts3[i,0], pts3[i,1], sym3, ms=2)

    ax = plt.gca()
    plt.xlabel('X [px]')
    plt.ylabel('Y [px]')
    # Draw the fitted circle
    circle1 = plt.Circle((xc, yc), Rm, color=ccolor, fill=False, linewidth=1)
    ax.add_artist(circle1)
    plt.axis([0, naxis1, 0, naxis2])

    plt.savefig(fname, dpi=180, bbox_inches='tight')

def export_limbfit_preview_quadrants(image, xc, yc, Rm, pts, pts3, fname):
    """
    Plot and save preview of the results of limb fitting. Show close-ups on each of the 4 quadrants
    Limb points from 1st and 3rd pass are shown, on top of the fitted circle.

    :param image: limb-fitted image
    :param xc: x-coordinate of disk center [pixels]
    :param yc: y-coordinate of disk center [pixels]
    :param Rm: disk radius
    :param pts: Series of points coordinates that were fitted. Dimensions: [# of points, 2]
    :param fname: full path of image file with extension (.png, .jpeg, ...)
    :return: figure written on disk
    """

    # Color of symbols in pyplot.plot sym1 for the 1st pass. sym2 for 3rd pass
    sym1 = 'co'
    sym3 = 'ro'
    # Color of limb-fitted circle
    ccolor = 'y'

    # Image dimension
    naxis1 = image.shape[1]
    naxis2 = image.shape[0]

    # Width and height [px] of the close-up of each quadrant
    width  = 500
    height = 200

    xgrid_tb = np.arange(width) + (naxis1/2 - width/2)
    ygrid_bot = np.arange(height)
    ygrid_top = ygrid_bot + (naxis2 - height)

    xgrid_east = np.arange(height)
    xgrid_west = xgrid_east + (naxis1 - height)
    ygrid_ew = np.arange(width) + (naxis2/2 - width/2)

    fig = plt.figure(1)
    fig.clear()
    # For full-screen in interactive mode:
    # manager = plt.get_current_fig_manager()
    # manager.window.showMaximized()

    # South limb
    # ax = plt.axes([x0, y0, width, height])
    ax = plt.axes([0.175, 0.06, 0.75, 0.3])
    plt.imshow(image, cmap='gray', origin='lower')
    plt.axis([xgrid_tb[0], xgrid_tb[-1], ygrid_bot[0], ygrid_bot[-1]])
    for i in range(0, pts.shape[0]):
       plt.plot(pts[i,0], pts[i,1], sym1, ms=2)
    for i in range(0, pts3.shape[0]):
       plt.plot(pts3[i,0], pts3[i,1], sym3, ms=2)

    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    # Draw the fitted circle
    circle1 = plt.Circle((xc, yc), Rm, color=ccolor, fill=False, linewidth=1)
    ax.add_artist(circle1)

    # North limb
    # ax = plt.axes([x0, y0, width, height])
    ax = plt.axes([0.175, 0.41, 0.75, 0.3])
    plt.imshow(image, cmap='gray', origin='lower')
    plt.axis([xgrid_tb[0], xgrid_tb[-1], ygrid_top[0], ygrid_top[-1]])
    for i in range(0, pts.shape[0]):
        plt.plot(pts[i, 0], pts[i, 1], sym1, ms=2)
    for i in range(0, pts3.shape[0]):
        plt.plot(pts3[i, 0], pts3[i, 1], sym3, ms=2)
    # plt.xlabel('X')
    # plt.ylabel('Y')
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    # Draw the fitted circle
    circle1 = plt.Circle((xc, yc), Rm, color=ccolor, fill=False, linewidth=1)
    ax.add_artist(circle1)

    # East limb
    # ax = plt.axes([x0, y0, width, height])
    ax = plt.axes([0, 0, 0.3, 0.75])
    plt.imshow(image, cmap = 'gray', origin='lower')
    plt.axis([xgrid_east[0], xgrid_east[-1], ygrid_ew[0], ygrid_ew[-1]])
    for i in range(0, pts.shape[0]):
       plt.plot(pts[i,0], pts[i,1], sym1, ms=2)
    for i in range(0, pts3.shape[0]):
       plt.plot(pts3[i,0], pts3[i,1], sym3, ms=2)
    # plt.xlabel('X')
    # plt.ylabel('Y')
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    # Draw the fitted circle
    circle1 = plt.Circle((xc, yc), Rm, color=ccolor, fill=False, linewidth=1)
    ax.add_artist(circle1)

    # West limb
    # ax = plt.axes([x0, y0, width, height])
    ax = plt.axes([0.8, 0, 0.3, 0.75])
    plt.imshow(image, cmap='gray', origin='lower')
    plt.axis([xgrid_west[0], xgrid_west[-1], ygrid_ew[0], ygrid_ew[-1]])
    for i in range(0, pts.shape[0]):
        plt.plot(pts[i, 0], pts[i, 1], sym1, ms=2)
    for i in range(0, pts3.shape[0]):
        plt.plot(pts3[i, 0], pts3[i, 1], sym3, ms=2)
    # plt.xlabel('X')
    # plt.ylabel('Y')
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    # Draw the fitted circle
    circle1 = plt.Circle((xc, yc), Rm, color=ccolor, fill=False, linewidth=1)
    ax.add_artist(circle1)

    plt.savefig(fname, dpi=180, bbox_inches='tight')

def circle_fit_by_Taubin(data_x, data_y):
    """
    Circle fit to a given set of data points (in 2D)

      This is an algebraic fit, due to Taubin, based on the journal article

      G. Taubin, "Estimation Of Planar Curves, Surfaces And Nonplanar
                  Space Curves Defined By Implicit Equations, With
                  Applications To Edge And Range Image Segmentation",
                  IEEE Trans. PAMI, Vol. 13, pages 1115-1138, (1991)

      Input:  data     - the class of data (contains the given points):

          data.n   - the number of data points
          data.X[] - the array of X-coordinates
          data.Y[] - the array of Y-coordinates

     Output:
               circle - parameters of the fitting circle:

           circle.a - the X-coordinate of the center of the fitting circle
           circle.b - the Y-coordinate of the center of the fitting circle
           circle.r - the radius of the fitting circle
           circle.s - the root mean square error (the estimate of sigma)
           circle.j - the total number of iterations

     The method is based on the minimization of the function

                  sum [(x-a)^2 + (y-b)^2 - R^2]^2
              F = -------------------------------
                      sum [(x-a)^2 + (y-b)^2]

     This method is more balanced than the simple Kasa fit.

     It works well whether data points are sampled along an entire circle or
     along a small arc.

     It still has a small bias and its statistical accuracy is slightly
     lower than that of the geometric fit (minimizing geometric distances),
     but slightly higher than that of the very similar Pratt fit.
     Besides, the Taubin fit is slightly simpler than the Pratt fit

     It provides a very good initial guess for a subsequent geometric fit.

       Nikolai Chernov  (September 2012)

    :return:
    """
    IterMax = 99

    meanX = np.mean(data_x)
    meanY = np.mean(data_y)
    n = data_x.shape[0]

    Mxx, Myy, Mxy, Mxz, Myz, Mzz = 0., 0., 0., 0., 0., 0.

    for i in range(0, n):

        Xi = data_x[i] - meanX
        Yi = data_y[i] - meanY
        Zi = Xi*Xi + Yi*Yi

        Mxy += Xi*Yi
        Mxx += Xi*Xi
        Myy += Yi*Yi
        Mxz += Xi*Zi
        Myz += Yi*Zi
        Mzz += Zi*Zi

    Mxx /= n
    Myy /= n
    Mxy /= n
    Mxz /= n
    Myz /= n
    Mzz /= n

    # computing coefficients of the characteristic polynomial

    Mz = Mxx + Myy
    Cov_xy = Mxx * Myy - Mxy * Mxy
    Var_z = Mzz - Mz * Mz
    A3 = 4.0 * Mz
    A2 = -3.0 * Mz * Mz - Mzz
    A1 = Var_z * Mz + 4.0 * Cov_xy * Mz - Mxz * Mxz - Myz * Myz
    A0 = Mxz * (Mxz * Myy - Myz * Mxy) + Myz * (Myz * Mxx - Mxz * Mxy) - Var_z * Cov_xy
    A22 = A2 + A2
    A33 = A3 + A3 + A3

    # finding the root of the characteristic polynomial
    # using Newton's method starting at x=0
    # (it is guaranteed to converge to the right root)

    x = 0
    y = A0
    n_iter = 0
    for iter in range(0, IterMax):
        n_iter = iter
        Dy = A1 + x * (A22 + A33 * x)
        xnew = x - y / Dy
        if (xnew == x) or np.isfinite(xnew):
            break

        ynew = A0 + xnew * (A1 + xnew * (A2 + xnew * A3))

        if (np.abs(ynew) >= np.abs(y)):
            break

        x = xnew
        y = ynew


    # Computing parameters of the fitting circle
    DET = x*x - x*Mz + Cov_xy
    Xcenter = (Mxz * (Myy - x) - Myz * Mxy) / DET / 2.0
    Ycenter = (Myz * (Mxx - x) - Mxz * Mxy) / DET / 2.0

    # Assemble the output
    a = Xcenter + meanX
    b = Ycenter + meanY
    r = np.sqrt(Xcenter**2.0 + Ycenter**2 + Mz)

    # Calculate standard deviation of the discrepancies
    dx = data_x - a
    dy = data_y - b
    dr = np.sqrt(dx**2 + dy**2) - r
    sigma = np.std(dr)

    return a, b, r, sigma, n_iter


def emil_circle_fit(x, y, debug, num_clipping_passes=1, clipping_sigma=2):
    """
    Fit a circle with data points (x,y) using a 3-pass least-square fit with sigma-clipping
    Uses method 2 from http://scipy-cookbook.readthedocs.io/items/Least_Squares_Circle.html

    :param x: numpy 1D array of the x-coordinates of points to fit
    :param y: numpy 1D array of the y-coordinates of points to fit
    :return: coordinates and radius of fitted circle for final pass

    Edit Emil: got rid of the Taubin method (it didn't seem to do any iterations ever, I think there is an implementation error; the breaks didn't make sense )
    Edit Emil: nicer sigma clipping
    """
    xc, yc, Rm, numDatapoints = np.nan, np.nan, np.nan, 0

    if (len(x) > 0) and (len(x) == len(y)):
        # coordinates of the barycenter used as initial estimate
        for clipping_pass in range(num_clipping_passes):
            numDatapoints = len(x)
            center_estimate = np.array([np.mean(x), np.mean(y)])
            center = optimize.leastsq(f_2, center_estimate, (x, y))[0]
            xc, yc = center
            R_dist = calc_R(x, y, xc, yc)
            Rm = np.mean(R_dist)
            if debug:
                print
                "pass %d Circle found!: %f %f radius %f (%d points)" % (clipping_pass, xc, yc, Rm, len(x))

            if clipping_pass < num_clipping_passes:
                Rmedian = np.median(R_dist)
                keep_mask = np.abs(R_dist - Rmedian) < R_dist.std() * clipping_sigma
                x = x[keep_mask]
                y = y[keep_mask]
                if len(keep_mask) == len(x):  # if we didn't throw away any datapoints, we can stop doing passes...
                    break

    return xc, yc, Rm, numDatapoints

def filter_hanning(image, w, filter_type):

    # 2D hanning filter, assumes a square image size
    # w: filter width
    n = image.shape[0]
    xgrid, ygrid = np.meshgrid(np.arange(n), np.arange(n))
    cutoff = 0

    # Get the fourier transform
    Fimage = np.fft.fftshift(np.fft.fftn(image))

    # grid of radial distances
    r = np.sqrt((xgrid - (n / 2 - 0.5)) ** 2 + (ygrid - (n / 2 - 0.5)) ** 2)
    han2D = 0.5 + 0.5 * np.cos(np.pi * (r - cutoff) / (2 * w))
    mask_r = r > cutoff + 2 * w
    han2D[mask_r] = 0

    # Apply 2D Hanning window to fourier transform
    windowed_Fimage = Fimage * han2D
    filtered_image = np.real(np.fft.ifftn(np.fft.ifftshift((windowed_Fimage))))

    return filtered_image

def correct_limb_darkening(image, center, radius):
    # Fit limb darkening.
    # This assume a square image in the polar transform

    # I(theta) = I(0) * ( 1 - ux)
    # I(theta) = I(0) * ( 1 - u(1-cos(theta)) )
    # cos(theta) = sqrt(radius**2 - r**2)/radius

    # Image dimensions
    naxis1 = image.shape[0]
    naxis2 = image.shape[1]
    # Coordinate of disk center
    xc = center[0]
    yc = center[1]

    imageP = cv2.linearPolar(image, center, naxis1, cv2.INTER_LANCZOS4 + cv2.WARP_FILL_OUTLIERS)
    # Take the median over theta (azimuthal average)
    Iaz_avg = np.median(imageP, 0)

    a = np.round(radius)
    r = np.arange(a - 10)

    Itheta = Iaz_avg[0:len(r)]

    x = 1 - np.sqrt(a ** 2 - r ** 2) / a
    y = Itheta / Itheta[0]

    fit = np.polyfit(x, y, 1)
    fit_fn = np.poly1d(fit)

    limb_factor = fit[0] * x + fit[1]

    # Define grid of cos(theta)
    xgrid1, ygrid1 = np.meshgrid(np.arange(naxis1), np.arange(naxis2))
    pixel_r = np.sqrt((xgrid1 - xc) ** 2 + (ygrid1 - yc) ** 2)

    cos_theta = np.sqrt(a ** 2 - pixel_r ** 2) / a

    x2 = 1 - cos_theta
    limb_factor2D = fit[0] * x2 + fit[1]

    mask_limb = pixel_r > a
    limb_factor2D[mask_limb] = 1

    image_limb = image / limb_factor2D

    return image_limb, fit, x, fit_fn(x)