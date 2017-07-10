import os
import glob
from astropy.io import fits
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import calibration.uset_calibration as uset
import sunpy.cm as cm

import cv2



# Disable interactive mode so figure goes directory to a file and does not show up on screen
plt.ioff()

# Set the directory where the FITS are.
data_dir    =  'D:\\USET\\DATA\\HAlpha\\UPH20161215_Short_Exp\\calibrated_level1.1'
# Get the list of files, change it according to where your data files are and how are they are named.
file_list   = glob.glob(os.path.join(data_dir, '*.fits'))
# Output directory of jpeg previews
outdir_jpeg  = os.path.join(data_dir, 'preview_limb_scaled_jpeg')
outdir_filtered = os.path.join(data_dir, 'filtered_fits')
if not os.path.isdir(outdir_jpeg):
    os.makedirs(outdir_jpeg)
if not os.path.isdir(outdir_filtered):
    os.makedirs(outdir_filtered)

# Image dimension
naxis1 = 2048
naxis2 = naxis1
nt     = 10
images = np.zeros([naxis1, naxis2, nt])
for i in range(0,nt):
    file = file_list[i]
    hdu = fits.open(file, ignore_missing_end=True)
    # Load header and image data from the 1st data unit: hdu[0]
    header = hdu[0].header
    images[:,:,i] = hdu[0].data
# Take median of image series
image = np.mean(images, 2)

# Coordinate of disk center
xc = naxis1/2 - 0.5
yc = naxis2/2 - 0.5

# Calculate azimuthal average
# Take the polar transform of the image; knowing the disk center it is rather straight forward
# rho = (image width / max_radius) * sqrt(x**2 + y**2) , phi = atan2(y/x)
center = (xc, yc)
# Set max_radius to image width = naxis1 so rho = 1.0 * sqrt(x**2 + y**2)
max_radius = naxis1

imageP = cv2.linearPolar(image, center, max_radius, cv2.INTER_LANCZOS4 + cv2.WARP_FILL_OUTLIERS)
# Invert the polar transform
# cv2.linearPolar(image, center, max_radius, cv2.INTER_LANCZOS4 + cv2.WARP_FILL_OUTLIERS)

# Take the median over theta (azimuthal average)
Iaz_avg = np.median(imageP, 0)


xgrid1, ygrid1 = np.meshgrid(np.arange(naxis1), np.arange(naxis2))
# Distance of each pixel to disk center
pixel_r = np.sqrt((xgrid1 - xc)**2 + (ygrid1 - yc)**2)



# apply a radius-based scaling to boost the limb
radius = header['SOLAR_R']
# Mask of off-limb pixels
pixels_off_limb = pixel_r > radius

image2 = np.array(image)
image2[pixels_off_limb] = image2[pixels_off_limb] * 5

# Get the maximum intensity for the rescaled image based on 99.97% percentile
new_max = uset.compute_intensity_high(image)


# 1D hanning
n = int(100)
x = np.arange(n)
cutoff = 50
a = 20
han = 0.5 + 0.5*np.cos(np.pi*(x - cutoff)/(2*a) )
xz0 = cutoff - 2*a
xz1 = cutoff + 2*a
han[0:xz0] = 0
han[xz1:] = 0

# 2D Hanning filter
for w in range(500, 1024, 4):
    filtered_image = uset.filter_hanning(image, w, 'blah')
    basename_fits = 'filtered_width_%d.fits' %(w)
    fname = os.path.join(outdir_filtered, basename_fits)
    print(basename_fits)
    uset.write_uset_fits(filtered_image, header, fname)

#
# filtered_image[pixels_off_limb] = filtered_image[pixels_off_limb] * 5
#
#
# fig1 = plt.figure(1)
# fig1.clear()
# # Plot image
# plt.subplot(1,2,1)
# plt.imshow(image, vmin=100, vmax=new_max, cmap = 'gray', origin='lower')
# plt.xlabel('X [px]')
# plt.ylabel('Y [px]')
#
# plt.subplot(1,2,2)
# # Plot image
# plt.imshow(image2, vmin=100, vmax=new_max, cmap = 'gray', origin='lower')
# plt.xlabel('X [px]')
# plt.ylabel('Y [px]')
#
#
# # Plot polar transform
# fig2 = plt.figure(2)
# plt.subplot(1,2,1)
# plt.imshow(np.transpose(imageP), vmin=100, vmax=new_max, cmap = 'gray', origin='lower')
# plt.xlabel('Theta (degrees)')
# plt.ylabel('Radius (px)')
# ax = plt.gca()
#
# # rho = (naxis1 / max_radius) * sqrt(x**2 + y**2) , phi = atan2(y/x)
#
# # ticks
# xticks = np.arange(0, naxis1, 200)
# yticks = np.arange(0, naxis2, 400)
#
# rho = np.round((naxis1 / max_radius) * np.sqrt(yticks**2))
# phi = np.linspace(0, 360, len(xticks))
#
# xticklabels = np.char.mod('%d', phi)
# yticklabels = np.char.mod('%d', rho)
#
# ax.set_xticks(xticks)
# #ax.set_yticks(yticks)
# ax.set_xticklabels(xticklabels)
# #ax.set_yticklabels(yticklabels)
#
# plt.subplot(1,2,2)
# plt.plot(Iaz_avg)
#
#
#
# fig5 = plt.figure(5)
# plt.subplot(1,2,1)
# plt.imshow(image2, vmin=100, vmax=new_max, cmap = 'gray', origin='lower')
# plt.axis([0, 500, 0, 500])
# plt.subplot(1,2,2)
# plt.imshow(filtered_image, vmin=100, vmax=new_max, cmap = 'gray', origin='lower')
# plt.axis([0, 500, 0, 500])
# plt.show()
#
