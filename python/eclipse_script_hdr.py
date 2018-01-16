import numpy as np
import cv2
import matplotlib
matplotlib.use('macosx')
from scipy.ndimage import generic_filter
#from skimage.filters.rank import entropy
#from skimage.morphology import disk
import matplotlib.pyplot as plt
from matplotlib import cm

def calc_entropy(array):
    # Apply this within a generic filter
    # https://docs.scipy.org/doc/scipy-0.19.0/reference/generated/scipy.ndimage.generic_filter.html
    hist, _ = np.histogram(array, bins = 100)
    hist = hist / array.size
    hist = hist[np.nonzero(hist)]
    log_hist = np.log(hist) / np.log(2)
    entrop = hist * log_hist
    entrop = -1 * entrop.sum()

    return entrop

def get_radius_array(radius):

    x = np.arange(0, radius*2)
    y = np.arange(0, radius*2)
    xv, yv = np.meshgrid(x, y)
    disk_center = (2*radius - 1) / 2
    radius_array = np.sqrt( (xv - disk_center)**2 + (yv - disk_center)**2 )

    return radius_array

def get_disk(radius):

    r = get_radius_array(radius)
    r[r < radius] = 1
    r[r > radius] = 0
    disk = r.astype(np.bool)

    return disk

def get_radius_array2(disk_center, nx, ny):

    x = np.arange(0, nx)
    y = np.arange(0, ny)
    xv, yv = np.meshgrid(x, y)
    radius_array = np.sqrt( (xv - disk_center[0])**2 + (yv - disk_center[1])**2 )

    return radius_array


# files = ['/Users/rattie/Data/Eclipse/Stack_mean_aligned_series/aligned/stack_mean_C_exp250.tiff',
#          '/Users/rattie/Data/Eclipse/Stack_mean_aligned_series/aligned/stack_mean_C_exp30.tiff',
#          '/Users/rattie/Data/Eclipse/Stack_mean_aligned_series/aligned/stack_mean_C_exp20.tiff',
#          '/Users/rattie/Data/Eclipse/Stack_mean_aligned_series/aligned/stack_mean_C_exp10.tiff',
#          '/Users/rattie/Data/Eclipse/Stack_mean_aligned_series/aligned/stack_mean_C_exp1.tiff']
# times = np.array([250, 30, 20, 10, 1])
files = ['/Users/rattie/Data/Eclipse/Stack_mean_aligned_series/aligned/stack_mean_C_exp250.tiff',
         '/Users/rattie/Data/Eclipse/Stack_mean_aligned_series/aligned/stack_mean_C_exp30.tiff',
         '/Users/rattie/Data/Eclipse/Stack_mean_aligned_series/aligned/stack_mean_C_exp1.tiff']
times_array = np.array([250, 30, 1])

# Exposure times
nexp = len(files)
sum_times = times_array.sum()
weights = np.reshape(times_array/sum_times, [nexp, 1, 1, 1])

max16 = 2**16-1
max14 = 2**14-1

img0 = cv2.imread(files[0], cv2.IMREAD_COLOR + cv2.IMREAD_ANYDEPTH)
img0RGB = np.fliplr(img0.reshape(-1,3)).reshape(img0.shape)
img0RGBN = img0RGB/ img0RGB.max()

nx = img0.shape[1]
ny = img0.shape[0]

images = np.zeros([nexp, ny, nx, 3])
for i in range(nexp):
    images[i, :, :, :] = cv2.imread(files[i], cv2.IMREAD_COLOR + cv2.IMREAD_ANYDEPTH)
    images[i, :, :, :] = np.fliplr(images[i,:,:,:].reshape(-1,3)).reshape(images[i,:,:,:].shape)
    #maxRGB = temp.reshape([nx*ny, 3]).max(axis=0)

images_RGBN = images/max14

reds     = np.squeeze(images[:,:,:,0])
sat_mask = reds >= max14

images[sat_mask, :] = max14

plt.figure(0)
plt.imshow(images[1,...]/max14, origin='lower')

# Normalizeto the exposure time, so we get the intensity per second
imagesN = images * np.reshape(times_array, [nexp, 1, 1, 1]) # * (2**16 -1) / (2**14 -1)
# RGB normalized
images_expN_RGB = np.zeros([nexp, ny, nx, 3])
for i in range(nexp):
    temp = np.fliplr(imagesN[i, :, :, :].reshape(-1, 3)).reshape(imagesN[i, :, :, :].shape)
    images_expN_RGB[i, :, :, :] = temp / max16



# Gray scaled image
gimages = imagesN.mean(axis=3)
# Coordinate of the center
center = (2120, 1361)
radius = 641/2

pimages = np.zeros([nexp, ny, nx])
for i in range(0,3):
    # Azimuthal average
    pimages[i,:,:] = cv2.linearPolar(gimages[i,:,:], center, gimages.shape[2], cv2.INTER_LANCZOS4 + cv2.WARP_FILL_OUTLIERS)
    # Take the median over theta (azimuthal average)

ymin = 1100
ymax = 1301
# Take the median over this selection on the y-axis.
az_avgs = np.median(pimages[:,ymin:ymax, :], 1)
zeroloc1 = np.where(pimages[2, ymin, :] == 0)[0].min()
zeroloc2 = np.where(pimages[2, ymax, :] == 0)[0].min()
zeroloc = np.min([zeroloc1, zeroloc2])
az_avgs[:, zeroloc:] = 0

ny, nx = gimages[0,...].shape


plt.figure(1, figsize=(19,10))
ax1 = plt.subplot(2,3,1)
ax1.imshow(gimages[0,:,:], cmap='gray', origin='lower')

ax2 = plt.subplot(2,3,2)
ax2.imshow(gimages[1,:,:], cmap='gray', origin='lower')

ax3 = plt.subplot(2,3,3)
ax3.imshow(gimages[2,:,:], cmap='gray', origin='lower')

ax4 = plt.subplot(2,3,4)
ax4.imshow(pimages[0,:,:], cmap='gray', origin='lower')
ax4.set_aspect('equal', 'box-forced')

ax41 = ax4.twinx()
ax41.plot(az_avgs[0,:])
ax41.set_xlim(ax4.get_xlim())

a = np.diff( ax41.get_ylim() )[0] / np.diff( ax4.get_xlim() ) * nx/ny

ax41.set_aspect(1./a, 'box-forced')

ax5 = plt.subplot(2,3,5)
ax5.imshow(pimages[1,:,:], cmap='gray', origin='lower')
ax5.set_aspect('equal', 'box-forced')

ax6 = plt.subplot(2,3,6)
ax6.imshow(pimages[2,:,:], cmap='gray', origin='lower')
ax6.set_aspect('equal', 'box-forced')


ax51 = ax5.twinx()
ax51.plot(az_avgs[1,:])
ax51.set_xlim(ax5.get_xlim())

a = np.diff( ax51.get_ylim() )[0] / np.diff( ax5.get_xlim() ) * nx/ny
ax51.set_aspect(1./a, 'box-forced')


ax61 = ax6.twinx()
ax61.plot(az_avgs[2,:])
ax61.set_xlim(ax6.get_xlim())

a = np.diff( ax61.get_ylim() )[0] / np.diff( ax6.get_xlim() ) * nx/ny
ax61.set_aspect(1./a, 'box-forced')

#plt.tight_layout()

plt.figure(2, figsize=(15,10))
plt.imshow(pimages[1,ymin:ymax,:], vmax=max14, cmap='gray', origin='lower')
plt.xlim(xmin= 0, xmax=2500)


# Fit from x=800 to x=2000 incl.
x1 = 700
x2 = zeroloc
data_to_fit = az_avgs[1, x1:x2].copy()
x = np.arange(x1, x2)
p = np.polyfit(x, data_to_fit, 6)
fit = np.poly1d(p)
# Draw the fit and extrapolate beyond the fitted data
xf = np.arange(radius, x2)

## Use the fit to reconstruct a background image to subtract.

# Make a radius map: at each pixel, the value is the distance to disc center
r = get_radius_array2(center, nx, ny)

imfit = fit(r)
imfit[r <= radius] = 0

# plot
plt.figure(3, figsize=(18,7))
plt.subplot(121)
plt.plot(az_avgs[0,:], 'k-')
plt.plot(az_avgs[1,:], 'g-')
plt.plot(az_avgs[2,:], 'r-')
plt.yscale('log')
plt.xlim(xmax=2500)
plt.ylim(ymin=1e3)
plt.grid(True, which='both')


plt.plot(xf, fit(xf), 'b--')
plt.subplot(122)
plt.imshow(imfit, origin='lower', cmap='gray')
plt.tight_layout()

imgN_float = imagesN/max16
imgN_float[sat_mask] = 1
imgRGBN = imagesN.copy()
for i in range(0, 3):
    # Swap blue <-> red channel
    fimgRGB = np.fliplr(imagesN[i,:,:,:].reshape(-1,3)).reshape(imagesN[i,:,:,:].shape)
    #fimgRGBN[i, :, :, :] = fimgRGB / imagesN[i, :, :, :].max()
    imgRGBN[i, :, :, :] = (fimgRGB - imfit[:,:,np.newaxis]) / max16


plt.figure(4, figsize=(18,10))
plt.subplot(121)
# Plot at 1/30s which was used for the fit
plt.imshow(imgN_float[1,...], origin='lower')

## Apply radial FILTERS


# plt.figure(0)
# plt.imshow(r, cmap='gray')

# Mask everything inside the disk
rr = r - radius
# Take care of singularity at r = 0
rr[rr <= 1] = 1
rrN = rr / rr.max()

rrs = np.tile(np.reshape(rr, [1, ny, nx, 1]), (3, 1, 1, 3))
rrsN = np.tile(np.reshape(rrN, [1, ny, nx, 1]), (3, 1, 1, 3))
rrs2 =  1000*rrsN**1.5 + 100/(0.1*rrs + 100)

#rrs2[np.logical_and(rrs >= 0, rrs <= 15) ] = rrs[np.logical_and(rrs >= 0, rrs <= 15)]**2

fimages = imagesN * rrs2

fimgRGBN = fimages.copy()
fimgRGBN2 = fimages.copy()

# When calculating the max for normalization, do not take into account the saturated pixels that are rescaled
fimgRGBN2[sat_mask] = 0


for i in range(0, 3):
    # Swap blue <-> red channel
    fimgRGB = np.fliplr(fimages[i,:,:,:].reshape(-1,3)).reshape(fimages[i,:,:,:].shape)
    #fimgRGBN[i, :, :, :] = fimgRGB / imagesN[i, :, :, :].max()
    fimgRGBN[i, :, :, :] = fimgRGB / (max16 * rrs2.max())


#fimgRGBN[sat_mask] = 1

scale = 8

plt.figure(2, figsize=(18,10))
plt.subplot(1,3,1)
plt.imshow(fimgRGBN[0, :, :, :]*scale, origin='lower')
plt.title('1/r compensated (1/250s)')

plt.subplot(1,3,2)
plt.imshow(fimgRGBN[1, :, :, :]*scale, origin='lower')
plt.title('1/r compensated (1/30s)')

plt.subplot(1,3,3)
plt.imshow(fimgRGBN[2, :, :, :]*scale, origin='lower')
plt.title('1/r compensated (1s)')

plt.tight_layout()



plt.figure(0, figsize=(20,10))
plt.subplot(1,2,1)
plt.imshow(images_RGBN[0, :, :, :]*2, origin='lower')
ax = plt.gca()
#plt.imshow(rrMask.astype(np.float32), origin='lower')
circle = plt.Circle(center, radius, color='g', fill=False)
ax.add_artist(circle)
plt.title('Original image')

plt.subplot(1,2,2)
plt.imshow(fimgRGBN[0, :, :, :]*scale, origin='lower')
plt.title('1/r compensated')

plt.tight_layout()

plt.figure(1, figsize=(20,10))
plt.subplot(1,2,1)
plt.imshow(images_RGBN[1, :, :, :], origin='lower')
ax = plt.gca()
#plt.imshow(rrMask.astype(np.float32), origin='lower')
circle = plt.Circle(center, radius, color='g', fill=False)
ax.add_artist(circle)
plt.title('Original image')

plt.subplot(1,2,2)
plt.imshow(fimgRGBN[1, :, :, :], origin='lower')
plt.title('1/r compensated')

plt.tight_layout()


plt.figure(2, figsize=(20,10))
plt.subplot(1,2,1)
plt.imshow(images_RGBN[2, :, :, :], origin='lower')
plt.title('Original image')

plt.subplot(1,2,2)
plt.imshow(fimgRGBN[2, :, :, :]*scale, origin='lower')
plt.title('1/r compensated (1s)')

plt.tight_layout()


# fimagesN16 = fimages.copy()
# for i in range(0, 3):
#     fimagesN16[i, :, :, :] = fimages[i, :, :, :] * scale * max16/ fimgRGBN2.max()
#
# cv2.imwrite('/Users/rattie/Data/Eclipse/Stack_mean_aligned_series/aligned/rf_img_250_r1p5.tiff', fimagesN16[0, :, :, :].astype(np.uint16))
# cv2.imwrite('/Users/rattie/Data/Eclipse/Stack_mean_aligned_series/aligned/rf_img_30_r1p5.tiff', fimagesN16[1, :, :, :].astype(np.uint16))
# cv2.imwrite('/Users/rattie/Data/Eclipse/Stack_mean_aligned_series/aligned/rf_img_1_r1p5.tiff', fimagesN16[2, :, :, :].astype(np.uint16))



### Using expoBlend

# images = np.zeros([nexp, ny, nx, 3])
# for i in range(nexp):
#     images[i, :, :, :] = cv2.imread(files[i], cv2.IMREAD_COLOR + cv2.IMREAD_ANYDEPTH)
#
# images_gray = images.mean(axis=3)
# reds = images[:, :, :, 2]
#
# imagesN     = images.astype(np.float) / images.max()
# images_grayN = images_gray / images_gray.max()
# redsN = reds / reds.max()
#
# Ilog = np.log(1 + imagesN)/np.log(2)
# Ilog_gray = np.log(1 + images_grayN)/np.log(2)
# Ilog_red = np.log(1 + redsN)/np.log(2)
#
# # Entropy transform
# radius = 10
# disk = get_disk(radius)
# Hklog = np.zeros([nexp, ny, nx], dtype = np.float)
# for i in range(nexp):
#     # Convert BGR to gray
#     imlog   = Ilog_red[i, :, :]
#     Hklog[i, :, :] = generic_filter(imlog, calc_entropy, footprint=disk)
#
#
# Hk_norm = Hklog / np.reshape(Hklog.sum(axis=0), [1, ny, nx])
#
# beta = 10
# Hk_nlnorm = np.exp(beta * Hk_norm)
#
# Hk_weights = Hk_nlnorm / np.reshape(Hk_nlnorm.sum(axis = 0), [1, ny, nx])
#
# log_images_w = Ilog * np.tile(Hk_weights[:, :, :, np.newaxis], (1, 1, 1, 3))
# summed_log_image = log_images_w.sum(axis = 0)
#
# final_blended = np.exp(summed_log_image) -1
# final_blended16 = final_blended * 65535 / final_blended.max()
#
# cv2.imwrite('/Users/rattie/Data/Eclipse/Stack_mean_aligned_series/aligned/expoBlend_beta10_red.tiff', final_blended16.astype(np.uint16))
