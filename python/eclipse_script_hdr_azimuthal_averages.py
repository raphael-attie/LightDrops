import matplotlib
matplotlib.use('TKagg')
import matplotlib.pyplot as plt
import numpy as np
import cv2
from scipy import optimize
import os
plt.ion()


def get_disk(radius):

    r = get_radius_array(radius)
    r[r < radius] = 1
    r[r > radius] = 0
    disk = r.astype(np.bool)

    return disk


def get_radius_array(disk_center, nx, ny):

    x = np.arange(0, nx)
    y = np.arange(0, ny)
    xv, yv = np.meshgrid(x, y)
    radius_array = np.sqrt( (xv - disk_center[0])**2 + (yv - disk_center[1])**2 )

    return radius_array


def fit_bkg(bkg_avg, xmin, xmax, radius_grid, rmin):
    # Define function for calculating a power law
    powerlaw = lambda x, amp, index: amp * (x ** index)
    # define our (line) fitting function
    fitfunc = lambda p, x: p[0] + p[1] * x
    errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err

    xdata = np.arange(xmin, xmax, 1)
    ydata = bkg_avg[xdata]
    yerr = 0.2 * ydata
    logx = np.log10(xdata)
    logy = np.log10(ydata)
    logyerr = 0.001 * np.ones(ydata.shape)

    pinit = [1.0, -1.0]
    out = optimize.leastsq(errfunc, pinit,
                           args=(logx, logy, logyerr), full_output=1)
    pfinal = out[0]
    covar = out[1]
    index = pfinal[1]
    amp = 10.0 ** pfinal[0]

    bkg_2d = powerlaw(radius_grid, amp, index)
    bkg_2d[rr <= rmin] = 0

    return bkg_2d, amp, index



files = ['/Users/rattie/Data/Eclipse/Stack_mean_aligned_series/aligned/stack_mean_C_exp250.tiff',
         '/Users/rattie/Data/Eclipse/Stack_mean_aligned_series/aligned/stack_mean_C_exp30.tiff',
         '/Users/rattie/Data/Eclipse/Stack_mean_aligned_series/aligned/stack_mean_C_exp1.tiff']

savedir = '/Users/rattie/Data/Eclipse/Background_removed/'

# Exposure times
nexp = 3
times_array = np.array([250, 30, 1])
sum_times = times_array.sum()
# Normalization factors
max16 = 2**16-1
max14 = 2**14-1
# Load a sample
img0 = cv2.imread(files[0], cv2.IMREAD_COLOR + cv2.IMREAD_ANYDEPTH)
nx = img0.shape[1]
ny = img0.shape[0]
# Load all images and convert from openCV BGR to RGB
images = np.zeros([nexp, ny, nx, 3])
for i in range(nexp):
    images[i, :, :, :] = cv2.imread(files[i], cv2.IMREAD_COLOR + cv2.IMREAD_ANYDEPTH)
    images[i, :, :, :] = np.fliplr(images[i,:,:,:].reshape(-1,3)).reshape(images[i,:,:,:].shape)


# Get a reference version that can be visualized as float, normalized from [0-1].
# Because of the flat fielding, saturation occurs above max14. So we choose to normalize to a meaningful value,
# maxval =  maximum near-saturating intensity in the prominence.
maxval = images[:, 1080:1120, 1900:1960, 0].max() # ~ 20000
imagesf = images/maxval
# Clip within 0, 1 - for visualization only. Processing is made on an unclipped version.
imagesf = np.clip(imagesf, 0, 1)
# Coordinate of the center
center = (2120, 1361)
radius = 323
# Plot one image and the circle to check the center and radius
#circle1 = plt.Circle((0, 0), 0.2, color='r')
plt.figure(0, figsize=(15,11))
ax1 = plt.gcf().add_subplot(111)
ax1.imshow(np.clip(imagesf[0,...]*5, 0, 1), origin='lower')
#ax1.imshow(np.clip(imagesf[0,...], 0, 1), origin='lower')
ax1.add_artist(plt.Circle(center, radius, color='green', fill=False, linewidth=1))
ax1.add_artist(plt.Circle(center, 330, color='yellow', fill=False))
ax1.axis([1000, 3000, 500, 2500])
plt.title('1/250s , brightness x 5')
plt.tight_layout()

# Normalize to the exposure time, so we get the same intensity per second + noise
tarray3 = np.reshape(times_array, [nexp, 1, 1, 1])
images_exp = images * tarray3 # * (2**16 -1) / (2**14 -1)
# Divide again to maxval to have a visualization equivalent to "imagesf" (consider that 1s should look the same).
images_expf = images_exp/maxval
# Clip saturated pixels in the exposure-normalized images
images_expf_clipped = np.clip(images_expf, 0, 1)
# Display the 3 images in RGB.
fig = plt.figure(1, figsize=(19,10))
plt.subplot(231)
plt.imshow(imagesf[0,...])
plt.title('1/250s')
plt.subplot(232)
plt.imshow(imagesf[1,...])
plt.title('1/30s')
plt.subplot(233)
plt.imshow(imagesf[2,...])
plt.title('1s')

plt.subplot(234)
plt.imshow(images_expf_clipped[0,...])
plt.title('x 250')
plt.subplot(235)
plt.imshow(images_expf_clipped[1,...])
plt.title('x 30')
plt.subplot(236)
plt.imshow(images_expf_clipped[2,...])
plt.title('x 1')
plt.tight_layout()
plt.savefig('/Users/rattie/Data/Eclipse/figures/imagesRGB_and_exp_normalized.png')


## Remap to polar coordinates
pimages_rgb = np.zeros([nexp, ny, nx, 3])
for exp in range(0,3):
    for ch in range(0,3):
        # Azimuthal average
        pimages_rgb[exp,: , :, ch] = cv2.linearPolar(images_exp[exp, :, :, ch], center, nx, cv2.INTER_LANCZOS4 + cv2.WARP_FILL_OUTLIERS)

pimages_rgbf = np.clip(pimages_rgb / maxval, 0, 1)

# Get a selection to take the median over the y-axis.
ymin = 990
ymax = 1100

fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(19, 8))
axes[0].imshow(np.clip(pimages_rgbf[0,...], 0, 1), origin='lower')
axes[0].set_title('1/250s')
axes[1].imshow(np.clip(pimages_rgbf[1,...], 0, 1), origin='lower')
axes[1].set_title('1/30s')
axes[2].imshow(np.clip(pimages_rgbf[2,...]*5, 0, 1), origin='lower')
axes[2].set_title('1s x 5')

axes[1].axhline(y=ymin, ls='--', linewidth=1, color='yellow')
axes[1].axhline(y=ymax, ls='--', linewidth=1, color='yellow')

plt.tight_layout()



# Extract and average background radial profile
bkg_az_avg1 = pimages_rgbf[1, ymin:ymax, :, :].mean(axis=0)
bkg_az_avg2 = pimages_rgbf[2, ymin:ymax, :, :].mean(axis=0)
rr = get_radius_array(center, nx, ny)

rmin = 1300
rmax = 1800
# Red
bkg_r, amp_r, index_r = fit_bkg(bkg_az_avg1[:,0], rmin, rmax, rr, radius)
# Green
bkg_g, amp_g, index_g = fit_bkg(bkg_az_avg1[:,1], rmin, rmax, rr, radius)
# Blue
bkg_b, amp_b, index_b = fit_bkg(bkg_az_avg1[:,2], rmin, rmax, rr, radius)

# 1s exposure background fit
rmin = 1300
rmax = 1800
# Red
bkg_r2, amp_r2, index_r2 = fit_bkg(bkg_az_avg2[:,0], rmin, rmax, rr, radius)
# Green
bkg_g2, amp_g2, index_g2 = fit_bkg(bkg_az_avg2[:,1], rmin, rmax, rr, radius)
# Blue
bkg_b2, amp_b2, index_b2 = fit_bkg(bkg_az_avg2[:,2], rmin, rmax, rr, radius)


# for plotting, get a longer array of radius values for extrapolating the profile to the center.
xdata2 = np.arange(radius, nx, 1)
powerlaw = lambda x, amp, index: amp * (x ** index)


fig = plt.figure(figsize=(15, 10))
plt.subplot(2,1,1)
plt.plot(bkg_az_avg1[:,0], 'r-')
plt.plot(bkg_az_avg1[:,1], 'g-')
plt.plot(bkg_az_avg1[:,2], 'b-')
plt.plot(xdata2, powerlaw(xdata2, amp_r, index_r), 'k--')
plt.plot(xdata2, powerlaw(xdata2, amp_g, index_g), 'k-.')
plt.plot(xdata2, powerlaw(xdata2, amp_b, index_b), 'k:')
plt.title('Background 1/30 x 30 azimuthal average over y=[{:d} ; {:d}]'.format(ymin, ymax))
plt.subplot(2,1,2)
plt.plot(bkg_az_avg2[:,0], 'r-')
plt.plot(bkg_az_avg2[:,1], 'g-')
plt.plot(bkg_az_avg2[:,2], 'b-')
plt.plot(xdata2, powerlaw(xdata2, amp_r2, index_r2), 'k--')
plt.plot(xdata2, powerlaw(xdata2, amp_g2, index_g2), 'k-.')
plt.plot(xdata2, powerlaw(xdata2, amp_b2, index_b2), 'k:')
plt.title('Background 1s azimuthal average over y=[{:d} ; {:d}]'.format(ymin, ymax))
plt.tight_layout()


vmin = 0
vmax = 0.5
fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(18, 6))
# Show fitted background image
axs[0].imshow(bkg_r, vmin=vmin, vmax=vmax, cmap='gray', origin='lower')
axs[1].imshow(bkg_g, vmin=vmin, vmax=vmax, cmap='gray', origin='lower')
axs[2].imshow(bkg_b, vmin=vmin, vmax=vmax, cmap='gray', origin='lower')
axs[0].set_title('Fitted red background')
axs[1].set_title('Fitted green background')
axs[2].set_title('Fitted blue background')
plt.tight_layout()
plt.savefig('/Users/rattie/Data/Eclipse/figures/fitted_backgrounds.png')


# Remove background from the images
bkg_stack = np.stack([bkg_r, bkg_g, bkg_b], axis=-1)
bkg_stack2 = np.stack([bkg_r2, bkg_g2, bkg_b2], axis=-1)
images_bkg_removed = np.clip(images_expf - bkg_stack, 0, 1)
images_bkg_removed2 = images_expf - bkg_stack2
images_bkg_removed2[images_bkg_removed2 < 0 ] = 0

images_bkg_removed_v = np.clip((images_expf - bkg_stack)*5, 0, 1)


# Display the 3 images in RGB with background removed
fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(19,11))
axs[0,0].imshow(images_expf_clipped[0,...])
axs[0,0].set_title('x 250')
axs[0,1].imshow(images_expf_clipped[1,...])
axs[0,1].set_title('x 30')
axs[0,2].imshow(images_expf_clipped[2,...])
axs[0,2].set_title('x 1')

axs[1,0].imshow(images_bkg_removed[0,...])
axs[1,0].set_title('bkg removed')
axs[1,1].imshow(images_bkg_removed[1,...])
axs[1,1].set_title('bkg removed')
axs[1,2].imshow(images_bkg_removed[2,...])
axs[1,2].set_title('bkg removed')

axs[2,0].imshow(images_bkg_removed2[0,...])
axs[2,0].set_title('bkg removed (2)')
axs[2,1].imshow(images_bkg_removed2[1,...])
axs[2,1].set_title('bkg removed (2)')
axs[2,2].imshow(images_bkg_removed2[2,...])
axs[2,2].set_title('bkg removed (2)')

plt.tight_layout()

plt.savefig('/Users/rattie/Data/Eclipse/figures/Images_background_removed.png')


# Compare the far end of the corona once the backgroun is removed => multiply brightness by some factor
images_bkg_removed2_v = np.clip((images_expf - bkg_stack2)*5, 0, 1)


fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(19,11))
axs[0,0].imshow(images_bkg_removed_v[0,...])
axs[0,0].set_title('bkg removed')
axs[0,1].imshow(images_bkg_removed_v[1,...])
axs[0,1].set_title('bkg removed')
axs[0,2].imshow(images_bkg_removed_v[2,...])
axs[0,2].set_title('bkg removed')

axs[1,0].imshow(images_bkg_removed2_v[0,...])
axs[1,0].set_title('bkg removed (2)')
axs[1,1].imshow(images_bkg_removed2_v[1,...])
axs[1,1].set_title('bkg removed (2)')
axs[1,2].imshow(images_bkg_removed2_v[2,...])
axs[1,2].set_title('bkg removed (2)')
plt.tight_layout()

plt.savefig('/Users/rattie/Data/Eclipse/figures/Images_background_removed_compare.png')

# Radial filter. Try with 1/r




# Exporting to Tiff files: must restore the equivalent data range by multiplying by the exposure and maxval

suffix = ['250th', '30th', '1']
exposures = [250, 30, 1]

# Save the images as numpy arrays
fname = os.path.join(savedir, 'img_bkg_removed.npz')
np.savez(fname, images_bkg_removed2)


#images = np.zeros([nexp, ny, nx, 3])
# "De-normalize" the data to restore the equivalent data range falling within 14-bit depth.
images = images_bkg_removed2 * maxval / tarray3

for exp, s in enumerate(suffix):
    # convert from openCV RGBB to BGR
    image_bgr = cv2.cvtColor(images[exp, ...].astype(np.uint16), cv2.COLOR_RGB2BGR)
    fname = os.path.join(savedir, 'bkg_removed_{:s}.tiff'.format(s))
    cv2.imwrite(fname, image_bgr)

#TODO: remove background independent from 30th and 1s exposure. The 30th is done already. Only the 1s needs to be done.
#TODO: No need to remove the background from the 250s
#TODO: The corona will have to be radially attenuated at the begnining for any overlay with prominence, and increased with radius.