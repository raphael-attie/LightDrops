import numpy as np
import cv2
import matplotlib
matplotlib.use('macosx')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from scipy.optimize import curve_fit
from scipy.signal import medfilt

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

def newline(p1, p2):
    ax = plt.gca()
    xmin, xmax = ax.get_xbound()

    if(p2[0] == p1[0]):
        xmin = xmax = p1[0]
        ymin, ymax = ax.get_ybound()
    else:
        ymax = p1[1]+(p2[1]-p1[1])/(p2[0]-p1[0])*(xmax-p1[0])
        ymin = p1[1]+(p2[1]-p1[1])/(p2[0]-p1[0])*(xmin-p1[0])

# def onclick(event):
#     print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
#           ('double' if event.dblclick else 'single', event.button,
#            event.x, event.y, event.xdata, event.ydata))
      ## Pos must be defined outside
#     pos.append([event.xdata, event.ydata])

def rational(x, q):
    """
    The general rational function description.
    p is a list with the polynomial coefficients in the numerator
    q is a list with the polynomial coefficients (except the first one)
    in the denominator
    The first coefficient of the denominator polynomial is fixed at 1.
    """
    return 1 / np.polyval(q, x)

def rational2(x, p, q):
    """
    The general rational function description.
    p is a list with the polynomial coefficients in the numerator
    q is a list with the polynomial coefficients (except the first one)
    in the denominator
    The first coefficient of the denominator polynomial is fixed at 1.
    """
    return np.polyval( [1] + p,x) / np.polyval(q, x)

def rational1_2(x, q0, q1):
    return rational(x, [q0, q1])

def rational1_3(x, q0, q1, q2):
    return rational(x, [q0, q1, q2])

def rational1_20(x, q0, q2):
    return rational(x, [q0, 0, q2])

def rational1_4(x, q0, q1, q2, q3):
    return rational(x, [q0, q1, q2, q3])

def rational2_3(x, p0, q0, q1, q2):
    return rational2(x, p0, [q0, q1, q2])

def rational_sum211(x, p0, p1, p2, q0, q1):
    return rational1_3(x, p0, p1, p2) + rational1_2(x, q0, q1)

# Does not converge
# def linear_rational_sum1211(x, n0, n1, p0, p1, p2, q0, q1):
#     return (n0*x + n1) * rational1_3(x, p0, p1, p2) + rational1_2(x, q0, q1)

def rational_sum21(x, p0, p2, q0, q1):
    return rational1_20(x, p0, p2) + rational1_2(x, q0, q1)

# def rational_sum121(x, n0, n1, p0, p2, q0, q1):
#     return (n0*x + n1)*rational1_20(x, p0, p2) + rational1_2(x, q0, q1)

def linear_rational_sum12(x, n0, n1, p0, p2):
    return (n0*x + n1) * rational1_20(x, p0, p2)


def fit01_bkg_rgb(ya, xa, xbkg, bkg_to_fit):

    b = 1
    a1 = (1 - ya * b) / (ya * xa)
    q12 = [a1, b]
    pbr12, _ = curve_fit(rational1_2, xbkg, bkg_to_fit[:, 0], p0=q12)
    pbg12, _ = curve_fit(rational1_2, xbkg, bkg_to_fit[:, 1], p0=q12)
    pbb12, _ = curve_fit(rational1_2, xbkg, bkg_to_fit[:, 2], p0=q12)

    return pbr12, pbg12, pbb12


def fit_bkg_rgb(ya, xa, xbkg, bkg_to_fit):

    pbackr, pbackg, pbackb = [], [], []

    b = 1
    a1 = (1 - ya * b) / (ya * xa)
    q12 = [a1, b]

    for exp in range(3):


        pbackr1, _ = curve_fit(rational1_2, xbkg[exp], bkg_to_fit[exp][:, 0], p0=q12)
        pbackg1, _ = curve_fit(rational1_2, xbkg[exp], bkg_to_fit[exp][:, 1], p0=q12)
        pbackb1, _ = curve_fit(rational1_2, xbkg[exp], bkg_to_fit[exp][:, 2], p0=q12)

        a2 = (1 - ya * b) / (ya * xa ** 2)

        q13r = [a2] + pbackr1.tolist()
        q13g = [a2] + pbackg1.tolist()
        q13b = [a2] + pbackb1.tolist()

        pbackr2, _ = curve_fit(rational1_3, xbkg[exp], bkg_to_fit[exp][:, 0], p0=q13r)
        pbackg2, _ = curve_fit(rational1_3, xbkg[exp], bkg_to_fit[exp][:, 1], p0=q13g)
        pbackb2, _ = curve_fit(rational1_3, xbkg[exp], bkg_to_fit[exp][:, 2], p0=q13b)

        pbackr.append(pbackr2)
        pbackg.append(pbackg2)
        pbackb.append(pbackb2)

    return pbackr, pbackg, pbackb

def fit_bkg_rgb211(ya, xa, xbkg, bkg_to_fit):

    p2 = 0.5
    p1 = (1 - ya * p2) / (ya * xa)
    p12 = [p1, p2]


    pbackr1, _ = curve_fit(rational1_2, xbkg, bkg_to_fit[:, 0], p0=p12)
    pbackg1, _ = curve_fit(rational1_2, xbkg, bkg_to_fit[:, 1], p0=p12)
    pbackb1, _ = curve_fit(rational1_2, xbkg, bkg_to_fit[:, 2], p0=p12)

    p0 = (1 - ya * p2) / (ya * xa ** 2)

    p13r = [p0] + pbackr1.tolist()
    p13g = [p0] + pbackg1.tolist()
    p13b = [p0] + pbackb1.tolist()

    pbackr2, _ = curve_fit(rational1_3, xbkg, bkg_to_fit[:, 0], p0=p13r)
    pbackg2, _ = curve_fit(rational1_3, xbkg, bkg_to_fit[:, 1], p0=p13g)
    pbackb2, _ = curve_fit(rational1_3, xbkg, bkg_to_fit[:, 2], p0=p13b)

    q0 = p1
    q1 = 0.5

    psum211r = pbackr2.tolist()
    psum211r.append(q0)
    psum211r.append(q1)

    psum211g = pbackg2.tolist()
    psum211g.append(q0)
    psum211g.append(q1)

    psum211b = pbackb2.tolist()
    psum211b.append(q0)
    psum211b.append(q1)

    pbackr211, _ = curve_fit(rational_sum211, xbkg, bkg_to_fit[:, 0], p0=psum211r)
    pbackg211, _ = curve_fit(rational_sum211, xbkg, bkg_to_fit[:, 1], p0=psum211g)
    pbackb211, _ = curve_fit(rational_sum211, xbkg, bkg_to_fit[:, 2], p0=psum211b)

    return pbackr211, pbackg211, pbackb211

# Does not converge
# def fit_bkg_rgb1211(ya, xa, xbkg, bkg_to_fit):
#
#     pr, pg, pb = fit_bkg_rgb211(ya, xa, xbkg, bkg_to_fit)
#
#     n0 = -0.01
#     n1 = 0.1
#
#     psum1211r = [n0] + [n1] + pr.tolist()
#     psum1211g = [n0] + [n1] + pg.tolist()
#     psum1211b = [n0] + [n1] + pb.tolist()
#
#     pbr1211, _ = curve_fit(rational_sum1211, xbkg, bkg_to_fit[:, 0], p0=psum1211r)
#     pbg1211, _ = curve_fit(rational_sum1211, xbkg, bkg_to_fit[:, 1], p0=psum1211g)
#     pbb1211, _ = curve_fit(rational_sum1211, xbkg, bkg_to_fit[:, 2], p0=psum1211b)
#
#     return pbr1211, pbg1211, pbb1211

def fit_bkg_120(ya, xa, xbkg, data_to_fit):
    p2 = 0.5
    p0a = (1 - ya * p2) / (ya * xa)
    p0 = [p0a, p2]

    p, _ = curve_fit(rational1_20, xbkg, data_to_fit, p0=p0)

    return p

def fit_bkg_L12(ya, xa, xbkg, data_to_fit):

    p00 = fit_bkg_120(ya, xa, xbkg, data_to_fit)

    n0 = -0.01
    n1 = 0.01

    p0 = [n0] + [n1] + p00.tolist()

    pL12, _ = curve_fit(linear_rational_sum12, xbkg, data_to_fit, p0=p0)

    p = np.insert(pL12, 3, 0)

    return p


def fit_bkg_rgb120(ya, xa, xbkg, data_to_fit):

    p2 = 0.5
    p0 = (1 - ya * p2) / (ya * xa)

    p13r = [p0, p2]
    p13g = [p0, p2]
    p13b = [p0, p2]

    pbr120, _ = curve_fit(rational1_20, xbkg, data_to_fit[:, 0], p0=p13r)
    pbg120, _ = curve_fit(rational1_20, xbkg, data_to_fit[:, 1], p0=p13g)
    pbb120, _ = curve_fit(rational1_20, xbkg, data_to_fit[:, 2], p0=p13b)

    return pbr120, pbg120, pbb120




def fit_bkg_rgb_L12(ya, xa, xbkg, bkg_to_fit):

    pr, pg, pb = zip(*[fit_bkg_120(ya, xa, xbkg, bkg_to_fit[:,rgb]) for rgb in range(3)])
    #pr, pg, pb = fit_bkg_rgb120(ya, xa, xbkg, bkg_to_fit)

    n0 = -0.01
    n1 = 0.01

    psum121r = [n0] + [n1] + pr.tolist()
    psum121g = [n0] + [n1] + pg.tolist()
    psum121b = [n0] + [n1] + pb.tolist()

    pbr, _ = curve_fit(linear_rational_sum12, xbkg, bkg_to_fit[:, 0], p0=psum121r)
    pbg, _ = curve_fit(linear_rational_sum12, xbkg, bkg_to_fit[:, 1], p0=psum121g)
    pbb, _ = curve_fit(linear_rational_sum12, xbkg, bkg_to_fit[:, 2], p0=psum121b)

    pbr = np.insert(pbr, 3, 0)
    pbg = np.insert(pbg, 3, 0)
    pbb = np.insert(pbb, 3, 0)

    return pbr, pbg, pbb


def fit_bkg_rgb21(ya, xa, xbkg, bkg_to_fit):

    p2 = 0.5
    p0 = (1 - ya * p2) / (ya * xa)

    p13r = [p0, p2]
    p13g = [p0, p2]
    p13b = [p0, p2]

    pbackr2, _ = curve_fit(rational1_20, xbkg, bkg_to_fit[:, 0], p0=p13r)
    pbackg2, _ = curve_fit(rational1_20, xbkg, bkg_to_fit[:, 1], p0=p13g)
    pbackb2, _ = curve_fit(rational1_20, xbkg, bkg_to_fit[:, 2], p0=p13b)

    q0 = p0
    q1 = 0.5

    psum32r = pbackr2.tolist()
    psum32r.append(q0)
    psum32r.append(q1)

    psum32g = pbackg2.tolist()
    psum32g.append(q0)
    psum32g.append(q1)

    psum32b = pbackb2.tolist()
    psum32b.append(q0)
    psum32b.append(q1)

    pbackr21, _ = curve_fit(rational_sum21, xbkg, bkg_to_fit[:, 0], p0=psum32r)
    pbackg21, _ = curve_fit(rational_sum21, xbkg, bkg_to_fit[:, 1], p0=psum32g)
    pbackb21, _ = curve_fit(rational_sum21, xbkg, bkg_to_fit[:, 2], p0=psum32b)

    return pbackr21, pbackg21, pbackb21


files = ['/Users/rattie/Data/Eclipse/Stack_mean_aligned_series/aligned/stack_mean_C_exp250.tiff',
         '/Users/rattie/Data/Eclipse/Stack_mean_aligned_series/aligned/stack_mean_C_exp30.tiff',
         '/Users/rattie/Data/Eclipse/Stack_mean_aligned_series/aligned/stack_mean_C_exp1.tiff']

# Exposure times
nexp = 3
times_array = np.array([250, 30, 1])
sum_times = times_array.sum()
weights = np.reshape(times_array/sum_times, [nexp, 1, 1, 1])
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
# Normalize to the exposure time, so we get the same intensity per second + noise
tarray3 = np.reshape(times_array, [nexp, 1, 1, 1])
images_exp = images * tarray3 # * (2**16 -1) / (2**14 -1)
images_expf = images_exp/maxval
# Clip saturated pixels in the exposure-normalized images
images_expf_masked = images_expf.copy()
images_expf_masked[images_expf > 1] = 1



## Remap to polar coordinates
# Coordinate of the center
center = (2120, 1361)
radius = 323
# Plot one image and the circle to check the center and radius
#circle1 = plt.Circle((0, 0), 0.2, color='r')
plt.figure(0, figsize=(15,11))
ax1 = plt.gcf().add_subplot(111)
ax1.imshow(imagesf[0,...])

ax1.add_artist(plt.Circle(center, radius, color='green', fill=False, linewidth=1))
ax1.add_artist(plt.Circle(center, 330, color='yellow', fill=False))
ax1.axis([1700,2500,1000,1700])
plt.title('1/250s')
plt.tight_layout()

plt.figure(1, figsize=(15,11))
ax1 = plt.gcf().add_subplot(111)
ax1.imshow(imagesf[0,...])

ax1.add_artist(plt.Circle(center, radius, color='green', fill=False, linewidth=1))
ax1.add_artist(plt.Circle(center, 330, color='yellow', fill=False))
ax1.axis([1700,2500,1000,1700])
plt.title('1/250s')
plt.tight_layout()



# Display the 3 images in RGB.
plt.figure(2, figsize=(19,10))
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
plt.imshow(images_expf_masked[0,...])
plt.subplot(235)
plt.imshow(images_expf_masked[1,...])
plt.subplot(236)
plt.imshow(images_expf_masked[2,...])
plt.tight_layout()

plt.savefig('/Users/rattie/Data/Eclipse/figures/imagesRGB_and_exp_normalized.png')


pimages0_rgb = np.zeros([nexp, ny, nx, 3])
pimages_rgb = np.zeros([nexp, ny, nx, 3])
for i in range(0,3):
    for k in range(0,3):
        # Azimuthal average
        pimages0_rgb[i, :, :, k] = cv2.linearPolar(imagesf[i, :, :, k], center, nx, cv2.INTER_LANCZOS4 + cv2.WARP_FILL_OUTLIERS)
        pimages_rgb[i,:,:,k] = cv2.linearPolar(images_exp[i,:,:,k], center, nx, cv2.INTER_LANCZOS4 + cv2.WARP_FILL_OUTLIERS)

pimages0_rgb = np.clip(pimages0_rgb, 0,1)
pimages_rgbf = np.clip(pimages_rgb / maxval, 0, 1)


# Take the median over a selection on the y-axis.
#ymin = 1100
ymin = 990
ymax = 1100 # 1301

ymin2 = 970
ymax2 = 1040

zeroloc1 = np.where(pimages_rgb[2, ymin, :, 0] == 0)[0].min()
zeroloc2 = np.where(pimages_rgb[2, ymax, :, 0] == 0)[0].min()
zeroloc = np.min([zeroloc1, zeroloc2])

minima = pimages_rgb[:,ymin:ymax, :, 0].argmin(axis=1)
minima2 = pimages_rgb[:,ymin2:ymax2, :, 0].argmin(axis=1)

az_avgs_rgb = np.zeros([3, nx, 3])
az_avgs_rgbM = np.zeros([nx, 3])

for i in range(minima.shape[1]):
    for exp in range(3):
        az_avgs_rgb[exp,i,:] = pimages_rgb[exp, minima[exp,i]+ymin, i, :]

# for i in range(430, minima.shape[1]):
#     for exp in range(3):
#         az_avgs_rgb[exp,i,:] = pimages_rgb[exp, minima2[exp,i]+ymin2, i, :]

for i in range(460, minima.shape[1]):
    for exp in range(3):
        for rgb in range(3):
            az_avgs_rgb[exp, i, rgb] = pimages_rgb[exp, ymin2:ymax2, i, rgb].mean()

az_avgs_gb = 0.5 * (az_avgs_rgb[:,:, 1] + az_avgs_rgb[:,:, 2])



## Alternative method: Take only the minimum values up to a given limit. Take mean beyond it.
for i in range(460):
    # From 0 to 430, take the minima of exp=0
    az_avgs_rgbM[i,:] = pimages_rgb[0, minima[0,i]+ymin, i, :]

# for i in range(460, minima.shape[1]):
#     az_avgs_rgbM[i, :] = pimages_rgb[1, minima2[1,i]+ymin2, i, :]

for i in range(460, minima.shape[1]):
        az_avgs_rgbM[i, 0] = pimages_rgb[1, ymin2:ymax2, i, 0].mean()
        az_avgs_rgbM[i, 1] = pimages_rgb[1, ymin2:ymax2, i, 1].mean()
        az_avgs_rgbM[i, 2] = pimages_rgb[1, ymin2:ymax2, i, 2].mean()


# Create a unique profile to fit, this will eventually be subtracted to all exposures.
az_avgs_rgb2 = az_avgs_rgbM.copy()
# # In [0 - 500], take profile of 1/250s
# az_avgs_rgb2[:, 0:360, :] = az_avgs_rgbM[0, 0:360, :][np.newaxis,...]
# # In [500 - end], take profile of 1/30s
# az_avgs_rgb2[:, 360:, :] = az_avgs_rgbM[1, 360:, :][np.newaxis,...]

# The maximum value to normalize the profiles should be taken from the 1/250s-normalized: the only one not clipped.
max_bkg = az_avgs_rgb[0,...].max()
max_bkg_gb = az_avgs_gb[0,:].max()
max_bkg2 = az_avgs_rgb[0,:,:].max(axis=0)
max_bkg2_loc = az_avgs_rgb[0,:,:].argmax(axis=0)

max_bkg22 = az_avgs_rgb2.max(axis=0)
max_bkg22_loc = az_avgs_rgb2.argmax(axis=0)

az_avgs_rgbN = az_avgs_rgb/max_bkg2[np.newaxis, np.newaxis,:]
az_avgs_gbN = 0.5 * (az_avgs_rgbN[..., 1] + az_avgs_rgbN[..., 2])


zoom1 = [200, 1800, 800, 1200]
zoom2 = [200, 1800, 1e-3, 1.1]
#1/250s
plt.figure(3, figsize=(19,10))
ax1 = plt.subplot(231)
ax1.imshow(np.clip(pimages0_rgb[0,...], 0, 1), origin='lower')
ax1.axhline(y=ymin, ls='--', linewidth=1, color='yellow')
ax1.axhline(y=ymax, ls='--', linewidth=1, color='yellow')
ax1.axhline(y=ymin2, ls='--', linewidth=1, color='white')
ax1.axhline(y=ymax2, ls='--', linewidth=1, color='white')
ax1.set_title('Normalized from 1/250s; Polar transform')
ax1.plot(minima[0,:]+ymin, marker='.',  color='red', ls='', ms=2, zorder=4)
ax1.plot(minima[1,:]+ymin, marker='.',  color='orange', ls='', ms=2, zorder=3)
ax1.axis(zoom1)

ax2 = plt.subplot(234)
ax2.semilogy(az_avgs_rgb[0,:,0]/max_bkg2[0], 'r-', linewidth=2, label='red 1/250s norm.')
ax2.semilogy(az_avgs_rgb[0,:,1]/max_bkg2[0], 'g--', linewidth=2, label='green 1/250s norm.')
ax2.semilogy(az_avgs_rgb[0,:,2]/max_bkg2[0], 'b-.', linewidth=2, label='blue 1/250s norm.')
ax2.axis(zoom2)
plt.legend()


ax3 = plt.subplot(232)
ax3.imshow(np.clip(pimages0_rgb[1,...], 0, 1), origin='lower')
ax3.axhline(y=ymin, ls='--', linewidth=1, color='yellow')
ax3.axhline(y=ymax, ls='--', linewidth=1, color='yellow')
ax3.axhline(y=ymin2, ls='--', linewidth=1, color='white')
ax3.axhline(y=ymax2, ls='--', linewidth=1, color='white')
ax3.set_title('Normalized from 1/30s; Polar transform')
ax3.axis(zoom1)

ax4 = plt.subplot(235)
ax4.semilogy(az_avgs_rgb[1,:,0]/max_bkg2[0], 'r-', linewidth=2.5, label='red 1/30s norm.')
ax4.semilogy(az_avgs_rgb[1,:,1]/max_bkg2[0], 'g--', linewidth=2, label='green 1/30s norm.')
ax4.semilogy(az_avgs_rgb[1,:,2]/max_bkg2[0], 'b-.', linewidth=1.5, label='blue 1/30s norm. ')
ax4.axis(zoom2)
plt.legend()


ax5 = plt.subplot(233)
ax5.imshow(np.clip(pimages0_rgb[2,...], 0, 1), origin='lower')
ax5.axhline(y=ymin, ls='--', linewidth=1, color='yellow')
ax5.axhline(y=ymax, ls='--', linewidth=1, color='yellow')
ax5.axhline(y=ymin2, ls='--', linewidth=1, color='white')
ax5.axhline(y=ymax2, ls='--', linewidth=1, color='white')
ax5.set_title('Normalized from 1/30s; Polar transform')
ax5.axis(zoom1)

ax6 = plt.subplot(236)
ax6.semilogy(az_avgs_rgb[2,:,0]/max_bkg2[0], 'r-', linewidth=2.5, label='red 1s norm.')
ax6.semilogy(az_avgs_rgb[2,:,1]/max_bkg2[0], 'g--', linewidth=2, label='green 1s norm.')
ax6.semilogy(az_avgs_rgb[2,:,2]/max_bkg2[0], 'b-.', linewidth=1.5, label='blue 1s norm. ')
ax6.axis(zoom2)
plt.legend()

plt.tight_layout()


# Fit from x1 to x2=zeroloc.
# Start with red channel
x0 = [400, 420, 600] #int(radius)
x1 = [700, 1700, 1700]#zeroloc
x00 = int(radius)

x12 = [800, 1000, 1000]
x13 = [900, 1000, 1000]

bkg_to_fit = [az_avgs_rgbN[exp, x0[exp]:x1[exp], :] for exp in range(3)]
bkg_to_fit_gb = [az_avgs_gbN[exp, x0[exp]:x1[exp]].copy() for exp in range(3)]
#bkg_to_fit2 = az_avgs_rgb2[x02:x12, :].copy()/max_bkg22[np.newaxis,:]
bkg_to_fit2 = [az_avgs_rgbN[exp, x0[exp]:x12[exp], :] for exp in range(3)]
bkg_to_fit3 = [az_avgs_rgbN[exp, x0[exp]:x13[exp], :] for exp in range(3)]


xbkg = [np.arange(x0[exp], x1[exp])-x00 for exp in range(3)]
xbkg22 = [np.arange(x0[exp], x12[exp])-x00 for exp in range(3)]
xbkg3 = [np.arange(x0[exp], x13[exp])-x00 for exp in range(3)]

# Fit the inverse polynomial and not the high-degree polynomial nonsense!!

xa = 380 - x00
ya = 0.15#0.1

# Fit 3 exposures with profile: 1/(ax2+bx+c)+ 1/(dx+e)
#pbr, pbg, pbb = zip(*[fit_bkg_rgb211(ya, xa, xbkg[exp], bkg_to_fit[exp]) for exp in range(3)])
pr, pg, pb = zip(*[fit_bkg_rgb_L12(ya, xa, xbkg[exp], bkg_to_fit[exp]) for exp in range(3)])
pgb = [fit_bkg_L12(ya, xa, xbkg[exp], bkg_to_fit_gb[exp]) for exp in range(3)]
pr2, pg2, pb2 = zip(*[fit_bkg_rgb_L12(ya, xa, xbkg22[exp], bkg_to_fit2[exp]) for exp in range(3)])
pr3, pg3, pb3 = zip(*[fit_bkg_rgb_L12(ya, xa, xbkg3[exp], bkg_to_fit3[exp]) for exp in range(3)])

xbkg2 = np.arange(300, zeroloc)- x00

# fit_bkg_red = [1/np.poly1d(pbr[exp][0:3])(xbkg2) + 1/np.poly1d(pbr[exp][3:])(xbkg2) for exp in range(3)]
# fit_bkg_green = [1/np.poly1d(pbg[exp][0:3])(xbkg2) + 1/np.poly1d(pbg[exp][3:])(xbkg2) for exp in range(3)]
# fit_bkg_blue = [1/np.poly1d(pbb[exp][0:3])(xbkg2) + 1/np.poly1d(pbb[exp][3:])(xbkg2) for exp in range(3)]

# fit_bkg_red2 = [np.poly1d(pbr2[exp][0:2])(xbkg2)/np.poly1d(pbr2[exp][2:5])(xbkg2) + 1/np.poly1d(pbr2[exp][5:])(xbkg2) for exp in range(3)]
# fit_bkg_green2 = [np.poly1d(pbg2[exp][0:2])(xbkg2)/np.poly1d(pbg2[exp][2:5])(xbkg2) + 1/np.poly1d(pbg2[exp][5:])(xbkg2) for exp in range(3)]
# fit_bkg_blue2 = [np.poly1d(pbb2[exp][0:2])(xbkg2)/np.poly1d(pbb2[exp][2:5])(xbkg2) + 1/np.poly1d(pbb2[exp][5:])(xbkg2) for exp in range(3)]

fit_bkg_red = [np.poly1d(pr[exp][0:2])(xbkg2) * 1/np.poly1d(pr[exp][2:5])(xbkg2)  for exp in range(3)]
fit_bkg_green = [np.poly1d(pg[exp][0:2])(xbkg2) * 1/np.poly1d(pg[exp][2:5])(xbkg2) for exp in range(3)]
fit_bkg_blue = [np.poly1d(pb[exp][0:2])(xbkg2) * 1/np.poly1d(pb[exp][2:5])(xbkg2) for exp in range(3)]

#fit_bkg_gb = [np.poly1d(pgb[exp][0:2])(xbkg2) * 1/np.poly1d(pgb[exp][2:5])(xbkg2) for exp in range(3)]

fit_bkg_red2 = [np.poly1d(pr2[exp][0:2])(xbkg2) * 1/np.poly1d(pr2[exp][2:5])(xbkg2)  for exp in range(3)]
fit_bkg_green2 = [np.poly1d(pg2[exp][0:2])(xbkg2) * 1/np.poly1d(pg2[exp][2:5])(xbkg2) for exp in range(3)]
fit_bkg_blue2 = [np.poly1d(pb2[exp][0:2])(xbkg2) * 1/np.poly1d(pb2[exp][2:5])(xbkg2) for exp in range(3)]

fit_bkg_red3 = [np.poly1d(pr3[exp][0:2])(xbkg2) * 1/np.poly1d(pr3[exp][2:5])(xbkg2)  for exp in range(3)]
fit_bkg_green3 = [np.poly1d(pg3[exp][0:2])(xbkg2) * 1/np.poly1d(pg3[exp][2:5])(xbkg2) for exp in range(3)]
fit_bkg_blue3 = [np.poly1d(pb3[exp][0:2])(xbkg2) * 1/np.poly1d(pb3[exp][2:5])(xbkg2) for exp in range(3)]


zoomin = [250, 1700, 1e-3, 1.1]
plot = plt.semilogy

# Test fitting ranges
plt.figure(4, figsize=(19,10))
ax1 = plt.subplot(131)

plot(az_avgs_rgb[0,:,0]/max_bkg2[0], 'r-', label='red 1/250s norm.', linewidth=2.5)
plot(az_avgs_rgb[0,:,1]/max_bkg2[1], 'g-', label='green 1/250s norm', linewidth=2, alpha=0.4)
plot(az_avgs_rgb[0,:,2]/max_bkg2[2], 'b-', label='blue 1/250s norm', linewidth=2, alpha=0.4)
plot(xbkg2+x00, fit_bkg_red[0], 'k-', label='red 1/250s linear rational', linewidth=1.5)
plot(xbkg2+x00, fit_bkg_red2[0], ls='--', color='gray', label='red 1/250s linear rational', linewidth=2)
plot(xbkg2+x00, fit_bkg_red2[0], 'k-.', label='red 1/250s linear rational', linewidth=3)


ax1.axvline(x=x0[0], ls='-', linewidth=1, color='black')
ax1.axvline(x=x1[0], ls='-', linewidth=1, color='black')
ax1.axvline(x=x12[0], ls='--', linewidth=1, color='gray')
ax1.axvline(x=x13[0], ls='-.', linewidth=1, color='black')

plt.axis(zoomin)
plt.legend()

ax2 = plt.subplot(132)
plot(az_avgs_rgb[0,:,1]/max_bkg2[1], 'g-', label='green 1/250s norm', linewidth=2)
plot(xbkg2+x00, fit_bkg_green[0], 'k-', label='green 1/250s linear rational', linewidth=1.5)
plot(xbkg2+x00, fit_bkg_green2[0], ls='--', color='gray', label='green 1/250s linear rational', linewidth=2)
plot(xbkg2+x00, fit_bkg_green3[0], 'k-.', label='green 1/250s linear rational', linewidth=3)

ax2.axvline(x=x0[0], ls='-', linewidth=1, color='black')
ax2.axvline(x=x1[0], ls='-', linewidth=1, color='black')
ax2.axvline(x=x12[0], ls='--', linewidth=1, color='gray')
ax2.axvline(x=x13[0], ls='-.', linewidth=1, color='black')


plt.axis(zoomin)
plt.legend()

ax3 = plt.subplot(133)
plot(az_avgs_rgb[0,:,2]/max_bkg2[2], 'b-', label='blue 1/250s norm', linewidth=2)
plot(xbkg2+x00, fit_bkg_blue[0], 'k-', label='blue 1/250s linear rational', linewidth=1.5)
plot(xbkg2+x00, fit_bkg_blue2[0], ls='--', color='gray', label='blue 1/250s linear rational', linewidth=2)
plot(xbkg2+x00, fit_bkg_blue3[0], 'k-.', label='blue 1/250s linear rational', linewidth=2.5)

ax3.axvline(x=x0[0], ls='-', linewidth=1, color='black')
ax3.axvline(x=x1[0], ls='-', linewidth=1, color='black')
ax3.axvline(x=x12[0], ls='--', linewidth=1, color='gray')
ax3.axvline(x=x13[0], ls='-.', linewidth=1, color='black')

plt.axis(zoomin)
plt.legend()

plt.tight_layout()


plt.figure(5, figsize=(19,10))
ax1 = plt.subplot(131)
#plot(az_avgs_gb[0,:]/max_bkg_gb, ls='--', color='gray', label='blue 1/250s original', linewidth=2)

#plot(az_avgs_rgb2[:,0]/max_bkg22[0], 'r--', label='red 1/250s norm. unique', linewidth=2.5)
plot(az_avgs_rgb[0,:,0]/max_bkg2[0], 'r-', label='red 1/250s norm.', linewidth=2.5)
plot(az_avgs_rgb[0,:,1]/max_bkg2[1], 'g--', label='green 1/250s norm', linewidth=2)
plot(az_avgs_rgb[0,:,2]/max_bkg2[2], 'b-.', label='blue 1/250s norm.', linewidth=1.5)

#plt.plot(xbkg2+x00, fit_bkg_red[0], 'k--', label='red 1/250s rational', linewidth=1)
plot(xbkg2+x00, fit_bkg_red[0], 'k-', label='red 1/250s linear rational', linewidth=1.5)
plot(xbkg2+x00, fit_bkg_green[0], 'k--', label='green 1/250s linear rational', linewidth=1.5)
plot(xbkg2+x00, fit_bkg_blue[0], 'k-.', label='blue 1/250s linear rational', linewidth=1.5)
#plot(xbkg2+x00, fit_bkg_gb[0], 'c:', label='average [g;b] 1/250s linear rational', linewidth=1.5)
#plot(xbkg2+x00, fit_bkg_red2, 'm--', label='red 1/250s unique linear rational', linewidth=1.5)

#plt.semi
ax1.axvline(x=x0[0], ls='-', linewidth=1, color='black')
ax1.axvline(x=x1[0], ls='-', linewidth=1, color='black')

plt.axis(zoomin)
plt.legend()

ax2 = plt.subplot(132)
plot(az_avgs_gb[1,:]/max_bkg_gb, ls='-.', color='gray', label='blue 1/30s original', linewidth=2)
plot(az_avgs_rgb[1,:,0]/max_bkg2[0], 'r-', label='red 1/30s norm.', linewidth=2.5)
plot(az_avgs_rgb[1,:,1]/max_bkg2[1], 'g--', label='green 1/30s norm.', linewidth=2)
plot(az_avgs_rgb[1,:,2]/max_bkg2[2], 'b-.', label='blue 1/30s norm.', linewidth=1.5)

plot(xbkg2+x00, fit_bkg_red[1], 'k-', label='red 1/30s norm. linear rational', linewidth=1.5)
# plot(xbkg2+x00, fit_bkg_blue[1], 'k--', label='blue 1/30s linear rational', linewidth=1.5)
# plot(xbkg2+x00, fit_bkg_gb[1], 'k:', label='average [g;b] 1/30s linear rational', linewidth=1.5)
plot(xbkg2+x00, fit_bkg_red2, 'k--', label='red 1/250s linear rational', linewidth=1.5)

plt.axis(zoomin)
plt.legend()

ax3 = plt.subplot(133)
plot(az_avgs_gb[2,:]/max_bkg_gb, ls='--', color='gray', label='blue 1s original', linewidth=2)
plot(az_avgs_rgb[2,:,0]/max_bkg2[0], 'r-', label='red 1s', linewidth=2)
plot(az_avgs_rgb[2,:,1]/max_bkg2[1], 'g-', label='green 1s', linewidth=2)
plot(az_avgs_rgb[2,:,2]/max_bkg2[2], 'b-', label='blue 1s', linewidth=2)

plot(xbkg2+x00, fit_bkg_red[2], 'k-', label='red 1s linear rational', linewidth=1.5)
# plot(xbkg2+x00, fit_bkg_blue[2], 'k--', label='blue 1s linear rational', linewidth=1.5)
# plot(xbkg2+x00, fit_bkg_gb[2], 'k:', label='average [g;b] 1s linear rational', linewidth=1.5)
plot(xbkg2+x00, fit_bkg_red2, 'k--', label='red 1/250s linear rational', linewidth=1.5)

plt.axis(zoomin)
plt.legend()

plt.tight_layout()



plt.savefig('/Users/rattie/Data/Eclipse/figures/polar_background_fit_RGB.png')


# From the plot above, we want the 2nd fit (pbr2, pbg2, pbb2) for the 1st exposure
# and the 2nd fit for the 2nd and 3rd exposure.
pbr3 = [pbr2[0], pbr[1], pbr[2]]
pbg3 = [pbg2[0], pbg[1], pbg[2]]
pbb3 = [pbb2[0], pbb[1], pbb[2]]
# Pack all that into another list for use in list comprehension. It makes the code more concise.
pbkg = [pbr3, pbg3, pbb3]


## Use the fit to reconstruct a background image to subtract.
# Make a radius map: at each pixel, the value is the distance to disc center
r = get_radius_array(center, nx, ny)
rr = r - x0[0]
#rr[rr < radius-x0] = 0

imfit_rgb = [[1/np.poly1d(pbkg[rgb][exp])(rr) for rgb in range(3)] for exp in range(3)]
imfit_rgb = [np.moveaxis(np.array(imfit_rgb[exp]), 0, -1) for exp in range(3)]
# The fitted background was originally from data normalized to max_bkg
#imfit_rgb2 = np.array(imfit_rgb) * max_bkg
max_bkg3 = max_bkg2[np.newaxis, np.newaxis, np.newaxis, :]
imfit_rgb2 = np.array(imfit_rgb) * max_bkg3

for exp in range(3):
    imfit_rgb[exp][r <= radius, :] =0
    imfit_rgb2[exp][r <= radius, :] = 0


r2 = r.copy()
r2[r > radius] = 0

images_exp_back = images_exp - imfit_rgb2
# Once the background is subtracted, need to define a new value to normalize the images.
maxval_exp_back = images_exp_back[0, 1080:1120, 1900:1960, 0].max()
images_exp_backf = np.clip(images_exp_back / maxval_exp_back, 0, 1)
images_exp_backf[0, r <= int(radius), :] = 0
# Original images, normalized and clipped
imagesf2 = np.clip(imagesf,0,1)

axis_zoom=[1000, 3200, 900, 1800]

plt.figure(6, figsize=(19,10))
ax1 = plt.gcf().add_subplot(221)
ax2 = plt.gcf().add_subplot(222)
ax3 = plt.gcf().add_subplot(223)
ax4 = plt.gcf().add_subplot(224)

ax1.imshow(np.clip(imagesf2[0,...]*2, 0, 1), origin='lower')
#ax1.add_artist(plt.Circle(center, radius, color='green', fill=False, linewidth=1, ls='--'))
ax1.axis(axis_zoom)
ax1.set_title('1/250s norm.')

ax2.imshow(np.clip(images_exp_backf[0,...]*2, 0, 1), origin='lower')
#ax2.add_artist(plt.Circle(center, radius, color='green', fill=False, linewidth=1, ls='--'))
ax2.axis(axis_zoom)
ax2.set_title('1/250s norm. backg removed')
ax3.imshow(np.clip(imagesf2[0,...]*256, 0, 1), origin='lower')
ax3.axis(axis_zoom)
ax3.set_title('1/250s norm. intensity 20x')
ax4.imshow(np.clip(images_exp_backf[0,...]*256, 0, 1), origin='lower')
ax4.axis(axis_zoom)
ax4.set_title('1/250s norm. backg removed, intensity 20x')
plt.tight_layout()

plt.savefig('/Users/rattie/Data/Eclipse/figures/images_30s_background_removed_RGB.png')

# Fit streamer rays
# remove background from polar image. Build a bacground directly in polar space instead of converting from cartesian to polar
#pbackground = np.array([cv2.linearPolar(imfit_rgb2[exp], center, nx, cv2.INTER_LANCZOS4 + cv2.WARP_FILL_OUTLIERS) for exp in range(3)])
r1d = np.arange(nx)[np.newaxis,:]
rpolar = np.repeat(r1d, ny, axis=0) - x0[0]
pbackground = [[1/np.poly1d(pbkg[rgb][exp])(rpolar) for rgb in range(3)] for exp in range(3)]
pbackground = np.array([np.moveaxis(np.array(pbackground[exp]), 0, -1) for exp in range(3)])* max_bkg3

# plt.figure(4, figsize=(10,10))
# ax1 = plt.gcf().add_subplot(111)
# ax1.imshow(np.clip(pbackground[0,...]/max_bkg3[0,...], 0, 1), origin='lower')
# plt.tight_layout()

pimages_rgb_back0 = pimages_rgb - pbackground
pimages_rgb_backf = np.clip(pimages_rgb_back0 / maxval_exp_back, 0 , 1)
pimages_rgb_backf[pimages_rgb==0] = 0

zoom5 = [0, ny, 0, ny]

plt.figure(7, figsize=(19,10))
plt.subplot(231)
plt.imshow(np.clip(pimages0_rgb[0,...]*5, 0, 1), origin='lower')
plt.axis(zoom5)
plt.title('1/250s intensity x5')

plt.subplot(232)
plt.imshow(np.clip(pimages0_rgb[1,...]*1, 0, 1), origin='lower')
plt.axis(zoom5)
plt.title('1/30s intensity x1')

plt.subplot(233)
plt.imshow(np.clip(pimages0_rgb[2,...], 0, 1), origin='lower')
plt.axis(zoom5)
plt.title('1s')


plt.subplot(234)
plt.imshow(np.clip(pimages_rgb_backf[0,...]*5, 0, 1), origin='lower')
plt.axis(zoom5)
plt.title('1/250s intensity x5 backg removed')

plt.subplot(235)
plt.imshow(np.clip(pimages_rgb_backf[1,...]*500, 0, 1), origin='lower')
plt.axis(zoom5)
plt.title('1/30s intensity x500 backg removed')

plt.subplot(236)
plt.imshow(np.clip(pimages_rgb_backf[2,...]*2000, 0, 1), origin='lower')
plt.axis(zoom5)
plt.title('1s intensity x2000 backg removed')

plt.tight_layout()

plt.savefig('/Users/rattie/Data/Eclipse/figures/polar_30s_background_removed_RGB.png')


plt.figure(6, figsize=(19,10))
ax1 = plt.subplot(131)
plt.imshow(np.clip(pimages_rgb_backf[0,...]*100, 0, 1), origin='lower')
plt.axis([0, ny, 0, ny])
plt.title('1/250s norm. bkg removed intensity x5')
ax2 = plt.subplot(132)
plt.imshow(np.clip(pimages_rgb_backf[1,...]*100, 0, 1), origin='lower')
plt.axis([0, ny, 0, ny])
plt.title('1/30s norm. bkg removed intensity x1000')
ax3= plt.subplot(133)
plt.imshow(np.clip(pimages_rgb_backf[2,...]*1000, 0, 1), origin='lower')
plt.axis([0, ny, 0, ny])
plt.title('1s norm. bkg removed intensity x2000')

plt.tight_layout()

# TODO: In streamer_to_fit_rgb, need to find proper max() since the values from pimages_rgbf at 1/30th are saturated and thus clipped.
# TODO: After doing the above, do proper fit with x0 = radius

ys = np.array([300, 1500, 1950])
hwidth = 100
ystreamers = [np.arange(y-hwidth, y+hwidth) for y in ys]
nstreamers = len(ystreamers)

for i in range(nstreamers):
    plt.axhline(y=ystreamers[i][0], ls='--', linewidth=1, color='yellow')
    plt.axhline(y=ystreamers[i][-1], ls='--', linewidth=1, color='yellow')

maxima = [pimages_rgb[:,ystreamers[s][0]:ystreamers[s][-1], :, 0].argmax(axis=1) for s in range(nstreamers)]

for s in range(3):
    ax2.plot(maxima[s][1,:]+ystreamers[s][0], marker='.',  color='cyan', ls='', ms=1)
    ax3.plot(maxima[s][2, :] + ystreamers[s][0], marker='.', color='cyan', ls='', ms=1)

# Make a list of 3 arrays
streamers_rgb = [np.zeros([3, nx, 3])]*3

for s in range(3):
    for i in range(minima.shape[1]):
        for j in range(3):
            streamers_rgb[s][j,i,:] = pimages_rgb[j, maxima[s][j,i]+ystreamers[s][0], i, :]

x0 = 324 #int(radius)
x1 = 1500
xstreamers = np.arange(x0,x1) -x0

# Extract streamer data from the normalized polar images in RGB, get the exposure at 30th.
# dimensions are [exp][streamer region over y-axis][x,y,rgb]
streamer_data_rgb = [[pimages_rgb[exp, ycoords, :, :].mean(axis=0) for ycoords in ystreamers] for exp in range(3)]
# Make an intermediate plot to see what the profiles look like

# Select over the x-axis, dimensions are [exp][streamer region over y-axis][x0:x1,y,rgb]
streamer_to_fit_rgb = [[data[x0:x1, :] for data in streamer_data_rgb[exp]] for exp in range(3)]

# Normalization by what: value at radius or close to it? or just the absolute max?
red_streamer = [streamer_to_fit_rgb[exp][0][:,0]/streamer_data_rgb[0][0][:,0].max() for exp in range(3)]

# plt.figure(6)
# plt.plot(xstreamers+x0, red_streamer[0])

#
# b = 1
# ya = 0.05
# xa = 460
# a1 = (1 - ya*b)/(ya*(xa-x0))
# q12=[a1, b]
# p_red0, p_cov_red0= curve_fit(rational1_2, xstreamers, red_streamer[0], p0=q12)
#
# a2 = (1 - ya*b)/(ya*(xa-x0)**2)
# q13= [a2] + p_red0.tolist()
# p_red1, p_cov_red0= curve_fit(rational1_3, xstreamers, red_streamer[0], p0=q13)
#
# fit_streamer0 = 1/np.poly1d(p_red0)(xstreamers)
# fit_streamer1 = 1/np.poly1d(p_red1)(xstreamers)
#
#
# fig = plt.figure(7, figsize=(19,10))
# ax1 = fig.add_subplot(121)
# ax2 = fig.add_subplot(122)
# ax1.imshow(pimage_rgb_back_scaled, origin='lower')
# ax1.set_title('1/30s norm. backg removed, intensity scaled up 5x')
# ax1.set_xlim(xmin=0, xmax=zeroloc+100)
#
# ax1.axhline(y=ystreamers[0][0], xmin=x0/ax1.get_xlim()[1], xmax=x1/ax1.get_xlim()[1], ls='--', linewidth=1, color='red')
# ax1.axhline(y=ystreamers[0][-1], xmin=x0/ax1.get_xlim()[1], xmax=x1/ax1.get_xlim()[1], ls='--', linewidth=1, color='red')
# ax1.axhline(y=ystreamers[1][0], xmin=x0/ax1.get_xlim()[1], xmax=x1/ax1.get_xlim()[1], ls='--', linewidth=1, color='cyan')
# ax1.axhline(y=ystreamers[1][-1], xmin=x0/ax1.get_xlim()[1], xmax=x1/ax1.get_xlim()[1], ls='--', linewidth=1, color='cyan')
# ax1.axhline(y=ystreamers[2][0], xmin=x0/ax1.get_xlim()[1], xmax=x1/ax1.get_xlim()[1], ls='--', linewidth=1, color='orange')
# ax1.axhline(y=ystreamers[2][-1], xmin=x0/ax1.get_xlim()[1], xmax=x1/ax1.get_xlim()[1], ls='--', linewidth=1, color='orange')
#
# ax2.plot(streamer_data_rgb[0][0][:,0]/streamer_data_rgb[0][0][:,0].max(), 'r-')
# ax2.plot(streamer_data_rgb[1][0][:,0]/streamer_data_rgb[0][0][:,0].max(), ls='--', color='orange')
# ax2.set_xlim(xmin=0, xmax=x1)
# ax2.set_ylim(ymin=0, ymax=1.5)
# # ax2.plot(xstreamers, streamer_to_fit_rgb[0][:,1], 'g-')
# # ax2.plot(xstreamers, streamer_to_fit_rgb[0][:,2], 'b-')
#
# #ax2.plot(streamer_data_rgb[1][:,0]/streamer_to_fit_rgb[1][:,0].max(), 'c-')
# # ax2.plot(xstreamers, streamer_to_fit_rgb[1][:,1], 'g--')
# # ax2.plot(xstreamers, streamer_to_fit_rgb[1][:,2], 'b--')
#
# #ax2.plot(streamer_data_rgb[2][:,0]/streamer_to_fit_rgb[2][:,0].max(), ls='-', color='orange')
# # ax2.plot(xstreamers, streamer_to_fit_rgb[2][:,1], 'g--')
# # ax2.plot(xstreamers, streamer_to_fit_rgb[2][:,2], 'b--')
# ax2.axvline(x=x0, ls='--', linewidth=1, color='black')
#
#
# ax2.plot(xstreamers + x0, fit_streamer0, 'b-')
# ax2.plot(xstreamers + x0, fit_streamer1, 'k--')
#
# plt.tight_layout()
#
#
# plt.savefig('/Users/rattie/Data/Eclipse/figures/streamer_profiles.png')

# # We defined the radius array above as
# # r = get_radius_array(center, nx, ny)
#
# # Get the images normalized to exposure and apply the radial filter based on what is fitted above
# rr = r - radius
# rr[rr < 0] = 0
# af = 0.3
# rfilter = 1/(af*rr+1)
#
# rescale = 1/rfilter[r==324].mean()
# # maxval was the max physical value in the 250th exposure image.
# # => maxval * 250 is the maximum physical value in the exposure-normalized image.
# # So I could normalize to this values and clip within [0,1] (for export, multiply by max16)
# # Due to rfilter that multiply values close to radius (rr = 4) by ~1.35, we first normalize by maxval*1.35
# images_exp_backf = np.clip(images_exp_back/(maxval_back * 250), 0,1)
#
# images_exp_backf_f = np.clip(images_exp_backf / rfilter[np.newaxis,:,:,np.newaxis], 0, 1)
#
# fig = plt.figure(8, figsize=(19,8))
# ax1=fig.add_subplot(231)
# ax1.imshow(images_exp_backf[0,...])
# ax2=fig.add_subplot(232)
# ax2.imshow(images_exp_backf[1,...])
# ax4=fig.add_subplot(234)
# ax4.imshow(images_exp_backf_f[0,...])
# ax4=fig.add_subplot(235)
# ax4.imshow(images_exp_backf_f[1,...])
# plt.tight_layout()








# imfit_rgbs = np.moveaxis(np.array([imfit_rgb/250, imfit_rgb / 30, imfit_rgb]), 0, -1)
# # HDR exposure blending with radial filter. The point is to apply to each exposure independently, and export 3 tiffs for blending in Photoshop
# # I still visualize the 3 exposures in 3 subplots
# rr = r - radius
# rr[rr < 0] = 0
# rrN = rr / rr.max()
# rrs = np.tile(np.reshape(rr, [1, ny, nx, 1]), (3, 1, 1, 3))
# rrsN = np.tile(np.reshape(rrN, [1, ny, nx, 1]), (3, 1, 1, 3))

