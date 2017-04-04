import os
import glob
from astropy.io import fits
import numpy as np
from scipy.optimize import curve_fit
from scipy.misc import factorial
from collections import namedtuple
import matplotlib.pyplot as plt



def load_cube(file_list):
    hdu = fits.open(file_list[0], ignore_missing_end=True)
    naxis1 = hdu[0].header['NAXIS1']
    naxis2 = hdu[0].header['NAXIS2']

    cube = np.zeros([naxis2, naxis1, len(file_list)], dtype=np.int16)
    for i in range(0, len(file_list)):
        file = file_list[i]
        hdu = fits.open(file, ignore_missing_end=True)
        # Load header and image data from the 1st data unit: hdu[0]
        h   = hdu[0].header
        cube[:,:,i] = hdu[0].data

    return cube, h

def make_master_frame(file_list):

    cube, h = load_cube(file_list)

    # Average
    master = np.median(cube, axis=2)
    mean = cube.mean()
    sigma = cube.std()
    return master, h, cube, mean, sigma


# poisson function, parameter lamb is the fit parameter
def poisson(k, lamb):
    return (lamb**k/factorial(k)) * np.exp(-lamb)


def fit_poisson(data):
    # the bins should be of integer width, because poisson is an integer distribution
    bin_edges = np.arange(0, 4095, 5)
    hist, bin_edges = np.histogram(data.ravel(), bins=bin_edges, density=True)
    # calculate bin middles
    bin_middles = 0.5*(bin_edges[1:] + bin_edges[:-1])
    # fit with curve_fit
    params, cov_matrix = curve_fit(poisson, bin_middles, hist)
    return params, hist, bin_edges


def master_statistics(directory, gains, master_name):
    # Study bias images
    ng = len(gains)
    statdict = {'masters':[None] * ng, 'headers': [None] * ng, 'cubes': [None] * ng, 'cube_mean':[None] * ng, 'cube_sigma':[None] * ng, 'params': [None] * ng,
                'hists': [None] * ng, 'bins': [None] * ng, 'xpoisson': [None] * ng, 'ypoisson': [None] * ng}

    for g in gains:
        i = gains.index(g)
        fits_files  = glob.glob(os.path.join(directory[i], '*.FTS'))
        statdict['masters'][i], statdict['headers'][i], statdict['cubes'][i], statdict['cube_mean'][i], statdict['cube_sigma'][i] = make_master_frame(fits_files)
        statdict['params'][i], statdict['hists'][i], statdict['bins'][i] = fit_poisson(np.round(statdict['cubes'][i]))
        # plot poisson-deviation with fitted parameter
        statdict['xpoisson'][i]  = np.arange(0, 1000)
        statdict['ypoisson'][i]   = poisson(statdict['xpoisson'][i] , *statdict['params'][i])
        # Write master to fits file
        rimage = np.rint(statdict['masters'][i])
        int_image = np.int16(rimage)
        fname = os.path.join(directory[gains.index(g)], 'master_%s_G%i.fits'%(master_name, g))
        fits.writeto(fname, int_image, header=statdict['headers'][i], output_verify='exception', overwrite=True)

    return statdict

# Set the directory where the FITS are and get the list of files
data_dir    =  '/Users/rattie/Data/USET/campaign/calibration/HALPHA'
# List of tested gain values
gains       = [1000, 2000]
gain_dirs   = [os.path.join(data_dir, 'Gain'+str(gain)) for gain in gains]
lights_dirs = [os.path.join(gdir, 'lights') for gdir in gain_dirs]
darks_dirs  = [os.path.join(gdir, 'darks') for gdir in gain_dirs]
bias_dirs   = [os.path.join(gdir, 'bias') for gdir in gain_dirs]

# Bias
bias_stat = master_statistics(bias_dirs, gains, 'bias')
# Darks
dark_stat = master_statistics(darks_dirs, gains, 'dark')

def plot_bar(stat):

    nbins = 500
    gain = 1000

    x = stat['bins'][gains.index(gain)][:-1]
    width = np.diff(stat['bins'][gains.index(gain)])
    npixels1000 = (2048 ** 2) * stat['cubes'][gains.index(gain)].shape[2]*width

    bar_label1000 = r"Gain=%i: Mean = %0.1f; $\sigma= %0.1f$" % (gain, stat['cube_mean'][gains.index(gain)], stat['cube_sigma'][gains.index(gain)])
    bars1000 = plt.bar(x, stat['hists'][gains.index(gain)]*npixels1000, width=width, color='blue')

    #bars1000 = plt.hist(stat['cubes'][gains.index(gain)].ravel(), nbins, color='blue')

    gain = 2000
    npixels2000 = (2048 ** 2) * stat['cubes'][gains.index(gain)].shape[2]*width
    x = stat['bins'][gains.index(gain)][:-1]
    bar_label2000 = r"Gain=%i: Mean = %0.1f; $\sigma= %0.1f$" % (gain, stat['cube_mean'][gains.index(gain)], stat['cube_sigma'][gains.index(gain)])
    bars2000 = plt.bar(x, stat['hists'][gains.index(gain)]*npixels2000, width=width, color='green', alpha=0.6)
    #bars2000 = plt.hist(stat['cubes'][gains.index(gain)].ravel(), nbins, color='green', alpha = 0.6)


    #plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Image intensity')
    plt.ylabel('Histogram')
    plt.grid('on')
    #plt.legend(loc='lower right')
    plt.xlim([0, 800])

    #plt.axis([0, nbins/2, 0, 1E7])
    #plt.xlim([0, nbins/2])
    # #plt.xticks(np.arange(0, nbins + 1))
    ax1 = plt.gca()
    ax2 = ax1.twinx()
    # cdf1000 = ax2.plot(x, np.cumsum(stat['hists'][gains.index(1000)]), 'k-', label='Gain=1000: cdf')
    # cdf2000 = ax2.plot(x, np.cumsum(stat['hists'][gains.index(2000)]), 'k--', label='Gain=2000: cdf')

    dx = stat['bins'][gains.index(1000)][1:] - stat['bins'][gains.index(1000)][:-1]
    cdf1000, = ax2.plot(x, np.cumsum(stat['hists'][gains.index(1000)]) * dx, 'k-')
    cdf2000, = ax2.plot(x, np.cumsum(stat['hists'][gains.index(2000)]) * dx, 'k--')
    #ax2.set_xlim([0, 1E3])
    ax2.set_ylim([0, 1.1])
    #ax2.set_xscale('log')
    ax2.set_ylabel('cdf')
    #ax2.set_ylim([0.9, 1.1])
    ax1.set_xlim([0, 100])
    plt.xlim([0, 100])


    #plt.legend(bbox_to_anchor=(0.4, 0.85), loc=3, frameon=False)
    #plt.legend(bbox_to_anchor=(0.4, 0.85), loc=3, frameon=False)

    plt.legend([bars1000, bars2000, cdf1000, cdf2000],
               [bar_label1000, bar_label2000, 'Gain=1000: cdf', 'Gain=2000: cdf'])



figdir = '/Users/rattie/Data/USET/campaign/calibration/HALPHA'

plt.figure(0, figsize=(15, 10))
plt.ion()

fname = os.path.join(figdir, 'master_frames.pdf')

x1, x2 = 1000, 1500
y1, y2 = 1000, 1500

# Plot master frames
plt.subplot(221)
plt.title('Master Bias; Gain 1000; Exp. time=%0.2f ms' %(bias_stat['headers'][0]['EXPTIME']*1000))
plt.imshow(bias_stat['masters'][gains.index(1000)][y1:y2, x1:x2], cmap='gray', vmin=0, vmax=10, extent=[x1, x2, y1, y2])
plt.colorbar()
plt.subplot(222)
plt.title('Master Bias; Gain 2000; Exp. time=%0.2f ms' %(bias_stat['headers'][0]['EXPTIME']*1000))
plt.imshow(bias_stat['masters'][gains.index(2000)][y1:y2, x1:x2], cmap='gray', vmin=0, vmax=10, extent=[x1, x2, y1, y2])
plt.colorbar()

plt.subplot(223)
plt.title('Master Dark; Gain 1000; Exp. time=%0.2f ms' %(dark_stat['headers'][0]['EXPTIME']*1000))
plt.imshow(dark_stat['masters'][gains.index(1000)][y1:y2, x1:x2], cmap='gray', vmin=0, vmax=10, extent=[x1, x2, y1, y2])
plt.colorbar()

plt.subplot(224)
plt.title('Master Dark; Gain 2000; Exp. time=%0.2f ms' %(dark_stat['headers'][0]['EXPTIME']*1000))
plt.imshow(dark_stat['masters'][gains.index(2000)][y1:y2, x1:x2], cmap='gray', vmin=0, vmax=10, extent=[x1, x2, y1, y2])
plt.colorbar()
plt.tight_layout()

plt.savefig(fname, dpi=300)


# Plot statistics of individual frames


plt.figure(1, figsize=(14, 8))
fname = os.path.join(figdir, 'statistics2.pdf')

plt.subplot(121)
plt.title('Bias; Gain 1000-2000; Exp. time=%0.2f ms' %(bias_stat['headers'][0]['EXPTIME']*1000))
plot_bar(bias_stat)
plt.xlim([0, 100])

plt.subplot(122)
plt.title('Dark; Gain 1000-2000; Exp. time=%0.2f ms' %(dark_stat['headers'][0]['EXPTIME']*1000))
plot_bar(dark_stat)

plt.tight_layout()

plt.savefig(fname, dpi=300)




print('done')


