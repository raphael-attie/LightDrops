import numpy as np
from astropy.io import fits
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import os

def fitsread(f, header=False):
    with fits.open(f) as hdul:
        data = hdul[0].data
        if header:
            h = hdul[0].header
            return data, h
        else:
            return data


def loop_wins(images, med, sigma, rejMask, verbose=False):
    n1 = 1
    qmask = np.ones(images.shape[0], dtype=np.bool)
    images[rejMask] = np.nan
    while n1 > 0:
        m0 = (med - 1.5 * sigma)[:, np.newaxis]
        m1 = (med + 1.5 * sigma)[:, np.newaxis]
        mask_min = images < m0
        mask_max = images > m1
        print(mask_min.sum())
        # mask_min[~qmask, :] = False
        # mask_max[~qmask, :] = False
        images[mask_min] = np.tile(m0, (1, images.shape[1]))[mask_min]
        images[mask_max] = np.tile(m1, (1, images.shape[1]))[mask_max]
        sigma0 = sigma.copy()
        # Apply pixel rejection mask before calculating new sigma
        med[qmask] = np.median(images, axis=1)[qmask]
        sigma[qmask] = 1.134 * np.std(images, axis=1)[qmask]

        qmask = np.abs(sigma - sigma0) / sigma0 > 0.0005

        n1 = qmask.sum()
        if n1 > 0:
            print(med.sum())
        if verbose: print('n1 = ', n1)
    return med, sigma


def loop_wins2(images, med, sigma, rejMask, verbose=False):
    rows0 = np.arange(images.shape[0])
    med_out = med.copy()
    sigma_out = sigma.copy()
    images[rejMask] = np.nan
    # 1st pass of winsorization using just the median as the replacement value
    mask_init = np.abs(images - med[:, np.newaxis])/sigma_out[:, np.newaxis] > 5
    images[mask_init] = np.tile(med[:, np.newaxis], (1, images.shape[1]))[mask_init]
    sigma0 = sigma.copy()
    # Apply pixel rejection mask before calculating new sigma
    if np.isnan(images.sum()):
        med = np.nanmedian(images, axis=1)
        sigma = 1.134 * np.nanstd(images, axis=1)
    else:
        med = np.median(images, axis=1)
        sigma = 1.134 * np.std(images, axis=1)
    med_out[rows0] = med
    sigma_out[rows0] = sigma

    rows = np.where(np.abs(sigma - sigma0) / sigma0 > 0.0005)[0]

    images = images[rows, :]
    med = med[rows]
    sigma = sigma[rows]

    rows0 = rows0[rows]
    n1 = len(rows)

    while n1 > 0:

        m0 = (med - 1.5 * sigma)[:, np.newaxis]
        m1 = (med + 1.5 * sigma)[:, np.newaxis]
        mask_min = images < m0
        mask_max = images > m1
        images[mask_min] = np.tile(m0, (1, images.shape[1]))[mask_min]
        images[mask_max] = np.tile(m1, (1, images.shape[1]))[mask_max]
        sigma0 = sigma.copy()
        # Apply pixel rejection mask before calculating new sigma
        if np.isnan(images.sum()):
            med = np.nanmedian(images, axis=1)
            sigma = 1.134 * np.nanstd(images, axis=1)
        else:
            med = np.median(images, axis=1)
            sigma = 1.134 * np.std(images, axis=1)
        med_out[rows0] = med
        sigma_out[rows0] = sigma

        rows = np.where(np.abs(sigma - sigma0) / sigma0 > 0.0005)[0]

        images = images[rows, :]
        med = med[rows]
        sigma = sigma[rows]

        rows0 = rows0[rows]

        n1 = len(rows)
        if n1 > 0:
            print(med_out.sum())
        if verbose: print('n1 = ', n1)
    med_out[rows0] = med
    sigma_out[rows0] = sigma

    return med_out, sigma_out


cal_dir = Path(os.environ['DATA'], 'DDS', 'Taka', 'Calibration')

bias_dir = Path(cal_dir, 'Bias_bin1_Mar_2018')
biasf = list(bias_dir.rglob('Bias*.fit'))
print(len(biasf))
d0, h0 = fitsread(biasf[0], header=True)
nx = h0['NAXIS1']
ny = h0['NAXIS2']
ny2 = int(ny/4)
print(nx, ny)
print(nx, ny2)
slice1 = np.s_[0:ny2, :]

datacube = np.moveaxis(np.array([fits.getdata(f) for f in biasf]), 0, 2)

flatc = datacube[0,...].astype(np.float32)
m = np.median(flatc, axis=1)
sigma = np.std(flatc, axis=1)
rej_mask = np.zeros(flatc.shape, dtype=np.bool)
m1, sigma1 = loop_wins(flatc.copy(), m, sigma, rej_mask, verbose=True)

m = np.median(flatc, axis=1)
sigma = np.std(flatc, axis=1)
rej_mask = np.zeros(flatc.shape, dtype=np.bool)
m2, sigma2 = loop_wins2(flatc.copy(), m, sigma, rej_mask, verbose=True)

print(np.array_equal(sigma1,sigma2))
