from astropy.io import fits
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import os
from multiprocessing import Pool
import multiprocessing as mp
from functools import partial
from multiprocessing import shared_memory
import time


def get_divisors(n):
    for i in range(1, int(n / 2) + 1):
        if n % i == 0:
            yield i
    yield n


def loop_wins(images, med, sigma, rejMask, verbose=False):

    rows0 = np.arange(images.shape[0])
    med_out = med.copy()
    sigma_out = sigma.copy()
    images[rejMask] = np.nan
    # # 1st pass of winsorization using just the median as the replacement value
    # mask_init = np.abs(images - med[:, np.newaxis])/sigma[:, np.newaxis] > 5
    # images[mask_init] = np.tile(med[:, np.newaxis], (1, images.shape[1]))[mask_init]
    # sigma0 = sigma.copy()
    # # Apply pixel rejection mask before calculating new sigma
    # if np.isnan(images.sum()):
    #     med = np.nanmedian(images, axis=1)
    #     sigma = 1.134 * np.nanstd(images, axis=1)
    # else:
    #     med = np.median(images, axis=1)
    #     sigma = 1.134 * np.std(images, axis=1)
    # med_out[rows0] = med
    # sigma_out[rows0] = sigma
    #
    # rows = np.where(np.abs(sigma - sigma0) / sigma0 > 0.0005)[0]
    #
    # images = images[rows, :]
    # med = med[rows]
    # sigma = sigma[rows]
    #
    # rows0 = rows0[rows]
    # n1 = len(rows)
    n1 = 1

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

    med_out[rows0] = med
    sigma_out[rows0] = sigma

    return med_out, sigma_out


def win_sigma_clip2(flatc, verbose=False):
    # Make a 1st pass while there's no need to consider NaN-flagged arrays
    # as median() and std() are faster than nanmedian() and nanstd().

    if verbose: print('flattened shape = ', flatc.shape)
    rej_mask = np.zeros(flatc.shape, dtype=np.bool)

    m = np.median(flatc, axis=1)
    sigma = np.std(flatc, axis=1)
    # Winsorized sigma
    m, sigma = loop_wins(flatc.copy(), m, sigma, rej_mask, verbose=verbose)

    mask = np.abs(flatc - m[:,np.newaxis]) > 5*sigma[:,np.newaxis]
    n = mask.sum()
    if n == 0:
        return rej_mask
    # Prepare new passes only on the rows that have outliers.
    clip_rows = np.where(np.any(mask, axis=1))[0]
    mask = mask[clip_rows, :]
    rej_mask[clip_rows, :] = mask
    if verbose:
        print('Total outliers n =', n)
        print('new shape = ', mask.shape)
    while n > 0:
        # Work only on the rows that have outliers. Flag them as NaN
        flatc = flatc[clip_rows, :]
        flatc[mask] = np.nan
        m = np.nanmedian(flatc, axis=1)
        sigma = np.nanstd(flatc, axis=1)
        # Winsorized sigma
        m, sigma = loop_wins(flatc.copy(), m, sigma, mask, verbose=verbose)

        mask = np.abs(flatc - m[:,np.newaxis]) > 5*sigma[:,np.newaxis]
        clip_rows0 = clip_rows.copy()
        clip_rows = np.where(np.any(mask, axis=1))[0]
        mask = mask[clip_rows, :]
        n = mask.sum()
        if verbose: print('total outliers n =', n)
        rej_mask[clip_rows0[clip_rows]] = rej_mask[clip_rows0[clip_rows]] | mask

    return rej_mask


def wins_wrapper(rows, shape, shared_name, verbose=False):
    existing_shm = shared_memory.SharedMemory(name=shared_name)
    shared_slice = np.ndarray(shape, dtype=np.float32, buffer=existing_shm.buf)
    data_slice = shared_slice[rows[0]:rows[1], ...]
    sz = data_slice.shape
    rejmask = np.array([win_sigma_clip2(data_slice[r, ...], verbose=verbose) for r in range(sz[0])]).reshape(sz)
    return rejmask


def wins_sigma_par(datacube):

    shm = shared_memory.SharedMemory(create=True, size=datacube.nbytes)
    sh_cube_slice = np.ndarray(datacube.shape, dtype=datacube.dtype, buffer=shm.buf)
    sh_cube_slice[:] = datacube[:]

    nrows = datacube.shape[0]
    nchunks = 151
    nr = nrows / nchunks
    r0 = np.arange(nchunks) * nr
    r1 = r0 + nr
    rows = [[int(a), int(b)] for a, b in zip(r0, r1)]

    wins_partial = partial(wins_wrapper, shape=datacube.shape, shared_name=shm.name, verbose=False)
    # rejmask_list = list(map(wins_partial, rows))
    # print('starting parallel pool')
    with Pool(processes=10) as p:
        rejmask_list = p.map(wins_partial, rows)

    rejmask = np.array(rejmask_list).reshape(datacube.shape)

    return rejmask


if __name__ == "__main__":

    cal_dir = Path(os.environ['DATA'], 'DDS', 'Taka', 'Calibration')
    bias_dir = Path(cal_dir, 'Bias_bin1_Mar_2018')
    biasf = list(bias_dir.rglob('Bias*.fit'))
    print(len(biasf))
    d0, h0 = fits.getdata(biasf[0], header=True)
    nx = h0['NAXIS1']
    ny = h0['NAXIS2']
    nframes = len(biasf)

    ny2 = int(ny/2) # = 906
    # list(get_divisors(906)) gives [1, 2, 3, 6, 151, 302, 453, 906]
    rejMasks = []
    start_time = time.time()
    datacube = np.moveaxis(np.array([fits.getdata(f) for f in biasf]), 0, 2)
    for i in range(2):
        yslice = np.s_[ny2*i:ny2*(i+1), :, :]
        print('loading slice data')
        # cube_slice = np.moveaxis(np.array([fits.getdata(f)[yslice] for f in biasf]), 0, 2)
        rej_mask = wins_sigma_par(datacube[yslice].astype(np.float32))
        rejMasks.append(rej_mask)
        print(f"--- {time.time() - start_time} seconds ---")
    rejMask = np.array(rejMasks).reshape([ny, nx, nframes])
    print(f"--- {time.time() - start_time} seconds ---")
    # 6 cores - nchuncks=6- 416.48 seconds for the entire image series
    # 10 cores - nchuncks=151 - 383 seconds for the entire image series
    # 6 cores - nchuncks=151 - xxx seconds for the entire image series
    # 10 codes - ny/2 - nchuncks=151 - --- 137.73691129684448 seconds ---
    # With new 1st pass: --- 216.1965844631195 seconds ---

    maskedBias = np.ma.masked_array(datacube, mask=rejMask)
    masterBias = maskedBias.mean(axis=2).filled()
    # masterBias_float = masterBias / 65535
    # fits.writeto(Path(cal_dir, 'mBias_bin1_wsc_Mars_2018.fits'), masterBias, overwrite=True)
    # fits.writeto(Path(cal_dir, 'mBias_bin1_wsc_Mars_2018_normalized_float.fits'), masterBias_float.astype(np.float32), overwrite=True)

    from astropy.stats import median_absolute_deviation as mad

    def stats(data):
        count_px = len(data.ravel())
        mean = data.mean()
        median = np.median(data)
        avgDev = np.abs(data - median).mean()
        med_abs_dev = mad(data)
        print(f'count_px = {count_px}')
        print(f'mean = {mean:0.3f}')
        print(f'median = {median:0.3f}')
        print(f'avgDev = {avgDev:0.3f}')
        print(f'MAD = {med_abs_dev:0.3f}')
        print(f'minimum = {data.min():0.3f}')
        print(f'maximum = {data.max():0.3f}')
        return None

    stats(masterBias)


