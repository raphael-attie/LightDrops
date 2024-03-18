from astropy.io import fits
from scipy.signal import convolve2d
from astropy.stats import median_absolute_deviation as mad
from pathlib import Path
import numpy as np
import os
from multiprocessing import Pool
from functools import partial
from multiprocessing import shared_memory
import time


def get_divisors(n):
    for i in range(1, int(n / 2) + 1):
        if n % i == 0:
            yield i
    yield n


def sigma_clip(images, rej_mask, slow=5, shigh=5):
    ma_images = np.ma.masked_array(images, mask=rej_mask.copy())
    median = np.ma.median(ma_images, axis=-1)
    sigma = ma_images.std(axis=-1)
    rej_mask1 = median - ma_images > slow * sigma
    rej_mask2 = ma_images - median > shigh * sigma
    rej_mask = rej_mask1 | rej_mask2
    return rej_mask


def percentile_clip(images, plow=0.01, phigh=0.01):
    print(f'plow = {plow:.3f} ; phigh = {phigh:.3f}')
    median = np.median(images, axis=0)
    rej_mask1 = median - images > plow * median
    rej_mask2 = images - median > phigh * median
    rej_mask = rej_mask1 | rej_mask2
    # for i, mask in enumerate(rej_mask):
    #     print(f'rejected pixels (i) : {mask.sum()}')
    return rej_mask, rej_mask1, rej_mask2


def rejection(images, rej_function, **kwargs):
    rej_mask = np.zeros(images.shape, dtype=np.bool)
    n = 1
    while n > 0:
        rej_mask0 = rej_mask.copy()
        rej_mask = rej_function(images, rej_mask, **kwargs)
        n = rej_mask.sum()
        print(f'Nb of rejected pixels : {n}')
        rej_mask = rej_mask0 | rej_mask.filled(False)

    return rej_mask


def loop_wins(images, med, sigma, rejMask, median_pass=False):

    rows0 = np.arange(images.shape[0])
    med_out = med.copy()
    sigma_out = sigma.copy()
    images[rejMask] = np.nan
    n1 = 1

    while n1 > 0:

        if median_pass:
            mask_med = np.abs(images - med[:, np.newaxis]) / sigma[:, np.newaxis] > 5
            images[mask_med] = np.tile(med[:, np.newaxis], (1, images.shape[1]))[mask_med]
        else:
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


def win_sigma_clip(flatc, slow, shigh, median_pass=False, verbose=False):
    # Make a 1st pass while there's no need to consider NaN-flagged arrays
    # as median() and std() are faster than nanmedian() and nanstd().

    if verbose: print('flattened shape = ', flatc.shape)
    rej_mask = np.zeros(flatc.shape, dtype=np.bool)

    m = np.median(flatc, axis=1)
    sigma = np.std(flatc, axis=1)
    # Winsorized sigma
    m, sigma = loop_wins(flatc.copy(), m, sigma, rej_mask, median_pass=median_pass)
    mask = np.abs(flatc - m[:, np.newaxis]) > 5*sigma[:, np.newaxis]
    n = mask.sum()
    if n == 0:
        return rej_mask
    # Prepare new passes only on the rows that have outliers.
    clip_rows = np.where(np.any(mask, axis=1))[0]
    clip_rows0 = clip_rows.copy()
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
        m, sigma = loop_wins(flatc.copy(), m, sigma, mask)

        # mask = np.abs(flatc - m[:, np.newaxis]) > 5*sigma[:, np.newaxis]
        mask1 = m[:, np.newaxis] - flatc > slow * sigma[:, np.newaxis]
        mask2 = flatc - m[:, np.newaxis] > shigh * sigma[:, np.newaxis]
        mask = mask1 | mask2
        clip_rows = np.where(np.any(mask, axis=1))[0]
        mask = mask[clip_rows, :]
        n = mask.sum()
        if verbose:
            print('total outliers n =', n)
        clip_rows0 = clip_rows0[clip_rows]
        rej_mask[clip_rows0] = rej_mask[clip_rows0] | mask

    return rej_mask


def wins_wrapper(rows, shape, shared_name, slow, shigh, median_pass=False, verbose=False):
    existing_shm = shared_memory.SharedMemory(name=shared_name)
    shared_slice = np.ndarray(shape, dtype=np.float32, buffer=existing_shm.buf)
    data_slice = shared_slice[rows[0]:rows[1], ...]
    sz = data_slice.shape
    rejmask = np.array([win_sigma_clip(data_slice[r, ...], shigh, slow, median_pass=median_pass, verbose=verbose) for r in range(sz[0])]).reshape(sz)
    return rejmask


def wins_sigma_par(datacube, slow, shigh, median_pass=False):

    shm = shared_memory.SharedMemory(create=True, size=datacube.nbytes)
    sh_cube_slice = np.ndarray(datacube.shape, dtype=datacube.dtype, buffer=shm.buf)
    sh_cube_slice[:] = datacube[:]

    nrows = datacube.shape[0]
    nchunks = 151
    nr = nrows / nchunks
    r0 = np.arange(nchunks) * nr
    r1 = r0 + nr
    rows = [[int(a), int(b)] for a, b in zip(r0, r1)]

    wins_partial = partial(wins_wrapper, slow=slow, shigh=shigh, shape=datacube.shape, shared_name=shm.name, median_pass=median_pass, verbose=False)
    # rejmask_list = list(map(wins_partial, rows))
    # print('starting parallel pool')
    with Pool(processes=10) as p:
        rejmask_list = p.map(wins_partial, rows)

    rejmask = np.array(rejmask_list).reshape(datacube.shape)

    return rejmask


def make_master(data_dir, input_stem, median_pass=False):

    dataf = list(data_dir.rglob(f'{input_stem}.fit'))
    print(len(dataf))
    d0, h0 = fits.getdata(dataf[0], header=True)
    nx = h0['NAXIS1']
    ny = h0['NAXIS2']
    nframes = len(dataf)
    print(f'Processing {nframes} frames...')

    ny2 = int(ny / 2)  # = 906
    # list(get_divisors(906)) gives [1, 2, 3, 6, 151, 302, 453, 906]
    rejMasks = []
    start_time = time.time()
    datacube = np.moveaxis(np.array([fits.getdata(f) for f in dataf]), 0, 2)
    for i in range(2):
        yslice = np.s_[ny2 * i:ny2 * (i + 1), :, :]
        print('loading slice data')
        rej_mask = wins_sigma_par(datacube[yslice].astype(np.float32), median_pass=median_pass)
        rejMasks.append(rej_mask)
        print(f"--- {time.time() - start_time} seconds ---")
    rejMask = np.array(rejMasks).reshape([ny, nx, nframes])
    print(f"--- {time.time() - start_time} seconds ---")

    masked_datacube = np.ma.masked_array(datacube, mask=rejMask)
    master_data = masked_datacube.mean(axis=2).filled()

    return master_data


def make_master_flat(datacube, rej_function, **kwargs):

    ny, nx, nframes = datacube.shape
    print(f'Processing {nframes} frames...')
    start_time = time.time()
    rej_mask = rejection(datacube, rej_function, **kwargs)
    print(f"--- {time.time() - start_time} seconds ---")

    masked_datacube = np.ma.masked_array(datacube, mask=rej_mask)
    master_data = masked_datacube.mean(axis=-1).filled()

    return master_data, rej_mask


def write_master_bias(data_dir, input_stem, out_file, median_pass=False):

    master_bias = make_master(data_dir, input_stem, median_pass=median_pass)
    fits.writeto(Path(data_dir.parent, f'{out_file}.fits'), master_bias.astype(np.float32), overwrite=True)
    # For compatibility with Pixinsight, normalize
    master_bias_n = master_bias / 65535
    fits.writeto(Path(data_dir.parent, f'{out_file}_norm.fits'), master_bias_n.astype(np.float32), overwrite=True)

    print('------master bias (clipped)-----')
    stats(master_bias, clip=True)
    print('------master bias (unclipped)-----')
    stats(master_bias)
    return master_bias


def write_master_dark(data_dir, input_stem, out_file, master_bias=None, median_pass=False):

    master_dark = make_master(data_dir, input_stem, median_pass=median_pass)
    if master_bias is not None:
        master_dark = master_dark - master_bias

    fits.writeto(Path(data_dir.parent, f'{out_file}.fits'), master_dark.astype(np.float32), overwrite=True)
    # For compatibility with Pixinsight, normalize
    master_dark_n = master_dark / 65535
    fits.writeto(Path(data_dir.parent, f'{out_file}_norm.fits'), master_dark_n.astype(np.float32), overwrite=True)

    print('------master dark (clipped)-----')
    stats(master_dark, clip=True)
    print('------master dark (unclipped)-----')
    stats(master_dark)
    return master_dark


def stats(data0, clip=False):
    data = data0.copy()
    data = data.ravel()
    if clip:
        data = data[data < 65535]
        data = data[data > 0]
    count_px = len(data)
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


def evaluate_noise(image):
    kernel_b3_spline = np.array([
        [1 / 256, 1 / 64, 3 / 128, 1 / 64, 1 / 256],
        [1 / 64, 1 / 16, 3 / 32, 1 / 16, 1 / 64],
        [3 / 128, 3 / 32, 9 / 64, 3 / 32, 3 / 128],
        [1 / 64, 1 / 16, 3 / 32, 1 / 16, 1 / 64],
        [1 / 256, 1 / 64, 3 / 128, 1 / 64, 1 / 256]
    ])

    c0 = image
    c1 = convolve2d(c0, kernel_b3_spline, mode='same', boundary='symm')
    w1 = c0 - c1
    # Work with the 1D-unwrapped array
    w1f = w1.ravel()
    sigma = w1f.std()
    rejmask = np.abs(w1f) < 3 * sigma
    w1f = w1f[rejmask]
    sigmac = sigma
    sigma = w1f.std()
    ds = abs(sigmac - sigma) / sigmac * 100
    while ds >= 1:
        rejmask = np.abs(w1f) < 3 * sigma
        w1f = w1f[rejmask]
        sigmac = sigma
        sigma = w1f.std()
        ds = abs(sigmac - sigma) / sigmac * 100
    return sigma/0.889


def same_sign(a, b):
    return (a if a < 0 else -a) if (b < 0) else (-a if a < 0 else a)


def bracket_dark_opt(image, dark):
    TEST_DARK = partial(evaluate_noise)
    GOLD = 1.618034
    GLIMIT = 10
    TINY = 1.0e-20
    ax = 0.5
    bx = 2.0
    fa = TEST_DARK(image - ax * dark)
    fb = TEST_DARK(image - bx * dark)
    if (fb > fa):
        #         print('swapping')
        ax, bx = bx, ax
        fa, fb = fb, fa
    cx = bx + GOLD * (bx - ax)
    fc = TEST_DARK(image - cx * dark)

    while fb > fc:
        r = (bx - ax) * (fb - fc)
        q = (bx - cx) * (fb - fa)
        u = bx - ((bx - cx) * q - (bx - ax) * r) / 2 / same_sign(max(abs(q - r), TINY), q - r)
        ulim = bx + GLIMIT * (cx - bx)
        if (bx - u) * (u - cx) > 0:
            fu = TEST_DARK(image - u * dark)
            if fu < fc:
                ax = bx
                bx = u
                break
            elif fu > fb:
                cx = u
                break
            u = cx + GOLD * (cx - bx)
            fu = TEST_DARK(image - u * dark)
        elif (cx - u) * (u - ulim) > 0:
            fu = TEST_DARK(image - u * dark)
            if fu < fc:
                du = GOLD * (u - cx)
                bx = cx
                cx = u
                u += du
                fb = fc
                fc = fu
                fu = TEST_DARK(image - u * dark)
        elif (u - ulim) * (ulim - cx) >= 0:
            u = ulim
            fu = TEST_DARK(image - u * dark)
        else:
            u = cx + GOLD * (cx - bx)
            fu = TEST_DARK(image - u * dark)

        ax = bx
        bx = cx
        cx = u
        fa = fb
        fb = fc
        fc = fu

    return ax, bx, cx


def dark_optimize(image, dark):
    TEST_DARK = evaluate_noise

    R = 0.61803399
    C = 1 - R

    ax, bx, cx = bracket_dark_opt(image, dark)
    # Define [x0,x3] as the total search interval.
    x0 = ax
    x3 = cx
    # Set and use [x1,x2] as the inner search interval
    if abs(bx - ax) < abs(cx - bx):
        x1 = bx
        x2 = bx + C * (cx - bx)
    else:
        x1 = bx - C * (bx - ax)
        x2 = bx

    f1 = TEST_DARK(image - x1 * dark)
    f2 = TEST_DARK(image - x2 * dark)

    while abs(x3 - x0) > 0.0005:

        if f2 < f1:
            x0 = x1
            x1 = x2
            x2 = R * x2 + C * x3
            f1 = f2
            f2 = TEST_DARK(image - x2 * dark)
        else:
            x3 = x2
            x2 = x1
            x1 = R * x1 + C * x0
            f2 = f1
            f1 = TEST_DARK(image - x1 * dark)

    return max(0, (x1 if f1 < f2 else x2))