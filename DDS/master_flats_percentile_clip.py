import astro
from pathlib import Path
import os
from astropy.io import fits
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches

plt.rcParams['font.size'] = '14'

if __name__ == "__main__":

    cal_dir = Path(os.environ['DATA'], 'DDS', 'Taka', 'Calibration')
    mbiasf = Path(cal_dir, 'master_bias_bin1_Mar_2018.fits')
    flats_dir = Path(cal_dir, 'Flats')
    flatsf = sorted(list(flats_dir.glob('Hflat_Mar_3_2019/HFlat*.fit')))
    nframes = len(flatsf)
    # load masters bias and dark
    mbias = fits.getdata(mbiasf)
    # Get reference flat image and calibrate with master bias
    # (for comparison, choose same as Pixinsight)
    ref = 6
    flath = fits.getheader(flatsf[ref])
    nx, ny = flath['NAXIS1'], flath['NAXIS2']
    npixels = nx * ny
    flatd = fits.getdata(flatsf[ref]) - mbias
    # Reference median for flux equalization prior to pixel rejection
    m0 = np.median(flatd)

    # Loop over the flat files and calibrate with optimized dark
    flat_cube = np.zeros([nframes, ny, nx], dtype=np.float32)
    print('Calibrating flat series...')
    for i, f in enumerate(flatsf):
        flat_header = fits.getheader(f)
        flat_cube[i] = fits.getdata(f) - mbias
        # At this point, images confirmed to be equivalent to Pixinsight processing
        # # Flux equalization factor
        flux_eq_scale = m0/np.median(flat_cube[i])
        flat_cube[i] *= flux_eq_scale

    print('Creating master flat with percentile clipping...')
    percentiles = np.arange(51)*0.001
    rej1_l = []
    rej2_l = []
    nlows = []
    nhighs = []

    for i, p in enumerate(percentiles):
        print(f'percentile: {p:.3f}')
        plow = p
        phigh = p
        rej_mask, rej1, rej2 = astro.percentile_clip(flat_cube, plow=plow, phigh=phigh)
        rej1_l.append(rej1)
        rej2_l.append(rej2)
        # flat_cube_ma = np.ma.masked_array(flat_cube, mask=rej_mask)
        # master_flat = flat_cube_ma.mean(axis=0).filled()
        # print('------master Flat (unclipped)-----')
        # astro.stats(master_flat)

        rej1s = rej1.sum(axis=0)
        rej2s = rej2.sum(axis=0)

        nlow = round(rej1s.sum() / nframes)
        nhigh = round(rej2s.sum() / nframes)
        nlows.append(nlow)
        nhighs.append(nhigh)

    low_rates = np.array(nlows) / npixels * 100
    high_rates = np.array(nhighs) / npixels * 100


    for i, p, in enumerate(percentiles[0:]):

        # i = 10
        rej1s = rej1_l[i].sum(axis=0)
        rej2s = rej2_l[i].sum(axis=0)
        roi = np.s_[100:3000, 800:4000]

        idx = 0
        vmin = 0
        vmax1 = max(rej1s.mean(), 1)
        vmax2 = max(rej2s.mean(), 1)

        plt.close('all')
        fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(22, 12))
        axs[0, 0].set_position([-0.05, 0.3, 0.65, 0.65])
        axs[0, 1].set_position([0.4, 0.3, 0.65, 0.65])

        axs[0, 0].imshow(rej1s, vmin=vmin, vmax=vmax1, cmap='gray', interpolation='kaiser')
        axs[0, 0].tick_params(labelbottom=False, labelleft=False)
        axs[0, 0].set_title(f'rejection mask ; percentile low : {p:.3f}')
        axs[0, 0].text(50, 300,
                       f'mean rejected pixels per image: {nlows[i]} ({low_rates[i]:.1f}%)'
                       f'\nmax frames rejected: {rej1s[roi].max()}'
                       f'\nmin frames rejected: {rej1s[roi].min()}',
                       color='black',
                       bbox=dict(facecolor='white', edgecolor='black', boxstyle='square', alpha=0.7))
        axs[0, 1].imshow(rej2s, vmin=vmin, vmax=vmax2, cmap='gray', interpolation='kaiser')
        axs[0, 1].tick_params(labelbottom=False, labelleft=False)
        axs[0, 1].set_title(f'rejection mask ; percentile high: {p:.3f}')
        axs[0, 1].text(50, 300,
                       f'mean rejected pixels per image: {nhighs[i]} ({high_rates[i]:.1f}%)'
                       f'\nmax frames rejected: {rej2s[roi].max()}'
                       f'\nmin frames rejected: {rej2s[roi].min()}',
                       color='black',
                       bbox=dict(facecolor='white', edgecolor='black', boxstyle='square', alpha=0.7))

        axs[1, 0].plot(low_rates)
        axs[1, 0].set_xlim([0, len(percentiles)])
        axs[1, 0].set_ylim([0, 45])
        axs[1, 0].axvline(x=i, color='black', ls='--', lw=1)
        axs[1, 0].text(i + 0.5, 20, f'{nhighs[i]} px', rotation=90)
        axs[1, 0].set_xlabel('Frame #')
        axs[1, 0].set_ylabel('mean rejection per frame [%]', fontsize=12)
        axs[1, 1].plot(high_rates)
        axs[1, 1].set_xlim([0, len(percentiles)])
        axs[1, 1].set_ylim([0, 45])
        axs[1, 1].axvline(x=i, color='black', ls='--', lw=1)
        axs[1, 1].text(i+0.5, 20, f'{nhighs[i]} px', rotation=90)
        axs[1, 1].set_xlabel('Frame #')
        axs[1, 1].set_ylabel('mean rejection per frame [%]', fontsize=12)

        axs[1, 0].set_position([0.08, 0.07, 0.3, 0.2])
        axs[1, 1].set_position([0.53, 0.07, 0.3, 0.2])

        roi = [570-50, 570+50, 1000, 1100]
        sub1 = rej1s[roi[2]:roi[3], roi[0]:roi[1]]
        sub2 = rej2s[roi[2]:roi[3], roi[0]:roi[1]]

        ax1 = fig.add_axes([0.34, 0.07, 0.2, 0.2])
        vmax1 = nframes-1
        vmax2 = nframes-1
        ax1.imshow(rej1s, vmin=vmin, vmax=vmax1, cmap='gray', interpolation='kaiser')
        ax1.axis(roi)
        ax1.tick_params(labelbottom=False, labelleft=False)
        ax1.text(roi[0]+2, roi[2]+85, f'max frames rejected: {sub1.max()}'
                                      f'\nmin frames rejected: {sub1.min()}',
                 color='white', fontsize=12)
                 # bbox=dict(facecolor='white', edgecolor='black', boxstyle='square', alpha=0.7))
        ax2 = fig.add_axes([0.79, 0.07, 0.2, 0.2])
        ax2.imshow(rej2s, vmin=vmin, vmax=vmax2, cmap='gray', interpolation='kaiser')
        ax2.axis(roi)
        ax2.tick_params(labelbottom=False, labelleft=False)
        ax2.text(roi[0]+2, roi[2]+85, f'max frames rejected: {sub2.max()}'
                                      f'\nmin frames rejected: {sub2.min()}',
                 color='white', fontsize=12)
                 # bbox=dict(facecolor='white', edgecolor='black', boxstyle='square', alpha=0.7))
        # Create a Rectangle patch
        rect1 = patches.Rectangle((roi[0], roi[2]), 100, 100, linewidth=1, edgecolor='red', facecolor='none')
        rect2 = patches.Rectangle((roi[0], roi[2]), 100, 100, linewidth=1, edgecolor='red', facecolor='none')
        # Add the patch to the Axes
        axs[0, 0].add_patch(rect1)
        axs[0, 1].add_patch(rect2)

        plt.savefig(Path(flats_dir, f'hflat_rejection_{i}.jpg'))

    # Check movie
    # decide on plow = phigh = 0.017
    rej_mask, rej1, rej2 = astro.percentile_clip(flat_cube, plow=0.017, phigh=0.017)
    flat_cube_ma = np.ma.masked_array(flat_cube, mask=rej_mask)
    master_flat = flat_cube_ma.mean(axis=0).filled()
    fits.writeto(Path(flats_dir.parent, f'master_hflat_m25_bin1_Mar_2019.fits'), master_flat.astype(np.float32), overwrite=True)
    master_flat_n = master_flat / 65535
    fits.writeto(Path(flats_dir.parent, f'master_hflat_m25_bin1_Mar_2019_norm.fits'), master_flat_n.astype(np.float32), overwrite=True)

