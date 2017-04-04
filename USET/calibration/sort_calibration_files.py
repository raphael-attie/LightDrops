"""
Move the calibration files to directories according to their nature: bias, dark or light
"""

import os
import glob
import shutil
from astropy.io import fits

# Set the directory where the FITS are and get the list of files
data_dir    =  '/Users/rattie/Data/USET/campaign/calibration/HALPHA'
file_list   = glob.glob(os.path.join(data_dir, '*.FTS'))
# List of tested gain values
gains       = [1000, 2000]

gain_dirs   = [os.path.join(data_dir, 'Gain'+str(gain)) for gain in gains]
lights_dirs = [os.path.join(gdir, 'lights') for gdir in gain_dirs]
darks_dirs  = [os.path.join(gdir, 'darks') for gdir in gain_dirs]
bias_dirs   = [os.path.join(gdir, 'bias') for gdir in gain_dirs]

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
    h   = hdu[0].header
    img = hdu[0].data
    # Get gain index of current image
    g           = gains.index(h['GAIN'])
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




print('done')