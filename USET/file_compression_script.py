import os
import glob
import warnings
from joblib import Parallel, delayed
import multiprocessing
from astropy.io import fits
import calibration.uset_calibration as uset

# Set the directory where the FITS are.
# with a more explicit path if it is not possible to set an environment variable.
data_dir1    =  '/Volumes/SDobo-A/Raphael/USET/campaign/20170327/HALPHA'
data_dir2    =  '/Volumes/SDobo-A/Raphael/USET/campaign/20170327/HALPHA_compressed'

# data_dir1    =  '/Users/rattie/Data/USET/campaign/temp1'
# data_dir2    =  '/Users/rattie/Data/USET/campaign/temp2'


# Get the list of files, change it according to where your data files are and how are they are named.
file_list   = glob.glob(os.path.join(data_dir1, '*.FTS'))

def compress_file(file1):
    try:
        hdu = fits.open(file1, ignore_missing_end=True)
        image = hdu[0].data
    except TypeError:
        hdu.close()
        print('Warning at file (ignored): %s ' %file1)
    else:
        # Load header and image data from the 1st data unit: hdu[0]
        header = hdu[0].header
        # New file name for the compressed file
        file2 = os.path.join(data_dir2, uset.get_basename(file1)+'.fits')
        uset.write_uset_fits(image, header, file2, compressed=True)
        hdu.close()

#num_cores = multiprocessing.cpu_count()
num_cores = 4

Parallel(n_jobs = num_cores)(delayed(compress_file)(file) for file in file_list)


# ignored_files_list = open('/Volumes/SDobo-A/Raphael/USET/campaign/20170327/corrupted_files.txt', 'w')
# for item in ignored_files:
#     ignored_files_list.write("%s\n" % item)
# ignored_files_list.close()
