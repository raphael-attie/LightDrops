import astro
from pathlib import Path
import os
import time


if __name__ == "__main__":

    bias_dir = Path(os.environ['DATA'], 'DDS', 'Taka', 'Calibration', 'Bias_bin1_Mar_2018')
    mbias = astro.write_master_bias(bias_dir, 'Bias*', 'master_bias_bin1_Mar_2018', median_pass=True)

    darks_dir = Path(os.environ['DATA'], 'DDS', 'Taka', 'Calibration', 'Darks_Jan_2019')
    mdark = astro.write_master_dark(darks_dir, '*1x1_600s*', 'master_dark_m25_600s_bin1_Jan_2019',
                                    master_bias=mbias, median_pass=True)




