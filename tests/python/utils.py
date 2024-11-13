import os
import shutil
import sys
import warnings
from subprocess import check_call, check_output

import numpy as np
from astropy.io import fits

# Append current directory to system path in order to import testconfig variables
sys.path.append(".")

# Import configuration variables as test configuration (tcf)
import config_vars as tcf


def assert_taql(command, expected_rows=0):
    assert shutil.which("taql") is not None, "taql executable not found!"
    result = check_output(["taql", "-noph", command]).decode().strip()
    assert result == f"select result of {expected_rows} rows", result


def basic_image_check(fits_file):
    """
    Checks that the fits file has no NaN or Inf values and that the number
    of zeroes is below a limit.
    """
    image = fits.open(fits_file)
    assert len(image) == 1
    data = image[0].data
    assert data.shape == (1, 1, 256, 256)
    for x in np.nditer(data):
        assert np.isfinite(x)

    # Test runs showed that 3.1 % of the values are zero. This regression test
    # has a higher limit, since results may vary due to floating point rounding.
    zeros = data.size - np.count_nonzero(data)
    ZERO_LIMIT = 0.035
    assert (zeros / data.size) <= ZERO_LIMIT


def check_and_remove_files(fpaths, remove=False):
    """
    Check whether the entries in the provided list of paths
    are files, and - optionally - remove these files.

    Parameters
    ----------
    fpaths : list
        List of file paths to check
    remove : bool, optional
        Remove files, by default False
    """
    for fpath in fpaths:
        assert os.path.isfile(fpath)

    if remove:
        [os.remove(fpath) for fpath in fpaths]


def compute_rms(fits_file):
    """
    Compute and return the root-mean-square of an input fits image.
    """
    data = fits.open(fits_file)[0].data[0, 0, ...]
    return np.sqrt(np.mean(data**2))


def compare_rms_fits(fits1, fits2, threshold):
    """
    Checks the root-mean square of the difference between
    two input fits files against a user-defined threshold.
    """
    image1 = fits.open(fits1)[0].data.flatten()
    image2 = fits.open(fits2)[0].data.flatten()
    dimage = image1 - image2
    rms1 = np.sqrt(image1.dot(image1) / image1.size)
    rms2 = np.sqrt(image2.dot(image2) / image2.size)
    drms = np.sqrt(dimage.dot(dimage) / dimage.size)
    print(
        f"compare_rms_fits():\n- RMS of {fits1}: {rms1}\n- RMS of {fits2}: {rms2}\n- Difference between {fits1} and {fits2} = {drms}; threshold: {threshold}"
    )
    assert drms <= threshold


def validate_call(cmdline):
    try:
        print("Running: " + " ".join(cmdline))
        check_call(cmdline)
    except:
        # To avoid having to work back what the command was, the command is reported:
        raise RuntimeError(
            "Command failed. To run manually, cd to "
            + tcf.WORKING_DIR
            + " and execute: "
            + " ".join(cmdline)
        )
