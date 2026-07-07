import os
import shutil
import sys
import warnings
from subprocess import check_call, check_output

import numpy as np
import pytest
from astropy.io import fits

# Append current directory to system path in order to import testconfig variables
sys.path.append(".")

# Import configuration variables as test configuration (tcf)
import config_vars as tcf


def assert_taql(command, expected_rows=0):
    assert shutil.which("taql") is not None, "taql executable not found!"
    print("Running: ", " ".join(["taql", "-noph", command]))
    result = check_output(["taql", "-noph", command]).decode().strip()
    assert result == f"select result of {expected_rows} rows", result


def check_image_values(fits_file, maximum_zero_fraction):
    """
    Checks that the fits file has no NaN or Inf values and that the number
    of zeroes is below a limit.
    """
    image = fits.open(fits_file)
    assert len(image) == 1
    data = image[0].data
    # Does the image have a limited nr of zeros?
    zeros = data.size - np.count_nonzero(data)
    assert (zeros / data.size) <= maximum_zero_fraction
    for x in np.nditer(data):
        assert np.isfinite(x)


def basic_mwa_image_check(fits_file):
    """
    Checks that the fits file has no NaN or Inf values and that the number
    of zeroes is below a limit.
    """
    image = fits.open(fits_file)
    assert len(image) == 1
    data = image[0].data
    assert data.shape[0] == 1
    assert data.shape[1] == 1
    assert data.shape[2] >= 256
    assert data.shape[3] >= 256

    check_image_values(fits_file, 0.05)


def basic_facet_image_check(fits_file):
    """
    Same as basic_mwa_image_check(), but more specific for facet image settings.
    """
    image = fits.open(fits_file)
    assert len(image) == 1
    data = image[0].data
    assert data.shape == (1, 1, 256, 256)

    # Does the image have a limited nr of zeros?
    # Test runs showed that 3.1 % of the values are zero. This regression test
    # has a higher limit, since results may vary due to floating point rounding.
    check_image_values(fits_file, 0.035)


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


def validate_failing_call(cmdline):
    got_exception = True
    try:
        print("Running: " + " ".join(cmdline))
        check_call(cmdline)
        got_exception = False
    except:
        pass

    if not got_exception:
        raise RuntimeError(
            "Command was expected to fail, but didn't. To run manually, cd to "
            + tcf.WORKING_DIR
            + " and execute: "
            + " ".join(cmdline)
        )


def require_idg():
    if "IDG" not in check_output([tcf.WSCLEAN, "--version"]).decode():
        pytest.skip("WSClean was not compiled with IDG.")


def require_everybeam():
    if "EveryBeam" not in check_output([tcf.WSCLEAN, "--version"]).decode():
        pytest.skip("WSClean was not compiled with EveryBeam.")
