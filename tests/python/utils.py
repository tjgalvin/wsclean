import os
import shutil
import sys
import warnings
from subprocess import check_call, check_output

import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import mad_std
from astropy.wcs import WCS

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


def compare_psf_coordinate_and_amplitude_accuracy(
    fits_baseline,
    fits_compare,
    position_threshold_arcsec=0.001,
    amplitude_threshold=0.00001,
):
    """
    Validates the physical celestial coordinate alignment and peak amplitude consistency
    of the synthesized Point Spread Function (PSF).

    Extracts the World Coordinate System (WCS) information from FITS headers to track
    visibility phase errors via true great-circle angular separation at the PSF peak.
    """
    with fits.open(fits_baseline) as h1, fits.open(fits_compare) as h2:
        wcs1 = WCS(
            h1[0].header if h1[0].data is not None else h1[1].header, naxis=2
        )
        wcs2 = WCS(
            h2[0].header if h2[0].data is not None else h2[1].header, naxis=2
        )
        img1_raw = h1[0].data if h1[0].data is not None else h1[1].data
        img2_raw = h2[0].data if h2[0].data is not None else h2[1].data
        img1 = np.nan_to_num(img1_raw).squeeze()
        img2 = np.nan_to_num(img2_raw).squeeze()

    y1, x1 = np.unravel_index(np.argmax(img1), img1.shape)
    y2, x2 = np.unravel_index(np.argmax(img2), img2.shape)

    coord1 = SkyCoord.from_pixel(x1, y1, wcs1.celestial)
    coord2 = SkyCoord.from_pixel(x2, y2, wcs2.celestial)

    angular_separation_arcsec = coord1.separation(coord2).arcsecond
    amplitude_diff = np.abs(img1[y1, x1] - img2[y2, x2])

    print(f"\nPSF Verification: {fits_baseline} vs {fits_compare}")
    print(
        f"- Baseline Peak:   RA={coord1.ra.deg:.6f}°, Dec={coord1.dec.deg:.6f}°"
    )
    print(
        f"- Comparison Peak: RA={coord2.ra.deg:.6f}°, Dec={coord2.dec.deg:.6f}°"
    )
    print(
        f"- Angular Shift:  {angular_separation_arcsec:.6e} arcsec (Target: < {position_threshold_arcsec})"
    )
    print(
        f"- Amplitude Diff: {amplitude_diff:.6e} (Target: < {amplitude_threshold})"
    )

    assert (
        angular_separation_arcsec < position_threshold_arcsec
    ), f"❌ PSF coordinate shift: {angular_separation_arcsec:.6e} arcsec exceeds threshold {position_threshold_arcsec}"
    assert (
        amplitude_diff < amplitude_threshold
    ), f"❌ PSF amplitude mismatch: {amplitude_diff:.6e} exceeds threshold {amplitude_threshold}"

    print("✅ PSF Verification Passed.")


def compare_cleaned_images_equivalence(
    fits_baseline, fits_compare, rms_target_pct, lenient_drms_factor
):
    """
    Validates that two generated images are scientifically equivalent.
    Accounts for cleaning loop divergence by combining global energy,
    spatial dRMS, and peak coordinate validation.
    """
    with fits.open(fits_baseline) as h1, fits.open(fits_compare) as h2:
        img1_raw = h1[0].data if h1[0].data is not None else h1[1].data
        img2_raw = h2[0].data if h2[0].data is not None else h2[1].data
        img1 = np.copy(np.nan_to_num(img1_raw)).flatten()
        img2 = np.copy(np.nan_to_num(img2_raw)).flatten()

    # TEST 1: Global RMS Percentage Difference
    rms1 = np.sqrt(np.mean(img1**2))
    rms2 = np.sqrt(np.mean(img2**2))
    rms_diff_pct = (abs(rms1 - rms2) / rms1) * 100

    # TEST 2: Spatial dRMS Threshold
    dimage = img1 - img2
    drms = np.sqrt(np.mean(dimage**2))
    bg_noise = mad_std(img1)
    lenient_drms_threshold = bg_noise * lenient_drms_factor
    precision_floor = 1e-6
    if lenient_drms_threshold < precision_floor:
        lenient_drms_threshold = precision_floor

    # TEST 3: Coordinate Shift Protection (Peak Alignment)
    peak_idx1 = np.argmax(img1)
    peak_idx2 = np.argmax(img2)

    # Diagnostic Printout
    print("\n--- Image Equivalence Test Report ---")
    print(f"Baseline:   {fits_baseline}")
    print(f"Comparison: {fits_compare}")
    print(
        f"- Global RMS Diff: {rms_diff_pct:.6f}% (Target: < {rms_target_pct}%)"
    )
    print(
        f"- Spatial dRMS:    {drms:.6e} (Threshold: < {lenient_drms_threshold:.6e})"
    )
    print(
        f"- Peak Alignment:  Baseline Index {peak_idx1} vs Comparison Index {peak_idx2}"
    )

    # Check 1: Global energy scaling
    assert (
        rms_diff_pct < rms_target_pct
    ), f"❌ Global energy mismatch: Diff is {rms_diff_pct:.6f}%, expecting < {rms_target_pct}%"
    # Check 2: Spatial structure fail-safe
    assert (
        drms <= lenient_drms_threshold
    ), f"❌ Spatial structure mismatch: dRMS ({drms:.6e}) exceeds threshold ({lenient_drms_threshold:.6e})"
    # Check 3: Mapping/Coordinate shifts
    assert (
        peak_idx1 == peak_idx2
    ), f"❌ Coordinate shift detected: Peak pixel moved from index {peak_idx1} to {peak_idx2}"
    print("✅ Test Passed: Images are within acceptable tolerances.")


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
