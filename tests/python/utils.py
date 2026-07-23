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


def get_wcs_and_2d_image(hdul):
    """Safely extracts a 2D celestial slice and its corresponding 2D WCS from an HDUList,

    automatically finding the first extension containing pixel data.
    """
    hdu = None
    for extension in hdul:
        if extension.data is not None:
            hdu = extension
            break

    if hdu is None:
        raise ValueError(
            "❌ No valid image data found in any FITS HDU extension."
        )

    header = hdu.header.copy()
    data = hdu.data

    full_wcs = WCS(header)
    celestial_wcs = full_wcs.celestial

    celestial_axes = [
        i
        for i, axis in enumerate(full_wcs.axis_type_names)
        if axis in ["RA", "DEC"]
    ]

    slices = tuple(
        slice(None) if i in celestial_axes else 0
        for i in reversed(range(data.ndim))
    )
    return celestial_wcs, np.nan_to_num(data[slices])


def compare_psf_spatial_and_astrometric_integrity(
    fits_baseline,
    fits_compare,
    position_threshold_arcsec=0.001,
    amplitude_threshold=0.00001,
    structural_threshold=0.95,  # Adjusted to a realistic, strict threshold for isolated structure
    correlation_window_size=128,
):
    """Compares the complete 2D physical and structural alignment of the Point

    Spread Function (PSF) using a localized structural mask.

    Isolates real sidelobe features from background noise to prevent false passes
    when comparing distinctly structured PSFs.
    """
    with fits.open(fits_baseline) as h1, fits.open(fits_compare) as h2:
        wcs1, img1 = get_wcs_and_2d_image(h1)
        wcs2, img2 = get_wcs_and_2d_image(h2)

    # 1. CORE SEARCH: Locate maximum pixel positions
    y1, x1 = np.unravel_index(np.argmax(img1), img1.shape)
    y2, x2 = np.unravel_index(np.argmax(img2), img2.shape)

    coord1 = SkyCoord.from_pixel(x1, y1, wcs1)
    coord2 = SkyCoord.from_pixel(x2, y2, wcs2)

    angular_separation_arcsec = coord1.separation(coord2).arcsecond
    amplitude_diff = np.abs(img1[y1, x1] - img2[y2, x2])

    # 2. STRUCTURAL ISOLATION: Crop the active core region
    half_w = correlation_window_size // 2
    y1_min, y1_max = max(0, y1 - half_w), min(img1.shape[0], y1 + half_w)
    x1_min, x1_max = max(0, x1 - half_w), min(img1.shape[1], x1 + half_w)
    y2_min, y2_max = max(0, y2 - half_w), min(img2.shape[0], y2 + half_w)
    x2_min, x2_max = max(0, x2 - half_w), min(img2.shape[1], x2 + half_w)

    sub_img1 = img1[y1_min:y1_max, x1_min:x1_max]
    sub_img2 = img2[y2_min:y2_max, x2_min:x2_max]

    if sub_img1.shape == sub_img2.shape:
        # Calculate a dynamic threshold to isolate the core & structural sidelobes.
        # Uses standard deviation of the patch perimeter as a rough background noise floor proxy.
        edge_pixels = np.concatenate(
            [sub_img1[0, :], sub_img1[-1, :], sub_img1[:, 0], sub_img1[:, -1]]
        )
        bg_std = np.std(edge_pixels) if np.std(edge_pixels) > 0 else 1e-6

        # Create a binary mask where baseline features stand significantly out of the noise (3-sigma)
        structural_mask = sub_img1 > (3 * bg_std)

        # Ensure we actually captured features (fallback to the core pixel if the patch is empty)
        if np.sum(structural_mask) < 5:
            structural_mask = sub_img1 > (0.01 * np.max(sub_img1))

        # Extract only the structural vectors
        feat1 = sub_img1[structural_mask]
        feat2 = sub_img2[structural_mask]

        # Compute Pearson R correlation strictly on the masked structural pixels
        diff1 = feat1 - np.mean(feat1)
        diff2 = feat2 - np.mean(feat2)
        std1 = np.sqrt(np.sum(diff1**2))
        std2 = np.sqrt(np.sum(diff2**2))

        if std1 > 0 and std2 > 0:
            structural_correlation = np.sum(diff1 * diff2) / (std1 * std2)
        else:
            structural_correlation = 0.0
    else:
        structural_correlation = 0.0

    # 3. DIAGNOSTIC PRINT OUTS
    print(
        f"\n=== Masked PSF Integrity Audit: {fits_baseline} vs {fits_compare} ==="
    )
    print(
        f"- Baseline Peak:     RA={coord1.ra.deg:.6f}°, Dec={coord1.dec.deg:.6f}°"
    )
    print(
        f"- Comparison Peak:   RA={coord2.ra.deg:.6f}°, Dec={coord2.dec.deg:.6f}°"
    )
    print(
        f"- Astrometric Shift: {angular_separation_arcsec:.6e} arcsec (Target: < {position_threshold_arcsec})"
    )
    print(
        f"- Intensity Diff:    {amplitude_diff:.6e} (Target: < {amplitude_threshold})"
    )
    print(
        f"- Masked Profile Match: {structural_correlation:.6f} (Target: >= {structural_threshold})"
    )
    print(f"- Active Mask Size:  {np.sum(structural_mask)} pixels")

    # 4. ASSERTION ENFORCEMENT
    assert (
        angular_separation_arcsec < position_threshold_arcsec
    ), f"❌ PSF coordinate shift: {angular_separation_arcsec:.6e} arcsec exceeds threshold {position_threshold_arcsec}"

    assert (
        amplitude_diff < amplitude_threshold
    ), f"❌ PSF amplitude mismatch: {amplitude_diff:.6e} exceeds threshold {amplitude_threshold}"

    assert (
        structural_correlation >= structural_threshold
    ), f"❌ PSF structural divergence: Masked correlation {structural_correlation:.5f} falls below threshold {structural_threshold}"

    print("✅ Complete PSF Structural and Astrometric Integrity Passed.")


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
