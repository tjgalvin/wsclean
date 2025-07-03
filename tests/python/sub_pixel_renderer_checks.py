import os
import sys
from pathlib import Path

import numpy as np
import pytest
from astropy import wcs
from astropy.io import fits
from utils import check_and_remove_files, compare_rms_fits, validate_call

# Append current directory to system path in order to import testconfig
sys.path.append(".")

# Import configuration variables as test configuration (tcf)
import config_vars as tcf


@pytest.fixture(scope="class")
def create_mock_sky_model():
    """
    Write a sky model (mock-model.txt) with 1 point source, 1 small Gaussian source and 1 larger Gaussian source.
    """
    mock_model = """FORMAT = Name, Type, Patch, Ra, Dec, I, SpectralIndex='[]', LogarithmicSI, ReferenceFrequency='150000000.0', MajorAxis, MinorAxis, Orientation

, , source_a, 08:21:00, -42.45.00
point, POINT, source_a, 08:21:00, -42.45.00, 2.0, [-0.05], true, 150000000.0, 0.0, 0.0, 0.0

, , source_b, 08:20:00, -42.45.00
big_gauss, GAUSSIAN, source_b, 08:20:00, -42.45.00, 50.0, [-0.05], true, 150000000.0, 250.0, 250.0, 90.0

, , source_c, 08:19:00, -42.45.00
small_gauss, GAUSSIAN, source_c, 08:19:00, -42.45.00, 1.0, [-0.05], true, 150000000.0, 50.0, 20.0, 45.0
"""

    with open("mock-model.txt", "w") as mock_file:
        mock_file.write(mock_model)


def generate_full_image_facet(image_filename, region_filename):
    """
    Creates a facet coverering the entire (square) image.

    Parameters
    ----------
    image_filename: str
        Name of the input FITS image
    region_filename: str
        Name of the output region file describing the facet

    Returns
    -------
    polygon: str
        DS9 compatible polygon describing the facet
    """
    hdu_list = fits.open(image_filename)
    w = wcs.WCS(hdu_list[0].header, hdu_list)

    image_size = hdu_list[0].data.shape[2]

    # Note: first 2 coordinates are x and y, the last two are frequency and stokes parameter; the latter two are not used!
    corner_pixels = np.array(
        [
            [0, 0, 0, 0],
            [0, image_size - 1, 0, 0],
            [image_size - 1, image_size - 1, 0, 0],
            [image_size - 1, 0, 0, 0],
        ]
    )
    world_coordinates = w.wcs_pix2world(corner_pixels, 1)

    polygon_vertices = []

    for ra, dec, freq, stokes in world_coordinates:
        polygon_vertices.append(ra)
        polygon_vertices.append(dec)

    polygon = "polygon({})".format(
        np.array2string(
            np.array(polygon_vertices), precision=3, separator=","
        )[1:-1]
    )

    with open(region_filename, "w") as region_file:
        print(polygon, file=region_file)

    return polygon


def run_dp3_wgridder_predict(ms_in, ms_out, model_prefix, n_terms):
    """
    Run DP3's WGridderPredict to generate model visibilities based
    on model images produced by WSClean's sub-pixel renderer.

    Parameters
    ----------
    ms_in: str
        Name of the input MeasurementSet
    ms_out: str
        Name of the output MeasurementSet
    model_prefix: str
        Model image prefix (i.e., before -term-0.fits)
    n_terms: int
        The number of term images WGridderPredict should use
    """
    model_images = [f"{model_prefix}-term-{i}.fits" for i in range(n_terms)]
    region_filename = f"{model_prefix}.reg"
    generate_full_image_facet(model_images[0], region_filename)
    term_images_string = ",".join([filename for filename in model_images])
    dp3_wgridder_predict_run = [
        "DP3",
        "checkparset=1",
        f"msin={ms_in}",
        f"msout={ms_out}",
        "msout.overwrite=True",
        "steps=[wgridderpredict]",
        f"wgridderpredict.images=[{term_images_string}]",
        f"wgridderpredict.regions={region_filename}",
        "wgridderpredict.sumfacets=True",
    ]
    validate_call(dp3_wgridder_predict_run)


@pytest.mark.usefixtures(
    "prepare_mock_ms", "create_mock_sky_model", "prepare_3c196_sky_model"
)
class TestSubPixelRenderer:
    def test_sub_pixel_renderer(self):
        # Execute the sub-pixel renderer with the -draw-model option
        s = f"{tcf.WSCLEAN} -draw-frequencies 150e6 10e6 {tcf.DIMS_LARGE} -draw-model mock-model.txt -name test {tcf.MWA_MOCK_MS}"
        validate_call(s.split())

        # Verify that the expected model image has been created as expected
        expected_model_image_filename = "test-term-0.fits"
        check_and_remove_files([expected_model_image_filename], remove=True)

    def test_direct_against_image_based_predict(self):
        # Create a model image using the -draw-model option
        n_terms = 2
        prefix = "model-image"
        sky_model = "mock-model.txt"
        s = f"{tcf.WSCLEAN} -draw-frequencies 150e6 10e6 {tcf.DIMS_LARGE} -draw-model {sky_model} -draw-spectral-terms {n_terms} -name {prefix} {tcf.MWA_MOCK_MS}"
        validate_call(s.split())

        # Run image-based predict with DP3
        image_based_predict_ms = "image-based-predict.ms"
        run_dp3_wgridder_predict(
            tcf.MWA_MOCK_MS, image_based_predict_ms, prefix, n_terms
        )

        # Run direct predict with DP3
        dp3_run = [
            "DP3",
            "checkparset=1",
            f"msin={tcf.MWA_MOCK_MS}",
            "msout=direct-predict.ms",
            "msout.overwrite=True",
            "steps=[predict]",
            f"predict.sourcedb={sky_model}",
        ]
        validate_call(dp3_run)

        # Image the predicted visibilities
        s = f"{tcf.WSCLEAN} -name direct-predict {tcf.DIMS_LARGE} direct-predict.ms"
        validate_call(s.split())

        s = f"{tcf.WSCLEAN} -name image-based-predict {tcf.DIMS_LARGE} image-based-predict.ms"
        validate_call(s.split())

        # Compare the RMS of the residual image
        compare_rms_fits(
            "direct-predict-dirty.fits",
            "image-based-predict-dirty.fits",
            threshold=1.0e-2,
        )

    def test_image_based_predict_at_different_resolution(self):
        n_terms = 3
        window_size = 256

        # Step 1: create a low-resolution model image, with 0.5 arcmin/pixel resolution
        low_res_prefix = "low-res-model-image"
        low_res_dims = "-size 512 512 -scale 0.5arcsec"
        s = f"{tcf.WSCLEAN} -draw-model {tcf.SKYMODEL_3C196} -draw-frequencies 150e6 10e6 {low_res_dims} -draw-spectral-terms {n_terms} -sinc-window-size {window_size} -name {low_res_prefix} {tcf.LOFAR_3C196_MS}"
        validate_call(s.split())

        # Run image-based predict with DP3
        low_res_ms = "low-res-predict.ms"
        run_dp3_wgridder_predict(
            tcf.LOFAR_3C196_MS, low_res_ms, low_res_prefix, n_terms
        )

        # Step 2: create a high-resolution model image, with 0.05 arcmin/pixel resolution
        high_res_prefix = "high-res-model-image"
        high_res_dims = "-size 5120 5120 -scale 0.05arcsec"
        s = f"{tcf.WSCLEAN} -draw-model {tcf.SKYMODEL_3C196} -draw-frequencies 150e6 10e6 {high_res_dims} -draw-spectral-terms {n_terms} -sinc-window-size {window_size} -name {high_res_prefix} {tcf.LOFAR_3C196_MS}"
        validate_call(s.split())

        # Run image-based predict with DP3
        high_res_ms = "high-res-predict.ms"
        run_dp3_wgridder_predict(
            tcf.LOFAR_3C196_MS, high_res_ms, high_res_prefix, n_terms
        )

        # Image the predicted visibilities
        s = f"{tcf.WSCLEAN} -name low-res-predict {low_res_dims} {low_res_ms}"
        validate_call(s.split())

        s = f"{tcf.WSCLEAN} -name high-res-predict {low_res_dims} {high_res_ms}"
        validate_call(s.split())

        # Step 3: Compare the RMS of the residual image
        compare_rms_fits(
            "low-res-predict-dirty.fits",
            "high-res-predict-dirty.fits",
            threshold=1.0e-3,
        )
