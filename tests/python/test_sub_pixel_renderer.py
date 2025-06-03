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


@pytest.mark.usefixtures("prepare_mock_ms", "create_mock_sky_model")
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
        model_images = [f"{prefix}-term-{i}.fits" for i in range(n_terms)]
        generate_full_image_facet(model_images[0], "mock.reg")
        term_images_string = ",".join([filename for filename in model_images])
        dp3_run = [
            "DP3",
            "checkparset=1",
            f"msin={tcf.MWA_MOCK_MS}",
            "msout=image-based-predict.ms",
            "msout.overwrite=True",
            "steps=[wgridderpredict]",
            f"wgridderpredict.images=[{term_images_string}]",
            f"wgridderpredict.regions=mock.reg",
            "wgridderpredict.sumfacets=True",
        ]
        validate_call(dp3_run)

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
