import numpy as np
import os
import pytest
from subprocess import check_call, check_output
import shutil
import sys
import warnings

# Append current directory to system path in order to import testconfig
sys.path.append(".")

import testconfig as tcf

MODEL_IMAGE = "point-source-model.fits"
MWA_MOCK_ARCHIVE = "MWA_ARCHIVE.tar.bz2"
MWA_MOCK_MS = "MWA_MOCK.ms"
MWA_MOCK_FULL = "MWA_MOCK_FULL.ms"
MWA_MOCK_FACET = "MWA_MOCK_FACET.ms"
MWA_COEFF_ARCHIVE = "mwa_full_embedded_element_pattern.tar.bz2"
EVERYBEAM_BASE_URL = "http://www.astron.nl/citt/EveryBeam/"
SIZE_SCALE = "-size 256 256 -scale 4amin"

@pytest.fixture(autouse=True)
def prepare():
    # Change to directory containing the data
    os.makedirs(tcf.WORKDIR, exist_ok=True)
    os.chdir(tcf.WORKDIR)

    if not os.path.isfile(MODEL_IMAGE):
        wget = f"wget -q http://www.astron.nl/citt/ci_data/wsclean/{MODEL_IMAGE}"
        check_call(wget.split())

    if not os.path.isfile(MWA_MOCK_ARCHIVE):
        wget = f"wget -q {EVERYBEAM_BASE_URL}MWA-single-timeslot.tar.bz2 -O {MWA_MOCK_ARCHIVE}"
        check_call(wget.split())

    if not os.path.isfile(MWA_COEFF_ARCHIVE):
        check_call(["wget", "-q", EVERYBEAM_BASE_URL + MWA_COEFF_ARCHIVE])
        check_call(["tar", "xf", MWA_COEFF_ARCHIVE])

    os.makedirs(MWA_MOCK_MS, exist_ok=True)
    check_call(
        f"tar -xf {MWA_MOCK_ARCHIVE}  -C {MWA_MOCK_MS} --strip-components=1".split()
    )

    # Not pythonic, but works fine and fast
    check_call(f"cp -r {MWA_MOCK_MS} {MWA_MOCK_FULL}".split())
    check_call(f"cp -r {MWA_MOCK_MS} {MWA_MOCK_FACET}".split())


@pytest.fixture(scope="session", autouse=True)
def cleanup(request):
    # Fixture is run at end of test session and only does
    # executes something if CLEANUP to be in your environment
    def remove_mwa():
        if "CLEANUP" in os.environ:
            os.chdir(tcf.WORKDIR)
            os.remove(MWA_MOCK_ARCHIVE)
            shutil.rmtree(os.environ["MWA_MOCK_FULL"])
            shutil.rmtree(os.environ["MWA_MOCK_FACET"])

    request.addfinalizer(remove_mwa)


def gridders():
    return {"wstacking": "", "wgridder": "-use-wgridder", "idg": "-use-idg"}


def assert_taql(command, expected_rows=0):
    assert shutil.which("taql") is not None, "taql executable not found!"
    result = check_output(["taql", "-noph", command]).decode().strip()
    assert result == f"select result of {expected_rows} rows", result


def predict_full_image(ms, gridder):
    # Predict full image
    s = f"{tcf.WSCLEAN} -predict {gridder} -name point-source {ms}"
    check_call(s.split())


def predict_facet_image(ms, gridder):
    # Predict facet based image
    s = f"{tcf.WSCLEAN} -predict {gridder} -facet-regions {tcf.FACETFILE_4FACETS} -name point-source {ms}"
    check_call(s.split())


def deconvolve_facets(ms, gridder, reorder, mpi, apply_beam=False):
    nthreads = 4
    mpi_cmd = f"mpirun -tag-output -np {nthreads} {tcf.WSCLEAN_MP}"
    thread_cmd = f"{tcf.WSCLEAN} -parallel-gridding {nthreads}"
    reorder_ms = "-reorder" if reorder else "-no-reorder"
    s = [
        mpi_cmd if mpi else thread_cmd,
        gridder,
        SIZE_SCALE,
        reorder_ms,
        "-niter 1000000 -auto-threshold 5 -mgain 0.8",
        f"-facet-regions {tcf.FACETFILE_4FACETS}",
        f"-name facet-imaging{reorder_ms}",
        "-mwa-path . -apply-facet-beam" if apply_beam else "",
        "-v",
        ms,
    ]
    print("WSClean cmd: " + " ".join(s))
    check_call(" ".join(s).split())


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

def basic_image_check(fits_filename):
    """
    Checks that the fits file has no NaN or Inf values and that the number
    of zeroes is below a limit.
    """
    try:
        from astropy.io import fits
    except:
        warnings.warn(
            UserWarning(
                "Could not import astropy, so fits image checks are skipped."
            )
        )
        return

    image = fits.open(fits_filename)
    assert len(image) == 1
    data = image[0].data
    assert data.shape == (1, 1, 256, 256)
    for x in np.nditer(data):
        assert np.isfinite(x)

    # Test runs showed that 3.1 % of the values are zero. This regression test
    # has a higher limit, since results may vary due to floating point rounding.
    zeroes = data.size - np.count_nonzero(data)
    ZERO_LIMIT = 0.035
    assert (zeroes / data.size) <= ZERO_LIMIT

def test_makepsfonly():
    """
    Test that wsclean with the -make-psf-only flag exits gracefully and
    that the psf passes basic checks.
    """
    s = [
        tcf.WSCLEAN,
        "-make-psf-only -name facet-psf-only",
        SIZE_SCALE,
        f"-facet-regions {tcf.FACETFILE_4FACETS} {MWA_MOCK_MS}"
    ]
    check_call(" ".join(s).split())

    basic_image_check("facet-psf-only-psf.fits")

# Test assumes that IDG and EveryBeam are installed
@pytest.mark.parametrize("gridder", gridders().items())
def test_stitching(gridder):
    """Test stitching of the facets"""
    prefix = f"facet-stitch-{gridder[0]}"
    s = [
       tcf.WSCLEAN,
       "-quiet",
       gridder[1],
       SIZE_SCALE,
       "" if (gridder[0] == "idg") else "-pol XX,YY",
       f"-facet-regions {tcf.FACETFILE_2FACETS}",
       f"-name {prefix}",
       MWA_MOCK_MS
    ]
    check_call(" ".join(s).split())
    fpaths = (
        [prefix + "-dirty.fits", prefix + "-image.fits"]
        if (gridder[0] == "idg")
        else [
            prefix + "-XX-dirty.fits",
            prefix + "-YY-dirty.fits",
            prefix + "-XX-image.fits",
            prefix + "-YY-image.fits",
        ]
    )
    check_and_remove_files(fpaths, remove=True)


# FIXME: we should test wstacking and wgridder here too
# but they fail on the taql assertion
@pytest.mark.parametrize("gridder", ["-use-wgridder"])
def test_predict(gridder):
    """
    Test predict only run

    Parameters
    ----------
    gridder : str
        wsclean compatible description of gridder to be used.
    """

    predict_full_image(MWA_MOCK_FULL, gridder)
    predict_facet_image(MWA_MOCK_FACET, gridder)
    taql_command = f"select from {MWA_MOCK_FULL} t1, {MWA_MOCK_FACET} t2 where not all(near(t1.MODEL_DATA,t2.MODEL_DATA,5e-3))"
    assert_taql(taql_command)


@pytest.mark.parametrize("gridder", ["-use-wgridder"])
@pytest.mark.parametrize("reorder", [False, True])
@pytest.mark.parametrize("mpi", [False, True])
def test_facetdeconvolution(gridder, reorder, mpi):
    """
    Test facet-based deconvolution

    Parameters
    ----------
    gridder : str
        wsclean compatible description of gridder to be used.
    reorder : bool
        Reorder MS?
    mpi : bool
        True: Use MPI for parallel gridding.
        False: Use multi-threading for parallel gridding.
    """
    # Parametrization causes some overhead in that predict of full image is run for
    # every parametrization
    predict_full_image(MWA_MOCK_FULL, gridder)

    # Make sure old versions of the facet mock ms are removed
    shutil.rmtree(MWA_MOCK_FACET)

    # Copy output to new MS, swap DATA column, and remove MODEL_DATA
    check_call(f"cp -r {MWA_MOCK_FULL} {MWA_MOCK_FACET}".split())
    assert shutil.which("taql") is not None, "taql executable not found!"

    check_call(["taql", "-noph", f"UPDATE {MWA_MOCK_FACET} SET DATA=MODEL_DATA"])
    check_call(
        ["taql", "-noph", f"ALTER TABLE {MWA_MOCK_FACET} DROP COLUMN MODEL_DATA"]
    )
    taql_command = f"select from {MWA_MOCK_FULL} t1, {MWA_MOCK_FACET} t2 where not all(near(t1.MODEL_DATA,t2.DATA, 4e-3))"
    assert_taql(taql_command)

    deconvolve_facets(MWA_MOCK_FACET, gridder, reorder, mpi)

    taql_command = (
        f"select from {MWA_MOCK_FACET} where not all(near(DATA,MODEL_DATA, 4e-3))"
    )
    assert_taql(taql_command)

@pytest.mark.parametrize("mpi", [False, True])
def test_facetbeamimages(mpi):
    """
    Basic checks of the generated images when using facet beams. For each image,
    test that the pixel values are valid (not NaN/Inf) and check the percentage
    of zero pixels.
    """

    deconvolve_facets(MWA_MOCK_FACET, "-use-wgridder", True, mpi, True)

    basic_image_check("facet-imaging-reorder-psf.fits")
    basic_image_check("facet-imaging-reorder-dirty.fits")
