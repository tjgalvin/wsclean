import pytest
import os
import warnings
from subprocess import check_call
import shutil

# Append current directory to system path in order to import testconfig
import sys

sys.path.append(".")

import testconfig as tcf

"""
Test script containing a collection of wsclean commands, tested on an MWA measurement set. Tests contained in this
file can be invoked via various routes:

- execute "make checksystemtests"  in your build directory
- execute "[python3 -m] pytest [OPTIONS] source/tests.py" in your build directory
- execute "./various-tests.sh <MWA_MS> [optional] <MWA_COEFFS_PATH>" from the scripts
  directory in the repository
"""

MWA_MS = "MWA-1052736496-averaged.ms"
MWA_COEFFS = "mwa_full_embedded_element_pattern"
TEST_RESULTS = "test_results"
CWD = os.getcwd()


def name(name: str):
    return os.path.join(TEST_RESULTS, name)


@pytest.fixture(autouse=True)
def startup():
    # Fixture is ran for every test function

    # Check if file exist in test_data
    os.makedirs(tcf.WORKDIR, exist_ok=True)
    os.makedirs(TEST_RESULTS, exist_ok=True)
    os.chdir(tcf.WORKDIR)

    if not "MWA_MS" in os.environ:
        if not os.path.isfile(f"{MWA_MS}/table.f1"):
            check_call(
                ["wget", "-q", f"www.astron.nl/citt/ci_data/EveryBeam/{MWA_MS}.tgz"]
            )
            check_call(["tar", "-xf", f"{MWA_MS}.tgz"])
            os.remove(f"{MWA_MS}.tgz")
        os.environ["MWA_MS"] = os.path.join(os.getcwd(), MWA_MS)

    if not "MWA_COEFFS_PATH" in os.environ:
        if not os.path.isfile(f"{MWA_COEFFS}.h5"):
            check_call(
                ["wget", "-q", f"www.astron.nl/citt/EveryBeam/{MWA_COEFFS}.tar.bz2"]
            )
            check_call(["tar", "-xf", f"{MWA_COEFFS}.tar.bz2"])
            os.remove(f"{MWA_COEFFS}.tar.bz2")
        os.environ["MWA_COEFFS_PATH"] = os.getcwd()

    # change to original working directory
    os.chdir(CWD)


@pytest.fixture(scope="session", autouse=True)
def cleanup(request):
    # Fixture is run only at end of test session
    # NOTE: to enable clean-up of the
    # MWA ms and the coefficients file, make sure
    # CLEANUP is in your environment variables

    # Remove measurement set
    def remove_mwa_ms():
        if "CLEANUP" in os.environ:
            shutil.rmtree(os.environ["MWA_MS"])
            os.remove(os.path.join(os.environ["MWA_COEFFS_PATH"], f"{MWA_COEFFS}.h5"))

    request.addfinalizer(remove_mwa_ms)


def test_dirty_image():
    # Make dirty image
    s = f"./wsclean -name {name('test-dirty')} {tcf.DIMS} {os.environ['MWA_MS']}"
    check_call(s.split())


def test_clean_rectangular_unpadded_image():
    # Clean a rectangular unpadded image
    s = f"./wsclean -name {name('clean-rectangular')} -padding 1 \
          -auto-threshold 5 -mgain 0.8 -niter 100000 {tcf.RECTDIMS} {os.environ['MWA_MS']}"
    check_call(s.split())


def test_automask_multiscale_clean():
    # Auto-masked multi-scale clean
    s = f"./wsclean -name {name('multiscale-automasked')} -auto-threshold 0.5 -auto-mask 3 \
          -mgain 0.8 -multiscale -niter 100000 {tcf.RECTDIMS} {os.environ['MWA_MS']}"
    check_call(s.split())


def test_multiple_intervals():
    # Multiple intervals
    s = f"./wsclean -name {name('intervals')} -intervals-out 3 {tcf.RECTDIMS} {os.environ['MWA_MS']}"
    check_call(s.split())


def test_multiple_intervals_and_channels():
    # Multiple intervals + multiple channels with some cleaning
    s = f"./wsclean -name {name('intervals-and-channels')} -intervals-out 3 \
        -channels-out 2 -niter 1000 -mgain 0.8 {tcf.DIMS} {os.environ['MWA_MS']}"
    check_call(s.split())


def test_multifrequency_hogbom():
    # Multi-frequency Högbom clean, no parallel gridding
    s = f"./wsclean -name {name('mfhogbom')} -channels-out 4 -join-channels -auto-threshold 3 \
        -mgain 0.8 -niter 1000000 {tcf.RECTDIMS} {os.environ['MWA_MS']}"
    check_call(s.split())


def test_multifrequency_hogbom_spectral_fit():
    # Multi-frequency Högbom clean with spectral fitting
    s = f"./wsclean -name {name('mfhogbom-fitted')} -channels-out 4 -join-channels -parallel-gridding 4 \
       -fit-spectral-pol 2 -auto-threshold 3 -mgain 0.8 \
           -niter 1000000 {tcf.RECTDIMS} {os.environ['MWA_MS']}"
    check_call(s.split())


def test_mutifrequency_multiscale_parallel():
    # Multi-frequency multi-scale clean with spectral fitting, pallel gridding & cleaning
    s = f"./wsclean -name {name('mfms-fitted')} -channels-out 4 -join-channels -parallel-gridding 4 \
         -parallel-deconvolution 1000 -fit-spectral-pol 2 -multiscale -auto-threshold 0.5 \
              -auto-mask 3 -mgain 0.8 -niter 1000000 {tcf.RECTDIMS} {os.environ['MWA_MS']}"
    check_call(s.split())


def test_save_components():
    # Save the list of components
    s = f"./wsclean -name {name('mfms-components')} -save-source-list -channels-out 4 \
        -join-channels -parallel-gridding 4 -fit-spectral-pol 2 \
            -auto-threshold 0.5 -auto-mask 3 -mgain 0.8 -niter 1000000 \
                -multiscale -parallel-deconvolution 1000 {tcf.DIMS} {os.environ['MWA_MS']}"
    check_call(s.split())


def test_linear_joined_polarizations():
    # Linear joined polarizations with 4 joined channels
    s = f"./wsclean -name {name('linearpol')} -niter 1000000 -auto-threshold 3.0 \
         -pol XX,YY,XY,YX -join-polarizations -join-channels -mgain 0.85 \
             -channels-out 4 -parallel-gridding 16 {tcf.DIMS} {os.environ['MWA_MS']}"
    check_call(s.split())


def test_two_timesteps():
    # Image two timesteps
    s = f"./wsclean -name {name('two-timesteps')} -niter 1000000 -auto-threshold 3.0 \
        -intervals-out 2 -interval 20 22 -mgain 0.85 {tcf.RECTDIMS} {os.environ['MWA_MS']}"
    check_call(s.split())


def test_stop_on_negative_components():
    # Stop on negative components
    s = f"./wsclean -name {name('stop-on-negatives')} -stop-negative -niter 100000 {tcf.RECTDIMS} {os.environ['MWA_MS']}"
    check_call(s.split())

def test_save_imaging_weights():
    s = f"./wsclean -name {name('store-imaging-weights')} -no-reorder -store-imaging-weights {tcf.RECTDIMS} {os.environ['MWA_MS']}"
    check_call(s.split())

@pytest.mark.parametrize(
    "gridder, test_name", (["", "shift-ws"], ["-use-wgridder", "shift-wg"])
)
def test_shift_image(gridder, test_name):
    # Shift the image with w-stacking and w-gridder gridder
    s = f"./wsclean {gridder} -name {name(test_name)} -mgain 0.8 -auto-threshold 5 -niter 1000000 -make-psf {tcf.RECTDIMS} -shift 08h09m20s -39d06m54s -no-update-model-required {os.environ['MWA_MS']}"
    check_call(s.split())


def test_two_facets():
    # Apply the facet to the image
    s = f"./wsclean -name {name('two-facets')} -facet-regions {tcf.FACETFILE_2FACETS} \
        {tcf.RECTDIMS} {os.environ['MWA_MS']}"
    check_call(s.split())


def test_nfacets_pol_xx_yy():
    # Request two polarizations on approximately 25 facets
    s = f"./wsclean -name {name('nfacets-XX_YY')} -pol XX,YY \
        -facet-regions {tcf.FACETFILE_NFACETS} {tcf.RECTDIMS} {os.environ['MWA_MS']}"
    check_call(s.split())


def test_facet_beam():
    # Test facet beam, using 4 polarizations
    if os.environ["MWA_COEFFS_PATH"]:
        s = f"./wsclean -name {name('nfacets-iquv-facet-beam')} -interval 10 14 -apply-facet-beam -pol iquv \
            -facet-regions {tcf.FACETFILE_NFACETS} {tcf.RECTDIMS} \
                -mwa-path {os.environ['MWA_COEFFS_PATH']} {os.environ['MWA_MS']}"
        check_call(s.split())
    else:
        warnings.warn(
            "MWA_PATH environment variable not set, tests that require beam corrections will be skipped."
        )


def test_mpi():
    # Test wsclean-mp command
    s = f"mpirun ./wsclean-mp -name {name('mpi')} {tcf.RECTDIMS} -scale 1amin -channels-out 2 -join-channels -niter 1000000 -mgain 0.8 -auto-threshold 5 -multiscale -no-update-model-required {os.environ['MWA_MS']}"
    check_call(s.split())
