import pytest
import os, glob
import warnings
from subprocess import check_call
import shutil
import sys

# Append current directory to system path in order to import testconfig
sys.path.append(".")

import testconfig as tcf

"""
Test script containing a collection of wsclean commands, tested on an MWA measurement set. Tests contained in this
file can be invoked via various routes:

- execute "make checksystemtests"  in your build directory
- execute "[python3 -m] pytest [OPTIONS] source/<test_name.py>" in your build/tests directory
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
    """Fixture is ran for every test function"""

    # Make directory for test results
    os.makedirs(os.path.join(CWD, TEST_RESULTS), exist_ok=True)

    # Check if file exist in test_data
    os.makedirs(tcf.WORKDIR, exist_ok=True)
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
        shutil.rmtree(os.environ["MWA_MS"])
        os.remove(os.path.join(os.environ["MWA_COEFFS_PATH"], f"{MWA_COEFFS}.h5"))

    def remove_fits_images():
        for f in glob.glob("*.fits"):
            os.remove(f)

    def remove_h5_files():
        for f in glob.glob("*.h5"):
            os.remove(f)

    def collect_cleanup():
        if "CLEANUP" in os.environ:
            remove_mwa_ms()
            remove_fits_images()
            remove_h5_files()

    request.addfinalizer(collect_cleanup)


def test_dirty_image():
    # Make dirty image
    s = f"{tcf.WSCLEAN} -name {name('test-dirty')} {tcf.DIMS} {os.environ['MWA_MS']}"
    check_call(s.split())


def test_clean_rectangular_unpadded_image():
    # Clean a rectangular unpadded image
    s = f"{tcf.WSCLEAN} -name {name('clean-rectangular')} -padding 1 \
          -auto-threshold 5 -mgain 0.8 -niter 100000 {tcf.RECTDIMS} {os.environ['MWA_MS']}"
    check_call(s.split())


def test_automask_multiscale_clean():
    # Auto-masked multi-scale clean
    s = f"{tcf.WSCLEAN} -name {name('multiscale-automasked')} -auto-threshold 0.5 -auto-mask 3 \
          -mgain 0.8 -multiscale -niter 100000 {tcf.RECTDIMS} {os.environ['MWA_MS']}"
    check_call(s.split())


def test_multiple_intervals():
    # Multiple intervals
    s = f"{tcf.WSCLEAN} -name {name('intervals')} -intervals-out 3 {tcf.RECTDIMS} {os.environ['MWA_MS']}"
    check_call(s.split())


def test_multiple_intervals_and_channels():
    # Multiple intervals + multiple channels with some cleaning
    s = f"{tcf.WSCLEAN} -name {name('intervals-and-channels')} -intervals-out 3 \
        -channels-out 2 -niter 1000 -mgain 0.8 {tcf.DIMS} {os.environ['MWA_MS']}"
    check_call(s.split())


def test_multifrequency_hogbom():
    # Multi-frequency Högbom clean, no parallel gridding
    s = f"{tcf.WSCLEAN} -name {name('mfhogbom')} -channels-out 4 -join-channels -auto-threshold 3 \
        -mgain 0.8 -niter 1000000 {tcf.RECTDIMS} {os.environ['MWA_MS']}"
    check_call(s.split())


def test_multifrequency_hogbom_spectral_fit():
    # Multi-frequency Högbom clean with spectral fitting
    s = f"{tcf.WSCLEAN} -name {name('mfhogbom-fitted')} -channels-out 4 -join-channels -parallel-gridding 4 \
       -fit-spectral-pol 2 -auto-threshold 3 -mgain 0.8 \
           -niter 1000000 {tcf.RECTDIMS} {os.environ['MWA_MS']}"
    check_call(s.split())


def test_mutifrequency_multiscale_parallel():
    # Multi-frequency multi-scale clean with spectral fitting, pallel gridding & cleaning
    s = f"{tcf.WSCLEAN} -name {name('mfms-fitted')} -channels-out 4 -join-channels -parallel-gridding 4 \
         -parallel-deconvolution 1000 -fit-spectral-pol 2 -multiscale -auto-threshold 0.5 \
              -auto-mask 3 -mgain 0.8 -niter 1000000 {tcf.RECTDIMS} {os.environ['MWA_MS']}"
    check_call(s.split())


def test_save_components():
    test_name = name("mfms-components")

    # Remove old component list if it exists
    component_file = test_name + "-sources.txt"
    if os.path.exists(component_file):
        os.remove(component_file)

    # Save the list of components
    s = f"{tcf.WSCLEAN} -name {test_name} -save-source-list -channels-out 4 \
        -join-channels -parallel-gridding 4 -fit-spectral-pol 2 \
            -auto-threshold 0.5 -auto-mask 3 -mgain 0.8 -niter 1000000 \
                -multiscale -parallel-deconvolution 1000 {tcf.DIMS} {os.environ['MWA_MS']}"
    check_call(s.split())

    # Check whether source files is generated
    assert os.path.isfile(component_file)


def test_linear_joined_polarizations():
    # Linear joined polarizations with 4 joined channels
    s = f"{tcf.WSCLEAN} -name {name('linearpol')} -niter 1000000 -auto-threshold 3.0 \
         -pol XX,YY,XY,YX -join-polarizations -join-channels -mgain 0.85 \
             -channels-out 4 -parallel-gridding 16 {tcf.DIMS} {os.environ['MWA_MS']}"
    check_call(s.split())


def test_two_timesteps():
    # Image two timesteps
    s = f"{tcf.WSCLEAN} -name {name('two-timesteps')} -niter 1000000 -auto-threshold 3.0 \
        -intervals-out 2 -interval 20 22 -mgain 0.85 {tcf.RECTDIMS} {os.environ['MWA_MS']}"
    check_call(s.split())


def test_stop_on_negative_components():
    # Stop on negative components
    s = f"{tcf.WSCLEAN} -name {name('stop-on-negatives')} -stop-negative -niter 100000 {tcf.RECTDIMS} {os.environ['MWA_MS']}"
    check_call(s.split())


def test_save_imaging_weights():
    s = f"{tcf.WSCLEAN} -name {name('store-imaging-weights')} -no-reorder -store-imaging-weights {tcf.RECTDIMS} {os.environ['MWA_MS']}"
    check_call(s.split())


@pytest.mark.parametrize(
    "gridder, test_name", (["", "shift-ws"], ["-use-wgridder", "shift-wg"])
)
def test_shift_image(gridder, test_name):
    # Shift the image with w-stacking and w-gridder gridder
    s = f"{tcf.WSCLEAN} {gridder} -name {name(test_name)} -mgain 0.8 -auto-threshold 5 -niter 1000000 -make-psf {tcf.RECTDIMS} -shift 08h09m20s -39d06m54s -no-update-model-required {os.environ['MWA_MS']}"
    check_call(s.split())


def test_missing_channels_in_deconvolution():
    # The test set has some missing MWA subbands. One MWA subband is 1/24 of the data (32/768 channels), so
    # by imaging with -channels-out 24, it is tested what happens when an output channel has no data.
    s = f"{tcf.WSCLEAN} -name {name('missing-channels-in-deconvolution')} -use-wgridder -size 1024 1024 -scale 1amin -baseline-averaging 2.0 -no-update-model-required -niter 150000 -auto-threshold 2.0 -auto-mask 5.0 -mgain 0.9 -channels-out 24 -join-channels -fit-spectral-pol 4 {os.environ['MWA_MS']}"
    check_call(s.split())


def test_grid_with_beam():
    """Requires that WSClean is compiled with IDG and EveryBeam
    """
    name = "idg-beam"
    if os.environ["MWA_COEFFS_PATH"]:
        # Remove existing component files if present
        for source_file in ["sources", "sources-pb"]:
            component_file = name + "-" + source_file + ".txt"
            if os.path.exists(component_file):
                os.remove(component_file)

        s = f"{tcf.WSCLEAN} -name {name} -use-idg -grid-with-beam -save-source-list -mgain 0.8 -auto-threshold 5 -niter 1000000 -interval 10 14 {tcf.DIMS} -mwa-path {os.environ['MWA_COEFFS_PATH']} {os.environ['MWA_MS']}"
        check_call(s.split())
        for image_type in [
            "psf",
            "beam",
            "dirty",
            "image",
            "image-pb",
            "model",
            "model-pb",
            "residual",
            "residual-pb",
        ]:
            image_name = name + "-" + image_type + ".fits"
            assert os.path.isfile(image_name)
        # Check whether source files are generated
        for source_file in ["sources", "sources-pb"]:
            assert os.path.isfile(name + "-" + source_file + ".txt")
    else:
        warnings.warn(
            "MWA_PATH environment variable not set, test_grid_with_beam test will be skipped."
        )


def test_two_facets():
    # Apply the facet to the image
    s = f"{tcf.WSCLEAN} -name {name('two-facets')} -facet-regions {tcf.FACETFILE_2FACETS} \
        {tcf.RECTDIMS} {os.environ['MWA_MS']}"
    check_call(s.split())


def test_nfacets_pol_xx_yy():
    # Request two polarizations on approximately 25 facets
    s = f"{tcf.WSCLEAN} -name {name('nfacets-XX_YY')} -pol XX,YY \
        -facet-regions {tcf.FACETFILE_NFACETS} {tcf.RECTDIMS} {os.environ['MWA_MS']}"
    check_call(s.split())


@pytest.mark.parametrize("npol", (2, 4))
def test_facet_h5solution(npol):
    # Test facet-based imaging and applying h5 solutions
    # where the polarization axis in the h5 file has size npol
    h5download = (
        f"wget -N -q www.astron.nl/citt/ci_data/wsclean/mock_soltab_{npol}pol.h5"
    )
    check_call(h5download.split())

    name = f"facet-h5-{npol}pol"
    s = f"{tcf.WSCLEAN} -use-wgridder -name {name} -apply-facet-solutions mock_soltab_{npol}pol.h5 ampl000,phase000 -pol xx,yy -facet-regions {tcf.FACETFILE_4FACETS} {tcf.DIMS} -join-polarizations -interval 10 14 -niter 1000000 -auto-threshold 5 -mgain 0.8 {os.environ['MWA_MS']}"
    check_call(s.split())

    # Check for output images
    assert os.path.isfile(f"{name}-psf.fits")
    for pol in ["XX", "YY"]:
        trunk = name + "-" + pol
        for image_type in [
            "image",
            "image-pb",
            "dirty",
            "model",
            "model-pb",
            "residual",
            "residual-pb",
        ]:
            image_name = trunk + "-" + image_type + ".fits"
            assert os.path.isfile(image_name)


def test_facet_beam():
    # Test facet beam, using 4 polarizations
    if os.environ["MWA_COEFFS_PATH"]:
        s = f"{tcf.WSCLEAN} -name {name('nfacets-iquv-facet-beam')} -interval 10 14 -apply-facet-beam -pol iquv \
            -facet-regions {tcf.FACETFILE_NFACETS} {tcf.RECTDIMS} \
                -mwa-path {os.environ['MWA_COEFFS_PATH']} {os.environ['MWA_MS']}"
        check_call(s.split())
    else:
        warnings.warn(
            "MWA_PATH environment variable not set, tests that require beam corrections will be skipped."
        )


# Test wsclean-mp command
def test_mpi_join_channels():
    s = f"mpirun {tcf.WSCLEAN_MP} -name {name('mpi-join')} {tcf.RECTDIMS} -scale 1amin -channels-out 2 -join-channels -niter 1000000 -mgain 0.8 -auto-threshold 5 -multiscale -no-update-model-required {os.environ['MWA_MS']}"
    check_call(s.split())


def test_mpi_split_channels():
    s = f"mpirun {tcf.WSCLEAN_MP} -name {name('mpi-split')} {tcf.RECTDIMS} -scale 1amin -channels-out 2 -niter 1000000 -mgain 0.8 -auto-threshold 5 -multiscale -no-update-model-required {os.environ['MWA_MS']}"
    check_call(s.split())


def test_idg_with_reuse_psf():
    # Test for issue #81: -reuse-psf gives segmentation fault in IDG
    # First make sure input files exist:
    s = f"{tcf.WSCLEAN} -name {name('idg-reuse-psf-A')} {tcf.DIMS} -use-idg -idg-mode cpu -grid-with-beam -interval 10 14 -mgain 0.8 -niter 1 -mwa-path {os.environ['MWA_COEFFS_PATH']} {os.environ['MWA_MS']}"
    check_call(s.split())
    os.rename(
        name("idg-reuse-psf-A") + "-model.fits", name("idg-reuse-psf-B") + "-model.fits"
    )
    os.rename(
        name("idg-reuse-psf-A") + "-beam.fits", name("idg-reuse-psf-B") + "-beam.fits"
    )
    # Now continue:
    s = f"{tcf.WSCLEAN} -name {name('idg-reuse-psf-B')} {tcf.DIMS} -use-idg -idg-mode cpu -grid-with-beam -interval 10 14 -mgain 0.8 -niter 1 -continue -reuse-psf {name('idg-reuse-psf-A')} -mwa-path {os.environ['MWA_COEFFS_PATH']} {os.environ['MWA_MS']}"
    check_call(s.split())


def test_masked_parallel_deconvolution():
    # Test for two issues:
    # - issue #96: Source edges in restored image after parallel deconvolution
    # - issue #31: Model images are masked in parallel cleaning
    # The result of this test should be a model image with an unmasked Gaussian and a
    # properly residual. Before the fix, the Gaussian was masked in the model, and
    # therefore only a single pixel was visible, and residual would only be updated
    # on the place of the pixel.

    # First create a mask image with one pixel set:
    s = f"{tcf.WSCLEAN} -name {name('masked-parallel-deconvolution-prepare')} -size 256 256 -scale 1amin -interval 10 14 -niter 1 {os.environ['MWA_MS']}"
    check_call(s.split())
    # Now use this as a mask, and force a Gaussian on the position
    s = f"{tcf.WSCLEAN} -name {name('masked-parallel-deconvolution')} -size 256 256 -scale 1amin -fits-mask {name('masked-parallel-deconvolution-prepare')}-model.fits -interval 10 14 -niter 10 -parallel-deconvolution 128 -multiscale -multiscale-scales 10 {os.environ['MWA_MS']}"
    check_call(s.split())
    for f in glob.glob(f"{name('masked-parallel-deconvolution-prepare')}*.fits"):
        os.remove(f)
