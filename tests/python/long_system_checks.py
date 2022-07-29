import pytest
import os, glob
import sys
from utils import validate_call

# Append current directory to system path in order to import testconfig
sys.path.append(".")

# Import configuration variables as test configuration (tcf)
import config_vars as tcf

"""
Test script containing a collection of wsclean commands, tested on an MWA
measurement set. Tests contained in this file can be invoked via various routes:

- execute "make longsystemcheck"  in your build directory
- execute "[python3 -m] pytest [OPTIONS] source/<test_name.py>" in your build/tests/python directory
"""


def name(name: str):
    return os.path.join(tcf.RESULTS_DIR, name)


@pytest.mark.usefixtures("prepare_large_ms")
class TestLongSystem:
    """
    Collection of long system tests.
    """

    def test_dirty_image(self):
        # Make dirty image
        s = f"{tcf.WSCLEAN} -name {name('test-dirty')} {tcf.DIMS_LARGE} {tcf.MWA_MS}"
        validate_call(s.split())

    def test_clean_rectangular_unpadded_image(self):
        # Clean a rectangular unpadded image
        s = f"{tcf.WSCLEAN} -name {name('clean-rectangular')} -padding 1 -local-rms \
              -auto-threshold 5 -mgain 0.8 -niter 100000 {tcf.DIMS_RECTANGULAR} {tcf.MWA_MS}"
        validate_call(s.split())

    def test_automask_multiscale_clean(self):
        # Auto-masked multi-scale clean
        s = f"{tcf.WSCLEAN} -name {name('multiscale-automasked')} -auto-threshold 0.5 -auto-mask 3 \
              -mgain 0.8 -multiscale -niter 100000 {tcf.DIMS_RECTANGULAR} {tcf.MWA_MS}"
        validate_call(s.split())

    def test_multiple_intervals(self):
        # Multiple intervals
        s = f"{tcf.WSCLEAN} -name {name('intervals')} -intervals-out 3 \
            {tcf.DIMS_RECTANGULAR} {tcf.MWA_MS}"
        validate_call(s.split())

    def test_multiple_intervals_and_channels(self):
        # Multiple intervals + multiple channels with some cleaning
        s = f"{tcf.WSCLEAN} -name {name('intervals-and-channels')} -intervals-out 3 \
            -channels-out 2 -niter 1000 -mgain 0.8 {tcf.DIMS_LARGE} {tcf.MWA_MS}"
        validate_call(s.split())

    def test_multiple_intervals_and_facets(self):
        # Multiple intervals + multiple facets with some cleaning
        s_base = f"{tcf.WSCLEAN} -name {name('intervals-and-facets')} -intervals-out 3 \
            -facet-regions {tcf.FACETFILE_4FACETS}"
        s = f"{s_base} -niter 1000 -mgain 0.8 {tcf.DIMS_LARGE} {tcf.MWA_MS}"
        validate_call(s.split())

        # Run predict, using the model generated above.
        s = f"{s_base} -predict {tcf.MWA_MS}"
        validate_call(s.split())

    def test_multifrequency_hogbom(self):
        # Multi-frequency Högbom clean, no parallel gridding
        s = f"{tcf.WSCLEAN} -name {name('mfhogbom')} -channels-out 4 -join-channels -auto-threshold 3 \
            -mgain 0.8 -niter 1000000 {tcf.DIMS_RECTANGULAR} {tcf.MWA_MS}"
        validate_call(s.split())

    def test_multifrequency_hogbom_spectral_fit(self):
        # Multi-frequency Högbom clean with spectral fitting
        s = f"{tcf.WSCLEAN} -name {name('mfhogbom-fitted')} -channels-out 4 -join-channels -parallel-gridding 4 \
           -fit-spectral-pol 2 -auto-threshold 3 -mgain 0.8 \
               -niter 1000000 {tcf.DIMS_RECTANGULAR} {tcf.MWA_MS}"
        validate_call(s.split())

    def test_mutifrequency_multiscale_parallel(self):
        # Multi-frequency multi-scale clean with spectral fitting, pallel gridding & cleaning
        s = f"{tcf.WSCLEAN} -name {name('mfms-fitted')} -channels-out 4 -join-channels -parallel-gridding 4 \
             -parallel-deconvolution 1000 -fit-spectral-pol 2 -multiscale -auto-threshold 0.5 \
                  -auto-mask 3 -mgain 0.8 -niter 1000000 {tcf.DIMS_RECTANGULAR} {tcf.MWA_MS}"
        validate_call(s.split())

    def test_save_components(self):
        test_name = name("mfms-components")

        # Remove old component list if it exists
        component_file = test_name + "-sources.txt"
        if os.path.exists(component_file):
            os.remove(component_file)

        # Save the list of components
        s = f"{tcf.WSCLEAN} -name {test_name} -save-source-list -channels-out 4 \
            -join-channels -parallel-gridding 4 -fit-spectral-pol 2 \
                -auto-threshold 0.5 -auto-mask 3 -mgain 0.8 -niter 1000000 \
                    -multiscale -parallel-deconvolution 1000 {tcf.DIMS_LARGE} {tcf.MWA_MS}"
        validate_call(s.split())

        # Check whether source files is generated
        assert os.path.isfile(component_file)

    def test_linear_joined_polarizations(self):
        # Linear joined polarizations with 4 joined channels
        s = f"{tcf.WSCLEAN} -name {name('linearpol')} -niter 1000000 -auto-threshold 3.0 \
             -pol XX,YY,XY,YX -join-polarizations -join-channels -mgain 0.85 \
                 -channels-out 4 -parallel-gridding 16 {tcf.DIMS_LARGE} {tcf.MWA_MS}"
        validate_call(s.split())

    def test_two_timesteps(self):
        # Image two timesteps
        s = f"{tcf.WSCLEAN} -name {name('two-timesteps')} -niter 1000000 -auto-threshold 3.0 \
            -intervals-out 2 -interval 20 22 -mgain 0.85 {tcf.DIMS_RECTANGULAR} {tcf.MWA_MS}"
        validate_call(s.split())

    def test_stop_on_negative_components(self):
        # Stop on negative components
        s = f"{tcf.WSCLEAN} -name {name('stop-on-negatives')} -stop-negative -niter 100000 {tcf.DIMS_RECTANGULAR} {tcf.MWA_MS}"
        validate_call(s.split())

    def test_save_imaging_weights(self):
        s = f"{tcf.WSCLEAN} -name {name('store-imaging-weights')} -no-reorder -store-imaging-weights {tcf.DIMS_RECTANGULAR} {tcf.MWA_MS}"
        validate_call(s.split())

    @pytest.mark.parametrize(
        "gridder, test_name", (["", "shift-ws"], ["-use-wgridder", "shift-wg"])
    )
    def test_shift_image(self, gridder, test_name):
        # Shift the image with w-stacking and w-gridder gridder
        s = f"{tcf.WSCLEAN} {gridder} -name {name(test_name)} -mgain 0.8 -auto-threshold 5 -niter 1000000 -make-psf {tcf.DIMS_RECTANGULAR} -shift 08h09m20s -39d06m54s -no-update-model-required {tcf.MWA_MS}"
        validate_call(s.split())

    def test_missing_channels_in_deconvolution(self):
        # The test set has some missing MWA subbands. One MWA subband is 1/24 of the data (32/768 channels), so
        # by imaging with -channels-out 24, it is tested what happens when an output channel has no data.
        s = f"{tcf.WSCLEAN} -name {name('missing-channels-in-deconvolution')} -use-wgridder {tcf.DIMS_LARGE} -baseline-averaging 2.0 -no-update-model-required -niter 150000 -auto-threshold 2.0 -auto-mask 5.0 -mgain 0.9 -channels-out 24 -join-channels -fit-spectral-pol 4 {tcf.MWA_MS}"
        validate_call(s.split())

    def test_grid_with_beam(self):
        """Requires that WSClean is compiled with IDG and EveryBeam
        """
        name = "idg-beam"

        # Remove existing component files if present
        for source_file in ["sources", "sources-pb"]:
            component_file = name + "-" + source_file + ".txt"
            if os.path.exists(component_file):
                os.remove(component_file)

        s = f"{tcf.WSCLEAN} -name {name} -use-idg -grid-with-beam -save-source-list -mgain 0.8 -auto-threshold 5 -niter 1000000 -interval 10 14 {tcf.DIMS_LARGE} -mwa-path . {tcf.MWA_MS}"
        validate_call(s.split())
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

    def test_two_facets(self):
        # Apply the facet to the image
        s = f"{tcf.WSCLEAN} -name {name('two-facets')} -facet-regions {tcf.FACETFILE_2FACETS} \
            {tcf.DIMS_RECTANGULAR} {tcf.MWA_MS}"
        validate_call(s.split())

    def test_nfacets_pol_xx_yy(self):
        # Request two polarizations on approximately 25 facets
        s = f"{tcf.WSCLEAN} -name {name('nfacets-XX_YY')} -pol XX,YY \
            -facet-regions {tcf.FACETFILE_NFACETS} {tcf.DIMS_RECTANGULAR} {tcf.MWA_MS}"
        validate_call(s.split())

    @pytest.mark.parametrize("npol", (2, 4))
    def test_facet_h5solution(self, npol):
        # Test facet-based imaging and applying h5 solutions
        # where the polarization axis in the h5 file has size npol
        h5download = f"wget -N -q {tcf.WSCLEAN_DATA_URL}/mock_soltab_{npol}pol.h5"
        validate_call(h5download.split())

        name = f"facet-h5-{npol}pol"
        s = f"{tcf.WSCLEAN} -use-wgridder -name {name} -apply-facet-solutions mock_soltab_{npol}pol.h5 ampl000,phase000 -pol xx,yy -facet-regions {tcf.FACETFILE_4FACETS} {tcf.DIMS_LARGE} -join-polarizations -interval 10 14 -niter 1000000 -auto-threshold 5 -mgain 0.8 {tcf.MWA_MS}"
        validate_call(s.split())

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

    def test_facet_beam(self):
        # Test facet beam, using 4 polarizations
        s = f"{tcf.WSCLEAN} -name {name('nfacets-iquv-facet-beam')} -interval 10 14 -apply-facet-beam -pol iquv \
            -facet-regions {tcf.FACETFILE_NFACETS} {tcf.DIMS_RECTANGULAR} \
                -mwa-path . {tcf.MWA_MS}"
        validate_call(s.split())

    def test_mpi_join_channels(self):
        # Test wsclean-mp command
        s = f"mpirun {tcf.WSCLEAN_MP} -name {name('mpi-join')} {tcf.DIMS_RECTANGULAR} -channels-out 2 -join-channels -niter 1000000 -mgain 0.8 -auto-threshold 5 -multiscale -no-update-model-required {tcf.MWA_MS}"
        validate_call(s.split())

    def test_mpi_split_channels(self):
        s = f"mpirun {tcf.WSCLEAN_MP} -name {name('mpi-split')} {tcf.DIMS_RECTANGULAR} -channels-out 2 -niter 1000000 -mgain 0.8 -auto-threshold 5 -multiscale -no-update-model-required {tcf.MWA_MS}"
        validate_call(s.split())

    def test_idg_with_reuse_psf(self):
        # Test for issue #81: -reuse-psf gives segmentation fault in IDG
        # First make sure input files exist:
        s = f"{tcf.WSCLEAN} -name {name('idg-reuse-psf-A')} {tcf.DIMS_LARGE} -use-idg -idg-mode cpu -grid-with-beam -interval 10 14 -mgain 0.8 -niter 1 -mwa-path . {tcf.MWA_MS}"
        validate_call(s.split())
        os.rename(
            name("idg-reuse-psf-A") + "-model.fits",
            name("idg-reuse-psf-B") + "-model.fits",
        )
        os.rename(
            name("idg-reuse-psf-A") + "-beam.fits",
            name("idg-reuse-psf-B") + "-beam.fits",
        )
        # Now continue:
        s = f"{tcf.WSCLEAN} -name {name('idg-reuse-psf-B')} {tcf.DIMS_LARGE} -use-idg -idg-mode cpu -grid-with-beam -interval 10 14 -mgain 0.8 -niter 1 -continue -reuse-psf {name('idg-reuse-psf-A')} -mwa-path . {tcf.MWA_MS}"
        validate_call(s.split())

    def test_idg_with_reuse_dirty(self):
        # Test for issue #80: -reuse-dirty option fails (#80)
        # First make sure input files exist:
        s = f"{tcf.WSCLEAN} -name {name('idg-reuse-dirty-A')} {tcf.DIMS_LARGE} -use-idg -idg-mode cpu -grid-with-beam -interval 10 14 -mgain 0.8 -niter 1 -mwa-path . {tcf.MWA_MS}"
        validate_call(s.split())
        os.rename(
            name("idg-reuse-dirty-A") + "-model.fits",
            name("idg-reuse-dirty-B") + "-model.fits",
        )
        # Now continue:
        s = f"{tcf.WSCLEAN} -name {name('idg-reuse-dirty-B')} {tcf.DIMS_LARGE} -use-idg -idg-mode cpu -grid-with-beam -interval 10 14 -mgain 0.8 -niter 1 -continue -reuse-dirty {name('idg-reuse-dirty-A')} -mwa-path . {tcf.MWA_MS}"
        validate_call(s.split())

    def test_idg_with_reuse_psf(self):
        # Test for issue #81: -reuse-psf gives segmentation fault in IDG
        # First make sure input files exist:
        s = f"{tcf.WSCLEAN} -name {name('idg-reuse-psf-A')} {tcf.DIMS_LARGE} -use-idg -idg-mode cpu -grid-with-beam -interval 10 14 -mgain 0.8 -niter 1 -mwa-path . {tcf.MWA_MS}"
        validate_call(s.split())
        # Model image A is copied to B-model-pb corrected image, to avoid
        # issues due to NaN values in the A-model-pb.fits file.
        # As such, this test is purely illlustrative.
        os.rename(
            name("idg-reuse-psf-A") + "-model.fits",
            name("idg-reuse-psf-B") + "-model-pb.fits",
        )
        os.rename(
            name("idg-reuse-psf-A") + "-beam.fits",
            name("idg-reuse-psf-B") + "-beam.fits",
        )
        # Now continue:
        s = f"{tcf.WSCLEAN} -name {name('idg-reuse-psf-B')} {tcf.DIMS_LARGE} -use-idg -idg-mode cpu -grid-with-beam -interval 10 14 -mgain 0.8 -niter 1 -continue -reuse-psf {name('idg-reuse-psf-A')} -mwa-path . {tcf.MWA_MS}"
        validate_call(s.split())

    def test_idg_with_reuse_dirty(self):
        # Test for issue #80: -reuse-dirty option fails (#80)
        # First make sure input files exist:
        s = f"{tcf.WSCLEAN} -name {name('idg-reuse-dirty-A')} {tcf.DIMS_LARGE} -use-idg -idg-mode cpu -grid-with-beam -interval 10 14 -mgain 0.8 -niter 1 -mwa-path . {tcf.MWA_MS}"
        validate_call(s.split())
        # Model image A is copied to B-model-pb corrected image, to avoid
        # issues due to NaN values in the A-model-pb.fits file.
        # As such, this test is purely illlustrative.
        os.rename(
            name("idg-reuse-dirty-A") + "-model.fits",
            name("idg-reuse-dirty-B") + "-model-pb.fits",
        )
        # Now continue:
        s = f"{tcf.WSCLEAN} -name {name('idg-reuse-dirty-B')} {tcf.DIMS_LARGE} -use-idg -idg-mode cpu -grid-with-beam -interval 10 14 -mgain 0.8 -niter 1 -continue -reuse-dirty {name('idg-reuse-dirty-A')} -mwa-path . {tcf.MWA_MS}"
        validate_call(s.split())

    def test_masked_parallel_deconvolution(self):
        # Test for two issues:
        # - issue #96: Source edges in restored image after parallel deconvolution
        # - issue #31: Model images are masked in parallel cleaning
        # The result of this test should be a model image with an unmasked Gaussian and a
        # properly residual. Before the fix, the Gaussian was masked in the model, and
        # therefore only a single pixel was visible, and residual would only be updated
        # on the place of the pixel.

        # First create a mask image with one pixel set:
        s = f"{tcf.WSCLEAN} -name {name('masked-parallel-deconvolution-prepare')} -size 256 256 -scale 1amin -interval 10 14 -niter 1 {tcf.MWA_MS}"
        validate_call(s.split())
        # Now use this as a mask, and force a Gaussian on the position
        s = f"{tcf.WSCLEAN} -name {name('masked-parallel-deconvolution')} -size 256 256 -scale 1amin -fits-mask {name('masked-parallel-deconvolution-prepare')}-model.fits -interval 10 14 -niter 10 -parallel-deconvolution 128 -multiscale -multiscale-scales 10 {tcf.MWA_MS}"
        validate_call(s.split())
        for f in glob.glob(f"{name('masked-parallel-deconvolution-prepare')}*.fits"):
            os.remove(f)

    @pytest.mark.parametrize("use_beam", (False, True))
    def test_idg_predict(self, use_beam):
        # Check whether primary beam corrected image is used in -predict
        # First make sure model images exist
        run_name = name("idg-predict")
        grid_with_beam = "-grid-with-beam" if use_beam else ""
        s0 = f"{tcf.WSCLEAN} -name {run_name} {tcf.DIMS_LARGE} -use-idg -idg-mode cpu {grid_with_beam} -interval 10 12 -mgain 0.8 -niter 1 -mwa-path . {tcf.MWA_MS}"
        validate_call(s0.split())

        # Remove the model image that shouldn't be needed for the predict
        if use_beam:
            # Move model.fits to model-pb.fits file. Formally, the model-pb.fits file
            # should be used directly, but as this file can contain NaN values, a predict run
            # can bail out on these NaN values.
            os.rename(f"{run_name}-model.fits", f"{run_name}-model-pb.fits")

        s1 = f"{tcf.WSCLEAN} -name {run_name} {tcf.DIMS_LARGE} -predict -use-idg -idg-mode cpu {grid_with_beam} -interval 10 12 -mwa-path . {tcf.MWA_MS}"
        validate_call(s1.split())

    def test_catch_invalid_channel_selection(self):
        # Invalid selection: people often forget the second value of -channel-range is an open interval end (i.e. excluded the value itself).
        s = f"{tcf.WSCLEAN} -name {name('test-caught-bad-selection')} -channels-out 256 -channel-range 0 255 {tcf.DIMS_LARGE} {tcf.MWA_MS}"
        with pytest.raises(Exception):
            validate_call(s.split())

    def test_catch_invalid_channel_selection_with_gaps(self):
        s = f"{tcf.WSCLEAN} -name {name('test-caught-bad-selection')} -gap-channel-division -channels-out 256 -channel-range 0 255 {tcf.DIMS_LARGE} {tcf.MWA_MS}"
        with pytest.raises(Exception):
            validate_call(s.split())

    def test_catch_invalid_channel_selection_with_division(self):
        s = f"{tcf.WSCLEAN} -name {name('test-caught-bad-selection')} -channel-division-frequencies 145e6 -channels-out 256 -channel-range 0 255 {tcf.DIMS_LARGE} {tcf.MWA_MS}"
        with pytest.raises(Exception):
            validate_call(s.split())

    def test_multiband_no_mf_weighting(self):
        # Tests issue #105: Segmentation fault (core dumped), when grouping spectral windows + no-mf-weighting Master Branch
        # The issue was caused by invalid indexing into the BandData object.
        s = f"{tcf.WSCLEAN} -name {name('vla-multiband-no-mf')} -size 768 768 -scale 0.05arcsec -pol QU -mgain 0.85 -niter 1000 -auto-threshold 3 -join-polarizations -squared-channel-joining -no-update-model-required -no-mf-weighting {tcf.JVLA_MS}"
        validate_call(s.split())
        for f in glob.glob(f"{name('vla-multiband-no-mf')}*.fits"):
            os.remove(f)

    def test_spectrally_fitted_with_joined_polarizations(self):
        s = f"{tcf.WSCLEAN} -name {name('iv-jointly-fitted')} {tcf.DIMS_LARGE} -parallel-gridding 4 -channels-out 4 -join-channels -fit-spectral-pol 2 -pol i,v -join-polarizations -niter 1000 -auto-threshold 5 -multiscale -mgain 0.8 {tcf.MWA_MS}"
        validate_call(s.split())

    def test_dd_psfs_call(self):
        s = f"{tcf.WSCLEAN} -size 200 200 -scale 2arcsec -make-psf -dd-psf-grid 5 5 {tcf.MWA_MS}"
        validate_call(s.split())
