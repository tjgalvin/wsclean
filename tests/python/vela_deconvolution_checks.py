import itertools
import os
import sys
from subprocess import check_output

import pytest
from utils import check_and_remove_files, compute_rms, validate_call

# Append current directory to system path in order to import testconfig variables
sys.path.append(".")

# Import configuration variables as test configuration (tcf)
import config_vars as tcf


@pytest.mark.usefixtures("prepare_mock_ms")
class TestVelaDeconvolution:
    @pytest.mark.parametrize("gridder", ["wstacking", "wgridder", "facet-idg"])
    def test_veladeconvolution(self, gridder):
        if (
            gridder == "wgridder"
            and "WGridder"
            not in check_output([tcf.WSCLEAN, "--version"]).decode()
        ):
            pytest.skip("WSClean was not compiled with WGridder.")

        # The parallel-gridding option is not compatible with IDG.
        parallel_gridding_option = (
            "" if "idg" in gridder else "-parallel-gridding 2"
        )

        nchannels = 8
        npixels = 1024
        name = "mwa_veladeconvolution_" + gridder
        s = f"{tcf.WSCLEAN} -quiet -gridder {gridder} -size {npixels} {npixels} -scale 1amin {parallel_gridding_option} -multiscale -parallel-deconvolution 512 -niter 1000000 -mgain 0.8 -channels-out {nchannels} -join-channels -deconvolution-channels 3 -fit-spectral-pol 2 -auto-threshold 1 -auto-mask 4 -name {name} {tcf.MWA_MOCK_MS}"
        validate_call(s.split())
        imagenames = ["dirty", "image", "model", "psf", "residual"]
        fpaths = [
            f"{name}-{chan:04d}-{image}.fits"
            for (chan, image) in itertools.product(
                list(range(nchannels)), imagenames
            )
        ]
        fpaths += [f"{name}-MFS-{image}.fits" for image in imagenames]
        check_and_remove_files(fpaths, remove=False)

        rms_dirty = compute_rms(f"{name}-MFS-dirty.fits")
        rms_residual = compute_rms(f"{name}-MFS-residual.fits")
        assert rms_residual < 0.15
        # RMS of residual should be considerably better than RMS of dirty image
        assert 2 * rms_residual < rms_dirty

        # Remove
        [os.remove(fpath) for fpath in fpaths]

    def test_vela_iuwt(self):
        npixels = 1024
        name = "mwa_vela_iuwt"
        s = f"{tcf.WSCLEAN} -quiet -size {npixels} {npixels} -scale 1amin -iuwt -niter 100 -gain 0.2 -mgain 0.8 -name {name} {tcf.MWA_MOCK_MS}"
        validate_call(s.split())
        imagenames = ["dirty", "image", "model", "psf", "residual"]
        fpaths = [f"{name}-{image}.fits" for image in imagenames]
        check_and_remove_files(fpaths, remove=False)

        rms_residual = compute_rms(f"{name}-residual.fits")
        assert rms_residual < 0.28

        # Remove
        [os.remove(fpath) for fpath in fpaths]
