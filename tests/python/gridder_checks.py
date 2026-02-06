import shutil
import sys

import pytest
from utils import check_output, compare_rms_fits, validate_call

# Append current directory to system path in order to import testconfig
sys.path.append(".")

# Import configuration variables as test configuration (tcf)
import config_vars as tcf


class TestGridders:
    def test_basic(self):
        if "W-Towers" not in check_output([tcf.WSCLEAN, "--version"]).decode():
            pytest.skip("WSClean was not compiled with W-Towers.")
        wgridder = (
            f"{tcf.WSCLEAN} -name wgridder_basic -gridder wgridder "
            "-no-update-model-required "
            "-mgain 0.95 -nmiter 1 -niter 6000 "
            f"-size 2000 2000 -scale 1asec {tcf.LOFAR_3C196_MS}"
        )
        validate_call(wgridder.split())
        wtowers = (
            f"{tcf.WSCLEAN} -name wtower_basic -gridder wtowers "
            "-no-update-model-required "
            "-mgain 0.95 -nmiter 1 -niter 6000 "
            f"-size 2000 2000 -scale 1asec {tcf.LOFAR_3C196_MS}"
        )
        validate_call(wtowers.split())
        compare_rms_fits(
            "wtower_basic-image.fits",
            "wgridder_basic-image.fits",
            threshold=1.0e-3,
        )

    def test_basic_facet_idg(self):
        if "IDG" not in check_output([tcf.WSCLEAN, "--version"]).decode():
            pytest.skip("WSClean was not compiled with IDG.")
        wgridder = (
            f"{tcf.WSCLEAN} -name wgridder_basic -gridder wgridder "
            "-no-update-model-required "
            "-mgain 0.95 -nmiter 1 -niter 6000 "
            f"-size 2000 2000 -scale 1asec {tcf.LOFAR_3C196_MS}"
        )
        validate_call(wgridder.split())
        facet_idg = (
            f"{tcf.WSCLEAN} -name facet_idg_basic -gridder facet-idg "
            "-no-update-model-required "
            "-mgain 0.95 -nmiter 1 -niter 6000 "
            f"-size 2000 2000 -scale 1asec {tcf.LOFAR_3C196_MS}"
        )
        validate_call(facet_idg.split())
        compare_rms_fits(
            "facet_idg_basic-image.fits",
            "wgridder_basic-image.fits",
            threshold=1.0e-3,
        )
