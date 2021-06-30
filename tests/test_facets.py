import os
import pytest
from subprocess import check_call, check_output
import shutil

# Append current directory to system path in order to import testconfig
import sys

sys.path.append(".")

import testconfig as tcf

# Prepend path with current working directory to make sure
# wsclean executable from the build directory
os.environ["PATH"] = f"{os.getcwd()}:{os.environ['PATH']}"

MWA_MOCK_ARCHIVE = "MWA_ARCHIVE.tar.bz2"
MWA_MOCK_MS = "MWA_MOCK.ms"
MWA_MOCK_FULL = "MWA_MOCK_FULL.ms"
MWA_MOCK_FACET = "MWA_MOCK_FACET.ms"
MODEL_IMAGE = "point-source-model.fits"


@pytest.fixture(autouse=True)
def prepare():
    # Change to directory containing the data
    os.makedirs(tcf.WORKDIR, exist_ok=True)
    os.chdir(tcf.WORKDIR)

    if not os.path.isfile(MODEL_IMAGE):
        wget = f"wget -q http://www.astron.nl/citt/ci_data/wsclean/{MODEL_IMAGE}"
        check_call(wget.split())

    if not os.path.isfile(MWA_MOCK_ARCHIVE):
        wget = f"wget -q www.astron.nl/citt/EveryBeam/MWA-single-timeslot.tar.bz2 -O {MWA_MOCK_ARCHIVE}"
        check_call(wget.split())

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


def assert_taql(command, expected_rows=0):
    assert shutil.which("taql") is not None, "taql executable not found!"
    result = check_output(["taql", "-noph", command]).decode().strip()
    assert result == f"select result of {expected_rows} rows", result


def predict_full_image(ms, gridder):
    # Predict full image
    s = f"wsclean -predict {gridder} -name point-source {ms}"
    check_call(s.split())


def predict_facet_image(ms, gridder):
    # Predict facet based image
    s = f"wsclean -predict {gridder} -facet-regions {tcf.FACETFILE_4FACETS} -name point-source {ms}"
    check_call(s.split())


def deconvolve_facets(
    ms,
    gridder,
    reorder,
    nthreads=4,
    image="-size 256 256 -scale 4amin",
    major_cycle="-niter 1000000 -auto-threshold 5 -mgain 0.8",
):
    reorder_ms = "-reorder" if reorder else "-no-reorder"
    s = f"wsclean -parallel-gridding {nthreads} {gridder} {reorder_ms} {major_cycle} -facet-regions {tcf.FACETFILE_4FACETS} -name facet-imaging{reorder_ms} {image} {ms}"
    check_call(s.split())


# FIXME: we should test wstacking and wgridder here too
# but the fail on the taql assertion
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
def test_facetdeconvolution(gridder, reorder):
    """
    Test facet-based deconvolution

    Parameters
    ----------
    gridder : str
        wsclean compatible description of gridder to be used.
    reorder : bool
        Reorder MS?
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

    deconvolve_facets(MWA_MOCK_FACET, gridder, reorder)

    taql_command = f"select from {MWA_MOCK_FACET} t1, {MWA_MOCK_FACET} t2 where not all(near(t1.DATA,t2.MODEL_DATA, 4e-3))"
    assert_taql(taql_command)
