import os
import shutil
import sys
from subprocess import check_call

import pytest

# Append current directory to system path in order to import testconfig
sys.path.append(".")

# Import configuration variables as test configuration (tcf)
import config_vars as tcf


@pytest.fixture(scope="session", autouse=True)
def prepare_workdir():
    # Changes to directory containing the data.
    # Optionally removes the entire working directory
    # at the end of a test session. Use with care!
    cwd = os.getcwd()
    os.makedirs(tcf.WORKING_DIR, exist_ok=True)
    os.chdir(tcf.WORKING_DIR)
    os.makedirs(tcf.RESULTS_DIR, exist_ok=True)

    # Download and untar beam pattern file
    if not os.path.isfile(tcf.MWA_COEFF_FILE):
        check_call(
            [
                "wget",
                "-q",
                os.path.join(tcf.EVERYBEAM_DATA_URL, tcf.MWA_COEFF_ARCHIVE),
            ]
        )
        check_call(["tar", "xf", tcf.MWA_COEFF_ARCHIVE])

    yield

    if (
        "CLEANUP_WSCLEAN_TESTS" in os.environ
        and int(os.environ["CLEANUP_WSCLEAN_TESTS"]) == 1
    ):
        os.chdir(cwd)
        shutil.rmtree(tcf.WORKING_DIR)


@pytest.fixture(scope="class")
def prepare_large_ms():
    if not os.path.isfile(f"{tcf.MWA_MS}/table.f1"):
        check_call(
            [
                "wget",
                "-q",
                os.path.join(tcf.IDG_DATA_URL, f"{tcf.MWA_MS}.tgz"),
            ]
        )
        check_call(["tar", "-xf", f"{tcf.MWA_MS}.tgz"])
        os.remove(f"{tcf.MWA_MS}.tgz")
    else:
        pass


@pytest.fixture(scope="class")
def prepare_mock_ms():
    os.makedirs(tcf.MWA_MOCK_MS, exist_ok=True)

    if not os.path.isfile(tcf.MWA_MOCK_ARCHIVE):
        url = os.path.join(
            tcf.EVERYBEAM_DATA_URL, "MWA-single-timeslot.tar.bz2"
        )
        check_call(f"wget -q {url} -O {tcf.MWA_MOCK_ARCHIVE}".split())
        check_call(
            f"tar -xf {tcf.MWA_MOCK_ARCHIVE} -C {tcf.MWA_MOCK_MS} --strip-components=1".split()
        )

    # From python 3.8 onwards, use copytree(..., dirs_exist_ok=True)
    if not os.path.isdir(tcf.MWA_MOCK_FULL):
        shutil.copytree(tcf.MWA_MOCK_MS, tcf.MWA_MOCK_FULL)

    if not os.path.isdir(tcf.MWA_MOCK_FACET):
        shutil.copytree(tcf.MWA_MOCK_MS, tcf.MWA_MOCK_FACET)


@pytest.fixture(scope="function")
def tmp_mwa_mock_facet(tmp_path):
    mwa_mock_facet_path = tmp_path / tcf.MWA_MOCK_FACET
    shutil.copytree(tcf.MWA_MOCK_FACET, mwa_mock_facet_path)
    return mwa_mock_facet_path


@pytest.fixture(scope="class")
def prepare_model_image():
    if not os.path.isfile(tcf.MODEL_IMAGE):
        wget = f"wget -q {os.path.join(tcf.WSCLEAN_DATA_URL, tcf.MODEL_IMAGE)}"
        check_call(wget.split())


@pytest.fixture(scope="class")
def prepare_mock_soltab():
    if not os.path.isfile(tcf.MOCK_SOLTAB_2POL):
        wget = f"wget -q {os.path.join(tcf.WSCLEAN_DATA_URL, tcf.MOCK_SOLTAB_2POL)}"
        check_call(wget.split())
