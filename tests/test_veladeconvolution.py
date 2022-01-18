from subprocess import check_call, check_output
import pytest
import os
import itertools
import sys
import warnings

# Append current directory to system path in order to import testconfig
sys.path.append(".")

import testconfig as tcf

MWA_MOCK_ARCHIVE = "MWA_ARCHIVE.tar.bz2"
MWA_MOCK_MS = "MWA_MOCK.ms"

def validate_call(cmdline):
    try:
      check_call(cmdline)
    except:
      # To avoid having to work back what the command was, the command is reported:
      raise RuntimeError('Command failed. To run manually, cd to ' + tcf.WORKDIR + ' and execute: ' + ' '.join(cmdline))
      
def gridders():
    return {"wstacking": "", "wgridder": "-use-wgridder"}


@pytest.fixture(autouse=True)
def prepare():
    # Change to directory containing the data
    os.makedirs(tcf.WORKDIR, exist_ok=True)
    os.chdir(tcf.WORKDIR)

    if not os.path.isfile(MWA_MOCK_ARCHIVE):
        wget = f"wget -q www.astron.nl/citt/EveryBeam/MWA-single-timeslot.tar.bz2 -O {MWA_MOCK_ARCHIVE}"
        check_call(wget.split())

    os.makedirs(MWA_MOCK_MS, exist_ok=True)
    check_call(
        f"tar -xf {MWA_MOCK_ARCHIVE}  -C {MWA_MOCK_MS} --strip-components=1".split()
    )


# TODO: should be moved to a common utils script, as it's being used by other tests too
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


@pytest.mark.parametrize("gridder", gridders().items())
def test_veladeconvolution(gridder):
    if "WGridder" not in check_output([tcf.WSCLEAN, "--version"]).decode():
        pytest.skip("WSClean was not compiled with WGridder.")

    nchannels = 8
    npixels = 1024
    name = "mwa_veladeconvolution_" + gridder[0]
    s = f"{tcf.WSCLEAN} -quiet {gridder[1]} -size {npixels} {npixels} -scale 1amin -parallel-gridding 2 -multiscale -parallel-deconvolution 512 -niter 1000000 -mgain 0.8 -channels-out {nchannels} -join-channels -deconvolution-channels 3 -fit-spectral-pol 2 -auto-threshold 1 -auto-mask 4 -name {name} {MWA_MOCK_MS}"
    validate_call(s.split())
    imagenames = ["dirty", "image", "model", "psf", "residual"]
    fpaths = [
        f"{name}-{chan:04d}-{image}.fits"
        for (chan, image) in itertools.product(list(range(nchannels)), imagenames)
    ]
    fpaths += [f"{name}-MFS-{image}.fits" for image in imagenames]
    check_and_remove_files(fpaths, remove=False)

    try:
        from astropy.io import fits
    except:
        warnings.warn(
            UserWarning(
                "Could not import astropy, so numerical checks in test_veladeconvolution are skipped."
            )
        )
        return

    import numpy as np

    hdul_dirty = fits.open(f"{name}-MFS-dirty.fits")
    hdul_residual = fits.open(f"{name}-MFS-residual.fits")
    data_dirty = hdul_dirty[0].data[0, 0, ...]
    data_residual = hdul_residual[0].data[0, 0, ...]
    assert data_residual.shape == data_dirty.shape == (npixels, npixels)

    rms_residual = np.sqrt(np.mean(data_residual ** 2))
    rms_dirty = np.sqrt(np.mean(data_dirty ** 2))

    assert rms_residual < 0.15
    # RMS of residual should be considerably better than RMS of dirty image
    assert 2 * rms_residual < rms_dirty

    # Remove
    [os.remove(fpath) for fpath in fpaths]

def test_vela_iuwt():
    npixels = 1024
    name = "mwa_vela_iuwt"
    s = f"{tcf.WSCLEAN} -quiet -size {npixels} {npixels} -scale 1amin -iuwt -niter 100 -gain 0.2 -mgain 0.8 -name {name} {MWA_MOCK_MS}"
    validate_call(s.split())
    imagenames = ["dirty", "image", "model", "psf", "residual"]
    fpaths = [
        f"{name}-{image}.fits"
        for image in imagenames
    ]
    check_and_remove_files(fpaths, remove=False)

    try:
        from astropy.io import fits
    except:
        warnings.warn(
            UserWarning(
                "Could not import astropy, so numerical checks in test_vela_iuwt are skipped."
            )
        )
        return

    import numpy as np

    hdul_residual = fits.open(f"{name}-residual.fits")
    data_residual = hdul_residual[0].data[0, 0, ...]
    assert data_residual.shape == (npixels, npixels)

    rms_residual = np.sqrt(np.mean(data_residual ** 2))
    assert rms_residual < 0.28

    # Remove
    [os.remove(fpath) for fpath in fpaths]
