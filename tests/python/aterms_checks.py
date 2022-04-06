import pytest
import sys
import numpy as np
from utils import check_and_remove_files, check_output, validate_call

# Append current directory to system path in order to import testconfig variables
sys.path.append(".")

# Import configuration variables as test configuration (tcf)
import config_vars as tcf

fits_filename = "screens.fits"
config_filename = "screen_config.cfg"


def create_dummy_fits_screen():
    """Create a dummy fits screen file compatible with tcf.MWA_MOCK_MS"""

    try:
        from astropy.io import fits
    except:
        pytest.skip("Could not import astropy, so aterms checks are skipped.")
    # For "diagonal" screen type, the fits file has shape [RA, DEC, MATRIX, ANTENNA, FREQ, TIME].T
    x_imsize = 32
    y_imsize = 32
    n_antennas = 128
    n_times = 1
    n_frequencies = 20
    n_polarizations = 4
    shape_out = [
        n_times,
        n_frequencies,
        n_antennas,
        n_polarizations,
        y_imsize,
        x_imsize,
    ]
    hdu = fits.PrimaryHDU(np.ones(shape_out, dtype=np.float32))

    hdulist = fits.HDUList([hdu])
    header = hdulist[0].header

    reference_ra_deg = 1
    reference_dec_deg = 1
    cellsize_deg = 3

    # Add RA, Dec info
    i = 1
    header[f"CRVAL{i}"] = reference_ra_deg
    header[f"CDELT{i}"] = -cellsize_deg
    header[f"CRPIX{i}"] = x_imsize / 2.0
    header[f"CUNIT{i}"] = "deg"
    header[f"CTYPE{i}"] = "RA---SIN"
    i += 1
    header[f"CRVAL{i}"] = reference_dec_deg
    header[f"CDELT{i}"] = cellsize_deg
    header[f"CRPIX{i}"] = y_imsize / 2.0
    header[f"CUNIT{i}"] = "deg"
    header[f"CTYPE{i}"] = "DEC--SIN"
    i += 1

    # add MATRIX info
    header[f"CRVAL{i}"] = 0.0
    header[f"CDELT{i}"] = 1.0
    header[f"CRPIX{i}"] = 1.0
    header[f"CUNIT{i}"] = ""
    header[f"CTYPE{i}"] = "MATRIX"
    i += 1

    # add ANTENNA info
    header[f"CRVAL{i}"] = 0.0
    header[f"CDELT{i}"] = 1.0
    header[f"CRPIX{i}"] = 1.0
    header[f"CUNIT{i}"] = ""
    header[f"CTYPE{i}"] = "ANTENNA"
    i += 1

    # Add frequency info
    ref_freq = 1e8
    del_freq = 1e8
    header["RESTFRQ"] = ref_freq
    header[f"CRVAL{i}"] = ref_freq
    header[f"CDELT{i}"] = del_freq
    header[f"CRPIX{i}"] = 1.0
    header[f"CUNIT{i}"] = "Hz"
    header[f"CTYPE{i}"] = "FREQ"
    i += 1

    # Add time info
    ref_time = 1.0
    del_time = 1.0
    header[f"CRVAL{i}"] = ref_time
    header[f"CDELT{i}"] = del_time
    header[f"CRPIX{i}"] = 1.0
    header[f"CUNIT{i}"] = "s"
    header[f"CTYPE{i}"] = "TIME"
    i += 1

    # Add equinox
    header["EQUINOX"] = 2000.0

    # Add telescope
    header["TELESCOP"] = "LOFAR"
    hdulist[0].header = header

    hdulist.writeto(fits_filename)


def create_aterm_config():
    """Create the configuration file needed by WSClean to apply the dummy aterms"""
    with open(config_filename, "w") as f:
        f.write("aterms = [diagonal]\r\n")
        f.write(f"diagonal.images = [{fits_filename}]\r\n")


@pytest.mark.usefixtures("prepare_mock_ms")
class TestAterms:
    def test_aterms(self):

        if "IDG" not in check_output([tcf.WSCLEAN, "--version"]).decode():
            pytest.skip("WSClean was not compiled with IDG.")

        if "EveryBeam" not in check_output([tcf.WSCLEAN, "--version"]).decode():
            pytest.skip("WSClean was not compiled with EveryBeam.")

        create_dummy_fits_screen()
        create_aterm_config()

        s = f"{tcf.WSCLEAN} -size 200 200 -scale 2arcsec -niter 1 -mgain 0.7 -log-time -nmiter 1 -threshold 0.001 -use-idg -aterm-config screen_config.cfg {tcf.MWA_MOCK_MS}"
        validate_call(s.split())

        # Remove
        check_and_remove_files([fits_filename, config_filename], remove=True)
