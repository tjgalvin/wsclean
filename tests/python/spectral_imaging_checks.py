import pytest
import os
import sys
from utils import validate_call
from astropy.io import fits

# Append current directory to system path in order to import testconfig
sys.path.append(".")

# Import configuration variables as test configuration (tcf)
import config_vars as tcf


def name(name: str):
    return os.path.join(tcf.RESULTS_DIR, name)


@pytest.mark.usefixtures("prepare_large_ms")
class TestSpectralImaging:
    def test_spectral_mask(self):
        width = 512
        height = 512

        source0_pos_x = 100
        source0_pos_y = 200
        source0_flux = 1

        source1_pos_x = 350
        source1_pos_y = 450
        source1_flux = 2

        # Create two "template" model images
        s = f"{tcf.WSCLEAN} -name {name('spectral-mask')} -channel-range 10 12 -interval 10 20 -channels-out 2 -niter 1 -size {width} {height} -scale 1amin {tcf.MWA_MS}"
        validate_call(s.split())

        def set_mask(filename, x, y, flux):
            with fits.open(filename) as hdu:
                hdu[0].data[:] = 0.0
                hdu[0].data[0, 0, y, x] = flux
                hdu[0].writeto(filename, overwrite=True)

        # Change the two model images for two different channels, where each channel has a source at a different position.
        set_mask(
            name("spectral-mask-0000-model.fits"),
            source0_pos_x,
            source0_pos_y,
            source0_flux,
        )
        set_mask(
            name("spectral-mask-0001-model.fits"),
            source1_pos_x,
            source1_pos_y,
            source1_flux,
        )

        # Predict visibilities using wsclean into the MODEL_DATA column
        s = f"{tcf.WSCLEAN} -name {name('spectral-mask')} -channel-range 10 12 -interval 10 20 -channels-out 2 -predict {tcf.MWA_MS}"
        validate_call(s.split())

        # Combine the two fits images into one to make a "mask cube"
        with fits.open(name("spectral-mask-0000-model.fits")) as hdu:
            import numpy

            # Dimensions of a fits cube are: STOKES, FREQ, Y, X
            hdu[0].data = numpy.zeros([1, 2, height, width])
            hdu[0].data[0, 0, source0_pos_y, source0_pos_x] = 1.0
            hdu[0].data[0, 1, source1_pos_y, source1_pos_x] = 1.0
            hdu[0].writeto(name("spectral-mask-cube.fits"), overwrite=True)

        # Perform the cleaning with the spectral mask
        s = f"{tcf.WSCLEAN} -v -name {name('spectral-mask')} -data-column MODEL_DATA -fits-mask {name('spectral-mask-cube.fits')} -channel-range 10 12 -interval 10 20 -channels-out 2 -niter 1000 -mgain 0.8 -size {width} {height} -scale 1amin {tcf.MWA_MS}"
        validate_call(s.split())

        # Validate the output
        def validate(filename, x, y, flux):
            with fits.open(filename) as hdu:
                # Validate flux of resulting source within 3% of original flux
                assert (
                    hdu[0].data[0, 0, y, x] > flux * 0.97
                    and hdu[0].data[0, 0, y, x] < flux * 1.03
                )
                # Check if other values are approximately zero
                hdu[0].data[0, 0, y, x] = 0.0
                assert (numpy.abs(hdu[0].data[:]) < 1e-6).all()

        validate(
            name("spectral-mask-0000-model.fits"),
            source0_pos_x,
            source0_pos_y,
            source0_flux,
        )
        validate(
            name("spectral-mask-0001-model.fits"),
            source1_pos_x,
            source1_pos_y,
            source1_flux,
        )
