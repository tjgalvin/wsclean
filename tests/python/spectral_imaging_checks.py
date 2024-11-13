import os
import sys

import numpy
import pytest
from astropy.io import fits
from utils import validate_call

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


@pytest.mark.usefixtures("prepare_large_ms")
class TestForcedSpectrum:
    def test_forced_spectrum(self):
        width = 512
        height = 512

        # Create a quick template image
        s = f"{tcf.WSCLEAN} -name {name('forced-spectrum-terms')} -channel-range 196 204 -interval 10 20 -size {width} {height} -scale 1amin {tcf.MWA_MS}"
        validate_call(s.split())
        # The PSF file is used as template
        spectral_terms_file = name("forced-spectrum-terms-dirty.fits")

        # Create a cube that contains a (constant) spectral index and curvature value to which the deconvolution is forced.
        with fits.open(spectral_terms_file) as hdu:
            hdu[0].data = numpy.ndarray([1, 2, height, width])
            hdu[0].data[:, 0, :, :] = -0.85
            hdu[0].data[:, 1, :, :] = 0.5
            hdu[0].writeto(spectral_terms_file, overwrite=True)

        s = f"{tcf.WSCLEAN} -name {name('forced-spectrum')} -force-spectrum {spectral_terms_file} -fit-spectral-log-pol 3 -join-channels -save-source-list -interval 10 20 -channels-out 16 -parallel-gridding 8 -niter 1000000 -auto-threshold 3 -mgain 0.8 -size {width} {height} -scale 1amin {tcf.MWA_MS}"
        validate_call(s.split())

        # Validate the output
        with open(name("forced-spectrum-sources.txt")) as source_file:
            lines = source_file.readlines()
            assert len(lines) > 100
            for line in lines[1:]:
                terms_start = line.split(",[")[1]
                assert len(terms_start) > 0
                terms_str = terms_start.split("],")[0]
                assert len(terms_str) > 0
                terms = terms_str.split(",")
                assert len(terms) == 2
                term0 = float(terms[0])
                assert -0.86 < term0 < -0.84
                term1 = float(terms[1])
                assert 0.49 < term1 < 0.51
