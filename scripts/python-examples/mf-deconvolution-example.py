#! /usr/bin/python

# This file demonstrates a deconvolution approach similar to simple-deconvolution-example.py,
# but supports multiple frequencies. It thereby also demonstrates fitting the components
# to a requested spectrum.

# run this e.g. with:
# wsclean \
#  -python-deconvolution mf-deconvolution-example.py \
#  -fit-spectral-pol 2 -channels-out 8 -parallel-gridding 8 -join-channels \
#  -niter 1000 -save-first-residual -auto-threshold 3 -mgain 0.8 \
#  -interval 10 11 -size 1024 1024 -scale 1amin \
#  1052736496-averaged.ms/

import numpy


def deconvolve(residual, model, psf, meta):
    nchan, npol, height, width = residual.shape
    print(
        "Python deconvolve() function was called for "
        + f"{width} x {height} x {npol} (npol) x {nchan} (chan) dataset"
    )

    # residual and model are numpy arrays with dimensions nchan x npol x height x width
    # psf is a numpy array with dimensions nchan x height x width.

    # This file demonstrates a very simple deconvolution strategy, which doesn't
    # support multiple polarizations:
    if npol != 1:
        raise NotImplementedError("npol must be one")
    residual = residual[:, 0, :, :]
    # If there are channels missing (flagged), they will be set to NaN
    # They're here set to zero to make it easy to calculate the integrated image
    residual[numpy.isnan(residual)] = 0
    model = model[:, 0, :, :]

    # find the largest peak in the integrated image
    integrated_residual = numpy.sum(residual, axis=0)
    peak_index = numpy.unravel_index(
        numpy.argmax(integrated_residual), integrated_residual.shape
    )
    peak_value = integrated_residual[peak_index]

    mgain_threshold = peak_value * (1.0 - meta.mgain)
    first_threshold = max(
        meta.major_iter_threshold, meta.final_threshold, mgain_threshold
    )

    print(
        f"Starting iteration {meta.iteration_number}, peak={peak_value}, first threshold={first_threshold}"
    )
    while (
        peak_value > first_threshold
        and meta.iteration_number < meta.max_iterations
    ):
        y = peak_index[0]
        x = peak_index[1]
        spectrum = residual[:, y, x]
        spectrum = meta.spectral_fitter.fit_and_evaluate(spectrum, y, x)
        model[:, y, x] += spectrum

        psf_shift = (y + height // 2, x + width // 2)
        residual = (
            residual
            - numpy.roll(psf, psf_shift, axis=(1, 2))
            * spectrum[:, numpy.newaxis, numpy.newaxis]
        )

        integrated_residual = numpy.sum(residual, axis=0)
        peak_index = numpy.unravel_index(
            numpy.argmax(integrated_residual), integrated_residual.shape
        )
        peak_value = integrated_residual[peak_index]

        meta.iteration_number = meta.iteration_number + 1

    print(
        f"Stopped after iteration {meta.iteration_number}, peak={peak_value}"
    )

    # Fill a dictionary with values that wsclean expects:
    result = dict()
    result["residual"] = residual
    result["model"] = model
    result["level"] = peak_value
    result["continue"] = (
        peak_value > meta.final_threshold
        and meta.iteration_number < meta.max_iterations
    )

    print("Finished deconvolve()")
    return result
