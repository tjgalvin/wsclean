"""Dummy deconvolution algorithm writing its inputs to disk

Usage:
    wsclean ... -python-deconvolution test_deconvolution_write_input.py ...

This is used by the test_direction_dependent_psfs() test in long_system_checks.py to
inspect dirty and psf images passed on by the parallel deconvolution algorithm
"""


import numpy as np


def deconvolve(residual, model, psf, meta):
    """Dummy deconvolution function to be called from the PythonDeconvolution class in Radler

    Input residual and psf images are written to disk as numpy (.npy) array
    The function keeps a counter to generate unique filenames per application run
    """

    # Only write out files for the first call per image
    # The first call is peak finding only, max_iterations is zero then
    if meta.max_iterations == 0:
        # deconvolve.counter is an attribute of the deconvolve() function (this function).
        # It keeps track of many times the function is called with meta.max_iteration = 0
        deconvolve.counter = getattr(deconvolve, "counter", 0)
        count = deconvolve.counter
        deconvolve.counter += 1

        # Write out dirty image
        filename_dirty = (
            f"test_results/test-deconvolution-write-input-dirty-{count}.npy"
        )
        np.save(filename_dirty, residual)

        # Write out psf
        filename_psf = (
            f"test_results/test-deconvolution-write-input-psf-{count}.npy"
        )
        np.save(filename_psf, psf)

    # Fill a dictionary with values that wsclean expects:
    result = dict()
    result["residual"] = residual
    result["model"] = model
    result["level"] = 0.0
    result["continue"] = False

    return result
