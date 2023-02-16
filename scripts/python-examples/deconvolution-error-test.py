#! /usr/bin/python


def deconvolve(residual, model, psf, meta):
    nchan, npol, height, width = residual.shape
    print(
        "Python deconvolve() function was called for "
        + f"{width} x {height} x {npol} (npol) x {nchan} (chan) dataset"
    )

    raise RuntimeError(
        "This is a test to see if WSClean handles a raise correctly"
    )
