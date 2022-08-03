Direction-dependent Point Spread Functions (PSFs)
======================

With large fields of view, the approximation of considering the dirty image as a convolution of the sky and the PSF yields a higher error.
To mitigate this error, WSClean can calculate multiple PSFs in different directions on a grid, using the parameter ``-psf-grid-size``: the user can specify the number of psfs which should be fit horizontally and vertically in the image.
The ``-psf-grid-size`` argument only works together with the ``-parallel-deconvolution`` argument. If the number of cells in the PSF grid exceeds the number of parallel sub-images in one or more dimensions, WSClean prints a warning and reduces the PSF grid size. When parallel deconvolution is disabled, WSClean resets to PSF grid to 1x1, which disables using multiple PSFs, too.

.. warning::

    This feature is still under developement as of June 2022. If triggered, WSClean will fall back to the existing single PSF implementation.
