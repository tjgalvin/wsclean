Direction-dependent Point Spread Functions (PSFs)
======================

With large fields of view, the approximation of considering the dirty image as a convolution of the sky and the PSF yields a higher error.
To mitigate this error, WSClean can calculate multiple PSFs in different directions on a grid, using the parameter ``-dd-psf-grid``: the user can specify the number of psfs which should be fit horizontally and vertically in the image.
The ``-dd-psf-grid`` argument only works together with the ``-make-psf-only`` or ``-parallel-deconvolution`` arguments:

- If ``-parallel-deconvolution`` is used, and the number of cells in the PSF grid exceeds the number of parallel sub-images in one or more dimensions, WSClean prints a warning and reduces the PSF grid size.
- If ``-make-psf-only`` is used without ``-parallel-deconvolution``, WSClean uses the PSF grid as supplied by ``-dd-psf-grid``.
- If both ``-make-psf-only`` and ``-parallel-deconvolution`` are not used, WSClean resets the PSF grid to 1x1 and thus disables using multiple PSFs, too.

.. warning::

    This feature is still under developement as of June 2022. If triggered, WSClean will fall back to the existing single PSF implementation.
