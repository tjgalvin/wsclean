Direction-dependent Point Spread Functions (PSFs)
======================

With large fields of view, the approximation of considering the dirty image as a convolution of the sky and the PSF yields a higher error.
To mitigate this error, WSClean can calculate multiple PSFs in different directions on a grid, using the parameter ``-psf-grid-size``: the user can specify the number of psfs which should be fit horizontally and vertically in the image.

.. warning::

    This feature is still under developement as of June 2022. If triggered, WSClean will fall back to the existing single PSF implementation.

