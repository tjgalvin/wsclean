Basic cleaning
==============

Image dimensions
----------------

The image size is set in pixels with the ``-size`` parameter, which takes a width and a height. An image is not required to be square shaped. It is however required to use even numbers for width and height. The ``-scale`` parameter takes an angle and sets the angular size of a single pixel. The pixel scale currently has to be square shaped. Together, the ``-size`` and ``-scale`` parameters set the angular size of the image. 

Threshold and maximum number of iterations
------------------------------------------

The basic required parameters for cleaning are ``threshold`` and ``niter``. Threshold defines where to stop cleaning: WSClean will continue cleaning until the peak residual flux is below the given threshold. It is given in Jy. The ``niter`` parameter sets the maximum number of minor iterations that are allowed to be used.

It is good practice to make sure cleaning has reached the threshold, and only use niter to make sure wsclean will not run for an excessively long time. One should also not clean deeper than the noise, unless a mask is used. A typical stopping criterium is 3 x stddev, where stddev is the standard deviation of the noise in the image. The easiest way of setting a stopping criterium based on the noise, is by using an automatic threshold. This is discussed next.

Automatic threshold
-------------------

To set the deconvolution threshold, one would need to know the level of the noise before running the clean. For automated processing, this is undesirable. Therefore, WSClean can automatically set a threshold based on the residual noise level. The option for this is ``-auto-threshold``. With this option, WSClean will calculate the standard deviation of the residual image before the start of every major deconvolution iteration, and clean up to the given factor times the found noise standard deviation. The standard deviation is calculated using the medium absolute deviation, which is a robust estimator that is not very sensitive to source structure still present in the image. When performing :doc:`wideband <wideband_deconvolution>` and/or :doc:`polarized deconvolution <polarimetric_deconvolution>`, the RMS is measured from the integrated image. An example:

.. code-block:: bash

    # Clean to a 3 sigma noise level
    wsclean -auto-threshold 3 -size 2048 2048 -scale 1amin \
      -mgain 0.8 -niter 50000 observation.ms

Note that the ``-mgain`` parameter is used, in order to enable the Cotton-Schwab style major iterations. While this is not necessary for the automatic threshold to work, if the image has a very high dynamic range, the initially computed standard deviation might not be a good estimate. By using Cotton-Schwab, the standard deviation is recalculated at the beginning of every major iteration, and this will be more accurate. The ``-mgain`` parameter is discussed in more detail in the next section.

One can also specify both an automatic threshold and a manual threshold. In this case, whenever one of the thresholds is reached, the cleaning stops.

Using Cotton-Schwab: the ``-mgain`` parameter
---------------------------------------------

The ``-mgain`` parameter sets the major iteration gain: during every major iteration, the peak is reduced by the given factor. With an ``mgain`` of 0.8, the peak is reduced by 80%. This is quite a common and safe option for cleaning. With a reasonable good PSF, using 0.9 is possible without loss of accuracy, and a bit faster. With a very bad PSF, it might be necessary to lower the ``mgain`` parameter.

When ``-mgain`` is not given, or when it is set to 1, WSClean will never go back to the visibilities. It will therefore perform a simple image-based Högbom clean. While this is fast, it limits the accuracy of the image and the dynamic range that can be reached. WSClean will also not write to the MODEL column, as would be required for self-calibration (see the chapter on :doc:`self-calibration <selfcal>`). To use Högbom clean and still fill the model column, an mgain of e.g. 0.9999 can be used. Alternatively, the final model :doc:`can be predicted <prediction>` in a second WSClean run.

Note that the ``-mgain`` parameter is not the same as the ``-gain``  parameter. The latter sets the minor loop cleaning gain, and is 0.1 by default. It is almost never required to change the ``gain`` parameter.

**Next chapter:** :doc:`prediction`
