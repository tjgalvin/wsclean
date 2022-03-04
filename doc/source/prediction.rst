Prediction
==========

WSClean can be used to fill the ``MODEL_DATA`` column with the visibilities corresponding to an image. This is called 'predicting' visibilities (CASA's corresponding task is 'ft' -- sometimes this is referred to as 'degridding').

Predicting can be performed by adding '-predict' to the command line. For this to work, the input image needs to be in the exact same projection as that WSClean would output it. If your image is in a different projection, you will have to regrid it first.

An example prediction run:

.. code-block:: bash

    # Predict visibilities from a Stokes I image
    wsclean -predict -name my-image obs.ms

Be aware that the specified name (here '``-name my-image``') still specifies the prefix of the file names in the same way as that it is normally used for imaging. That means that in this run, prediction will look for an image called 'my-image-model.fits'. If you predict for multiple polarizations, e.g. with "``-pol xx,xy,yx,yy``", the files should be named accordingly, so the XX name is "``my-image-XX-model.fits``", etc. The same applies when predicting with multiple frequency intervals with '``-channels-out``'.

.. note::
    For predicts that include a beam while gridding (facet-based or using IDG), the beam-corrected model image is read. These images have a '-pb.fits' suffix, e.g. 'my-image-model-pb.fits'.

The normal use-case for using prediction is for self-calibration or subtraction. For these cases one should use a "model" image, which contains the clean components. If the input image is in absolute flux, it might be necessary to first apply the beam to the model image.

Prediction is supported since :doc:`WSClean version 1.2 <changelogs/v1.2>`. Since :doc:`WSClean 2.1 <changelogs/v2.1>`, it is no longer necessary to manually specify the image dimensions and pixel size. It is still allowed though, in which case the specified dimensions are checked against the image dimensions.

Frequency information during prediction
---------------------------------------

Some care should be taken when predicting at a different frequency than the model image or when predicting from image cubes made with WSClean (described in :doc:`making image cubes <making_image_cubes>`), because this works slightly different than CASA. At this point, WSClean does not look at the frequency of the FITS file, and therefore won't interpolate or extrapolate the model to the right frequency. For example, if an image is made with

.. code-block:: bash

    wsclean [...] ms-at-100MHz.ms

and the generated model is then predicted with

.. code-block:: bash

    wsclean -predict [...] ms-at-200MHz.ms

the model is not extrapolated to the right frequency, so one will end up with the flux levels of the 100 MHz model.

WSClean will also not look at the frequency of the FITS files when you make use of the :doc:`wide-band modes <wideband_deconvolution>` and create multiple model images. E.g. when you image with:

.. code-block:: bash

    wsclean -channels-out 4 [...] ms-at-100MHz.ms ms-at-110MHz.ms ms-at-120MHz.ms ms-at-130MHz.ms

and predict those four models into only one MS with:

.. code-block:: bash

    wsclean -channels-out 4 -predict [...] ms-at-100MHz.ms

this will split the channels in the 100MHz MS into four groups and predict the four channel images into these groups. So one will get different results in this MS compared to what one would get with:

.. code-block:: bash

    wsclean -channels-out 4 -predict [...] ms-at-100MHz.ms ms-at-110MHz.ms ms-at-120MHz.ms ms-at-130MHz.ms

Which will split the full available bandwidth into four groups and thus predict the first channel image into the first MS and so on.

Some MWA specifics
------------------

Applying the beam to an MWA image is a bit tricky, because the feeds are not orthogonal for anything but zenith. You can use the 'pbcorrect' tool (in my MWA repository) to apply a beam to an image. The basic syntax is:

.. code-block:: bash

    pbcorrect -uncorrect <image-prefix> <image-postfix> <beam-prefix> <stokes-prefix>

The input are absolute Stokes images and the output are apparent flux images with linear polarizations. The ``-uncorrect`` parameter specifies it should do the opposite of its normal operation, as it would normally make Stokes images out of wsclean's output. For example, if the ``<stokes-prefix>`` is "stokes", then ``pbcorrect`` will look for stokes-I.fits, stokes-Q.fits, stokes-U.fits and stokes-V.fits. If any of these is not present, it will be assumed zero (and a warning is issued). The 'beam' files are 8 files containing all real/imaginary components for the four linear polarizations. These can be created with the 'beam' tool in my MWA repository.

If you use ``pbcorrect`` to prepare an image for wsclean prediction, you should set ``<image-postfix>`` to "model.fits", and the image prefix is the same prefix you will specify to wsclean.

**Next chapter:** :doc:`WSClean and self-cal <selfcal>`
