Making spectral image cubes
===========================

Imaging can be done at several frequencies to create a 'cube' of images, by dividing the visibilities into several frequency ranges before imaging. To make image cubes with WSClean, the '``-channels-out``' parameter is used. This parameter specifies the number of different frequency ranges that the data is divided in and which are separately imaged.

As an example, we assume we have a measurement set 'myobservation.ms' with 100 channels. To image each channel separately, '``-channels-out 100``' is added to the command line, for example:

.. code-block:: bash

    wsclean -scale 1amin -size 3072 3072 -niter 1000 -auto-threshold 5 \
      -channels-out 100 myobservation.ms

This will output 100 images, named wsclean-0000-image.fits, wsclean-0001-image.fits, ..., wsclean-0099-image.fits (and similarly for the PSF, residual, and dirty). Additionally, images named like wsclean-MFS-image.fits will be output, which are the weighted sum over all 100 images.

WSClean does not output these images in a normal "image cube" like CASA does, i.e., a single fits file with several images in it. For now I've decided not to implement this (one of the reasons for this is that information about the synthesized beam is not properly stored in a multi-frequency fits file). Users can, of course, combine the outputs manually, e.g. with a simple Python script.

When you are interested in only a partial range of channels, the '``-channel-range``' parameter can be used to select that range. For example:

.. code-block:: bash

    wsclean -scale 1amin -size 3072 3072 -niter 1000 -auto-threshold 5 \
      -channel-range 60 70 -channels-out 10 myobservation.ms

This would image channels 60 to 70 and output 10 separate images, one for each channel.

The '``-channels-out``' parameter can be smaller than the available number of input channels. For example, a measurement set with 100 channels can be imaged with '``-channels-out 2``', which will divide the bandwidth in 2 parts and output 2 images. Furthermore, you can also specify multiple measurement sets at different frequencies. In such a case, the total bandwidth will be divided into the requested number of channels. For example:

.. code-block:: bash

    wsclean -channels-out 3 -niter 10000 -mgain 0.8 -auto-threshold 5 \
      band-100-MHz.ms band-110-MHz.ms band-145-MHz.ms \
      band-155-MHz.ms band-190-MHz.ms band-200-MHz.ms

This will output images at 105 MHz, 150 MHz, and 195 MHz.

Some more notes:

* The '``-channels-out``' parameter can be used in combination with '``-join-channels``' to improve the deconvolution of MFS imaging in wideband scenarios, or to get better SNR during the deconvolution of spectral imaging. This is explained in the :doc:`wideband deconvolution chapter <wideband_deconvolution>`.
* WSClean normally weights each output channel separately. This is equivalent to how other imagers do this. However, this is not always desired, especially not with modern correlators with large number of channels. Therefore, a special weighting mode called 'MF weighting' has been implemented in WSClean to improve weighting of spectral cubes. This is further explained in the :doc:`MF weighting chapter <mf_weighting>`.

