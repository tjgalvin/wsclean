Masks and auto-masking
======================

Masks are often used to be able to clean up to/below the noise level. Deeper cleaning leaves fewer undeconvolved residuals behind, hence the image quality is higher, and it decreases the self-cal bias when used inside a self-cal loop.

Providing a mask
----------------

WSClean accepts masks in CASA format and in FITS file format. A mask is a normal, single frequency, single polarization image file, where all zero values are interpreted as being not masked, and all non-zero values are interpreted as masked.

A basic example:

.. code-block:: bash

    wsclean -fits-mask mymask.fits -niter 1000 -mgain 0.8 \
      -size 1024 1024 -scale 10asec myobservation.ms

A model image can be used as a mask file. In standard clean (i.e., not multi-scale clean) this can be used to limit the clean components to significant features. For example, one can first clean to 3 sigma (with the :doc:`-auto-threshold option <basic_cleaning>`), then use the model as fitsmask and clean to e.g. 0.3 sigma. This does require two runs though, which is not necessary. Also, while this same method can be used for multi-scale cleaning, in general a multi-scale clean will make large parts of the model image non-zero, and cleaning will therefore not be very well constrained. To overcome both issues, WSClean implements a technique called auto-masking.

Auto-masking
------------

Auto-masking allows automated deep cleaning and solves the two problems mentioned above:

 * Only one run of wsclean is required;
 * It maintains scale-dependent masks, which improves multi-scale cleaning.

Auto-masking works in all modes since :doc:`WSClean 2.2 <changelogs/v2.2>`. The general syntax is as follows:

.. code-block:: bash

    wsclean -multiscale -auto-mask 3 -auto-threshold 0.3 \
      -niter 1000000 -mgain 0.8 -size 8000 8000 -scale 2asec \
      obs.ms
    
This will start multi-scale cleaning first towards a threshold of 3 sigma. During multi-scale cleaning, the positions and scale of each component is recorded and put in a scale-dependent mask. Only the center pixel of the kernel is placed in the mask. Once the 3 sigma level is reached, cleaning will continue towards the final threshold of 0.3 sigma. During that stage, the scale-dependent masks are used to constrain the cleaning.

A scale-dependent mask makes sure that when a certain scale-kernel size was cleaned, only that same scale is allowed at that position when the mask is used.

The combination ``-auto-mask 3 -auto-threshold 0.3`` seems like a good general setting, which normally leaves almost no residuals behind.

In cases where the RMS varies strongly over the field of view, for example because of calibration artefacts, it might be useful to use [local RMS thresholding](LocalRMSThresholding), instead of thresholding relative to the global RMS.

Combining auto-masking with a regular mask
------------------------------------------

Auto-masking works in combination with normal masking, so that WSClean can be forced to to find components only in a certain area, and after reaching the auto-masking threshold it would continue with the automask. Since the automask will only contain components within the fitsmask, no components are found outside the fitsmask. Auto-masking and regular masking can be combined by specifying both ``-fits-mask <file>`` and ``auto-mask <threshold>`` on the command line.

Availability
------------

Auto-masking is available in multi-scale since :doc:`WSClean version 2.1 <changelogs/v2.1>`, and is available in all modes since :doc:`WSClean 2.2 <changelogs/v2.2>`.

References
----------
The auto-masking algorithm is described and demonstrated in `Offringa and Smirnov (2017)  <http://arxiv.org/abs/1706.06786>`_.
