Rotation-measure (RM) synthesis
===============================

WSClean has some options which aid rotation-measure synthesis. The option ``-squared-channel-joining`` is in particular aimed at performing RM synthesis when combined with ``-join-channels`` and ``-join-polarizations``. This is explained below.

Before reading this chapter, it is useful to understand the following concepts:

 * :doc:`Making image cubes <making_image_cubes>`
 * :doc:`Wideband deconvolution <wideband_deconvolution>`
 * :doc:`Polarimetric deconvolution <polarimetric_deconvolution>`

Basic RM-synthesis imaging 
--------------------------

For RM-synthesis, one images the Stokes Q and Stokes U polarizations of a measurement set and produces images at many frequencies. After that, this QU-frequency cube can be Fourier transformed with an external tool to produce an RM cube.

Of course, making those QU images can be done with WSClean as expected, but there are some extra options which can be useful for RM-synthesis. The basic approach to make the cube is described in the manual page :doc:`Making image cubes <making_image_cubes>`, where it is described how to use the ``-channels-out`` option to partition the bandwidth. The ``-pol`` option can be added to do this for Q and U:

.. code-block:: bash

    wsclean -pol QU -channels-out 100 \
      -scale 1amin -size 1024 1024 observation.ms
    
This produces an *uncleaned* image cube for the Q and U polarizations. Standard cleaning parameters can be added to clean each polarization and each channel individually. However, by doing so one is limited by the noise in a single channel. This can be slightly improved by using the ``-join-polarizations`` option, and perform peak finding in QU or IQUV space. It is however beneficial to use the full bandwidth for cleaning. This is described below.

Using entire bandwidth for cleaning QU cubes
--------------------------------------------

WSClean has multi-frequency deconvolution options as described on the :doc:`multi-frequency deconvolution chapter <wideband_deconvolution>`. This mode allows using all images together for the peak finding, while still subtracting the peaks from the individual images. The default is to combine the images by summing over the channels. For example, if one is imaging Stokes Q over 100 channels and ``-join-channels`` is used, peak finding is performed on the integrated bandwidth in Stokes Q. This is undesirable when expecting flux with non-zero RM values, because the flux will average out over the bandwidth. Therefore, the option ``-squared-channel-joining`` was added. When added, WSClean will perform peak finding in the sum of squares of the images. Values with high RM values will thus not average out. This option can be combined with ``-join-polarizations`` to take the sum over :math:`Q^2 + U^2`. This is therefore the recommended way to perform RM-synthesis:

.. code-block:: bash

    wsclean -pol QU -join-polarization \
      -join-channels -squared-channel-joining -channels-out <nr> \
      [-niter/-mgain/-scale/... etc.] observation.ms

This mode allows one to clean deeper and more accurately compared to per-channel QU cleaning.
 
When using ``-squared-joining``, an MFS image will be outputted as normal in addition to the channel images, and the MFS image will still have the normal average value (not the sum of squares). A cleaning threshold can be given as normal in Jy, and cleaning stops when the sqrt of the average of squares is below that value. The statistics of this value is thus slightly different as normal, and in general one will start cleaning the noise quicker. Some experimentation with the threshold might be required.

The squared channel joining also works together with the multi-scale mode. I have noticed though that the multi-scale algorithm regularly gets stuck when it is asked to search on squared values. This is because the sum of squares might show a large structure, while it is actually just a composition of small structures in the individual channels/polarizations, so that the fit to each individual channel will not remove anything, and WSClean will continue to find that large structure. If this happens, I would first suggest to clean the Q and U polarizations independently: the squared sum over one polarization is much less likely to behave this way. Otherwise, you could stop cleaning before the problem occurs, or turn multi-scale off altogether. There is a somewhat trivial fix for wsclean to git rid of this issue (namely to first convolve each channel and then integrate instead of integrating and then convolving), which I hope to implement in the future.
 
What is calculated?
-------------------

Here are a few examples, along with a description of how peak finding is performed with the given settings:

.. code-block:: bash

    wsclean -pol QU -channels-out 100 ...

Cleaning is performed indepently for each polarization, peaks are found in each image.

.. code-block:: bash

    wsclean -pol QU -join-polarizations -channels-out 100 ...

Cleaning is performed in Q^2 + U^2, but independently for each channel.

.. code-block:: bash

    wsclean -pol Q -join-channels -channels-out 100 

Cleaning is performed in sum over channels of Q_ch: Pixels with non-zero RM values will average out and will not be cleaned (``-squared-channel-joining`` should be added).

.. code-block:: bash

    wsclean -pol QU -join-channels -channels-out 100 

Cleaning is performed in sum over channels of Q_ch and separately for Q and U. Again, pixels with non-zero RM values will average out and will not be cleaned  (``-squared-channel-joining`` should be added).

.. code-block:: bash

    wsclean -pol QU -join-polarizations -join-channels -squared-channel-joining -channels-out 100 ...

Cleaning is performed in sum over channels of :math:`Q_{ch}^2 + U_{ch}^2`. **When doing RM-synthesis, this is the most sensible option.**

.. code-block:: bash

    wsclean -pol QU -join-polarizations -join-channels -channels-out 100 ...

Cleaning is performed in sum over channels of :math:`\sqrt{Q_{ch}^2 + U_{ch}^2}`. Note that in this mode, flux with non-zero RM-values also does not get averaged out, hence squaring is not stricly necessary. The only difference between this example and the above example including ``-squared-channel-joining`` is the noise properties during peak finding: the square root makes the noise behave slightly worse, hence squaring is preferred (albeit that the difference is probably minor).

Note that these examples only differ in how cleaning is performed, they do not affect the output images otherwise.

Availability
------------

``-squared-channel-joining`` is available since :doc:`WSClean 1.12 <changelogs/v1.12>`.
