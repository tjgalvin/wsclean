Rotation-measure (RM) synthesis
===============================

WSClean has some options which aid rotation-measure synthesis. The options ``-squared-channel-joining`` and ``-fit-rm`` are in particular aimed at performing RM synthesis when combined with ``-join-channels`` and ``-join-polarizations``. This is explained below.

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

.. note::
    RM-synthesis often involves deconvolving a large number of channels. During deconvolution, all these channels need to be loaded in memory, and it is therefore quite easy to run out of memory. It may help to start with small images and keep an eye on memory usage. Deconvolving big images with a large number of channels might only be possible on very large-memory machines.

Using the full bandwidth for cleaning QU cubes
----------------------------------------------

WSClean has multi-frequency deconvolution options as described in the :doc:`multi-frequency deconvolution chapter <wideband_deconvolution>`. Multi-frequency deconvolutions makes it possible to use all images together for the peak finding, while still subtracting the peaks from the individual images. When joining channels, the default is to do the peak finding in the sum over the channels. For example, if 100 channels of Stokes Q are imaged and ``-join-channels`` is used, peak finding is performed on the integrated bandwidth in Stokes Q. This is undesirable when expecting signals with non-zero RM values, because these signals will average out over the bandwidth. Therefore, the option ``-squared-channel-joining`` was added. When added, WSClean will perform peak finding in the sum of squares of the images. Values with high RM values will thus not average out. This option can be combined with ``-join-polarizations`` to take the sum over :math:`Q^2 + U^2`.

Additionally, it is possible to fit each component to a rotation measure. This will model each component by a sinusodial signal (combining the value from stokes :math:`Q` and :math:`U` in the fit) in squared wavelength space, thereby performing rotation measure synthesis within the cleaning procedure. Because only one Faraday depth is deconvolved in each iteration, it may take more iterations compared to not using ``-fit-rm``, but the benefit is that the spectra of a component is more constrained, and therefore less noisier. This method is described in `Offringa & Smirnov (2017) <https://arxiv.org/abs/1706.06786>`_.

Combining this all together, this is an example of how to perform RM-synthesis:

.. code-block:: bash

    wsclean -pol QU -fit-rm -join-polarizations \
      -join-channels -squared-channel-joining -channels-out <nr> \
      [-niter/-mgain/-scale/... etc.] observation.ms

This mode (often) allows one to clean deeper and more accurately compared to per-channel QU cleaning.
 
When using ``-squared-joining``, an MFS image will be stored in addition to the channel images, as is normal in multi-channel imaging, and the MFS image will still have the normal average value (not the sum of squares).

A cleaning threshold can be given as normal in Jy, and cleaning stops when the square root of the average of squares is below that value. The statistics of this value is thus slightly different as normal, and in general one will start cleaning the noise quicker. Some experimentation with the threshold might be required.

These options work together with auto-masking, which can be effective in removing all components up to the noise.

The squared channel joining also works together with the multi-scale mode. I have noticed though that the multi-scale algorithm regularly gets stuck when it is asked to search on squared values. This is because the sum of squares might show a large structure, while it is actually just a composition of small structures in the individual channels/polarizations, so that the fit to each individual channel will not remove anything, and WSClean will continue to find that large structure. If this happens, I would first suggest to clean the Q and U polarizations independently: the squared sum over one polarization is much less likely to behave this way. Otherwise, you could stop cleaning before the problem occurs, or turn multi-scale off altogether. There is a somewhat trivial fix for wsclean to git rid of this issue (namely to first convolve each channel and then integrate instead of integrating and then convolving), which I hope to implement in the future.
 
What is calculated?
-------------------

Here are a few examples, along with a description of how peak finding is performed with the given settings:

.. code-block:: bash

    wsclean -pol QU -channels-out 100 ...

Cleaning is performed independently for each polarization, peaks are found in each individual image.

.. code-block:: bash

    wsclean -pol QU -join-polarizations -channels-out 100 ...

Peak finding is performed in :math:`Q^2 + U^2`, but independently for each channel.

.. code-block:: bash

    wsclean -pol Q -join-channels -channels-out 100 

Peak finding is performed in the sum over channels of one polarization, :math:`\sum\limits_{ch} Q_{ch}`. Pixels with non-zero RM values will average out and will not be cleaned (adding ``-squared-channel-joining`` would remedy this).

.. code-block:: bash

    wsclean -pol QU -join-channels -channels-out 100 

Peak finding is performed in the sum over channels, :math:`\sum\limits_{ch} Q_{ch}`, and separately for :math:`Q` and :math:`U`. Again, pixels with non-zero RM values will average out and will not be cleaned  (adding ``-squared-channel-joining`` would remedy this).

.. code-block:: bash

    wsclean -pol QU -join-polarizations -join-channels -squared-channel-joining -channels-out 100 ...

Peak finding is performed in the :math:`Q` and :math:`U` squared sum over channels: :math:`\sum\limits_{ch} Q_{ch}^2 + U_{ch}^2`. The spectra of each component is unconstrained.

.. code-block:: bash

    wsclean -pol QU -join-polarizations -join-channels -channels-out 100 ...

Peak finding is performed in :math:`\sum\limits_{ch} \sqrt{Q_{ch}^2 + U_{ch}^2}`. Note that in this mode, flux with non-zero RM-values also does not get averaged out, hence squaring is not stricly necessary. The only difference between this example and the above example including ``-squared-channel-joining`` is the noise properties during peak finding: the square root makes the noise behave slightly worse, hence squaring is preferred (albeit that the difference is probably minor).

.. code-block:: bash

    wsclean -pol QU -fit-rm -join-polarizations -join-channels -squared-channel-joining -channels-out 100 ...

This is the most advanced method, where peak-finding is performed on the sum over channels of :math:`Q_{ch}^2 + U_{ch}^2`, and each component is modelled by the signal from a single Faraday depth. **When doing RM-synthesis, this is often the most sensible option.**

Note that these examples only differ in how cleaning is performed, they do not affect the output images otherwise.

Availability
------------

* ``-squared-channel-joining`` is available since :doc:`WSClean 1.12 <changelogs/v1.12>`.
* ``-fit-rm`` is available since :doc:`WSClean 3.7 <changelogs/v3.7>`.
