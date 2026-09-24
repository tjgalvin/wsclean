Rotation measure (RM) synthesis
===============================

WSClean provides several options that facilitate rotation measure synthesis. The options ``-rm-deconvolution``, ``-squared-channel-joining`` and ``-fit-rm`` are specifically intended for RM synthesis when combined with ``-join-channels`` and ``-join-polarizations``. These are explained below.

Before reading this chapter, it is useful to understand the following concepts:

 * :doc:`Making image cubes <making_image_cubes>`
 * :doc:`Wideband deconvolution <wideband_deconvolution>`
 * :doc:`Polarimetric deconvolution <polarimetric_deconvolution>`

Basic RM synthesis imaging
--------------------------

RM synthesis requires imaging the Stokes Q and U parameters across many frequency channels. After that, this QU-frequency cube can be Fourier transformed with an external tool to produce an RM cube.

Apart from the normal Stokes I imaging settings, several additional options are particularly useful to create QU image cubes intended for RM synthesis. The general procedure for creating image cubes is explained in :doc:`Making image cubes <making_image_cubes>`, including the use of ``-channels-out`` option to partition the bandwidth. The ``-pol`` option can be added to do this for Q and U:

.. code-block:: bash

    wsclean -pol QU -channels-out 100 \
      -scale 1amin -size 1024 1024 observation.ms
    
This produces an *uncleaned* image cube for the Q and U polarizations. Standard deconvolution parameters can be added to clean each polarization and frequency channel independently. However, the results will be limited by the noise in a single channel. This can be improved by using ``-join-polarizations`` to perform peak finding in QU or IQUV space. A better approach is to use the full bandwidth for cleaning. This is described below.

.. note::
    RM synthesis often involves deconvolving a large number of channels. During deconvolution, all channels need to be loaded into memory, and it is therefore quite easy to run out of memory. It may help to start with small images, while keeping an eye on memory usage. Deconvolving big images with a large number of channels may require systems with large amounts of memory.

Deconvolving in Faraday-depth space
-----------------------------------

WSClean can deconvolve in Faraday-depth space (RM space) by adding ``-rm-deconvolution`` to the command line. This selects a different deconvolution algorithm that first transforms the images to RM space. If multiscale is disabled, a Clark-like approach is used, similar to FSCLEAN described by `Bell and Enßlin (2012) <https://arxiv.org/abs/1112.4175>`_. If ``-multiscale`` is additionally specified, the MSRM algorithm is selected, which performs multiscale deconvolution in Faraday-depth space. This is described by Offringa, Mahatma, de Jong, and Edler (in preparation).

Compared to the other methods described below, RM-space deconvolution often produces the cleanest RM cubes, deconvolving many of the faint sources that are only visible in RM space. However, depending on the approach, these methods can be slow and are sensitive to calibration errors and RFI. They also do not directly optimize the residuals in the aggregated total-intensity image, :math:`<P^2>`.

A few more notes on the combination of options for RM-space deconvolution:

- Options ``pol qu -join-polarizations`` should be used to image :math:`Q` and :math:`U` and combine these in deconvolution.
- ``-join-channels -channels-out <N>`` to perform multi-frequency deconvolution with the selected number of channels.
- Because this uses a specialized algorithm, spectral fitting options are not available.
- Option ``-squared-channel-joining`` should not be specified.
- Threshold options are applied in RM space. Because peak finding is performed in squared (:math:`Q^2 + U^2`) space, the threshold is offset by the mean of the image, to accomodate for the non-zero mean in squared space.

Faraday-depth deconvolution has been available in WSClean since :doc:`version 3.8 <changelogs/v3.8>`.

Using the full bandwidth for cleaning QU cubes
----------------------------------------------

WSClean has multi-frequency deconvolution options as described in the :doc:`multi-frequency deconvolution chapter <wideband_deconvolution>`. Multi-frequency deconvolution makes it possible to use all images together during peak finding, while still subtracting the corresponding model components from the individual images. When joining channels, the default is to do the peak finding in the sum over the channels. For example, if 100 channels of Stokes Q are imaged and ``-join-channels`` is used, peak finding is performed on the integrated bandwidth in Stokes Q. This is undesirable when expecting signals with non-zero RM values, because these signals will average out over the bandwidth. Therefore, the option ``-squared-channel-joining`` was added. When used, WSClean will perform peak finding in the sum of squares of the images. Components with large RM values will thus not average out. This option can be combined with ``-join-polarizations`` to take the sum over :math:`Q^2 + U^2`.

Additionally, each CLEAN component can be fitted with a rotation measure. This will model each component by a sinusoidal signal (combining the values from Stokes :math:`Q` and :math:`U` in the fit) in squared wavelength space, effectively performing RM synthesis during deconvolution. Because only one Faraday depth is deconvolved in each iteration, it may take more iterations compared to not using ``-fit-rm``, but the benefit is that the spectra of components are more constrained, and therefore less noisy. This method is described in `Offringa & Smirnov (2017) <https://arxiv.org/abs/1706.06786>`_.

Putting everything together, this is an example of how to perform RM synthesis:

.. code-block:: bash

    wsclean -pol QU -fit-rm -join-polarizations \
      -join-channels -squared-channel-joining -channels-out <nr> \
      [-niter/-mgain/-scale/-multiscale/... etc.] observation.ms

This mode often enables deeper and more accurate cleaning than per-channel QU imaging.
 
When using ``-squared-channel-joining``, an MFS image will be stored in addition to the channel images, as is normal in multi-channel imaging, and the MFS image will still have the normal average value (not the sum of squares).

A cleaning threshold can be specified as normal, and cleaning stops when the square root of the average of squares is below that value. The statistical properties of this quantity therefore differ from the standard case, and deconvolution will generally begin cleaning noise sooner. Some experimentation with the threshold might be required.

These options work together with auto-masking, which can be effective in removing all components up to the noise.

Polarimetric multiscale deconvolution
-------------------------------------

The multi-scale mode works together with squared channel joining and RM fitting. This did not work well before version 3.8 of WSClean, where it could regularly become stuck by repeatedly selecting the same structure [#f1]_, but :doc:`WSClean 3.8 <changelogs/v3.8>` introduced a new algorithm. Be aware though that this algorithm is considerably slower than the normal multi-scale algorithm: some of the optimizations that are possible when adding channels linearly are not possible in squared-channel mode. With this algorithm, the scale is determined by correlating the individual images with the scale kernel before taking the sum of the squared images.

.. [#f1] In versions before 3.8, WSClean determined the scale and position from the integrated image, i.e. the image formed by squaring and summing the channels and polarizations. With this approach, the sum of squares may show a large structure that is actually a composition of smaller structures in the individual channels or polarizations. As a result, the subsequent fit to each individual image may fail to remove significant emission, causing WSClean to select the same large structure again in the next iteration and become stuck in this state.

Effect of different option combinations
---------------------------------------

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

Peak finding is performed independently for Stokes Q and Stokes U, using the summed image over frequency: :math:`\sum\limits_{ch} Q_{ch}` and similarly for U. Again, pixels with non-zero RM values will average out and will not be cleaned  (adding ``-squared-channel-joining`` would remedy this).

.. code-block:: bash

    wsclean -pol QU -join-polarizations -join-channels -squared-channel-joining -channels-out 100 ...

Peak finding is performed in the :math:`Q` and :math:`U` squared sum over channels: :math:`\sum\limits_{ch} Q_{ch}^2 + U_{ch}^2`. Each component spectrum is unconstrained.

.. code-block:: bash

    wsclean -pol QU -join-polarizations -join-channels -channels-out 100 ...

Peak finding is performed in :math:`\sum\limits_{ch} \sqrt{Q_{ch}^2 + U_{ch}^2}`. Note that in this mode, emission with non-zero RM values also does not get averaged out, hence squaring is not strictly necessary. One difference between this example and the above example including ``-squared-channel-joining`` is the noise properties during peak finding: taking the square root results in slightly less favourable noise statistics, hence squaring is preferred (although the difference is likely to be small). Another difference is that this example will not trigger the polarimetric multiscale algorithm, whereas the example before this one would use it.

.. code-block:: bash

    wsclean -pol QU -fit-rm -join-polarizations -join-channels -squared-channel-joining -channels-out 100 ...

This is the most advanced method, where peak finding is performed on the sum over channels of :math:`Q_{ch}^2 + U_{ch}^2`, and each component is modelled by the signal from a single Faraday depth. **When doing RM synthesis, this is often the most sensible option.** The `-multiscale` option can be added for cleaning diffuse emission.

Note that these examples differ only in how deconvolution is performed; they do not otherwise affect the output images.

Availability
------------

* ``-squared-channel-joining`` is available since :doc:`WSClean 1.12 <changelogs/v1.12>`.
* ``-fit-rm`` is available since :doc:`WSClean 3.7 <changelogs/v3.7>`.
* The polarimetric multiscale algorithm is available since :doc:`WSClean 3.8 <changelogs/v3.8>`.
