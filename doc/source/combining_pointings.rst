Combining pointings
===================

WSClean makes it possible to image and clean a combination of pointings and/or a drift-scan, also known as *joint deconvolution*. It is advisable to use :doc:`version 2.10 <changelogs/v2.10>` -- older versions support some of the features, but in particular the 'measurementset multi-field' support was added only in 2.10.

What is discussed in this page, is to grid multiple pointings together on one grid, and applying the primary beam of each observation separately. It also allows to combine observations with different antenna responses (e.g. LOFAR international stations + LOFAR core stations; heterogenous VLBI arrays; interferometers with PAF like Apertif, etc.). It is somewhat similar to mosaicking a number of pointings together with the beam in image space, except that it additionally allows deconvolution on the full, combined image. This makes it possible to deconvolve images deeper.

Setting the phase centre
------------------------

For an observation, the pointing centre can be different from the phase centre: the pointing centre is to what direction the dish is pointed and affects the primary beam. The phase centre direction is the direction at which the phases are zero, and is more or less a reference direction that affects what is actually imaged. To make it possible to combine different pointings, all pointings should be set to have the same phase centre. This can be done using the ``chgcentre`` tool that is shipped with wsclean. The phase centre of an observation can be set as follows:

.. code-block:: bash

    chgcentre observation.ms 10h01m30s 02d26m00s

This implies your combined image will be centered on 10h01m30s 02d26m00s. If the observation contains multiple fields, all fields will be set to have this phase centre. See :doc:`chgcentre <chgcentre>` for more info.

Gridding with the beam
----------------------
To combine observations weighted with the beam, each visibility should be gridded with the corresponding beam. This is only possible with the IDG gridder, hence the IDG gridder has to be used (add ``-use-idg`` and optionally ``-idg-mode hybrid`` to use GPUs: see :doc:`image domain gridding <image_domain_gridding>`).

WSClean needs to know the beam shape. A few arrays are supported 'out of the box', such as LOFAR, MWA and VLA. In that case, simply adding ``-grid-with-beam`` is enough. Otherwise, a FITS file can be specified with the appropriate beams; see :doc:`a-term correction <a_term_correction>`.

WSClean needs the beam to have a certain level of smoothness, which results in a certain support of the kernel. In case the beam has unsmooth parts (e.g. near the horizon, or near the estimated response), the default settings might not be enough. Option ``-aterm-kernel-size`` can be used to change the default. Most beams are quite smooth, and don't require a kernel larger than 5..16 pixels. Hence, a safe start is to use a size of 16 and tweak later on for performance.

It is possible to apply other direction-dependent corrections during the imaging. To do so, the list of a-terms needs to be defined using a separate config file as explained on :doc:`the a-term correction page <a_term_correction>`. Such a config file can also specify the beam gridding settings, which means that in that case the parameters like ``-grid-with-beam`` do not have to be set on the command line.

PAF settings
------------

WSClean has some specific features for joint convolution of multi-beam system, such as phased-array feeds
(e.g. APERTIF or ASKAP).
For such telescopes, it is often desirable to provide response images per beam, for which the individual beam images
have also been centred on the pointing direction of each beam. This is normally simpler than
providing response images that all have been gridded to the full field centre. What this implies is that
each beam image only has to capture the part of the beam with significant gain.

To perform such 'paf' imaging, it is expected that the data for each beam is stored in a separate
measurement set. Each beam and/or antenna may provide a separate beam fits file. The fits file may also have
a frequency axis. Additionally, the beam can be corrected for the frequency of  the visibilities.

This information can be provided by an a-term config file as described in the chapter :doc:`a-term correction <a_term_correction>`.
Here is an example:

.. code-block:: text

    aterms = [ paf ]

    paf.antenna_map = [ RT0 RT1 RT2 RT3 RT4 RT5 RT6 RT7 RT8 RT9 RT10 RT11 ]
    paf.beam_map = [ 00 01 ]
    paf.beam_pointings = [ -09h25m00.0s 55d49m59.0s -09h35m33.0s 54d47m29.0s ]
    paf.file_template = observation_cbeam-$ANT-$BEAM.fits
    paf.window = raised-hann
    paf.reference_frequency = 1355e6

To allow some flexibility, the filenames of the beams are specified with a file template. This
template may use the special variables ``$ANT`` and ``$BEAM``, which will be replaced by
values from the ``antenna_map`` and ``beam_map``, respectively. If ``$ANT`` and/or ``$BEAM`` is not used,
the map still needs to be defined, but its values are not used. The ``beam_map`` maps
each measurement set to the beam string that can be used in the file_template.
The ``antenna_map`` links the antennas in the measurement set to an antenna string that can be used in
the file_template. With the above example, visibilities from the first measurement set that need to
be corrected for antenna 0 are corrected by the fits file ``observation_cbeam-RT0-00.fits``.

The ``beam_pointings`` variable is required to hold the centres of the beam images
(ra and dec values, separated by spaces). The ``reference_frequency`` is used as the central
"unit scale" frequency of the fits file: for higher frequencies, the image size is shrunk
by their ratios, and similarly stretched for lower frequencies. If this kind of
frequency correction is not desired, the ``reference_frequency`` can be set to 0 or left out.

Other image settings
--------------------

* The image dimensions should be set such that all pointings are included. For example, if you have 2x2 pointings that connect at the fwhm, and the fwhm corresponds to 1000 pixels, the phase centre should be set to the centre of the 4 pointings and by making an image of 2000x2000, all pointings are included up to their fwhm.
* When using 'multiple field' observations, that have been phase centered to the same direction, option ``-field all`` can be used to combine all fields, or a comma-separated numbered list can be specified if only a part of the fields are to be imaged, e.g. ``-field 0,1,2`` . Note that these option are only supported when all fields have the same phase center.
* Start with a low value for ``-mgain``, like 0.6, or even lower when cleaning diverges or the uv-coverage is not so well to begin with. Cleaning a mosaick of pointings with different responses is not as stable as cleaning a homogenous observation.
* Because it is likely that the images are going to be very big when combining pointings, it is probably advisable to use the :doc:`parallel deconvolution option <parallel_deconvolution>`.
* Because the noise probably changes over the image, it is likely useful to use the ``-local-rms`` option. See the :doc:`local RMS cleaning page <local_rms_thresholding>` for more info.

Here is an example run:

.. code-block:: bash

    wsclean \
        -use-idg -grid-with-beam -aterm-kernel-size 16 -multiscale \
        -field all \
        -mem 10 -temp-dir ~ -name fullfield -weight briggs 0 \
        -size 8192 8192 -scale 500masec \
        -niter 1000000 -nmiter 10 -mgain 0.5 -auto-threshold 1 \
        -auto-mask 5 -channels-out 4 -join-channels \
        -local-rms -parallel-deconvolution 4000 \
        vla-observations.ms
        
