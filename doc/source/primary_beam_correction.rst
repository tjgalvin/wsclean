Primary beam correction
=======================

WSClean can calculate and apply the primary beams of several instruments (LOFAR and MWA since :doc:`version 2.7 <changelogs/v2.7>`, and later versions supported more). It uses the EveryBeam library for this. There are two ways to correct the beam: during gridding (using IDG) and as an image plane correction. This page discusses image plane correction, whereas beam correction during gridding is discussed on the following page: :doc:`Gridding with the beam using IDG <image_domain_gridding>`. Image plane correction is faster and less complex, but will not correct for time variation of the primary beam, and won't weight the data optimally. Beam correction during gridding is therefore superior to image plane correction, but is not always required to get good results.

Image-based primary beam correction is enabled with the option ``-apply-primary-beam``. Primary beam correction can either be applied on Stokes I or on the four Stokes polarizations. In Stokes I-only mode, it is assumed that the field does not contain polarized flux (i.e., only a scalar correction is performed). For image-based correction, the beam's Jones matrix are summed as Mueller matrices. These are then inverted and applied:

    corrected = *B*\ :sup:`-1` vec(*V*)

Where *B* is the complex 4x4 Mueller matrix of the average beam, and *V* is the visibility matrix.

Supported telescopes
--------------------

Primary beams are calculated using the EveryBeam library. This library, which was originally part of WSClean, supports numerous telescopes. See the `EveryBeam README <https://git.astron.nl/RD/EveryBeam/-/blob/master/README.md>`_ for a list of supported instruments. The EveryBeam will automatically determine what instrument is associated with the measurement set, and use the right beam (if available).


Beam sensitivity limit
----------------------

By default, WSClean will clip the primary-beam-corrected image for which the response is smaller than 0.5%. This can be modified with the ``-primary-beam-limit`` parameter.

The reason for clearing pixels with low response is that the primary beam may have nulls or almost-nulls in them, and because the image is divided by those value, they may cause numerical issues. The non-primary-beam-corrected images will not be clipped, and can be inspected to analyze the full imaged area. The clipping can be turned off by setting ``-primary-beam-limit 0``.

LOFAR specifics
---------------

It is common in LOFAR observations to correct the visibilities for the beam at the phase-centre. WSClean can also correct for this. This is referred to as "differential beam correction", while the normal situation is referred to as "full beam correction". Using a differential beam leads to better polarization correction.

Full beam correction
~~~~~~~~~~~~~~~~~~~~

An example to perform Stokes I imaging and full beam correction:

.. code-block:: bash

    wsclean -apply-primary-beam -size 1024 1024 -scale 20asec \
      observation.ms

WSClean outputs the normal uncorrected image (``wsclean-image.fits``), the 8 components of the beam Jones matrix (``wsclean-beam-XX.fits``, ``wsclean-beam-XXi.fits``, ``wsclean-beam-XY.fits``, ...) and the primary beam corrected image (``wsclean-image-pb.fits``). Note that the dirty, residual and model images are not corrected.

Example to make beam-corrected Stokes I, Q, U and V images:

.. code-block:: bash

    wsclean \
      -apply-primary-beam \
      -pol iquv -size 1024 1024 -scale 20asec \
      observation.ms

WSClean outputs the 4 normal uncorrected images, the 8 components of the beam Jones matrix and the 4 primary beam images (which have ``pb`` in them).

Other combinations of polarizations, such as ``xx,yy`` or ``iq`` are currently not supported and will result in an error. If you need other modes let me know.

The beam correction can also be applied together with :doc:`multi-frequency output <making_image_cubes>`, :doc:`joined channel mode <wideband_deconvolution>` and/or the :doc:`snapshot mode <snapshot_imaging>`. In those cases, a ``pb`` image is saved for every output image. The channel-integrated "MFS" image is not corrected.

Beam correction works together with :doc:`baseline-dependent averaging <baseline_dependent_averaging>`, but only since :doc:`WSClean 2.5 <changelogs/v2.5>`. Before that, WSClean would crash or give incorrect results when combining primary beam correction with baseline-dependent averaging.

Differential beam
~~~~~~~~~~~~~~~~~

To correct an image for the beam when the phase-centre beam has been applied to the visibilities, the option '``-use-differential-lofar-beam``' can be added. ("``-apply-primary-beam``" also still needs to be given).

In summary: use "``-apply-primary-beam``" when no beam has been applied yet, and use "``-apply-primary-beam -use-differential-lofar-beam``" to apply the differential beam.

The ``REFERENCE_DIR`` column is used for determining what phase centre the beam has been applied to. Mathematically, WSClean then applies the differential beam Di as derived below. The data *V* being imaged have been premultiplied with the central beam *C* for baseline *ij*, and we want to
return a matrix that corrects the data for the full beam *B*. Given our data *R*:

.. math::

    V_{ij} = C_i^{-1} R_{ij} C_j^{-H}

we want to multiple data with a differential beam matrix *D* such that

.. math::

    D_i^{-1} V_{ij} D_j^{-H} = B_i^{-1} R_{ij} B_j^{-H}
  
With *B* the full beam matrix. We can solve for *D*:sub:`i`\ :

.. math::

    D_i^{-1} C_i^{-1} &= B_i^{-1} \\
    D_i^{-1} &= B_i^{-1} C_i \\
    D_i &= C_i^{-1} B_i \\
    
(The same could be achieved by solving for the *D*:sub:`j` term in :math:`C_j^{-H} D_j^{-H} = B_j^{-H}`).

MWA specifics
-------------
          
:doc:`Version 2.7 <changelogs/v2.7>` and upwards can directly apply the MWA beam during imaging. This avoids having to separately image XX and YY if only Stokes I is needed.

As for the other telescopes, the option to make this happen is ``-apply-primary-beam``. WSClean will determine from the telescope name stored in the measurement set that this is an MWA observation, and uses the MWA specific keywords that describe the pointing (antenna delays) of the tiles.

Usage of the MWA beam requires having installed the HDF5 file that is installed as part of the MWA repository, which will be searched at runtime. See also https://github.com/MWATelescope/mwa_pb.

Time-varying beams
------------------

When using image plane beam correction, WSClean calculates the time-integrated beam by summing snapshot beams; a beam is calculated for every 30 min and every *output* channel. Be aware that the beam correction is a single correction, and is not time-dependent. Hence, if the beam changes over time, information might smear out over the polarizations. This is less of an issue when the beam was taken out in the visibilities.

Installation information
------------------------

LOFAR beam correction is available since :doc:`WSClean version 1.11 <changelogs/v1.11>`, AARTFAAC beam correction since :doc:`WSClean version 2.6 <changelogs/v2.6>`. To use either beam, you need to have compiled WSClean with the EveryBeam library. CMake reports whether it has found the library. If WSClean has been compiled without the library, and you ask to correct for the primary beam, WSClean will report an error and stop.

