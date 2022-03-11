General usage
=============

WSClean is a command line program. The "``wsclean``" executable accepts many parameters. The "``wsclean``" program can be run without parameters to get a list of allowed parameters together with a short description. The general syntax of wsclean is as follows:

.. code-block:: bash

    wsclean [-options] <obs1.ms> [<obs2.ms> ..]

Here's an example of a typical run for an MWA observation:

.. code-block:: bash

    # Make a Stokes I image of size 3072 x 3072 with 0.7' pixels,
    # do 10000 iterations, use CS cleaning and stop when the peak
    # flux has reached 3σ level.
    wsclean -size 3072 3072 -scale 0.7amin -niter 10000 \
      -mgain 0.8 -auto-threshold 3 obs.ms

This performs a Cotton-Schwab clean. The Cotton-Schwab algorithm is enabled with the "``-mgain 0.8``" parameter, which means that the peak flux is reduced by 80% until a new major iteration is started. Multiple measurement sets can be specified on the command line to image the integration of those observations.

For fast imaging with less accuracy, you can perform a Högbom clean and disable padding:

.. code-block:: bash

    # Similar to above statement, but now only Högbom cleaning
    # and no padding
    wsclean -size 3072 3072 -scale 0.7amin -niter 10000 \
      -auto-threshold 3 -padding 1 obs.ms

Because the '``mgain``' parameter was left out, WSClean will not iteratively go back to the visibilities. This also implies that the MODEL_DATA column will not be filled (see :doc:`self-calibration with WSClean <selfcal>`) for more info).

A description of the basic cleaning parameters is given in the :doc:`basic cleaning chapter <basic_cleaning>`.

Multiple measurement sets
~~~~~~~~~~~~~~~~~~~~~~~~~

When specifying multiple measurement sets on the command line, WSClean images the combined data from all measurement sets together. For example,

.. code-block:: bash

    wsclean -size 1024 1024 -scale 1asec \
      -niter 1000000 -mgain 0.8 \
      -multiscale -auto-threshold 1 -auto-mask 5 \
      lofar-100mhz.ms vla-l-band.ms vla-c-band.ms

will make deconvolved images from the combination of the 3 specified measurement sets.
      
The measurement sets may be observations at different frequencies, at different times and/or from different instruments. They are required to have the same phase centre though. To image measurement sets with different phase centres, the measurement sets should be :doc:`phase rotated to the same direction <chgcentre>` before imaging. If the measurement sets have different pointings and it is desirable to correct for the primary beam, it is possible to :doc:`combine pointings by gridding with the beam <combining_pointings>`. 

When splitting the data into so-called output channels or when selecting a subset of channels using the ``-channel-range`` option, the full combined bandwidth (over all measurement sets) is split / divided. See the :doc:`making image cubes <making_image_cubes>` and :doc:`wideband deconvolution <wideband_deconvolution>` chapters for more information. This is different from selecting sub-intervals using ``-interval`` or ``-intervals-out``, for which timestep indices are determined from the first measurement set only, and the resulting timestep index selection is then applied to all the measurement sets. More information about time selections can be found in the :doc:`snapshot imaging chapter <snapshot_imaging>`.

An advanced MWA example
~~~~~~~~~~~~~~~~~~~~~~~

As a more enhanced example, here is a command to clean `MWA GLEAM <https://www.mwatelescope.org/gleam>`_ data:

.. code-block:: bash

    wsclean -name obs-1068210256 \
      -size 4000 4000 -niter 1000000 -mgain 0.95 \
      -weight briggs -1.0 -scale 0.75amin \
      -auto-threshold 1 -auto-mask 5 -multiscale \
      -channels-out 4 -join-channels \
      -pol xx,yy -join-polarizations \
      1068210256.ms

The explanation of this command:

* Briggs' weighting with robustness of -1 is used. For the MWA, this decreases the noise in single snapshots slightly. See :doc:`image weighting <image_weighting>` for more info on supported weightings.
* The (instrumental) XX and YY polarizations are imaged separately and cleaned together. This allows more accurate primary beam correction for the MWA. See :doc:`polarimetric deconvolution <polarimetric_deconvolution>` for more info.
* A high ``mgain`` value is used because the MWA synthesized beam is well behaved.
* A larger image is made because GLEAM includes lower frequencies, at which the primary beam is larger.
* For GLEAM, the W-snapshot algorithm is used by executing the ``chgcentre`` command prior to imaging, as described on the :doc:`w_snapshot_algorithm`.

**Next chapter:** :doc:`basic_cleaning`
