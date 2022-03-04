Continuing a deconvolution
==========================

WSClean uses the option '``-continue``' to continue a previous deconvolution run, and ``-reuse-psf`` and ``-reuse-dirty`` to skip imaging of those two and directly clean existing images.

There are multiple scenarios why continuing can be useful, for example:

 * To finish a previously unfinished deconvolution run when the initial deconvolution was prematurely stopped because of the stopping criteria (*not* because of a crash).
 * To use different imaging settings at different stages; e.g. clean up to 3 sigma, then construct a mask manually and continue a deeper deconvolution with a mask. (note that :doc:`auto-masking <masking>` is often a better way to do this)
 * If you want different image weighting settings with the same deconvolution. For example, one can run the deconvolution with uniform weights, and then reimage the set with the same deconvolution results but with different weighting.

There might be other use-cases for ``-continue``. However, not that what ``-continue`` does is a bit complex, so you should understand the specific steps that ``-continue`` performs to make effective use of it.

Syntax
------

The syntax to continue a run is to simply add '``-continue``'  to the command line:

First run:

.. code-block:: bash

    wsclean -size 1024 1024 -scale 1amin -niter 1000 \
      -auto-threshold 10 -mgain 0.8                   \
      observations.ms

Continued run:

.. code-block:: bash

    wsclean -continue -size 1024 1024 -scale 1amin -niter 1000 \
      -auto-threshold 3 -mgain 0.8                             \
      observations.ms

There are few considerations:

 * The constructed model visibilities need to be stored in the measurement set. This implies that ``-mgain`` should be used (see the :doc:`Selfcal instructions <selfcal>` for more info on mgain) in the first run, or you have to manually predict the model image from the first run, before continuing. WSClean does not verify whether this assumption holds.
 * As a result, it is not directly possible to use ``-no-update-model-required`` in the first run, because the second run requires the ``MODEL_DATA`` column to be filled. If ``-no-update-model-required`` was enabled, or only a model image without the corresponding predicted visibilities is available, it is still possible to continue the run by first :doc:`predicting <prediction>` the model data (using ``wsclean -predict ...``) from the model image before continuing.
 * In a continued run, WSClean expects the model image to be there. It also expects it to have its normal name, i.e. something like `prefix-model.fits`. For predicts that include a beam while gridding (facet-based or using IDG), the beam-corrected model image is read, with a name such as `prefix-model-pb.fits`.
 * One cannot change the image dimensions, phase center or pixel scale between the first and a continued run.
 * The second run can not include more visibilities than the first run. E.g., if selection parameters such as ``-maxuv-l`` are used, the first run should not be more restrictive than the second run, because the ``MODEL_DATA`` column would not be set for those visibilities.

This option is available since :doc:`WSClean version 2.0 <changelogs/v2.0>`.

Example: change weighting
-------------------------

A common situation for using continued imaging is to image diffuse structure in the presence of point sources. In the first run, the point sources are subtracted, while in the second run, the weighting is changed to focus on the diffuse structure.

The following sequence of statements would do this:

.. code-block:: bash

    # Run with uniform weighting to subtract point sources:
    wsclean -size 2048 2048 -scale 1amin      \
      -niter 10000 -threshold 0.1 -mgain 0.8  \
      -weight uniform                         \
      -name diffuse-field observation.ms

    # Since the images will be overwritten in the second run,
    # here I copy the images so they can be inspected for
    # debugging purposes.
    cp diffuse-field-model.fits uniform-diffuse-field-model.fits
    cp diffuse-field-image.fits uniform-diffuse-sources-image.fits

    # Finally, the run is continued with different weighting
    # settings to highlight the residual diffuse structure.
    wsclean -size 2048 2048 -scale 1amin      \
      -continue                               \
      -weight natural -taper-gaussian 3amin   \
      -name diffuse-field observations.ms

See the :doc:`chapter on weighting <image_weighting>` for more info on weighting.

What ``-continue`` really does
------------------------------

By adding ``-continue`` to the command line, WSClean will do the following things:

 * The previous model image will be read.
 * During the first inversion, WSClean will image the PSF (even thought it might already exist -- since e.g. the weights might have been changed).
 * During the first inversion, WSClean will immediately image the residual data (data - model data). This image is still named 'dirty image', even though it is actually the residual image of the first deconvolution.
 * Any new components found during cleaning will be added to the previous model image.
 * The previous model image will be overwritten.
 * The final image will be restored with the full model.

The implication of this strategy, is that if you continue a first run with a second run with the same settings, both using ``-niter N``, the resulting restored and model images are the same as when you would have run the full deconvolution with ``-niter 2N`` at once. However, during a continued run the PSF and residual image will be reimaged (to allow changing the weighting), so using the two runs will be slower. The bottom line is that one should only use ``-continue`` if there are good reasons for it.

Difference with ``-subtract-model``
-----------------------------------

The ``-subtract-model`` option only makes WSClean subtract the model column from the data column during the first imaging iteration. The useful use-case for this is to make it directly image the residual without any extra cleaning. Any previous sources won't be restored. Also, if you enable cleaning something rather confusing happens: during the first iteration, the "current" residual is imaged. This image is cleaned and a new model is formed with a few residual sources. The models of the residual sources is written to the MODEL_DATA column, and in the next iteration, the new MODEL_DATA is subtracted from the DATA, making all sources that were previously subtracted appear. The bottom line is that this option should only be used to quickly reimage a residual image without extra cleaning. If you do want cleaning, you have to manually subtract the MODEL_DATA from the DATA, which is easy e.g. with this Taql statement:

.. code-block:: bash

    taql update obs.ms set DATA=DATA-MODEL_DATA

Reusing PSF / dirty image
-------------------------

Existing PSF or dirty images can be reused to:
 - run the deconvolution algorithms of wsclean without doing the inversion (i.e. without ever going back to the visibilities)
 - speed up a second run of imaging when the PSF/dirty already exist, and no change in imaging settings (pixel scale, size, weighting, etc.) is made.

A first "regular" run to make the PSF and dirty image:

.. code-block:: bash

    wsclean -make-psf -size 1024 1024 -scale 30asec -channels-out 4 \
      obs.ms

A run that reuses the PSF and dirty images from the previous run:

.. code-block:: bash

    wsclean -reuse-psf wsclean -reuse-dirty wsclean \
      -no-reorder -name secondrun \
      -size 1024 1024 -scale 30asec -channels-out 4 \
      -niter 1000 -auto-threshold 5 obs.ms

The use of '``-no-reorder``' skips the reordering of visibilities by wsclean, which is useful when wsclean would never go back to the visibilities, as then the reordering is just overhead. It is allowed to keep the name between the runs the same (so to remove '``-name secondrun``' from the second run).
