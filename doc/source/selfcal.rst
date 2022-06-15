Self-calibration
================

WSClean can be used to perform self-cal. There are two main approaches, supporting a range of scenarios:

- Use WSClean to fill the ``MODEL_DATA`` with predicted visibilities from an imaging run. When using ``-mgain`` with a value less than 1, WSClean will fill the column. This column can then be used during calibration (e.g. with DP3). The advantage of this is that it is relatively fast and easy.
- Make WSClean output a source component list (see :doc:`component_list`) and use this list during the calibration (DP3 supports this format). The advantage of this approach is that it is easier to apply the beam (and other effects) on a source list, it is very accurate, and it is possible to prune/edit the source list. A disadvantage can be that if the number of sources is very large, calibration will be very slow.


Image-based self-calibration
----------------------------

One scenario where image-based self-cal can be useful is when combined with the ability to do (almost-)full-sky imaging, since in some situations MWA's field of view (for example) might require including highly off-axis sources in the calibration model.

To perform self-cal, you need to use major iterations, as these will fill/update the ``MODEL_DATA`` column. When using Casa's tasks for calibration, the calibration tasks will use the ``MODEL_DATA`` and calibrate the ``CORRECTED_DATA`` using this column. Other calibration utilities, like MWA's ``mitchcal`` or LOFAR's ``DP3`` can similarly be instructed to calibrate using the ``MODEL_DATA`` column.

Filling the ``MODEL_DATA`` requires a setting of mgain < 1, e.g. an mgain of 0.9. As long as mgain is not 1, WSClean will end with a major iteration, and the ``MODEL_DATA`` column will be set to the "best" model from the cleaning model.

Self-calibration from existing image
------------------------------------

The '``-predict``' option can be used to fill the ``MODEL_DATA`` column with a prediction from a pre-existing image (see :doc:`prediction <prediction>`). After having predicted model visibilities, these visibilities can be used to calibrate the data (e.g. with `DP3 <https://dp3.readthedocs.io/>`_).

Polarized imaging & calibration
-------------------------------

You can self-cal on Stokes I or on multiple polarizations, where the latter is more accurate (at least in the case of the MWA), but takes more time. If you run WSClean on the desired polarizations one by one, e.g. on XX and then on YY, or joinedly clean them, the ``MODEL_DATA`` column will have all imaged polarizations correctly filled in. The first run will create the ``MODEL_DATA`` column, set all values to zero and then fill the XX column, the second run will notice the ``MODEL_DATA`` already exists, and only update the YY column. B. McKinley has used this method and did a few self-cal loops to create very deep and well-calibrated Fornax A images. You can run all polarizations at once with ``-pol xx,xy,yx,yy``, ``-pol iquv``, or ``-pol rr,rl,lr,ll``. For info on polarimetric deconvolution settings, see :doc:`polarimetry deconvolution <polarimetric_deconvolution>`.

CASA on-the-fly mode
--------------------

Certain CASA commands (e.g. ft) will put keywords in a measurement set that turn on the "on-the-fly" (otf) mode. In OTF mode, CASA will ignore the ``MODEL_DATA`` column and use other keywords to determine the model data. To make use of the ``MODEL_DATA`` afterwards, you can use the delmod CASA command to disable OTF mode:

    delmod(vis='myobs.ms',otf=True,scr=False)

WSClean will never use or change OTF keywords in the measurement set.

**Next chapter:** :doc:`image_weighting`
