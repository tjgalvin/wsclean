Self-calibration
================

Using an imaging step to make a self-calibration model is quite a common scenario in radio astronomy. WSClean can be used to perform this imaging step. Imaging is often used in self-calibration when it is not sufficient to just transfer the solutions from a calibration observation, for example because the ionosphere needs to be calibrated for in the direction of the target observation. There are three main approaches, supporting a range of scenarios:

- Perform self-calibration from model visibilities. When performing a normal imaging run with WSClean, with an ``-mgain`` value less than 1 (which will enable major iterations in the imaging run), WSClean will fill the ``MODEL_DATA`` column. This column can then be used during calibration (e.g. with `DP3 <https://dp3.readthedocs.io/>`_). The advantage of this is that it is relatively fast and easy.

  Tools like ``DP3`` can be instructed to calibrate from the ``MODEL_DATA`` column. When using Casa's tasks for calibration, the calibration tasks will also use the ``MODEL_DATA`` and calibrate the ``CORRECTED_DATA`` using this column.
  
  Filling the ``MODEL_DATA`` requires a setting of ``mgain`` < 1. A typical value of ``mgain`` is ``0.8``. As long as ``mgain`` is not 1, WSClean will end with a major iteration, and the ``MODEL_DATA`` column will be calculated for the final model that was constructed during the deconvolution.
  
  This approach does not allow direction-dependent calibration. For that, use one of the following methods.

- Let WSClean output a source component list (see :doc:`component_list`) and use this list during the calibration. Tools like DP3 support this format directly. The advantages of this approach are i) that it is easier to apply the beam (and other effects) on a source list; ii) it is very accurate; iii) it is possible to prune/edit the source list. A disadvantage can be that if the number of sources is very large, the direct prediction of the component list (inside calibration) will be very slow.

  Direction-dependent calibration can be performed by clustering the components, and using a direction-dependent calibration tool such as DP3. This is the approach used by the `Rapthor <https://rapthor.readthedocs.io/>`_ pipeline.

  When performing an imaging run with the purpose of creating a component list, some performance can be gained by using the ``-skip-final-iteration`` option. This option skips the prediction-imaging round after finishing the deconvolution. The final images may have slightly lower quality, but the model (images/list) is complete.

- Perform self-calibration from an already existing model image. The ``-predict`` option can be used to fill the ``MODEL_DATA`` column with a prediction from a pre-existing image (see :doc:`prediction <prediction>`). After having predicted model visibilities, these visibilities can be used to calibrate the data. Direction-dependent calibration can be performed by predicting into different columns (using the ``-model-column`` option).

Polarized imaging & calibration
-------------------------------

It is possible to self-calibrate on Stokes I or on multiple polarizations. If you run WSClean on the desired polarizations one by one, e.g. on ``XX`` and then on ``YY``, or jointly clean them, the ``MODEL_DATA`` column will have all imaged polarizations correctly filled in. Assuming the column does not exist yet, the first run will create the ``MODEL_DATA`` column, set all values to zero and then fill the ``XX`` column. The second run will notice the ``MODEL_DATA`` already exists, and only update the ``YY`` column. B. McKinley has used this method and did a few self-cal loops to create very deep and well-calibrated `Fornax A images <https://arxiv.org/abs/1411.1487>`_. All polarizations can be selected by using the options ``-pol xx,xy,yx,yy``, ``-pol iquv``, or ``-pol rr,rl,lr,ll``. For info on polarimetric deconvolution settings, see :doc:`polarimetry deconvolution <polarimetric_deconvolution>`.

CASA on-the-fly mode
--------------------

Certain CASA commands (e.g. ft) will put keywords in a measurement set that turn on the "on-the-fly" (otf) mode. In OTF mode, CASA will ignore the ``MODEL_DATA`` column and use other keywords to determine the model data. To make use of the ``MODEL_DATA`` afterwards, you can use the delmod CASA command to disable OTF mode:

```
delmod(vis='myobs.ms',otf=True,scr=False)
```

WSClean will never use or change OTF keywords in the measurement set.

**Next chapter:** :doc:`image_weighting`
