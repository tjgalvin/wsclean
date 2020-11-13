Multi-scale cleaning
====================

Multi-scale cleaning is useful for accurate deconvolution of resolved sources.

Introduction
------------

WSClean supports multi-scale deconvolution with a new algorithm. Like Casa, WSClean's multi-scale deconvolution selects the highest peak and subtracts it with the best fitting 'scale', although internally it works somewhat different.

In WSClean, the scales to fit for do not have to be specified; WSClean will automatically use as many scales as necessary. The delta scale is always present. The next scale is calculated relative to the synthesized beam. Further scales are added by continuing to multiply by two until the scale is larger than the image size.

Multi-scale deconvolution can be turned on by added ``-multiscale`` to the command line. The other clean parameters that are also used for "normal" clean runs, work in multisscale mode in the same way. Here's an MWA example of a typical multi-scale deconvolution:

.. code-block:: bash

    wsclean -multiscale -mgain 0.8 -niter 50000 -auto-threshold 5 \
      -size 2048 2048 -scale 0.8amin vela.ms
    
This is a simple cleaning run, using Cotton-Schwab iterations subtracting 80% of the flux in each major iteration. The output is a normal dirty, residual, etc. image, but the model will be composed of the combination of several scales.

Bias parameter
--------------

There's one specific multi-scale parameter to tweak the results: "``-multiscale-scale-bias``". This parameter balances between how sensitive the algorithm is towards large scales compared to smaller scales. Lower values will clean larger scales earlier and deeper. It's default is 0.6, which means something like "if a peak is 0.6 times larger at a 2x larger scale, select the larger scale". (Peaks are normalized by their scale size, so the actual selection is a bit more complex). This implies that a smaller value favours selection of larger scales. The value 0.6 seems to generally work well, and was also shown to have favourable properties (Offringa & Smirnov 2017). If you use Briggs or Natural weighting to accentuate larger scales, it might be necessary/better to increase the value. In a LOFAR imaging with Briggs weighting and robustness of 0.5 (somewhat towards natural), I noticed slightly better results with a value of 0.7. With 0.6, the delta scale was almost never used. The following command runs WSClean with a scale bias of 0.7:

.. code-block:: bash

   wsclean -multiscale -multiscale-scale-bias 0.7 -mgain 0.8 \
      -niter 50000 -auto-threshold 5 -size 2048 2048 -scale 0.8amin \
      vela.ms

CASA's MSMFS and Moresane determine automatic multiscale bias values based on the strength of the PSF and/or the SNR at a given scale. This however causes image weighting to no longer have any effect on the deconvolution, which makes it impossible to tweak and/or focus on larger scales first. Therefore, WSClean uses the bias value.

Multi-scale with multi-frequency or polarimetric cleaning
---------------------------------------------------------

The multi-scale algorithm works well in combination with the parameters "``-join-channels -channels-out``". This woul perform peak finding and scale selection on the integrated image, and fit the found scale to each output frequency. This mode is explain in the :doc:`wide-band imaging <wideband_deconvolution>` chapter. In this mode, corrections are made for a wide band with a changing flux over frequency. Therefore, this mode is somewhat similar to CASA's MSMFS mode, but it is a different algorithm. The following statement would split the total bandwidth in four parts, and perform multi-frequency deconvolution.

.. code-block:: bash

    wsclean -multiscale -channels-out 4 -join-channels -mgain 0.8 \
      -niter 50000 -auto-threshold 5 -size 2048 2048 -scale 0.8amin \
      vela.ms

Multi-scale clean also works with :doc:`polarimetrc cleaning <polarimetric_deconvolution>` in a similar way. The combination of multi-frequency and multi-polarization imaging is also possible. For example, the following statement makes eight images; four images of both the xx and yy polarization, and cleaning is performed on the full integrated image:

.. code-block:: bash

     wsclean -multiscale -pol xx,yy -join-polarizations \
       -channels-out 4 -join-channels -mgain 0.8 \
       -niter 50000 -auto-threshold 5 -size 2048 2048 -scale 0.8amin \
       vela.ms
      
Masks
-----

Multi-scale cleaning works in combination with :doc:`masks <masking>`. If a mask is used, WSClean will only consider peaks that are in selected areas. For non-delta-scales, this implies that the centre of the structure has to lie in a selected area.

Auto-masking is also possible in combination with multi-scale cleaning, and is in fact one of WSClean's major improvements over other multi-scale algorithms, in terms of deconvolution quality (as was shown in Offringa & Smirnov 2017), and does even better than compressive sensing algorithms algorithms in terms of residual RMS (but not always in terms of model image quality). Because thresholds are calculated automatically, a cleaning mode such as the following command works often well without further tweaking:

.. code-block:: bash

    wsclean -multiscale -auto-threshold 1 -auto-mask 5 \
      -niter 1000000 -mgain 0.8 \
      -scale 1amin -size 4096 4096 obs.ms

For more info about masking, see the chapter on :doc:`masks and auto-masking <masking>`.

References
----------

WSClean uses a multi-scale algorithm that is an optimized version of `Cornwell (2018) <http://ieeexplore.ieee.org/document/4703304/>`_. The full algorithm is described in `Offringa and Smirnov (2017) <https://arxiv.org/abs/1706.06786>`_.

History
-------
This section will describe the multi-scale algorithm introduced in :doc:`WSClean 1.9 <changelogs/v1.9>`. Previous WSClean versions had a different multi-scale implementation, which is now deprecated.
