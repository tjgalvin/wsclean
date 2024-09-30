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

Scale kernel shape
------------------

WSClean supports two different kernel shapes for multi-scale cleaning: a tapered quadratic function and a Gaussian. The tapered quadratic function is the function that was originally introduced by `Cornwell (2008) <https://arxiv.org/abs/0806.2228>`_. It is smooth and goes to zero quite quickly, which can have advantages in some cases. A Gaussian function is preferred when the model needs to be described in analytical components. Almost all software supports Gaussians as sky model components, because it has a analytical Fourier transform. WSClean uses by default the tapered quadratic function, unless a :doc:`component list <component_list>` is requested, in which case Gaussians are used. The difference between these two is rarely noticeable.

Be aware that WSClean lists scale sizes on the command line in terms of the quadratic kernel size :math:`\alpha` instead of the Gaussian FWHM. If Gaussians are requested, the Gaussian width parameter :math:`\sigma` is calculated using :math:`\sigma = \frac{3}{16}\alpha`, which approximately matches the width of the tapered quadratic function of size :math:`\alpha`. This implies that the full-width half-maximum (FWHM) of the Gaussian is given by :math:`\textrm{FWMH} = 2 \frac{3}{16} \alpha \sqrt{2 \ln 2} \approx 0.45 \alpha`. When a component list output is requested, the kernel scales are corrected for this factor, such that the output sky model describes the Gaussian FWHM size.

Bias parameter
--------------

There's one specific multi-scale parameter to tweak the results: "``-multiscale-scale-bias``". This parameter balances between how sensitive the algorithm is towards large scales compared to smaller scales. Lower values will clean larger scales earlier and deeper. Its default is 0.6, which means something like "if a peak is 0.6 times larger at a 2x larger scale, select the larger scale". (Peaks are normalized by their scale size, so the actual selection is a bit more complex.) This implies that a smaller value favours selection of larger scales. The value 0.6 seems to generally work well, and was also shown to have favourable properties (`Offringa & Smirnov 2017 <https://arxiv.org/abs/1706.06786>`_). If you use Briggs or Natural weighting to accentuate larger scales, it might be necessary/better to increase the value. In LOFAR imaging with Briggs weighting and robustness of 0.5 (somewhat towards natural), I noticed slightly better results with a value of 0.7. With 0.6, the delta scale was almost never used. The following command runs WSClean with a scale bias of 0.7:

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

Multi-scale clean also works with :doc:`polarimetric cleaning <polarimetric_deconvolution>` in a similar way. The combination of multi-frequency and multi-polarization imaging is also possible. For example, the following statement makes eight images; four images of both the xx and yy polarization, and cleaning is performed on the full integrated image:

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

Adaptive-scale pixel algorithm
------------------------------
WSClean supports a variant of multi-scale called adaptive-scale pixel (ASP), originally described in `Bhatnagar & Cornwell (2004) <https://www.aanda.org/articles/aa/abs/2004/41/aa0354-04/aa0354-04.html>`_. Similar to the multi-scale algorithm, ASP starts by finding the brightest scale from a set of pre-selected Gaussian scales. Once found, ASP performs an additional fitting procedure to fit the position, Gaussian axes and strength of the found shape (taking into account the shape of the point-spread function).

ASP can be simple enabled by adding `-asp` to the command line, e.g.:

.. code-block:: bash

    wsclean -asp -auto-threshold 5 -niter 1000000 -mgain 0.8 \
      -scale 1amin -size 1024 1024 obs.ms

The advantage of ASP is that it can more accurately fit certain source components, particularly elliptical shaped ones. This freedom can also be a downside, causing over-fitting or mis-fitting of the Gaussian. For example, in complex shapes with sharp edges, ASP can produce very elliptical, thin components that "overshoot" at the feature's edges.

WSClean supports :doc:`multi-frequency <wideband_deconvolution>` and :doc:`multi-polarization <polarimetric_deconvolution>` ASP, and it allows the use of the :doc:`parallel deconvolution option <parallel_deconvolution>`. These are extensions of the original ASP algorithm. The ASP algorithm is *very* slow; it is not uncommon to take more than 100x to 1000x more time than WSClean's multi-scale algorithm. This is in part because no real scenarios were found in which it outperformed multi-scale, hence not a lot of effort was spent on making the algorithm faster. There may nevertheless be some cases, particular with very smooth sources, in which ASP may perform better than multi-scale in terms of accuracy.

The ASP algorithm in WSClean does not support auto-masking. One reason for this is that with ASP, every found component will be a unique component. Hence, to limit the deconvolution to the set of already found components does not have as much benefit as it has for multi-scale. Auto-masking is an important improvement to the accuracy of multi-scale, and the lack of auto-masking in ASP is another reason why it is hard to get better results with ASP, compared to multi-scale clean.

References
----------

A first multi-scale algorithm was described by `Cornwell (2008) <https://arxiv.org/abs/0806.2228>`_. WSClean implements an alternative to this algorithm, that is significantly faster and has support for multi-frequency deconvolution. It is described in `Offringa and Smirnov (2017) <https://arxiv.org/abs/1706.06786>`_.

History
-------
This section will describe the multi-scale algorithm introduced in :doc:`WSClean 1.9 <changelogs/v1.9>`. Previous WSClean versions had a different multi-scale implementation, which is now deprecated.
