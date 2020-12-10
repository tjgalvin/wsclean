Image weighting
===============

WSClean supports uniform, natural and Briggs' weighting. These work exactly like in Casa and other imagers, and are further described below.

Uniform weighting
-----------------

Uniform weighting will give each scale the same weight. If multiple visibilities fall in one uv cell, they are downweighted by the number of visibilities in that cell. In this weighting mode, the image will have the highest resolution, but the system noise will not be optimal. This gridmode is the default, but can be explicitly selected with ``-weight uniform``.

Natural weighting
-----------------

In natural weighting, the visibilities are not weighted before gridding, causing scales that have more baselines to dominate the image. Arrays have more small baselines, which will mean that large scales will be more dominantly present. This mode will make sure that the system noise has the least effect in the image, but resolution and sidelobe noise are worse. This mode is activated with option ``-weight natural``.

Briggs' weighting
-----------------

Briggs' weighting is a compromise between uniform and natural weighting. A parameter, called the robustness parameter, selects the desired level between uniform and natural. Typical values are between -1 and 1. Lower values give more uniform-like weights, higher values give more natural-like weights.  Note that values of -1 / 1 do not *equal* uniform / natural weighting. Briggs weighting mode is selected with ``-weight briggs <robustness>``, for example ``-weight briggs 0.5``.

Super-uniform and super-Briggs weighting
----------------------------------------

In some cases, slightly better results can be obtained by using super weighting. Super weighting refers to performing the counting of visibilities on larger grid cells, such that neighbouring uv-samples share their weight. This can be selected with ``-super-weight <value>``. The value represents how much larger the weight-counting cell should be. CASA's super uniform weighting mode is equivalent to ``-super-weight 3 -weight uniform``. Additionally, it is also possible to perform sub-uniform or sub-Briggs weighting, by giving a value less than one. The super-weighting can also be used in combination with Briggs' weighting. Super-weighting was implemented because many astronomers are familiar with the "super uniform" weighting option in CASA. In many cases, the weighting rank filter described below is a better option.

Super-pixel weighting can be nice on very large (full-sky) images, since the immense uv resolution might otherwise start to make your beam become natural-ish weighted. `Dan Briggs' thesis <http://www.aoc.nrao.edu/dissertations/dbriggs/>`_ describes sub/super-pixel weighting very well.

Weighting rank filter
---------------------

The weighting rank filter works in combination with one of the other weighting modes, and suppresses visibilities that would receive a much larger weight than its uv-neighbouring visibilities. This is basically done by smoothing the weights -- see the :doc:`weight rank filter page <weight_rank_filter>` for a detailed description. In uniform weighting, some visibilities can receive order of magnitudes more weight, causing these visibilities to dominate the image noise. By smoothing the weights, this issue is mitigated. The weighting rank filter does the smoothing in such a way that the PSF is not visibly affected.

Relation to visibility weights
------------------------------

In WSClean, image weighting normally works 'on top of' the visibility weights stored in the measurement set  (the ``WEIGHT_SPECTRUM`` column). This in particular means that uniform weighting will weight all scales (baselines) the same independent of their visibility weight. One therefore cannot increase the weight on long baselines by using uniform weighting and increasing the visibility weights of those baselines. The visibility weights nevertheless still have a function: they determine the weight in their particular uv-cell with respect to other visibilities *in that same uv-cell*. All in all, this means that WSClean assumes that the visibility weights specify the inverse variance weight of the visibility. Up or downweighting of scales has to be done by the image weighting (which can include :doc:`tapering <tapering>`). I've had reports that this behaviour is not consistent with how CASA treats the combination of visibility weights and imaging weights.

This behaviour can be altered by using the ``-use-weights-as-taper`` option. When specified, imaging weights will be determined without taking the visibility weights into account. Uniform weighting, for example, will make the visibilities have uniform weights based on the visibility *count* per uv-cell, without accounting for the visibility weights. The visibility weights are applied afterwards. This mode allows the visibility weights to be used to increase/decrease the weights of certain scales, just like an image taper.

An example: consider three visibilities, A, B and C. Visibilities A and B fall into the same uv-cell, and visibility C falls in a different UV cell. Assume the visibility weights of A, B and C are given by 1, 10 and 100. With normal uniform weighting, this implies that:

- Visibilities A and B are added together with weights, giving visibility B 10x more weight than visibility A.
- Visibility C is in a different cell, and since there are no other visibilities in this cell, its weight of 100 is irrelevant. In the image, the uv-cell of A and B receives the same weight as the uv-cell of C, as this is what uniform imaging weights implies.

Even though visibility C has a 10x higher weight than B, it will not be upweighted in the final image. Now, if ``-use-weights-as-taper`` would be specified:

- Again, visibilities A and B are added together with weights, giving visibility B 10x more weight than visibility A.
- The UV cell of A and B receives a weight of 11, whereas the UV cell of C receives a weight of 100.

Other weighting-related settings
--------------------------------

* :doc:`MF weighting <mf_weighting>` -- If the bandwidth is split in multiple subbands for imaging (e.g. with ``-channels-out``), it is good to know about MF weighting.
* :doc:`tapers <tapering>` -- can be used to further shape the synthesized beam.
* :doc:`The weight rank filter <weight_rank_filter>` -- used to decrease noise by avoiding giving too much weight to a single visibility.
* :doc:`Continuing a deconvolution <continue_deconvolution>` -- can be used to deconvolve a model with one weighting setting, and image the final result with different weighting settings.

**Next chapter:** :doc:`Tapering <tapering>`
