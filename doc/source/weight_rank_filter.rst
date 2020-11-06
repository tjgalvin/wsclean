Weight rank filter
==================

The weight rank filter can be used to decrease noise by avoiding giving too much weight to a single visibility.

To determine the :doc:`image weights <image_weighting>`, the weights are gridded on a regular grid. With uniform or Briggs' weighting, the total weight of a visibility is determined by the number of visibilities that fall into that cell: the fewer visibilities fall into the same cell, the more weight this visilibity receives. Because the visibilities form a track in uv-space, it can happen that a particular visibility is the only visibility gridded in a particular uv-cell, while other uv-cells have on average hundred, thousands, ... of visibilities in them. This implies that this single visibility will receive a lot of weight.

The weight rank filter in WSClean avoids this by calculating a local weight RMS of the grid, and truncating all weights that are some factor above the local RMS. This effectively decreases the weight of those visibilities that would receive excessive weight otherwise. This results in very little change in the PSF, but can effect the noise in the image significantly.

Effect & improvement
--------------------

The exact effect of this option changes very much on the type of observation, image size and uv-coverage. I've seen improvements of 20-30% reduction in noise in certain VLA observations, while I hardly saw any improvement in certain LOFAR and MWA observation.

To observe the effect of this filter, it can be useful to save the uv-weights with ``-save-weights``, and perform an imaging run with and without the filter.

Syntax
------

The parameter to enable the weight rank filter is ``-weighting-rank-filter``, which takes the factor relative to the RMS above which the values are truncated. For example:

.. code-block:: bash

    wsclean -weighting-rank-filter 3 [other options] observation.ms

The default size over which the RMS is calculated, is 16 uv-pixels. During this calculation, only uv-cells with coverage are used in the calculation. If the visibility coverage is very sparse, the size of 16 uv-pixels might be too small. It can be changed with the ``-weighting-rank-filter-size`` option, for example:

.. code-block:: bash

    wsclean -weighting-rank-filter 3 -weighting-rank-filter-size 64 \
        [other options] observation.ms

Since :doc:`WSClean 2.2 <changelogs/v2.2>`, the default in WSClean is to apply a weighting rank filter with a value of 3. It is recommended to always apply a weak rank filter, because it improves the noise with hardly any change on the PSF. In case there are no visibilities with excessive weight, this filter will do nothing.

The weight rank filter is available since :doc:`WSClean 1.8 <changelogs/v1.8>`. In versions before :doc:`WSClean 2.2 <changelogs/v2.2>`, the default setting was to not apply a weight rank filter, and since 2.2 it was turned on with a value of 3. As far as I'm aware, there's no weight rank filter option in CASA.

Related topics
--------------

* :doc:`Image weighting <image_weighting>`, which discusses natural, uniform and Briggs' weighting, super-uniform weighting, etc.
