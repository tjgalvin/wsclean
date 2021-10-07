Computational performance
=========================

This chapter explains some of the ideas to speed up imaging. If you are looking to speed up imaging, you should first determine if your imaging is dominated by gridding or by deconvolution. WSClean outputs the speed of gridding (prediction + inversion) and deconvolution at the end of a run.

Gridder
-------

The biggest choice to make that typically affects computational performance the most is what gridder to use. The default gridder, the *w*-stacking gridder, is very fast for relatively small images (up to about 5k), in particular when it can be combined with the ``-parallel-gridding`` option. When parallel gridding is not possible, or when image sizes are reasonably large (say 5k-10k), the :doc:`w-gridder <wgridding>` might be a better choice. That gridder has the additional advantage that it is (much!) more accurate compared to the default gridder, although the level of accuracy of the *w*-stacking gridder is good enough for most science cases to not be a concern. Finally, for big images the :doc:`idg <image_domain_gridding>` can be even faster, in particular when gpus are available. 

Each gridder performs different for different imaging setups: image size, number of visibilities, size of *w*-terms, resolution, these all impact the performance of the various gridders in different ways. Therefore, to really know which gridder is most efficient I recommend timing each of them.

General improvements
--------------------

A few generic tips to speed up processing:

* Use :doc:`parallel deconvolution <parallel_deconvolution>`.
* If the deconvolution is a bottleneck, don't use :doc:`multi-scale <multiscale_cleaning>` if you don't need it.
* Use the parallel reordering option.
* When using *w*-stacking, use the parallel gridding option.
* Or even better, irregardless of the gridder, parallellize over multiple nodes using :doc:`wsclean-mp <distributed_imaging>`.
* When using *a*-terms, make sure to keep your *a*-term kernel as small as possible (see :doc:`a-term correction <a_term_correction>`).
* Do not include baselines that provide a resolution higher than the imaging resolution. Despite that these baselines fall outside of the *uv*-plane, and therefore are discarded and don't affect performance, at lower elevations these baselines might still fall inside the *uv* plane and might have very large *w*-terms. Yet they don't contribute to imaging quality or the resolution. The best way to filter these baselines is using a ``-maxuvw-m`` limit. One use-case where this is important is when the international LOFAR baselines are present in the observation, but imaging is done at lower resolution (e.g. up to 2'').
* Don't split your data in too many measurement sets. For example, keeping each LOFAR subband in a separate measurement set and imaging the full bandwidth is not efficient. Better is to use DP3 to concatenate the frequencies into only a few measurement sets and image those.

*w*-stacking
------------

The *w*-stacking gridder of WSClean has very different performance characteristics to for example the *w*-projection algorithm. When wsclean is executed on a machine that does not have enough memory to store all *w*-layers at once in memory, the program will make several passes over your measurement set; in the first pass it will grid and FFT only the layers with lowest w terms, in the second pass it will process the data with second-lowest *w*-terms, etc. A few passes typically do not slow down the imaging much, but if on the order of ten or more passes are executed, the performance might no longer be acceptable. In that case you can:

* Decrease the number of *w*-layers (possibly by change the phase centre to zenith); or
* Decrease the size of your image; or
* Use a machine with more memory. 

WSClean makes by default use of the `small inversion optimalization <small_inversion>`.

In the *w*-projection algorithm, the number of *w*-layers hardly affects the speed of the algorithm. However, in *w*-stacking it is almost linear in the number of *w*-layers. Hence, you should not make this number larger than necessary. If the observation has large *w*-values, because it is far off-zenith, a large number of *w*-layers might be required. Note that the *w*-projection algorithm will in such a case have a very large *w*-kernel as well, and becomes also extremely slow. In such cases it might be better to change the phase centre before imaging: see the :doc:`chgcentre <chgcentre>` documentation.
