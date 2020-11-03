Weighting and gridding
======================

I got a few questions about weighting and gridding, so I decided to write down a technical description of the way the *w*-stacking algorithm of WSClean performs weighting and gridding. Note that this does not hold for the :doc:`w-gridding <wgridding>` and :doc:`idg <image_domain_gridding>` gridders.

Weighting
---------

Uniform, natural and Briggs' weighting work exactly like in Casa and other imagers:

* Uniform weighting will grid the weights on one "uv" plane (not accounting for w), and each sample will be weighted by this value before being gridded on one of the w-layers. Because it needs to grid the weights prior to filling the wlayers, an extra pass through the measurement set is required. Therefore, this option is slightly slower (~a percent) than the natural and mwa weighting modes, but normally looks the best.
* Natural weighting won't weight the samples before gridding, causing scales that have more baselines to dominate the image. For the case of the MWA, it means that because of the many small baselines, the large scales will dominate the image.
* I've also implemented Briggs' weighting (use with e.g. "``-weight briggs 0.5``") and super-pixel weighting (e.g. "``-super-weight 3``"). Briggs weighting works exactly like Casa's Briggs weighting, although I've heard that the same WSClean Briggs robust value is not exactly identical to the CASA robust value.

Super-pixel weighting is similar to Casa's "super-uniform" weighting, but is more flexible and can also be applied to other weighting schemes, for example to get "super-Briggs'" weighting. You can also set the super-pixel weighting factor to less than one to get sub-pixel weighting, which changes the image more towards naturally weighted. Super-pixel weighting might be nice on very large (full-sky) images, since the immense uv resolution might otherwise start to make your beam become natural-ish weighted (`Dan Briggs' thesis <http://www.aoc.nrao.edu/dissertations/dbriggs/>`_) describes sub/super-pixel weighting very well). I think a superweight of "image resolution / 2.5k" should be good to keep the beam uniform/Briggs' weighted. 

The weighting scheme is multiplicative to the normal sample weight, i.e., prior to the weighting scheme, samples will be weighted by their weight as specified in the measurement set, and receive zero weight if they are flagged.

Gridding
--------

WSClean can use several Gridding functions. Options include:

* Nearest-neighbour gridding
* A truncated sinc (rectangular function)
* The Kaisser-Bessel function
* A sinc function windowed by a Kaisser-Bessel function
* A Gaussian
* The Blackman-Nuttall window
* The Blackman-Harish window

The KB window is the default, and for most purposes accurate enough. By default, the function is 1023x supersampled with a kernel of size 7. With NN, a sample is placed at the nearest uv-pixel. Especially for small images, that might create aliasing and inaccurate fluxes. The supersampling places samples more accurately on the uv-grid, giving slightly more accurate fluxes. A windowed low-pass filter attenuates aliased sources. Sources that are just outside the field of view might otherwise be aliased back into the image. A Kaiser-Bessel function is very similar, but is considered slightly less optimal compared to the optimized prolate spheroidal function.

A 7 pixel low-pass filter is not very strong, and strong off-axis sources might show up as aliases ("ghost sources"). The best method to deal with those is to increase the field of view; that will also allow proper cleaning of those. Image accuracy is not so sensitive to the size or type of the windowing function, but is much more controlled by the oversampling rate. WSClean uses by default a quite large oversampling rate, because this is trivial to do with the w-stacking algorithm. This is also the reason why WSClean's output images are more accurate than when using a w-projection implementation (e.g. as implemented in CASA).

When making an image cube, it is possible to either grid the weights for each channel separately, or all on the same grid. The latter is called 'MFS weighting'; this is explained in more detail on the :doc:`MF weighting <mf_weighting>` page.
