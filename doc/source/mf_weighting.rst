Multi-frequency weighting
=========================

WSClean has a special weighting mode called 'MF weighting', which is turned on with the '``-mf-weighting``' or '``-join-channels``' parameter. This is a weighting mode that changes the way in which spectral imaging (i.e., making image cubes) is performed. MF weighting is therefore only relevant when using a '``-channels-out``' parameter larger than 1. This parameter is explained on the manual page :doc:`Making image cubes <making_image_cubes>` page.

MF weighting improves the frequency-integrated (multi-frequency, "MF") images. If the individual (non-integrated) spectral images are to be used for science, **MF weighting is discouraged**. This is because MF weighting adjusts the weight of individual frequencies to make the integrated image look best. This however can create artificial spectral structures when the individual frequency images are inspected / used in science. It is fine to use ``-join-channels`` for those science cases, but make sure MF weighting is turned off with ``-no-mf-weighting``.

If MF weighting is turned on, the weights of all imaged channels will be gridded on one grid. When turned off, the weights of each imaged channel are gridded on a separate grid.

Turning MF weighting on can be useful, for example, to get uniformly weighted image cubes with the MWA. When not using MF weighting, a single MWA channel creates such a fine track in uv space, that different tracks almost never share the same uv pixel, and hence all the weights tend to be closer to '1' in non-MF, but uniform, weighting mode. The effect is that each image is weighted more towards natural weighting compared to MFS imaging. In this case, when all images are summed later on, the integrated image does not look like a uniform image, because each channel has been weighted separately.

So, this is resolved by using MF weighting. In MF weighting, all weights are gridded on the same grid, and hence the sum of the image cube equals a MF image with the same weighting.

Note that MF weighting is off by default, but using the '``-join-channels``' automatically enables MFS weighting.

Example of using MF weighting to make a briggs -1 image cube:

.. code-block:: bash

    wsclean -channels-out 768 -mf-weighting -weight briggs -1 myobservation.ms

A related option is the :doc:`weight rank filter <weight_rank_filter>` that is enabled with the '``-weighting-rank-filter``' option. The weight rank filter prevents outlier weights that might increase the noise of the image.

**Next chapter:** :doc:`Weight rank filter <weight_rank_filter>`
