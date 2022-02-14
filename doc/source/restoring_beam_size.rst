Restoring beam size
===================

This chapter describes WSClean's algorithm for determining the synthesized beam major and minor axes and the position angle (PA).

The input to the fitting algorithm is the PSF size and an estimate of the PSF FWHM that is determined from the PSF size. Given this first estimate of the beamsize, a box width and height is calculated to be :math:`10 \times` the estimated FWHM of the PSF around the central pixel. A 3-parameter least-squares fit is then performed to find the beam's major, minor and PA values inside that box, to minimize the sum-of-squares difference between the fitted model and the real PSF over all pixels in that box of :math:`10 \times` the size of the PSF. So this includes values out to a few lobes.

Using the full image to do the fit is too expensive -- the factor of 10 is to make the fitting go fast. The fit is generally equal to a fit using a larger area.

There are some added heuristical things when fits are very different from the initial estimated PSF: WSClean will iteratively recalculate the fitting box and refit.

Highly elliptical restoring beams
---------------------------------

When the fit is highly elliptical and when this is undesirable, the option ``-circular-beam`` can be added. This constrains the fit to be circular. When forcing a circular beam, WSClean makes still sure that the resulting flux scale of restored sources is correct. This may improve the esthetics of produced images, but it might of course give a false sense of a circular beam.

History
-------

* In :doc:`WSClean 2.6 <changelogs/v2.6>`, the option ``-beam-fitting-size`` was added. It sets the box size to be used during the fit. If the beam that is fitted is different from what is expected, changing this value might help. It seems that in particular lowering the value to 1-3 can solve certain (rare) issues.
* In :doc:`WSClean 2.5 <changelogs/v2.5>`, an improved fit for small and forced circular beams was added.
