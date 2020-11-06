Restoring beam size
===================

This chapter describes WSClean's algorithm for determining the synthesized beam major and minor axes and the position angle (PA).

The input to the fitting algorithm is the PSF size and an estimate of the PSF FWHM that is determined from the PSF size. Given this first estimate of the beamsize, a box width and height is calculated to be 10 x the estimated FWHM of the PSF around the central pixel. A 3-parameter least-squares fit is then performed to find the beam's major, minor and PA values inside that box, to minimize the sum-of-squares difference over all pixels in that 10xpsf sized box between the fitted model and the real PSF. So this includes values out to a few lobes.

Using the full image to do the fit is too expensive -- the factor of 10 PSF box is to make the fitting go fast. It doesn't change the fit much whether more is used.

There's some added heuristical things when fits are very different from the initial estimated psf, it will recalculate the fitting box and refit.

When the fit is highly elliptical and when this is undesirable, the option ``-circular-beam`` can be added, which constrains the fit to come up with a circular beam. This generally improves the esthetics of produced images and does not change the flux scale, but it might of course give a false sense of a circular beam.

History
-------

* In :doc:`WSClean 2.6 <changelogs/v2.6>`, the option ``-beam-fitting-size`` was added. It sets the box size to be used during the fit. If the beam that is fitted is not what you can expect, you can try changing this value. It seems that in particular lowering the value to 1-3 can sometimes solve certain issues.
* In :doc:`WSClean 2.5 <changelogs/v2.5>`, an improved fit for small and forced circular beams was added.
