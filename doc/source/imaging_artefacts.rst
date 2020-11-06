Imaging artefacts
=================

Making a more accurate image is often more expensive, hence the default settings of WSClean are balanced between speed and accuracy, such that average imaging scenarios show accurate fluxes and no visible artefacts. Nevertheless, in high dynamic range scenarios or strong off-axis sources near the edge of the image, artefacts may appear. In those cases, you might want to increase some of the accuracy settings. Here is a short description of the most common problems.

*W*-term correction
-------------------

WSClean uses by default the *w*-stacking quite a low number of *w*-layers. This is because *w*-layers can increase memory consumption and speed quite significantly, and in most cases the default doesn't introduce any artefacts. The default settings that WSClean use correspond to a decorrelation over 1 radian at the image edge. This translates to quite a small error on sources not near the image edge. However, in the case of strong off-axis sources, those sources might show errors. In the case of snapshot imaging (typical for the MWA), faint "aliased" sources might show. In most cases, these average down when integrating more data.

If you suspect your image is affected by *w*-term corrections, you could:
* Use the :doc:`w-gridding gridder <wgridding>`;
* Increase the number of *w*-layers. The easiest way to do that is by setting the `-nwlayers-factor` option. The default is 1. With a value of 3, the number of *w*-layers will be 3 times higher, implying that the decorrelation between *w*-layers is over 1/3 radian. (The `-nwlayers-factor` option is available since :doc:`WSClean version 2.7 <changelogs/v2.7>`
