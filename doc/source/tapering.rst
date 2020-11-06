Tapering
========

Tapering can be used to shape the synthesized beam, by tapering off the weight of visibilities in the uv plane. This page lists the tapering options supported by WSClean.

Tapering works cumulatively with weights set by the imaging weights (uniform, Briggs, natural, see :doc:`the image weighting chapter <image_weighting>`) and the visibility weight. Tapering is available in WSClean since :doc:`version 1.10 <changelogs/v1.10>`. 

Gaussian taper
--------------

A Gaussian taper multiplies the uv-weights with a Gaussian function, which therefore makes the synthesized beam approach a Gaussian function. A Gaussian taper is selected with `-taper-gaussian <beamsize>`. The beamsize is by default in arcseconds, but can be given with different units, e.g. "2amin". A Gaussian taper in uv space will be calculated such that the FWHM of the Gaussian in real space has the given beam size.

Tukey taper
-----------

Tukey tapers can be used to smooth edges with a *Tukey window*, also known as a *tapered cosine window*. WSClean allows tapering three edges:
    * ``-taper-tukey``, a circular taper that smooths the outer edge set by ``-maxuv-l``
    * ``-taper-inner-tukey``, a circular taper that smooths the inner edge set by ``-minuv-l``
    * ``-taper-edge-tukey``, a square-shaped taper that smooths the edges as set by the uv grid and ``-taper-edge``.

The following plot shows how the Tukey-taper parameters influence the weights:

.. image:: images/WeightTapers.svg
    :alt: LOFAR RFI spectrum
    :width: 100%

**Next chapter:** :doc:`MF image weighting <mf_weighting>`
