Introduction
============

WSClean stands for "w-stacking clean". W-stacking is an alternative to the w-projection algorithm, in which the uv-samples are not convolved with a w-term correcting kernel before the FFT, but are instead corrected with a multiplication after the FFT. Multiplying each pixel is generally much faster than convolving an image, but the downside is that samples can not be added together on one grid until the w-term has been corrected. This means that in w-stacking, each w-layer needs to be stored and FFT-ed individually. The result is that w-stacking requires more memory and spends relative more time doing the FFTs, instead of gridding+convolving.

It turns out that for many telescopes, the WSClean approach is a very good trade-off, and WSClean is typically an order of magnitude faster than Casa's w-projection on MWA data. It can also handle full-sky 10k x 10k images on which Casa runs out of memory. WSClean with the same number of layers is at least as accurate as w-projection. In tests it seems even to be slightly more accurate, probably because there is no need to trim down a convolution kernel. 

:doc:`The "usage" chapter <usage>` gives an overview of how to call WSClean from the commandline.

WSClean Article
---------------

The WSClean implementation has been described and tested in the following paper:

`WSClean: an implementation of a fast, generic wide-field imager for radio astronomy <http://arxiv.org/abs/1407.1943>`_ (Offringa et al., 2014, Oct, MNRAS 444 (1): 606-619)

Abstract
~~~~~~~~

Astronomical widefield imaging of interferometric radio data is computationally expensive, especially for the large data volumes created by modern non-coplanar many-element arrays. We present a new widefield interferometric imager that uses the *w*-stacking algorithm and can make use of the *w*-snapshot algorithm. The performance dependencies of CASA's *w*-projection and our new imager are analysed and analytical functions are derived that describe the required computing cost for both imagers. On data from the Murchison Widefield Array, we find our new method to be an order of magnitude faster than *w*-projection, as well as being capable of full-sky imaging at full resolution and with correct polarisation correction. We predict the computing costs for several other arrays and estimate that our imager is a factor of 2-12 faster, depending on the array configuration. We estimate the computing cost for imaging the low-frequency Square-Kilometre Array observations to be 60 PetaFLOPS with current techniques. We find that combining *w*-stacking with the *w*-snapshot algorithm does not significantly improve computing requirements over pure *w*-stacking. The source code of our new imager is publicly released. 

**Next chapter:** :doc:`usage`
