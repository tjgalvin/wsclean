Introduction
============

WSClean stands for "w-stacking clean". W-stacking is an alternative to the w-projection algorithm, in which the uv-samples are not convolved with a w-term correcting kernel before the FFT, but are instead corrected with a multiplication after the FFT. Multiplying each pixel is generally much faster than convolving an image, but the downside is that samples can not be added together on one grid until the w-term has been corrected. This means that in w-stacking, each w-layer needs to be stored and FFT-ed individually. The result is that w-stacking requires more memory and spends relative more time doing the FFTs, instead of gridding+convolving.

It turns out that for many telescopes, the WSClean approach is a very good trade-off, and WSClean is typically an order of magnitude faster than Casa's w-projection on MWA data. It can also handle full-sky 10k x 10k images on which Casa runs out of memory. WSClean with the same number of layers is at least as accurate as w-projection. In tests it seems even to be slightly more accurate, probably because there is no need to trim down a convolution kernel. See the :doc:`further information page <further_information>` for more references.

**Next chapter:** :doc:`usage`
