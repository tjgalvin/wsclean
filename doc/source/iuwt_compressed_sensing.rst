IUWT compressed sensing
=======================

IUWT stands for isotropic undecimated wavelet transform, which is a wavelet transform suitable for astronomical imaging. The IUWT is a way to decompose an image into different scales, similar to multi-scale clean. The IUWT approach is combined with a conjugate gradient deconvolution step. This is fundamentally different to the CLEAN approach that removes a single peak at a time.

The IUWT and :doc:`MORESANE <moresane_deconvolution>` methods are somewhat experimental algorithms. They have shown good results on well calibrated data, but at the same time have shown rather bad results for data that (still) contain calibration errors. In the latter case, normal cleaning performs better. Also, with the advancement of :doc:`multi-scale cleaning <multiscale_cleaning>`, compressed sensing methods have become relevant only to corner cases.

Relation to MORESANE
--------------------

WSClean's IUWT deconvolution method is an extension/rewrite of the [MORESANE](MORESANE) deconvolution algorithm. PyMORESANE also uses IUWT and the conjugate gradient method, and is originally implemented by Jonathan Kenyon. It has shown great results on some obervations. It is implemented in Python, and can be used inside WSClean as described on the [MORESANE](MORESANE) page. The native Python implementation is quite slow, although it was specialized for GPUs, and when these could be used, the algorithm was faster. Reasons to implement a new algorithm are: to have a fast implementation in native C++ even without GPUs; to support multi-frequency and joined polarization deconvolution; and to support masks.

Usage
-----

The IUWT deconvolution method is enabled with option "``-iuwt``". This method accepts all of the normal parameters. IUWT supports masking since :doc:`WSClean 1.12 <changelogs/v1.12>`. Like the ``-multiscale`` parameter (see :doc:`multiscale <multiscale_cleaning>`), it is enabled with a single option and has reasonably optimized default settings. However, iterations are much slower. At the same time, it subtract much more flux per iteration, hence also needs fewer iterations. Therefore, ``-niter`` should normally be much lower than for cleaning; typical values are ~200 for complex images. Also, "``-gain``" can be somewhat higher, because each minor iteration performs a subminor iteration to deconvolve selected scales. A typical value for ``-gain`` is 0.2. Altogether, a typical IUWT run looks like this:

.. code-block:: bash

    wsclean -iuwt -gain 0.2 -mgain 0.6 -niter 200 \
      -size 1024 1024 -scale 2.4amin cyga.ms
    
:doc:`Wideband cleaning <wideband_deconvolution>` and :doc:`polarimetric deconvolution <polarimetric_deconvolution>` are supported, and implemented in the same way as for :doc:`multi-scale cleaning <multiscale_cleaning>`.
