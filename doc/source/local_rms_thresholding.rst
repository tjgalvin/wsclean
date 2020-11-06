Local RMS thresholding
======================

Since :doc:`WSClean 2.3 <changelogs/v2.3>`, WSClean supports setting thresholding levels that are relative to the local RMS. This is useful when the RMS varies strongly over the field of view. A typical use-case for local RMS thresholding is for cleaning images that contain strong calibration artefacts. In this case, strong sources might generate strong artefacts. Those artefacts might be much stronger than faint sources elsewhere in the field of view. To be able to clean those faint sources, but not clean the artefacts, local RMS thresholding is required.

There are two ways of using local RMS thresholding:

 * WSClean can automatically generate an RMS map; or
 * An pre-made RMS map can be supplied to WSClean.

The two cases will be described in the sections below.

Note that when using local RMS thresholding, values reported during the minor cycles are no longer absolute fluxes, but have been scaled in some way to account for the local RMS. 

Auto-generated RMS map
----------------------

When WSClean is asked to use a local RMS without supplying an RMS image, WSClean will generate an RMS map. This is done every major iteration, before the clean loop is called. A typical run looks like this:

.. code-block:: bash

    wsclean -size 1024 1024 scale 10asec -local-rms \
      -auto-threshold 3 -mgain 0.8 -niter 1000000 \
      [...] observation.ms
    
The local RMS thresholding is typically to be used together with the ``-auto-threshold`` parameter. While it is possible to combine it with a normal ``-threshold`` parameter, which will scale the threshold relatively to the minimum in the RMS map, this is not recommended. Local RMS thresholding can also be used together with [automatic masking](Masking) using the ``-auto-mask`` parameter. A typical run to do this looks like this:

.. code-block:: bash

    wsclean -size 1024 1024 scale 10asec -local-rms \
      -auto-threshold 0.3 -auto-mask 3 -mgain 0.8 -niter 1000000 \
      -multiscale [...] observation.ms
    
Note that the ``-local-rms`` parameter was called ``-rms-background`` before :doc:`WSClean 2.4 <changelogs/v2.4>`.

The default configuration is to calculate the RMS over a Gaussian kernel with dimensions of 25 times the PSF. This value can be changed with the ``-local-rms-window`` option, for example:

.. code-block:: bash

    wsclean -size 1024 1024 scale 10asec -local-rms \
      -local-rms-window 10 -auto-threshold 3 -mgain 0.8 \
      -niter 1000000 [...] \
      observation.ms
    
Instead of just using the RMS, it is also possible to use a combination of the local RMS and the negated minimum pixel. The formula used for this is max(window_rms, -1.5/5 x window_min) . This can be selected by adding ``-local-rms-method rms-with-min``. We have not seen much difference between using the normal local RMS or this combined quantity.

Using a pre-existing RMS map
----------------------------

A pre-existing RMS map can be supplied with the ``-rms-background-image`` parameter, like this:

.. code-block:: bash

    wsclean -size 1024 1024 -scale 10asec -local-rms-image rmsmap.fits \
      -auto-threshold 3 -mgain 0.8 -niter 1000000 [...] observation.ms
      
The input map (``rmsmap.fits`` in this case) will have to have the size of the (trimmed) output image.

A typical use-case for this is to supply an RMS map created by a source detector. Generally, to use this method, it is required to image the field twice. Also, when an automatically-created RMS map is used instead of a fixed pre-made map, WSClean will adapt the RMS map each major iteration, and in most cases the automatic RMS map will therefore be more accurate compared to a partially cleaned RMS map created by (for example) a source detector. The pre-existing RMS map should therefore normally not be the first choice, but there are probably cases where it *is* useful.
