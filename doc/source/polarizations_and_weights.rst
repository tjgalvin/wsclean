Polarizations and weights
=========================

Calculating Stokes I values
---------------------------

Forming Stokes I from linear/circular polarizations (XX/YY or LL/RR) is not straightforward when visibilities are weighted. When imaging Stokes I, WSClean combines XX and YY (or LL and RR). If either of the polarizations is flagged, the visibility is completely flagged. This means for example that an observation that has all LL visibilities flagged, will result in an empty Stokes I image. Some imagers "solve" this by being able to image "psuedo Stokes I", which is defined as the average over the unflagged polarizations. In case sources in the target image are polarized, this obviously leads to incorrect Stokes I values. This is often not such a concern, although with instrumental polarization leakage caused by fixed dipoles (LOFAR, MWA, etc.) this clearly is undesirable. WSClean therefore does not implement this directly, but can achieve the same effect: by imaging both XX and YY (or LL and RR) polarizations, and using the ``-join-polarizations`` option, you will be able to image and clean on the combination of the two polarizations. An example command:

.. code-block:: bash

    wsclean -join-polarizations -pol xx,yy -niter 1000000 \
        -auto-threshold 3 -mgain 0.8 myobservation.ms

Visibility weights
------------------

Visibility weights are the weights that are read from the ``WEIGHT_SPECTRUM`` column: each visibility has an associated weight value (do not confuse these with :doc:`image weights <image_weighting>`). When imaging one instrumental polarization, such as XX, just the XX weights are used, so no weights are averaged. When imaging Stokes I, the minimum value of XX and YY (or LL and RR) weights are used. A weight of 0 therefore prevents imaging the average I value altogether, as it should (because a missing XX value means that Stokes I has a completely undetermined value). It's not the statistically optimal, which would be:
    
.. math::

    w_I = \frac{2}{\frac{1}{w_{XX}} + \frac{1}{w_{YY}}},
    
but it is easier/faster, almost as good and has some implementation advantages.
