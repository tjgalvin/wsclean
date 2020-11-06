Storing imaging weights
=======================

WSClean offers the option ``-store-imaging-weights`` (since :doc:`WSClean 2.5 <changelogs/v2.5>`) to store the imaging weights in a column named ``IMAGING_WEIGHT_SPECTRUM``. Its use is very similar to the ``IMAGING_WEIGHT`` column that CASA uses, but WSClean stores the imaging weights per polarization and channel, instead of only per channel as CASA does. The weights are written only when this parameter is given, and are not written by default (for reasons of performance).

Meaning of weight values
------------------------

The values in ``IMAGING_WEIGHT_SPECTRUM`` are relative unnormalized weights, and exclude the visibility weights that are stored in ``WEIGHT_SPECTRUM``. Hence, the desired weighting is achieved by multiplying the values in ``IMAGING_WEIGHT_SPECTRUM`` by values in ``WEIGHT_SPECTRUM``. To convert a sum of weighted values back to proper units, the weighted sum should be divided by the sum of weights. The ``IMAGING_WEIGHT_SPECTRUM`` values will be set to unity when requesting natural weighting (without tapers). Any kind of weighting, tapering and weighting rank-filter is taken into account in the weights.

Polarizations and writing
-------------------------

The values are written/updated during the first inversion. Currently, **imaging weight values are only updated when reordering is disabled** with ``-no-reorder``. This is hopefully fixed in a later version. Because no-reordering is often slower, if you need imaging weights the best approach is to only make a dirty image with the appropriate weighting.

Weights can only be stored when either imaging a one-to-one polarization in the measurement set, or when imaging Stokes I in a circular or linearly polarized measurement set. For example, in the case the measurement set contains XX/XY/YX/YY polarizations, one can request weights when imaging any of those polarizations, or when imaging Stokes I, but not when imaging stokes Q/U/V or LL/LR/.. polarizations. Polarizations that are not imaged will not be changed. Imaging Stokes I in a set with LL/RR or XX/YY polarizations will set the weights of the two polarizations (i.e. both LL/RR or both XX/YY, respectively).
