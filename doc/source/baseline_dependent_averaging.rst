Baseline-dependent averaging
============================

This page explains baseline-dependent averaging, a feature that can increase the performance of gridding. This feature is available in WSClean since :doc:`version 2.0 <changelogs/v2.0>`.

Baseline-dependent averaging is a technique that allows more efficient averaging of data. In normal averaging, long baseline are averaged at the same level as small baselines. However, since small baselines can be averaged much more, this is inefficient. In baseline-dependent averaging, short baselines are averaged more than the long baselines.

By decreasing the number of visibilities, the imaging will become faster. By chosing the proper averaging factor, this can be quite significant with negligible effects on the image quality. The total time it will save is dependent on many factors ; the image size, the size of the w-terms, the number of visibilities, etc. If the gridding is the dominant cost, then baseline-dependent averaging is probably a significant improvement. In cases where the image size is very large and/or the w-terms are very large, it might be that the Fourier transforms are the dominant cost, instead of the gridding, and the savings might be less.

In principle, the baseline-dependent averaging can be done both in time and in frequency. However, WSClean can only do this in the time direction so far.

Baseline-dependent averaging in WSClean
---------------------------------------

WSClean can perform baseline-dependent time averaging during its "reordering" step. In this mode, it reads in a normal measurement set and averages into an internal format. The baseline-dependent averaging happens thus all internally, and other programs don't handle the averaged data. There's currently no good definition for how to store baseline-averaged data in Measurement Sets, so one cannot e.g. perform a self-calibration loop that images the baseline-averaged data, and then calibrate the baseline-averaged data.

Setting the averaging factor
----------------------------

Baseline-dependent average is enabled with the option "``-baseline-averaging <nwavelengths>``", which takes a parameter that is the number of wavelengths it is allowed to average over. This value can be related to the averaging factor as follows:

nwavelengths = max baseline in nwavelengths * 2pi * integration time in seconds / (24*60*60)

E.g. if a baseline of 20 kLambda can be averaged to 12 seconds without decorrelation, a corresponding baseline-dependent averaging factor can be calculated with:

20000 * 2pi * 12 / (24*60*60) = 17 wavelengths.

Using ``-baseline-averaging 17`` would therefore average any baseline of 20 kLambda to 12 seconds. Smaller baseline will be averaged more, and longer baselines will be averaged less.

Of course, the number of wavelengths to average over can also be calculated from first principles by determining how much phase-change in one averaging cell is acceptable.

For LOFAR, we've seen speed ups of a factor of 4-8 and for MWA data likewise.

Further options
---------------

Cleaning and major iterations etc. are all possible like normal.

Currently, option '``-baseline-averaging``' requires option '``-no-update-model-required``', because no un-averager (extrapolating) method to write the model data back at full resolution is available. This means this option can not be trivally used inside a :doc:`self-cal loop <selfcal>` where the model data is required. This can be overcome by performing the imaging with baseline averaging, and :doc:`predicting <prediction>` the final model image back into the measurement set with a separate WSClean run, that excludes baseline averaging.

This is a simple command to show the syntax of baseline-dependent averaging:

.. code-block:: bash

    wsclean -size 1024 1024 -scale 5asec -mgain 0.8 \
       -niter 10000 -threshold 0.002 \
       -baseline-averaging 16 -no-update-model-required \
       -maxuvw-m 100000 observation.ms

Note that baseline-dependent averaging currently only works together with :doc:`primary-beam correction <primary_beam_correction>` since :doc:`version 2.5 <changelogs/v2.5>`.

Some LOFAR specifics
--------------------

If a measurement set includes the international baselines, it can be helpful to remove these by a row-selection maximum uvw-m value. If one doesn't do this, the reordered file can still be very large because the international baselines won't be averaged at all. The ``-maxuvw-m`` option works before reordering (whereas ``-maxuv-l`` does not), and can thus remove those. Hence, adding a maxuvw-m that is just smaller than the first international (or unused) baseline can improve the speed. This even holds when those baselines are flagged.

Of course, even better is to explicitly remove the international stations from the measurement set when they are not used, e.g. using `DP3 <https://www.astron.nl/citt/DP3>`_.
