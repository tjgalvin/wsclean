w-snapshot algorithm
====================

WSClean provides an algorithm that has approximately the same performance as the W-snapshot algorithm suggested by `Cornwell et al. (2012) <https://arxiv.org/abs/1207.5861>`_, although the algorithm in WSClean is slightly different. In WSClean, it consists of phase rotating the visibilities to zenith, and then recentring the image during w-stacking. Mathematical details are explained in `the WSClean paper <http://arxiv.org/abs/1407.1943>`_. It seems to be worthwhile for MWA snapshots of a few minutes with images of 3072 x 3072 pixels at zenith angles > 20 degree, and provides a speed-up of about a factor of 3 at zenith angles > 45 degrees. For larger images the speed-up will be greater.

To use this method, start by making a copy of your measurement set. This measurement set should be phased up to the phase centre that you want to image (as it normally would be). Then, run :doc:`the chgcentre tool <chgcentre>` with the following parameters:

.. code-block:: bash

    chgcentre -minw -shiftback copy.ms

This will calculate the optimal projection direction, perform required phase rotations and shifts, and add some keywords to the measurement set. WSClean will recognize the presence of those keywords and perform the recentring, so WSClean can be run in the normal way. Since the visibilities are rewritten by chgcentre in a non-standard projection, you cannot use that measurement set with most other tools anymore (that's why it should be copied).

:doc:`The chgcentre page <chgcentre>` provides further information about the chgcentre tool.
