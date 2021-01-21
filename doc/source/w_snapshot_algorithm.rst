*w*-snapshot algorithm
======================

WSClean provides an algorithm that has approximately the same performance as the W-snapshot algorithm suggested by `Cornwell et al. (2012) <https://arxiv.org/abs/1207.5861>`_, although the algorithm in WSClean is slightly different. In WSClean, it consists of phase rotating the visibilities to zenith, and then recentring the image during w-stacking. Mathematical details are explained in `the WSClean paper <http://arxiv.org/abs/1407.1943>`_. It seems to be worthwhile for MWA snapshots of a few minutes with images of 3072 x 3072 pixels at zenith angles > 20 degree, and provides a speed-up of about a factor of 3 at zenith angles > 45 degrees. For larger images the speed-up will be greater. It is not useful for long integrations.

To use this method, the measurement set should be phase rotated to the phase centre that you want to image (as it normally would be). This chapter assumes this is ``08h20m00.0s -42d45m00s`` (which is Puppis A). To prepare for the *w*-snapshotting approach, first run :doc:`the chgcentre tool <chgcentre>` with the following parameters:

.. code-block:: bash

    chgcentre -minw copy.ms

This will calculate the optimal projection direction and perform the required phase rotations. Finally, during imaging the phases should be *shifted* to the original target, e.g.:

.. code-block:: bash

    wsclean -shift 08h20m00.0s -42d45m00s ...

The ``chgcentre`` command *rotates* the coordinate system whereas WSClean's ``-shift`` parameter *shifts* the coordinate system along the tangent plane. The result is that imaging is done in a different projection that is still centered on the target location, but that has smaller *w*-terms. The :doc:`chgcentre page <chgcentre>` provides further information about the ``chgcentre`` tool. 

The ``-shift`` parameter is new in  :doc:`WSClean 2.11 <changelogs/v2.11>`. Before 2.11, this could be achieved by using the ``-shiftback`` option of ``chgcentre``. That approach is no longer supported in WSClean.

The :doc:`idg <image_domain_gridding>`, :doc:`w-gridding <wgridding>` gridder and *w*-stacking gridders all support this approach. Support for the :doc:`w-gridding <wgridding>` gridder was added in WSClean 2.11.
