Facet-based imaging
===================

WSClean provides experimental support for facet-based imaging, and applying direction-dependent effects (DDEs) per facet.
A facet should be understood as a (polygonal shaped) subdomain of the full image, where both convex and concave polygons are supported.
A first advantage of facet-based imaging is that DDEs can be applied per facet. Parallelization and reduction of the memory footprint
can be other other advantages of facet-based imaging, because visibility gridding and predicting can largely be done for each facet independently, rather than requiring the full image.

Availability
------------
Facetting is available in WSClean :doc:`version 3.0 <changelogs/v3.0>` and later.

Command-line options
--------------------

To enable facet-based imaging in WSClean, a file containing the facet definitions should be provided via the :code:`-facet-regions` option on the command line:

.. code-block:: text

    -facet-regions <MY_REGIONS_FILE>

in which :code:`[MY_REGIONS_FILE]` is a file the contains the facet definitions using the DS9 region file format.

Enabling the facet beam correction can be done with the option

.. code-block:: text

    -apply-facet-beam

The facet beam update interval (in seconds) can be defined by specifying:

.. code-block:: text

    -facet-beam-update <seconds>

The default value for the update interval is 120s.

Direction-dependent corrections per facet can also be read and applied from an H5Parm file - which in essence is a HDF5 file with some prescribed lay-out.
This is done via the command-line option:

.. code-block:: text

    -apply-facet-solutions <path-to-h5parm> <name1[,name2]>

where :code:`<path-to-h5parm>` is the path to the H5Parm file and :code:`<name1[,name2]>`
is a comma-separated list of strings specifying which "soltabs" from the provided H5Parm file are used.
Acceptable names are :code:`ampl000` and/or :code:`phase000`.

**NOTE:** the ordering of the directions in the H5Parm file is currently **assumed** equal to the ordering of the facets in the DS9 region file.

Example command
---------------
An example facet-based imaging command in WSClean, applying both a facet-based beam correction as well as a gain correction from an H5Parm file could be:

.. code-block:: bash

    wsclean \
    -apply-facet-solutions mock_soltab_2pol.h5 ampl000,phase000 \
    -facet-regions ds9.reg \
    -apply-facet-beam \
    -facet-beam-update 120 \
    -niter 1000000 -auto-threshold 5 -mgain 0.8 \
    -size 1024 1024 -scale 1amin \
    ${ms}


Caveats
-------

Facet-based imaging is currently an experimental feature, and therefore should be used with care.
A (probably non-exhaustive) list of caveats is presented below:

- Parallel processing can be enabled with the :code:`-parallel-gridding` option (multi-threaded) or the :code:`wsclean-mp` (using MPI). Parallelization over facets is however barely tested, and may return unexpected errors or results, in particular when applying DDEs.
- Facet-based imaging in conjunction with the Image Domain Gridder (IDG) is only possible without applying DDEs.
- When applying solutions in WSClean for facetted imaging, only scalar solutions are currently applicable.
- Be aware of the (direction) ordering restriction when applying DD gains from an H5Parm file as mentioned above.
