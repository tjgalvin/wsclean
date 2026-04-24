Facet-based imaging
===================

WSClean supports facet-based imaging with correction of direction-dependent effects.
A facet is a (polygonal-shaped) subsection of the full image. 

Facet-based imaging supports the correction of each facet with a direction-dependent gain from a h5parm solution file. Additionally, the facet can be corrected for the primary beam of the instrument, which is calculated using EveryBeam. In some cases, facet-based imaging can also speed up imaging or reduce the memory footprint.

Enabling faceting mode
----------------------

To enable facet-based imaging in WSClean, a file containing the facet definitions should be provided via the :code:`-facet-regions` option on the command line:

.. code-block:: text

    -facet-regions <MY_REGIONS_FILE>

in which :code:`[MY_REGIONS_FILE]` is a file containing the facet definitions in the DS9 region file format.
For a detailed explanation on the expected file format, see the explanation on :doc:`ds9_facet_file`. WSClean supports both convex and concave polygons in the region file.

Enabling shared reads
---------------------

When working with facets it is possible for data reads and some computation to be done once only per group of facets instead of once for each individual facets gridder.
The individual gridders can then share a single large chunk of visibility memory instead of each having their own, allowing for more visibilities to be chunked at a time when parallel gridding.

To enable shared reads faceting mode should be enabled and then the :code:`-shared-facet-reads` option set on the command line:
Shared reads are currently only supported for the ``wgridder`` gridder.

When working without corrections and with parallel gridding shared reads are generally faster than repeated reads. This is as a consequence of being able to fit more visibilities in memory as well as to some extent the reduced IO and computation.
When working with corrections, in order to facilitate the use of a shared memory buffer, it becomes necessary for the corrections of each facet to be applied "on the fly" during the gridding as the gridder accesses the visibilities.
This can lead to increased time spent computing corrections. These can partially, or even completely, negate the performance gains and lead to slower performance instead.
Whether shared reads will harm or aid performance in such scenarios is complex and depends on many factors. The two main factors that have been identified are data vs memory size and how dominant FFT computation is in the gridding performance.
The larger the ratio of `visibilities/memory` the higher the chances that shared reads will outperform repeated reads and the more dominant FFT computation is the more this becomes true.

When using shared reads it is important to carefully pick your other parallel gridding settings:
1. Ensure `-parallel-gridding` is a factor of `-j` to avoid unused cores that will be idle for the entire execution.
2. Consider the number of facets. `-parallel-gridding` should generally be equal or less to the number of facets, as otherwise a lot of cores that are assigned to parallel gridders may sit idle with no work to do.
3. `-parallel-gridding` should optimally also be a factor of the number of facets, or if not should have as large a remainder as possible when dividing, as otherwise there can be a slight increase in idle cores when dealing with the remaining facets/

Beam correction
~~~~~~~~~~~~~~~

Enabling the facet beam correction can be done with the option

.. code-block:: text

    -apply-facet-beam

The facet beam update interval (in seconds) can be defined by specifying:

.. code-block:: text

    -facet-beam-update <seconds>

The default value for the update interval is 120s.

Correction with gain solutions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Direction-dependent corrections per facet can be read and applied from an ``h5parm`` file, which is a HDF5 file with gain solutions in a particular format. This format is for example supported by tools like `DP3 <https://dp3.readthedocs.io/>`_ and `losoto <https://github.com/revoltek/losoto>`_. WSClean supports correction using various gain solution configurations, including scalar, diagonal and full-Jones solutions with amplitude, phase or combined terms. 

The ``h5parm`` file is specified via the command-line option:

.. code-block:: text

    -apply-facet-solutions <path-to-h5parm> <name1[,name2]>

where :code:`<path-to-h5parm>` is the path to the ``h5parm`` solution file and :code:`<name1[,name2]>` is a comma-separated list of strings specifying which "soltabs" from the provided ``h5parm`` solution file are used. Example names are :code:`amplitude000` and/or :code:`phase000`. 

In case multiple measurement sets are specified, it is possible to either specify one ``h5parm`` solution file, or a separate ``h5parm`` solution file per measurement set. The correction that should be applied (:code:`ampl000`, :code:`phase000`, or both) is required to be identical for all ``h5parm`` solution files. As an illustration, assume that :code:`N` measurement sets are passed to WSClean, with corresponding solution files :code:`h5parm1.h5, h5parm2.h5, ..., h5parmN.h5` containing a scalar amplitude correction. The syntax for applying the facet solution files on its corresponding measurement set thus becomes:

.. code-block:: text

    -apply-facet-solutions h5parm1.h5,h5parm2.h5,...,h5parmN.h5 ampl000

.. note::
    To find the matching direction in the solution file for the specified facets,
    the (RA, Dec) pointing of each facet is matched against the direction with
    the smallest (Euclidean) distance in the solution file.
    For further information on the (RA, Dec) pointing of a facet, see :doc:`ds9_facet_file`.

Order of beam and solution gains
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
WSClean supports applying direction-dependent full-Jones solutions. It is important though that if this is combined with applying the beam, it is applied after the beam (on the model / in the forward direction), whereas e.g. Faraday rotation is practically an effect that happens before the beam. This issue is minimized by applying the pointing centre beam on the observed data, but that will not be perfect. The rationale for this is a combination of things:

- For direction *independent* solving, the correct order of solutions can be maintained by applying an effect (like differential Faraday rotation) on the observed data when it needs to applied before the beam, and applying an effect (like clock or cable delay) to the model when it needs to be applied after the beam.
- Current solvers also solve after applying the beam on the model (because it is a lot easier), so even if WSClean would first apply the full-Jones solutions and then the beam, there would directly be a way to solve for this order. This effectively means that when solving for and applying direction-dependent differential Faraday rotation, what basically is solved for is more-or-less a solution projected by the beam, so not really Faraday rotation.
- Another use-case of full Jones solutions is to correct for beam leakage modeling errors, and for those it wouldn't matter much if you do it before or after the beam.
    
Beam output file and ``-pb.fits`` files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When either beam or gain solutions are applied, WSClean will output a "beam" fits file (per output channel), and each output image and model file will be accompanied by a ``-pb.fits`` file. The beam fits file represents the Stokes I response. It is a combination of the average beam and average gain solution corrections. It can therefore be used for weighting when mosaicking (optimal weighting should square the values). The ``-pb.fits`` file holds the fully corrected images that hold correct flux values, whereas the normal (non-``-pb.fits``) files contain "flat noise" images.

The average beam is corrected "smoothly", which (in perfect situations) means that the facet edge is not visible. In reality, the facet edge may be still visible in the ``-pb.fits`` and beam images because of the gain solutions, because the average correction is corrected discretely per facet and not smoothly.

Selecting visibilities to read / reorder
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When applying solutions or the beam in facet mode, WSClean will by default reorder and go through the four (instrumental) polarizations, e.g. XX/XY/YX/YY, for every requested output polarization. This is necessary to correct for leakage terms that the beam or solutions might have. This will obviously cause more reading and writing, while it might not always be necessary to use all visibilities. For the case where no leakage terms are expected, the option ``-diagonal-visibilities`` can be used. This will in the case of linear polarizations *only* use the diagonal visibilities (e.g. XX/YY) and assume the off-diagonal visibilities (e.g. XY/YX) are zero. When differences between the two diagonal correlations can also be ignored, option ``-scalar-visibilities`` can be used, which will reduce polarizations to the single requested polarization (e.g. Stokes I) and therefore limit I/O further.

.. note::
    The effect of using ``-diagonal-visibilities`` changed between WSClean versions :doc:`v3.4 <changelogs/v3.4>` and :doc:`v3.5 <changelogs/v3.5>`: the text above describes the new behaviour. An older name for ``-diagonal-visibilities`` was ``-diagonal-solutions``; these are synonyms, but the second is deprecated now.

Examples
--------
This is an example facet-based imaging command that applies both a facet-based beam correction and a scalar gain correction from an ``h5parm`` file:

.. code-block:: bash

    wsclean \
    -apply-facet-solutions mock_soltab_2pol.h5 ampl000,phase000 \
    -facet-regions ds9.reg \
    -apply-facet-beam \
    -facet-beam-update 120 \
    -niter 1000000 -auto-threshold 5 -mgain 0.8 \
    -size 1024 1024 -scale 1amin \
    ${ms}

Availability
------------
Initial support for faceting is made available in WSClean :doc:`version 3.0 <changelogs/v3.0>`. In subsequent versions,
several bugs were fixed and support for different solution types was added. WSClean :doc:`version 3.4 <changelogs/v3.4>`
has support for scalar and diagonal solutions, and is considered stable.

Facet-based imaging in conjunction with the Image Domain Gridder (IDG) is only possible without applying DDEs.
