Polarimetric deconvolution
==========================

WSClean can perform deconvolution by combining polarizations. Similar to :doc:`joining channels <wideband_deconvolution>`, this implies that peaks (and component shape) are searched in the integrated values over the selected polarizations, but the components strengths are determined from each polarization once its location is selected.

This method improves the quality of polarized images. It can also be used to perform self-cal to keep the polarized flux accurate, although it requires some flux to be present in the polarized images (instrumental leakage). In general, Stokes I quality is not improved by using this method. The method increases computational and memory cost.

Here are two examples for using joined polarization cleaning:

.. code-block:: bash

    wsclean -pol xx,yy -join-polarizations -threshold 1 -mgain 0.8 -niter 10000 \
    -scale 15amin -size 1024 1024 -weight briggs 0 galaxy.ms

    wsclean -pol iquv -join-polarizations -channels-out 8 -join-channels \
    -fits-mask 3c196-mask.fits -threshold 1 -mgain 0.8 -niter 10000 \
    -scale 5asec -size 6000 6000 3c196.ms

Joined-polarization cleaning works in combination with :doc:`multi-scale cleaning <multiscale_cleaning>`, :doc:`multi-frequency deconvolution <wideband_deconvolution>` and :doc:`masked cleaning <masking>`.

If you are looking to do rotation measure synthesis with WSClean, have a look at the :doc:`RM-synthesis manual page <rm_synthesis>` that describes some particulars for this.

Linking polarizations
---------------------

Since :doc:`WSClean 2.6 <changelogs/v2.6>` it is possible to "link" polarizations. This has the effect of selecting a subset of polarizations for cleaning, but then subtracting components from all imaged polarizations. This results in that the subset of polarizations are cleaned as if they are the only polarizations being imaged, while the others polarizations are still saved, and sources that appear in the subset of polarizations are also cleaned in the other polarizations. Effectively, structure that appears mostly in the linked polarization(s) will be cleaned better. Hence, this is particularly useful for cleaning Stokes I or the XX,YY polarizations.

Note that this is not the recommended setting for RM synthesis; see above.

Linking polarizations is turned on with the ``-link-polarizations`` option, after which the set of polarizations (or single polarization) should be specified similar to the syntax for ``-pol``. An example:

.. code-block:: bash

    wsclean -pol xx,xy,yx,yy -link-polarizations xx,yy -auto-threshold 3 -mgain 0.8 \
    -niter 10000 -scale 1amin -size 4096 4096 observation.ms

