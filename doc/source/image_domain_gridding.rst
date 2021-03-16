Image domain gridding
=====================

The `image-domain gridder (IDG) <https://gitlab.com/astron-idg/idg>`_ is a new, fast gridder that makes *w*-term correction and *a*-term correction computationally very cheap. It performs extremely well on gpus.

To use IDG, the `IDG library <https://gitlab.com/astron-idg/idg>`_ needs to be available during compilation. :doc:`WSClean 2.5 <changelogs/v2.5>` is compatible with IDG version 0.2. See :doc:`installation <installation>` for more information. 

Considerations
--------------

IDG supports the same cleaning and data selections options that WSClean offers in normal mode (without IDG -- so e.g. Cotton-Schwab, multi-frequency multi-scale cleaning, auto-masking, etc.). However there are a few points which require a bit of care to use IDG:

 * :doc:`Baseline-dependent averaging <baseline_dependent_averaging>` does not work together with IDG and should not be used together.
 * Imaging multiple polarizations independently (without joining polarizations) is not possible with IDG, see the discussion on polarization settings below.
 * IDG can only be used on measurement sets that have four polarizations (e.g. XX/XY/YX/YY).

Command-line options
--------------------

There are several IDG related WSClean parameters:

.. code-block:: text

    -use-idg
    
This option turns on the IDG gridder. If IDG is not available, WSClean will show an error.  

.. code-block:: text

    -idg-mode [cpu/gpu/hybrid]
    
IDG mode to use. The default is CPU, which works on any computer. GPU and Hybrid both make use of the GPU and require the availability of GPU libraries when compiling IDG (see the IDG manual). The "pure" GPU mode is however **not recommended** to be used, as it does not have all the features that the hybrid mode has, hence **hybrid is the recommended method when GPUs are available**.

Gridding with the beam
----------------------

IDG allows gridding with a time-variable beam -- currently (as of 2.6), this includes the LOFAR, AARTFAAC and MWA beam. This would increase the accuracy of with which sources are deconvolved (in theory, but it depends of course on the accuracy of the beam model). The option to do this is ``-grid-with-beam``, and optionally, if the differential beam needs to be applied (which is common; see the LOFAR beam page for more info) one would also add ``-use-differential-lofar-beam``.

In summary; IDG with full beam correction:

.. code-block:: text

    wsclean -use-idg -grid-with-beam ...
    
IDG with differential beam correction:    

.. code-block:: text

    wsclean -use-idg -grid-with-beam -use-differential-lofar-beam ...

For more information about the differential lofar beam, see :doc:`primary beam correction <primary_beam_correction>`. Note that you do not have to add ``-apply-primary-beam``: if you don't specify this parameter, WSClean will grid with the beam and output both the "flat noise" image and the pb corrected image.

NB: In earlier versions, specifying both ``-apply-primary-beam`` and ``-use-idg`` would apply the wrong beam!

Polarization settings for IDG
-----------------------------

IDG always images all Stokes parameters. Therefore, it makes sense to set the polarizations to iquv (with ``-pol iquv``) when using the IDG, as otherwise the other polarizations are imaged and then thrown away. IDG does however require joined or linked polarization cleaning: it is not possible to clean the polarizations individually with IDG in a single run. If this is really desired it can be done by running wsclean 4 times, imaging one polarization at a time, however most science cases will want to use one of the other options described below.

There are several distinct cases of interest:

* Total intensity science, i.e., not interested in QUV: only image I with ``-pol I``.
* Total intensity science, but "nice to have QUV" for e.g. sensivity analysis (V is useful for that): image all polarizations with ``-pol IQUV`` and *link* the deconvolution on polarization I with ``-link-polarization i``.
* Interested in rotation measure synthesis, not directly in I and V but nice to have: image all polarizations with ``-pol IQUV`` and link the deconvolution on polarization Q and U with ``-link-polarization qu`` (possibly with squared deconvolution etc., see [polarized cleaning page](PolarizedCleaning)).
* Interested in all stokes parameter, cleaning each polarization in a joined way: ``-pol IQUV -join-polarizations``.

Some further explanation on this: IDG in combination with ``-pol iquv -join-polarizations`` can be the best choice in some cases, but if one is only interested in Stokes I, it has the downside that QUV are involved in the cleaning of I, and in case these are mostly empty, the net effect is that the noise in the clean-component finding step is increased slightly. Since [WSClean 2.6](Changelog-2.6) it is possible to "link" polarizations. Its use is explained on the [polarized cleaning page](PolarizedCleaning). To overcome the downside described above, the recommended setting for imaging with IDG is to use polarization linking on I, similar to:

.. code-block:: bash

    wsclean -link-polarization i -use-idg [..]

This will clean Stokes I fully, and clean the components found in I also from the other polarizations. However, it will *not* clean structure from QUV that is not found in I.

Performance
-----------

On my 4-core home machine, the CPU version is not as fast as the wstacking gridder. However, it can apply *a*-terms and uses less memory, which can be advantageous in some cases. In contrast to the CPU gridder, the GPU gridder *is* considerably faster (up to an order of magnitude), but only for large images (>6k or so).

The exact performance benefits depend heavily on image size, bandwidth and number of cores available. Hence if performance is important, I recommend to make a careful comparison for the particular test case you are interested in.

IDG has been tested and shown to work on MWA data. It performs well, but does require a lot of memory, caused by the wide field-of-view of the MWA and therefore high w-terms.

Advanced *a*-term corrections
-----------------------------

WSClean+IDG allows a combination of several direction-dependent corrections to be applied, including TEC screens, diagonal gain correction and position shifts ("dldm screens"). These are discussed on the [a-term correction page](ATermCorrection).

Information for older IDG versions
----------------------------------

Before :doc:`WSClean version 2.9 <changelogs/v2.9>`, IDG's memory usage was highly dependant on the number of channels in the set. IDG can be made to use considerably less memory by splitting the bandwidth using ``-channels-out`` and wide-band deconvolution (see :doc:`making image cube <making_image_cubes>` and :doc:`wideband deconvolution <wideband_deconvolution>`). For example, splitting the bandwidth in 4 output channels has allowed imaging one of the MWA sets on a 32 GB machine:

  wsclean -grid-with-beam -beam-aterm-update 10 -channels-out 4 -join-channels -link-polarizations i -use-idg -size 1536 1536 -scale 1amin -niter 1000000 -auto-threshold 0.5 -auto-mask 4 -multiscale -mgain 0.8 observation.ms
  
This should no longer be necessary for WSClean 2.9 and later. In those versions, IDG should honour the requested memory settings and available memory. If you do expect memory issues, you can tweak the memory usage using the ``-mem`` and ``-absmem`` parameters of WSClean.
