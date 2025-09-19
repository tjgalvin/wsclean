Sky-model rendering with sub-pixel accuracy
===========================================

WSClean can be used to render sky models with sub-pixel accuracy, enabled by the ``-draw-model <sky model>`` option.

The main use-case is enabling image-based visibility prediction from a sky model; see the :doc:`chapter on prediction <prediction>`.

Sub-pixel rendering
-------------------

The sub-pixel rendering approach implemented in WSClean is a combination of two methods, one for rendering point sources and another for Gaussian sources. Point sources are convolved with a sinc function, which allows them to be highly accurately placed on non-integer pixel positions. Gaussian sources are drawn by sampling the sky model in the Fourier domain, whereby a phase-shift places the source on a non-integer pixel position.

Using the sub-pixel renderer
----------------------------

To use the sub-pixel renderer in WSClean, a file containing the component list should be provided via the :code:`-draw-model` option on the command line:

.. code-block:: bash

    -draw-model <MY_SKY_MODEL>

in which :code:`<MY_SKY_MODEL>` is a file containing the model components in DP3's sky model format. For a detailed explanation on the expected file format, see the :doc:`Component Lists <component_list>` page.

Examples
--------
This is an example of a sub-pixel renderer command that draws model images of a sky model ``skymodel.txt`` at 150 MHz:

.. code-block:: bash

    wsclean -draw-model skymodel.txt -draw-frequencies 150e6 1e6 \
      -draw-spectral-terms 2 -sinc-window-size 256 \
      -size 1024 1024 -scale 1arcsec \
      -name skymodel observation.ms

This command creates two images, each representing one spectral term, i.e., the coefficients of the polynomial representing the spectral energy density function. In this example, we create images for 2 spectral terms,,the first representing the Stokes I flux and the other being a higher-order term.

In the example above, the sky model is drawn at the phase-centre read from the MeasurementSet. However, it is also possible to render a model by manually providing coordinates:

.. code-block:: bash

    wsclean -draw-model skymodel.txt -draw-frequencies 150e6 1e6 \
      -draw-centre 19h59m28.35s -40d44m2.09s -draw-spectral-terms 2 \
      -sinc-window-size 256 -size 1024 1024 -scale 1arcsec -name skymodel

Using DP3 to predict the rendered sky models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

At the moment of writing, there are two `DP3 <https://dp3.readthedocs.io/>`_ steps that predict visibility based on input images: `wgridderpredict <https://dp3.readthedocs.io/en/latest/steps/WGridderPredict.html>`_ and `idgpredict <https://dp3.readthedocs.io/en/latest/steps/IDGPredict.html>`_; we use the former in this example.

Provided that we have drawn the two term images using (one of) the aforementioned WSClean commands, we should have two FITS files: ``skymodel-term-0.fits`` and ``skymodel-term-1.fits``. These can be used to run image-based prediction with DP3, for example:

.. code-block:: bash

    DP3 msin=observation.ms msout=. msout.datacolumn=MODEL_DATA steps=[wgridderpredict] \
        wgridderpredict.images=[skymodel-term-0.fits, skymodel-term-1.fits] \
        wgridderpredict.regions=facets.reg wgridderpredict.sumfacets=True

In this example, the visibilities for all facet in the prodived region file (``facets.reg``) are summed and stored in the ``MODEL_DATA`` column.

It is worth noting that the predicted model is frequency dependent, provided that more than one term image has been drawn. This is distinct from WSClean's `-predict` option, which only predicts a single term, thus requires multiple ``-draw-model`` calls to obtain a frequency dependent model.
