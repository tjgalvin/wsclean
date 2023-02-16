Implementing custom strategies in Python
========================================

WSClean has an option to call a Python script instead of using one of its built-in approaches to perform the deconvolution. The option is called ``-python-deconvolution`` and takes as parameter the filename of the Python script, e.g.:

.. code-block:: bash

    wsclean -python-deconvolution my_deconvolution_script.py \
      -niter 1000 -auto-threshold 5 -mgain 0.8 input.ms

The Python script should declare a function with the name and signature ``deconvolve(residual, model, psf, meta)``. The first three parameters are numpy arrays. Both ``residual`` and ``model`` are 4-dimensional arrays with dimensions ``n_channels * n_polarizations * height * width``. The ``psf`` is a 3-dimensional cube with dimensions ``n_channels * height * width``; it does not have the polarization dimension because the PSF is assumed to be the same for all polarizations. The ``meta`` parameter is a dictionary with various meta information about the data or run:

* ``meta.major_iter_threshold``: Current requested major iteration threshold.
* ``meta.final_threshold``: Requested absolute final threshold. If a relative stopping threshold was specified, it is converted to an absolute value before calling this function.
* ``meta.mgain``: Requested major loop gain.
* ``meta.iteration_number``: Total number of iterations performed. The deconvolution function is responsible for updating this value.
* ``meta.max_iterations``: Requested maximum number of iterations.
* ``meta.channels``: An array, where each element has the properties ``frequency`` and ``weights``.
* ``meta.square_joined_channels``: Set to ``True`` when the user has requested squared channels during joining.
* ``meta.spectral_fitter``: An object that can apply the user-requested spectral fitting (see :doc:`wideband_deconvolution`). It has two methods that each take a spectrum, which is a list of ``n_channels`` flux density values.
  - ``fit(spectrum)``: returns the fitted coefficients
  - ``fit_and_evaluate(spectrum)``: returns the fitted spectrum

There are a few example Python files in the ``wsclean/scripts/python-examples`` directory of the repository, notably:
* `simple-deconvolution-example.py <https://gitlab.com/aroffringa/wsclean/-/blob/master/scripts/python-examples/simple-deconvolution-example.py>`_, which is a single-frequency, single-polarization algorithm that is similar to a basic clean method. It demonstrates the interface (note that it is slow and not aimed at performance comparisons).
* `mf-deconvolution-example.py <https://gitlab.com/aroffringa/wsclean/-/blob/master/scripts/python-examples/mf-deconvolution-example.py>`_, a similar algorithm, but this one makes use of the spectral fitter to demonstrate multi-frequency deconvolution.
