*w*-stacking interface
======================

WSClean's gridder is inside a single C++ unit file and is written to be modular. Therefore, it is not hard to integrate and reuse it in other software. The gridder can do both prediction and inversion.

The gridder is only a small part of the WSClean imager, and handles the low-level placing of visibilities on the uv-grid and the involved FFTs (inversion), or the other way around; predicting visibilities from a model image. There are also higher-level interfaces to WSClean, which support calling the whole WSClean program, including cleaning and handling measurement sets etc. If you need cleaning and/or your data is in the measurement set format, it might be useful to use other interfaces instead. However, if you would like to use the low-level interface, the following information might be useful.

Since WSClean 1.9, the gridder is contained in the unit file ``wsclean/wstackinggridder.cpp`` and its header ``wsclean/wstackinggridder.h``. The header has extensive code documentation that should give enough information about how to call the gridder. You can create the Doxygen documentation for this class with a "make doc" within the build dir (which places the info in build/doc/html/index.html, etc). The code documentation can also be found here:

* `WStackingGridder class API <http://www.andreoffringa.org/wsclean/doxygen/classWStackingGridder.html>`_

The gridder can be compiled with two external libraries: FFTW and Boost. To avoid a Casacore dependency, you need to define AVOID_CASACORE while compiling the gridder. There's a prediction example in wsclean/examples called "`wspredictionexample.cpp <https://gitlab.com/aroffringa/wsclean/-/blob/development/wsclean/examples/wspredictionexample.cpp>`_)", which can be compiled with:

.. code-block:: bash
    
    g++ -o wspredictionexample -O3 -march=native -pthread \
      -I ../../external/aocommon/include/ \
      -DAVOID_CASACORE \
      wspredictionexample.cpp \
      ../wstackinggridder.cpp \
      ../../io/logger.cpp \
      ../../structures/image.cpp \
      -lfftw3f -lfftw3 -lfftw3f_threads -lfftw3_threads \
      -lboost_date_time -lboost_system

This assumes you are in the ``wsclean/examples`` subdirectory. There's a separate makefile for this: `wsclean/examples/Makefile`.

The prediction example shows how to initialize the gridder and get visibility samples from it and includes inline documentation, so should be a good starting point for using the gridder.
