.. toctree::
   :maxdepth: 2
   :hidden:

   changelogs/list.rst

Installation instructions
=========================

Getting the source code
-----------------------

Unless you need specific new features, it is recommended to use a stable (tagged) version of WSClean. These can be downloaded from https://gitlab.com/aroffringa/wsclean/-/releases.

To retrieve the (experimental) master branch of WSClean, use git:

.. code-block:: bash

    git clone -b master git@gitlab.com:aroffringa/wsclean.git

This will retrieve the *master* branch of WSClean.

Manual compilation
------------------

After downloading the source code, one will need to compile WSClean. WSClean requires:

* `Casacore <https://github.com/casacore/casacore>`_, for opening measurement sets. Version >=2.0 is required, not lower. Casacore is required even if you already have Casa installed. Casacore needs to be compiled with C++11 support, which is the default for the latest version.
* `FFTW <http://www.fftw.org/>`_ version 3.3.5 or newer, used to perform Fourier transformations.
* `Boost <http://www.boost.org/>`_, used for date and time calculations and some other general functionalities.
* `CFITSIO <http://heasarc.nasa.gov/fitsio/>`_, for reading and writing FITS files.
* `GSL <https://www.gnu.org/software/gsl/>`_, the GNU Scientific Library, used for certain computations.
* `HDF5 <https://www.hdfgroup.org/solutions/hdf5/>`_, for reading solution files
* Python 3 libraries
* BLAS and LAPACK libraries

WSClean uses C++17 features. Because of this, building WSClean with GCC requires at least GCC version 8.

To use the :doc:`image-domain gridder <image_domain_gridding>` (a fast GPU gridder), you will need to install the IDG libraries from https://gitlab.com/astron-idg/idg. To apply primary beams, the EveryBeam package is required from https://git.astron.nl/RD/EveryBeam. To use the distributed mode, `OpenMPI <https://www.open-mpi.org/>`_ is required.

After installing these dependencies, compile WSClean with the following commands:

.. code-block:: bash

    mkdir -p build
    cd build
    cmake ../
    make -j 4
    sudo make install

If cmake reports errors, you might have to install or specify your libraries in the call to cmake if they are not in standard directories, e.g.:

.. code-block:: bash

    cmake ../ -DCMAKE_PREFIX_PATH=/opt/cep/casacore/

to add that directory to the search path. To add multiple directories to the search path, put the paths between double quotes and separate them with a semicolon:

.. code-block:: bash

    cmake ../ -DCMAKE_PREFIX_PATH="/path/to/casacore;/path/to/cfitsio"

On Ubuntu and Debian
--------------------

Binary packages are available on Ubuntu and Debian, and can be installed with ``sudo apt install wsclean``. In case that version is new enough for your purpose, you're all done. Be aware that this provides a WSClean version that is not optimized for your platform (see next section). If you want to compile WSClean from source, the following packages need to be installed:

.. code-block:: bash

    apt install  \
      casacore-dev cmake g++ \
      git pkg-config libblas-dev \
      libboost-filesystem-dev libboost-program-options-dev \
      libcfitsio-dev libfftw3-dev libgsl-dev \
      libhdf5-dev liblapack-dev libopenmpi-dev \
      libpython3-dev pkg-config

The LOFAR beam and IDG libraries are optional, but need to be installed manually from source if they are required (see elsewhere on this page).

Compiling platform independently / portability
----------------------------------------------

By default, cmake will create a binary that is optimized for the machine that it is compiled on, and will only work on machines that contain the same instruction set as the compiling machine. In other words, the same binary might not work on other machines, because it can use advanced instructions such as AVX-512 instructions that may not be available. It has been reported that this can e.g. lead to an "Illegal instruction" error. A common use-case where you can run into this issue is when working with Docker containers: when compiling the container on a different machine as where it runs, this can happen.

To make a binary that can be used on multiple platforms, there are two solutions. The recommended method is to set the ``TARGET_CPU`` variable in cmake to a reasonable base architecture that is supported by all machines that you want to run WSClean on. It is for example generally safe to set it to ``haswell``, available since 2013, as almost all machines currently available support this architecture. Moreover, the ``haswell`` target includes AVX, AVX2 and FMA extensions, which can have a significant performance improvement (see e.g. `Wikipedia's AVX page <https://en.wikipedia.org/wiki/Advanced_Vector_Extensions>`_). A slightly newer target would be ``skylake``, which is available since 2015. A full list of supported machines can be found in `the gcc manual <https://gcc.gnu.org/onlinedocs/gcc/Submodel-Options.html>`_. To compile WSClean with skylake, you would run cmake like this:

.. code-block:: bash

    cmake ../ -DTARGET_CPU=skylake

A more rigorous option is to set ``-DPORTABLE=True`` in cmake. This disables any extra instruction sets, and by that it basically uses only basic instructions that have been available since the beginning of the architecture. This may make the binary considerably slower, so should not be used unless strictly necessary (e.g. for compiling generic packages).

On Red Hat
----------

The following document lists some instructions that can be helpful for instaling WSClean on Red Hat and CentOS: :doc:`Installing WSClean on Red Hat and CentOS <installation-rhel>`.

Using the LOFAR/MWA/... beam
----------------------------

The latest versions of WSClean require the `EveryBeam package <https://git.astron.nl/RD/EveryBeam>`_ to apply the LOFAR beam and other known beams (VLA, MWA, LWA, ATCA, ...). Older versions used the LOFAR beam from the LOFAR repository or (since :doc:`version 2.6 <changelogs/v2.6>`) the `"LOFARBeam" package from Github <https://github.com/lofar-astron/LOFARBeam>`_ instead.

To have these primary beams available in WSClean, CMake needs to find the EveryBeam installation on your computer. If you have installed EveryBeam in a custom directory, you can add it to your search path. For example, if EveryBeam been installed to ``~/Software/EveryBeam-install``, this cmake command will use it:

.. code-block:: bash

    cmake ../ -DCMAKE_PREFIX_PATH="~/Software/EveryBeam-install/"

CMake will tell whether the LOFAR tree was found:

.. code-block:: bash

    cmake ../ -DCMAKE_PREFIX_PATH="~/Software/EveryBeam-install/"
    [..]
    EveryBeam beam library found.

Extra paths can be added to e.g. also include IDG. Paths can be separated with a ; (semicolon).

The MWA beam
------------

The MWA beam requires that a ``.h5`` file with beam coefficients is present in your path, which you can find in the MWA software tools. This file needs to be in your Python search path. WSClean will iterate over your Python search paths by executing the following Python program:

.. code-block:: python

    from __future__ import print_function
    import sys
    for a in sys.path:
        print(a)

.. toctree::
   :maxdepth: 2
   :hidden:

   installation-rhel

Installed files
---------------

After succesfully running ``make install``, the following programs will be installed:

``wsclean``
    The main executable

``wsclean-mp``
    The main executable for distributed runs

``chgcentre``
    A program to change the phase centre of measurement sets, see :doc:`chgcentre`.
