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
* `Boost <http://www.boost.org/>`_, used for threading, date and time calculations and some other general functionalities.
* `CFITSIO <http://heasarc.nasa.gov/fitsio/>`_, for reading and writing FITS files.
* `GSL <https://www.gnu.org/software/gsl/>`_, the GNU Scientific Library, used for certain computations.

WSClean uses C++11 features. Because of this, building WSClean with GCC requires at least GCC version 4.8.1.

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

Installed files
---------------

After succesfully running ``make install``, the following programs will be installed:

``wsclean``
    The main executable

``wsclean-mp``
    The main executable for distributed runs

``chgcentre``
    A program to change the phase centre of measurement sets, see :doc:`chgcentre`.
    
On Ubuntu and Debian
--------------------

Binary packages are available on Ubuntu and Debian, and can be installed with ``sudo apt-get install wsclean``. In case that version is new enough for your purpose, you're all done. If you want to compile WSClean from source, the following packages need to be installed:

.. code-block:: bash

    apt-get install  \
      casacore-dev libgsl-dev libhdf5-dev \
      libfftw3-dev libboost-dev \
      libboost-date-time-dev libboost-filesystem-dev \
      libboost-program-options-dev libboost-system-dev \
      libcfitsio-dev cmake g++

The LOFAR beam and IDG libraries are optional, but need to be installed manually from source if they are required (see elsewhere on this page).

On Red Hat
----------

The following document lists some instructions that can be helpful for instaling WSClean on Red Hat: [Installing WSClean on rhel 7.6](Installation instructions RHEL 7.6) (Thanks to Leonardo Saavedra from NRAO).

Linking errors
--------------

If you get undefined reference errors in casacore-code when you compile WSClean, e.g. similar to "undefined reference to ``casa::ArrayColumn<float>::get(unsigned int) const``", it probably means that you did not compile Casacore with C++11 support. When you compile Casacore, you need to turn this on explicitly: ::

    anoko@leopard:~/casacore/build$ cmake ../ -DCXX11="ON"

after which cmake should respond with a list of enabled features, including: ::

    -- C++11 support ......... = ON

Whether it is necessary to switch on C++11 depends on the version of the compiler -- with newer compilers it is turned on by default.

Compiling platform independently / portability
----------------------------------------------

By default, cmake will create a binary that is optimized for the machine that it is compiled on, and will only work on machines that contain the same instruction set as the compiling machine. In other words, the same binary might not work on other machines, because it might use advanced instructions such as AVX instructions that might not be available on another machine. It has been reported that this can e.g. lead to an "Illegal instruction" error (see [ticket 50](tickets:#50)).

If you want to make a binary that can be used on different platforms, you can add -DPORTABLE=True to cmake:

.. code-block:: bash

    cmake ../ -DPORTABLE=True
    
You should note that this makes the binary in general slower, so you should not use it unless you have to.

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

