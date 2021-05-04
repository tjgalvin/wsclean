Installation of WSClean on Red Hat-based systems
================================================

This chapter contains some instructions for Red Hat based systems.

CentOS 8
~~~~~~~~

The source repository contains a Docker file `scripts/docker/Centos8 <https://gitlab.com/aroffringa/wsclean/-/blob/master/scripts/docker/Centos8>`_. This is a script to build a Docker container based on CentOS 8. These same steps can be used as a reference when building the latest WSClean on a (non-virtual) CentOS system. Other Red Hat systems should be very similar.

WSClean 2.7 on Red Hat 7.6
~~~~~~~~~~~~~~~~~~~~~~~~~~

The text below was written by Leonardo Saavedra from NRAO.

.. note::
   Be aware that these instructions are not for the latest WSClean version.

Some WSClean dependencies are provided by RHEL 7.6, but you have to install the latest fftw and casacore.

* http://www.fftw.org/
* https://sourceforge.net/p/wsclean
* https://github.com/casacore/casacore

Download the packages
---------------------

.. code-block:: bash

  mkdir ~/wsclean
  cd wsclean
  wget -c http://www.fftw.org/fftw-3.3.8.tar.gz
  wget -c https://github.com/casacore/casacore/archive/v3.1.1.tar.gz
  wget -c https://sourceforge.net/projects/wsclean/files/wsclean-2.7/wsclean-2.7.tar.bz2

.. note::
    Be aware that these instructions refer to the old SourceForge version of WSClean.
    It is highly recommended to use newer versions from GitLab.

For this example I am going to install under ``/export/local``

Install fftw
------------

.. code-block:: bash

  cd ~/wsclean
  tar xzvf fftw-3.3.8.tar.gz
  cd fftw-3.3.8/
  ./configure --prefix=/export/local  --enable-threads --enable-openmp  --enable-shared
  make
  make install

Install Casacore
----------------

.. code-block:: bash

  export LD_LIBRARY_PATH=/export/local/lib:$LD_LIBRARY_PATH
  cd ~/wsclean
  tar xzvf v3.1.1.tar.gz
  cd casacore-3.1.1/
  mkdir build
  cd build
  cmake ../ -DCMAKE_PREFIX_PATH=/export/local/
  make -j `nproc`
  vim cmake_install.cmake <-- modified CMAKE_INSTALL_PREFIX
  make install

Install WSClean 2.7
-------------------

.. code-block:: bash

  cd ~/wsclean
  tar xvfj wsclean-2.7.tar.bz2
  cd wsclean-2.7/
  mkdir build
  cd build/
  cmake ../ -DCMAKE_PREFIX_PATH=/export/local/
  make -j `nproc`
  vim cmake_install.cmake <-- modified CMAKE_INSTALL_PREFIX
  make install

Check WSClean
-------------

.. code-block:: bash

  pwd
  /export/local/bin
  ./wsclean -version

  WSClean version 2.7.0 (2019-04-19)
  This software package is released under the GPL version 3.
  Author: André Offringa (offringa@gmail.com).

