FROM kernsuite/base:4

RUN export DEBIAN_FRONTEND=noninteractive && apt-get update && \
    apt-get install -y \
    cmake \
    build-essential \
    pkg-config \
    casacore-data casacore-dev \
    libblas-dev liblapack-dev \
    libboost-date-time-dev libboost-numpy-dev libboost-python-dev libboost-test-dev \
      libboost-system-dev libboost-filesystem-dev libboost-program-options-dev \
    libcfitsio-dev \
    libfftw3-dev \
    libgsl-dev \
    libhdf5-dev \
    libopenmpi-dev \
    libpython-dev \

ADD . /src
WORKDIR /src

RUN mkdir /build && cd /build && cmake ../src

RUN cd /build && make -j24 && make install && make check -j24
