Distributed imaging
===================

Multiple nodes can be used to speed up the imaging by running wsclean in distributed mode, which makes use of openmp to parallelize the imaging over multiple nodes. To use this, you need a recent WSClean version that provides the executable ``wsclean-mp``.

A typical run:

.. code-block:: bash

    mpirun --hostfile host_file -np 8 wsclean-mp -size 10000 10000 ...
    
And a host-file could look like this:

.. code-block:: text

    node100 slots=1
    node116 slots=1
    node118 slots=1
    node119 slots=1
    node122 slots=1
    node123 slots=1
    node124 slots=1
    node125 slots=1
    node126 slots=1
    node127 slots=1
    node128 slots=1
    node129 slots=1
    node130 slots=1
    
The host file should specify one slot per host, otherwise multiple wsclean's are executed on the same host, and that has (normally) no benefit, and would actually make those processes compete for memory and cpu. If you use wsclean-mp, all paths should be absolute for the hosts participating, and the paths should be reachable by all nodes (so ``-name`` should specify an absolute path name, and the input files should have an absolute path name).

``wsclean-mp`` will distribute the different channels to different nodes. This implies that if you don't use ``-channels-out`` there's no benefit, whereas using ``-channels-out 8`` with ``-np 8`` gives you a speed-up of 8. If multiple output channels are not necessary for your science goal, one can use ``-fit-spectral-pol 1 -deconvolution-channels 1``.
