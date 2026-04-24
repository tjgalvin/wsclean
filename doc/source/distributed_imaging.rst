Distributed imaging
===================

Multiple nodes can be used to speed up the imaging by running wsclean in
distributed mode, which makes use of MPI to parallelize the imaging over
multiple nodes. If MPI is found during compilation, the compiler will produce
a ``wsclean-mp`` binary which can exploit parallelism across multiple compute
nodes.  This executable supports the same command-line options as
the regular ``wsclean`` executable, but will instead perform
parallel gridding using multiple MPI processes.

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
    
The host file should specify one slot per host, otherwise multiple wsclean's are executed on the same host, and that has (normally) no benefit, and would actually make those processes compete for memory and cpu. If you use ``wsclean-mp``, all paths should be absolute for the hosts participating, and the paths should be reachable by all nodes (so ``-name`` should specify an absolute path name, and the input files should have an absolute path name).

``wsclean-mp`` will distribute the different channels to different nodes. This implies that if you don't use ``-channels-out`` there's no benefit, whereas using ``-channels-out 8`` with ``-np 8`` gives you a speed-up of 8. If multiple output channels are not necessary for your science goal, one can use ``-fit-spectral-pol 1 -deconvolution-channels 1``.

The multi-threading arguments (as described in the
:doc:`parallelization section <parallelization>`)
apply to each MPI process.
For example, when using ``-j 4``, each MPI process will use 4 threads.

Note that WSClean only uses MPI during gridding. Other parts, such as
deconvolution and reordering, only use the main MPI process and therefore only
use multi-threading.

When using a single compute node, using MPI is discouraged, since multi-threaded
gridding using the ``-parallel-gridding`` is more efficient. When using multiple
nodes, the best performance is normally achieved using one MPI process with
multiple parallel gridders per node. The total number of threads per node
(``-j`` option) should equal to the number of cores per node,
which is the default setting.

Multi-node processing with facets
---------------------------------

When the MPI method described above is applied to facetted imaging, the
distribution over nodes is done both over facets and over output channels. This
is often not efficient for facetted imaging.
To improve this, it is recommended to use the "shared reads and writes" option.
This is described in the
:doc:`facet-based imaging section <facet_based_imaging>`. There's also an
option ``-compound-tasks`` that makes sure to send the facets from the same
channel to the same node. Generally, using shared reads and writes is
the faster approach, though there may be circumstances where the use of
compound tasks is faster.

Pinning work to nodes
---------------------

The ``-channel-to-node`` argument specifies the mapping of output channels to
node. For example, ``-channel-to-node 0,0,1,2`` schedules the tasks for
output channel indices 0, 1, 2, 3 at compute nodes 0, 0, 1, 2, respectively.
The length of the list must equal the number of output channels
(``-channels-out`` argument).
If ``-no-work-on-master`` is specified, the list may not contain ``0``.

By default, WSClean uses a round-robin distribution of output channels to nodes,
because later channels with a higher frequency are often more expensive to grid.
For example, with 10 output channels, 4 nodes, and an active main node, the
mapping is 0,1,2,3,0,1,2,3,0,1.
