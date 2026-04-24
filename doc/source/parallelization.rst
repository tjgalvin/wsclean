Parallization options
=====================

WSClean can exploit parallelism by using multiple threads and by using multiple
MPI processes. The MPI processes may run at different compute nodes, which
enables distributed processing.

Multi-threading
---------------
WSClean has several command-line options for configuring multi-threading:

* Number of cores to use: ``-j <N>``

  Tells WSClean to use ``N`` threads. By default, WSClean uses as many threads
  as there are cores in the system.
  WSClean uses this setting everywhere possible, including
  when WSClean uses external libraries, such as IDG and Radler.

* Do :doc:`parallel deconvolution<parallel_deconvolution>`:
  ``-parallel-deconvolution <maxsize>``

  Tells WSClean to do parallel deconvolution where ``maxsize`` is the maximum
  subimage size. For more information, see
  :doc:`parallel deconvolution<parallel_deconvolution>`.

* Reduce memory usage for deconvolution: ``-deconvolution-threads <N>``

  On machines with a large number of cores, using less deconvolution threads
  helps if WSClean runs out of memory during deconvolution. WSClean allocates
  memory for each thread, e.g., for storing the sub-image a number of times.

  This option tells WSClean to use a maximum of ``N`` threads during
  deconvolution. The default value is the number of threads as specified using
  the ``-j`` argument.

* Enable parallel reordering: ``-parallel-reordering <N>``

  Tells WSClean to use parallelism when reordering
  the input measurement sets on disk.
  WSClean uses one reordering task for each input measurements set, and
  executes up to ``N`` tasks in parallel.

  Reordering is bound by disk speed when the number of cores is high.
  Using 4 reordering threads is generally a good balance between disk and CPU
  speed. By default, WSClean therefore uses 4 reordering threads, even if ``-j``
  is specified. The optimal number of reordering threads depends on disk speed.

* Enable parallel gridding: ``-parallel-gridding <N>``

  Tells WSClean to use ``N`` threads during gridding and degridding.
  By default, parallel gridding is disabled, even if ``-j`` is specified.
  Parallel gridding does not work with the ``idg`` gridder, or when using MPI
  (see below).

  The gridders themselves also exploit parallelism internally, however,
  gridders can't always scale well to all cores, in particular when using beam
  or h5parm corrections. In these cases, using parallel gridding might yield
  better performance.

  All gridders get an amount of threads equal to the total number of threads,
  divided by the number of parallel gridders, rounded up. For example,
  ``-parallel-gridding 3 -j 16`` yields 3 parallel gridders that use
  6 threads each.


Distributed imaging
-------------------

WSClean allows the use of multiple nodes for imaging. This is described in
the :doc:`distributed imaging section <distributed_imaging>`.
