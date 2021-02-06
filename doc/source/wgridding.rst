W-gridding
==========

W-gridding is a wide-field gridding algorithm similar to W-stacking, but instead
of assigning every visibility to its nearest w-plane, it extends the concept
of uv-gridding to the w-direction and therefore grids each visibility to a small
range of w-planes, weighted with a kernel function.
The theoretical foundations of the algorithm are described in
`the thesis of H. Ye (2019) <https://www.repository.cam.ac.uk/handle/1810/292298>`_
(pp. 139); a technical
description of the implementation is given in
`Arras, Reinecke, Westermann & Enßlin (2020) <https://arxiv.org/abs/2010.10122>`_.

W-gridding is enabled by the command-line flag ``-use-wgridder``,
and its accuracy can be controlled via the parameter ``-wgridder-accuracy``,
which is set to ``1e-4`` by default and can be varied in the range from ``1e-2``
(very coarse, but fast gridding) down to ``1e-6`` (most accurate gridding
achievable with single-precision floating point). This value specifies
the maximum acceptable rms error of the result when compared to a direct Fourier
transform. The algorithm will select
appropriate parameters (like amount of padding, kernel shape and kernel support)
automatically to reach the requested accuracy in the least amount of time.
Therefore, many parameters accepted by ``wsclean``'s w-stacking gridder (e.g.
``-padding``, ``-nwlayers*``, ``-grid-mode``, ``-kernel-size`` and ``-oversampling``
and ``-parallel-gridding``) will be ignored in that mode.

The algorithm has a very small memory footprint: it only requires storage for
a single complex w-plane, a copy of the dirty image and some indexing data
structures, which are typically much smaller than the visibility data.

Both gridding and degridding directions are available and support shared-memory
parallelization that can be controlled using the ``-j`` parameter.

The w-gridder is available since :doc:`WSClean version 2.9 <changelogs/v2.9>`,
and was further improved in speed in versions
:doc:`2.10 <changelogs/v2.10>` and :doc:`2.11 <changelogs/v2.11>`.
