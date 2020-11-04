W-gridding
==========

W-gridding is a wide-field gridding algorithm similar to W-stacking, but instead
of assigning every visibility to its nearest w-plane, it extends the concept
of uv-gridding to the w-direction and therefore grids each visibility to a small
range of w-planes, weighted with a kernel function.
The theoretical foundations of the algorithm are described in
<https://www.repository.cam.ac.uk/handle/1810/292298> (pp. 139); a technical
description of the implementation is given in
<https://arxiv.org/abs/2010.10122>.

W-gridding is enabled via the command-line ``-use-wgridder``. The algorithm will
select appropriate parameters (like amount of padding, kernel shape and kernel
support) automatically to reach the requested accuracy (currently hard-wired to
1e-4) in the least amount of time. Therefore, many paremeters accepted by
``wsclean``'s w-stacking gridder (e.g. ``-padding``, ``-nwlayers*``, ``-grid-mode``,
``-kernel-size``, ``-oversampling`` and ``-parallel-gridding``) will be ignored in
that mode.

The algorithm has a very small memory footprint: it only requires storage for
a single complex w-plane, a copy of the dirty image and some indexing data
structures, which are typically much smaller than the visibility data.

Both gridding and degridding directions are available and support shared-memory
parallelization that can be controlled using the ``-j`` parameter.
