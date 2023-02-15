W-gridding
==========

W-gridding is a wide-field gridding algorithm similar to W-stacking, but instead
of assigning every visibility to its nearest w-plane, it extends the concept
of uv-gridding to the w-direction and therefore grids each visibility to a small
range of w-planes, weighted with a kernel function.
The theoretical foundations of the algorithm are described in
`Ye, Gull, Tan & Nikolic (2021) <https://arxiv.org/abs/2101.11172>`_; a technical
description of the implementation is given in
`Arras, Reinecke, Westermann & Enßlin (2020) <https://arxiv.org/abs/2010.10122>`_.

W-gridding is enabled by the command-line flag ``-gridder wgridder``,
and its accuracy can be controlled via the parameter ``-wgridder-accuracy``,
which is set to ``1e-4`` by default and can be varied in the range from ``1e-2``
(very coarse, but fast gridding) down to ``1e-12``.
This value specifies the maximum acceptable rms error of the result when compared
to a direct Fourier transform.
A value of ``1e-5`` or smaller will enable the double-precision floating
point mode of the wgridder, and larger values will use single-precision calculations.
Be aware that an accuracy as low as ``1e-12`` is not actually reachable, because
the final image is converted to single precision. 

The algorithm will select
appropriate parameters (like amount of padding, kernel shape and kernel support)
automatically to reach the requested accuracy in the least amount of time.
Therefore, many parameters accepted by ``wsclean``'s w-stacking gridder (e.g.
``-nwlayers*``, ``-grid-mode``, ``-kernel-size`` and ``-oversampling``)
will be ignored in that mode.
One parameter that can further be used to influence the gridding accuracy
is ``[-no]-small-inversion``. The default is to enable this optimization,
which speeds up imaging considerably, but be aware that it reduces the
accuracy. The requested accuracy is therefore no longer guaranteed, especially if
there is emission at the image edges.
The ``-padding`` parameter does not directly
influence the padding of the wgridder and should best be left to the default
value.

The algorithm has a very small memory footprint: it only requires storage for
a single complex w-plane, a copy of the dirty image and some indexing data
structures, which are typically much smaller than the visibility data.

Both gridding and degridding directions are available and support shared-memory
parallelization that can be controlled using the ``-j`` parameter.

The w-gridder is available since :doc:`WSClean version 2.9 <changelogs/v2.9>`,
and was further improved in speed in versions
:doc:`2.10 <changelogs/v2.10>` and :doc:`3.0 <changelogs/v3.0>`.
