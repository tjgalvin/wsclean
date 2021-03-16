Snapshot imaging
================

WSClean accepts a convenience parameter to make it easy to separate a measurement set in small snapshots and image each one separately. The parameter to do this is:

.. code-block:: text

    -intervals-out <count>

This tells WSClean to split the measurement set in the given number of intervals, and image each interval separately. The file name of all the image products of an interval will have the interval number preceded by a 't'. For example:

.. code-block:: bash

    wsclean -intervals-out 100 -scale 1asec -mgain 0.8 -threshold 1 \
      -niter 10000 -name cyga MyCygASet.ms
    
This will output images ``cyga-t0000-dirty.fits``, ``cyga-t0000-image.fits``, and all the other products for t0000, and so for t0001...t0100.

To image only part of the full time span, the ``-interval`` parameter can be used to select the range of timesteps that is to be imaged. So, with "``-intervals-out 100 -interval 200 300``", the first image will be made timestep 200, the next one from timestep 201 and further until timestep 299. Units are in "number of timesteps" (so the integration time). 

Note that the selection of timesteps with ``-interval`` is done before splitting into snapshots: if a measurement set has 10 timesteps, giving ``-intervals-out 10 -interval 2 4`` is malformed, as it would first select the interval from timestep 2 to 4, and then split those 2 timesteps into 10 intervals.

You can use all other :doc:`cleaning parameters <basic_cleaning>` together with ``-intervals-out`` (``-mgain``, ``-auto-threshold``, ``-niter``... etc.). 
