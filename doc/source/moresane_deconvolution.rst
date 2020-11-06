MORESANE deconvolution
======================

WSClean supports deconvolution with the MORESANE algorithm, which is a compressed sensing method. The MORESANE algorithm is described by `Dabbech et al. (2015) <http://arxiv.org/abs/1412.5387>`_.

WSClean also implements a newer :doc:`IUWT deconvolution algorithm <iuwt_compressed_sensing>` that supersedes MORESANE in most cases. The new IUWT mode supports more features, such as multi-frequency and multi-polarization deconvolution, and the CPU version is faster. A full proper comparison between MORESANE and IUWT has yet to be performed.

Using MORESANE
--------------

WSClean uses the Python implementation of MORESANE, written by Jonathan Kenyon. To use the MORESANE algorithm, you have to manually install the PyMORESANE source code. See https://github.com/ratt-ru/PyMORESANE.

Once installed, you can run wsclean with the -moresane-ext parameter. Directly after the -moresane-ext parameter, the full path of the '``runsane``' (older versions: '``pymoresane.py``') file needs to be specified, e.g.:

.. code-block:: bash

    wsclean -moresane-ext /usr/local/bin/runsane [other parameters] \
        obs.ms

Most normal Clean parameters do not have effect on the moresane deconvolution. The following parameters influence MORESANE deconvolution:

 * The ``-mgain`` parameter turns on a major loop, thereby correcting w-terms in the model. In general, w-terms cause the PSF not to be invariant and moresane does by itself not correct for this. At this point, its exact value does not matter, except that setting it to 1 (the default) implies no major iterations, while setting it lower than one turns on major iterations. Setting it thus lower than 1 is recommended for good imaging.

 * The ``-niter`` parameter sets the number of major (!) iterations to be performed. Setting it to 1 will run MoreSane once, and (with ``mgain``!=1) it will subsequently immediately finish after a prediction-imaging round. Setting it higher will call MoreSane more often. In later MoreSane rounds, WSClean will provide MoreSane with the residual image to which the model convolved with the invariant static PSF has been added. This improves accuracy somewhat, because intuitively the image now consists of a PSF that is "less varying". In images where w-values are relevant, I notice significant improvement in the second and third iteration, after which the image seems to have been converged and further iterations do not noticably change the image.

 * Specifying '``-no-negative``' will constrain the model to not have negative values. This turns on the '``-ep``' parameter of PyMORESANE. So far, I've not seen this make the model better; in general, false negative components are replaced by other false positive components.

 * Parameter "``-moresane-arg``" allows forwarding extra arguments to Moresane.
 
 * Masks are also forwarded to Moresane.

 * Parameter "-moresane-sl" sets the sigma depth level for different iterations.

NB: According to Jonathan, PyMORESANE only works with image sizes that are a power of two. Apparently this is solved in newer versions, but only in CPU mode.

NB2: The ``-niter`` has a significantly different meaning when using MoreSane deconvolution. Specifying hundreds of iterations will not improve quality and will take forever!

I have tested the MORESANE algorithm on images with the (highly resolved) Vela and Puppis A supernova remnants and got very good deconvolution results. My best results were imaged with these settings:

.. code-block:: bash

    wsclean -niter 3 -mgain 0.9 -moresane `which runsane` \
      -size 2048 2048 -scale 1.2amin \
      -name vela-and-puppis vela-and-puppis.ms

Cleaning large images can be very expensive; 1-2k pixel images is doable; 1k takes tens of minutes, 2k takes a few hours. It goes up by more than the square, so you probably won't want to try this on 16k LOFAR images, but rather on small images centred on the diffuse sources.

Finally, the MORESANE algorithm does not support joined deconvolutions, so it will not work with the ``-join-polarizations`` or ``-join-channels`` modi. The :doc:`IUWT deconvolution mode <iuwt_compressed_sensing>` can be used for this.

History
-------
WSClean supports deconvolution with the MORESANE algorithm since :doc:`version 1.7 <changelogs/v1.7>`. :doc:`WSClean 1.8 <changelogs/v1.8>` added parameters "``-moresane-arg``" and added support for passing masks to Moresane. :doc:`Version 1.9 <changelogs/v1.9>` added parameter "``-moresane-sl``".
