Wideband deconvolution
======================

Since :doc:`version 1.1 <changelogs/v1.1>`, WSClean has a 'wideband' multi-frequency deconvolution mode, which allows cleaning channels joinedly. This means that peak finding is performed in the sum of all channels, allowing deep cleaning, and the psf is subtracted from each channel & polarization individually, scaled to the value of the peak in that image, which takes care of spectral variation.

Usage
-----

WSClean's wideband mode can work together with joined polarization cleaning or for imaging a single polarization (from :doc:`version 1.5 <changelogs/v1.5>`).

A typical run in multi-frequency deconvolution would look like:

.. code-block:: bash

    wsclean -join-channels -channels-out 4 [other parameters] \
      <measurement set>

This outputs 4 deconvolved images at the different frequencies and the weighted average of those 4. It is of course expensive: the deconvolution performance of the above statement is about four times more expensive than a bandwidth-integrated clean, and especially the minor clean iterations become much slower.  This mode can be combined with Cotton-Schwab imaging (with ``-mgain ...``) which would fill the ``MODEL_DATA`` column with frequency dependent model info, so that it is possible to perform self-cal with the proper frequency information. It can also be combined with -joinpolarizations to clean polarizations and channels both joinedly.

Something to be aware of is that the ``-join-channels`` parameter turns on :doc:`MF weighting <mf_weighting>`.

Multiple measurement sets and subbands
--------------------------------------

WSClean calculates the full available bandwidth (over all measurement sets) and splits that into chunks with the same amount of input channels (or 1 off if the nr of input channels is not divisable by the nr of output channels). For example, if you have four measurement sets where each contains data from a different subband, e.g. 100 MHz, 110 MHz, 200 MHz and 210 MHz, and you specify "``-join-channels -channels-out 2``", the first images will be centred at 105 MHz and contain data from the first two MSes, and the other images will be centred at 205 MHz.

Example:

.. code-block:: bash

    wsclean -join-channels -channels-out 3 -niter 10000 \
      -mgain 0.8 -threshold 0.01 \
      band-100-MHz.ms band-110-MHz.ms \
      band-145-MHz.ms band-155-MHz.ms \
      band-190-MHz.ms band-200-MHz.ms

This will perform multi-frequency deconvolution and output images at 105 MHz, 150 MHz and 195 MHz. This works since :doc:`WSClean version 1.9 <changelogs/v1.9>`. 

In case the measurement sets have a different amout of channels *and* have gaps between them, the normal division of WSClean might divide the channels somewhat undesirably. For example, the data of the MSSS survey observations are stored in 10 measurement sets, each containing one 'band', but one of those 10 has fewer channels. The default WSClean division divides those "missing" channels out over the bandwidth, producing this imaging table:

.. code-block:: text

    .      # Pol Ch JG ²G In Freq(MHz) ch
    | Independent group:
    +-+-J- 0  I   0  0  0  0  119-121 (39)
      +-J- 1  I   1  0  1  0  121-126 (40)
      +-J- 2  I   2  0  2  0  126-130 (39)
      +-J- 3  I   3  0  3  0  130-136 (40)
      +-J- 4  I   4  0  4  0  136-146 (39)
      +-J- 5  I   5  0  5  0  146-150 (40)
      +-J- 6  I   6  0  6  0  150-152 (39)
      +-J- 7  I   7  0  7  0  156-158 (40)

There are two ways to remedy this. The first, most automated way is to use an option called ``-gap-channel-division`` (available since :doc:`WSClean 2.6 <changelogs/v2.6>`), which calculates the gaps between channels, and splits the input channels into output channels by splitting the largest gap until the number of output channels has been reached. For the above situation, this gives the following table: 

.. code-block:: text

    .      # Pol Ch JG ²G In Freq(MHz) ch
    | Independent group:
    +-+-J- 0  I   0  0  0  0  119-121 (40)
      +-J- 1  I   1  0  1  0  124-126 (40)
      +-J- 2  I   2  0  2  0  128-130 (40)
      +-J- 3  I   3  0  3  0  134-136 (40)
      +-J- 4  I   4  0  4  0  142-144 (36)
      +-J- 5  I   5  0  5  0  146-148 (40)
      +-J- 6  I   6  0  6  0  150-152 (40)
      +-J- 7  I   7  0  7  0  156-158 (40)

Notice the different bandwidth per channel and output channel index 4 which now contains 36 channels.

The second option is to manually add splits by specifying frequencies. For this, the option ``-channel-division-frequencies`` can be used. The list doesn't necessarily contain all splits necessary. If fewer splits are given than required, WSClean will further subdivide the remaining parts as necessary. Be aware that splits are in Hz, but can be in scientific notation (e.g. 1.4e9 is 1.4 GHz). An example:

.. code-block:: text

    wsclean -size 512 512 -scale 1amin -channels-out 8 \
      -channel-division-frequencies 140e6,141e6,142e6 \
      observation.ms
    [..]
    === IMAGING TABLE ===
          # Pol Ch JG ²G In Freq(MHz)
    +-+-J- 0  I   0  0  0  0  134-137 (39)
    +-+-J- 1  I   1  1  1  0  137-140 (39)
    +-+-J- 2  I   2  2  2  0  140-140 (6)
    +-+-J- 3  I   3  3  3  0  140-141 (7)
    +-+-J- 4  I   4  4  4  0  141-142 (6)
    +-+-J- 5  I   5  5  5  0  142-142 (6)
    +-+-J- 6  I   6  6  6  0  142-153 (140)
    +-+-J- 7  I   7  7  7  0  153-164 (141)

Fitting smooth spectra
----------------------

The joined channel deconvolution method discussed above does not enforce a smooth spectra; each channel gets a separate solution. If one wants to image many spectral channels separately, while it is known that the sources have smooth behaviour, it is possible to enforce this during cleaning. Currently, WSClean supports fitting a polynomial and fitting a double-logarithmic polynomial. The command line parameters for this are ``-fit-spectral-pol`` and ``-fit-spectral-log-pol``. Both require an extra parameter specifying the number of terms (degrees of freedom). For example, ``-fit-spectral-log-pol 2`` will fit a power law through all the output channels.

A simple -- but very slow -- example to perform cleaning with spectral fitting:

.. code-block:: bash

    wsclean -multiscale -join-channels -channels-out 64 -niter 10000 \
      -mgain 0.8 -auto-threshold 3 -fit-spectral-pol 4 observations.ms

This will fit a polynomial with 4 terms (i.e., a third-order polynomial). This would be similar to CASA multi-term deconvolution with ``nterms=4``. During each minor clean cycle, the 64 images at different frequencies will be added together, the pixel with the highest summed brightness is selected, the brightness for that pixel is found for each image, a 3rd order polynomial is fitted through those measurements and the smoothed "model" component is added to the model, as well as convolved with the PSF and subtracted from the residual dirty image. Spectral fitting works in all joined-channel modes (i.e., hogbom, multi-scale, iuwt, moresane).

As you might imagine, doing a clean with 64 images in memory is expensive, both in terms of memory and computing. Since it is not necessary to have that many images in memory when fitting only a few terms, it is also possible to decrease the number of output channels just during deconvolution. This is done with the ``-deconvolution-channels`` parameter, for example:

.. code-block:: bash

    wsclean -join-channels -channels-out 64 -niter 10000 \
      -mgain 0.8 -auto-threshold 3 -fit-spectral-pol 4 \
      -deconvolution-channels 8 observations.ms

This will decrease the number of images from 64 to 8 before starting the deconvolution by averaging 8 groups together. Cleaning is then performed with just 8 images. After cleaning, the requested function (3rd order polynomial in this case) is fitted to the model, and the model is interpolated using that function.
This is much faster than the previous command, and equally precise. Setting the deconvolution channels is supported in all modes since :doc:`WSClean 2.2 <changelogs/v2.2>`. The spectral-fitting features were added in :doc:`WSClean version 1.11 <changelogs/v1.11>`.

Forced spectral indices
-----------------------

Spectral fitting is useful to reduce the degrees of freedom in deconvolution. It does not always produce accurate, physical spectral indices. This can either be because the bandwidth is too small, sources are very complex causing degeneracies, or both. WSClean version 2.12 has therefore a method called 'forced spectrum fitting'. In this mode, a pre-existing spectral index map is used during the deconvolution, and the resulting spectra of the model components are forced onto this spectral index map.

The spectral index map may be the result of earlier runs or from fitting between data from other bandwidths / telescopes. An added advantage is that those maps can be smoothed to limit the effect of noise and imaging artefacts. The mode is enabled by combining ``-force-spectrum <fits filename>`` with ``-fit-log-pol 2``. The fits file should have the same dimensions and coordinate system as used for the imaging, and each pixel should hold a spectral index value. This mode should currently not be used together with the ``-deconvolution-channels`` option.

Together with :doc:`multiscale cleaning <multiscale_cleaning>` and :doc:`source list output <component_list>`, this mode allows building (text) models of sources with accurate spectral index information. 

This method is currently being written up into an article -- to be submitted somewhere in 2021.

Fit normal or logarithmic polynomials?
--------------------------------------

In general, fitting polynomials (``fit-spectral-pol``) works better than fitting logarithmic functions (``fit-spectral-log-pol``).  The latter option fits every component to a power law or higher-termed logarithmic function. This works fine on strong sources, but when cleaning picks up some partly-negative artefacts, it can create bad results and fail. As an example, when one spectral value is negative and the others are positive, a power-law fit that minimizes the squared error, might create extremely high values at the edge channels.

As a side note, I believe that CASA also doesn't fit logarithmic polynomials, but rather fits normal polynomials (I'm not 100% sure from the papers about the technique, but from inspecting the CASA code, it seems to use normal polynomials). Hence, getting CASA's "nterms" behaviour is closest to ``-fit-spectral-pol``.

Avoiding "steps" in the model data visibilities
-----------------------------------------------

When splitting up bandwidth with the ``-channels-out`` option, the output ``MODEL_DATA`` visibilities will have a step function. For example, when splitting up the bandwidth of a 256 channel set in 4 parts (``-channels-out 4``), there will be a step each 64 channels. When the ``MODEL_DATA`` is important, for example when removing the continuuum data from a HI set, these steps might be undesirable. Using ``-fit-spectral-pol`` does not remove these steps; that option only forces the individual images to obey a polynomial. To remove the steps, one would have to use ``-fit-spectral-pol`` and specify as many channels to ``-channels-out`` as that are present in the measurement set. This is of course very computationally and memory expensive, in particular during deconvolution. To make that possible in some cases, the above option ``-deconvolution-channels`` can be used in combination with the other options. 

Another option to solve these step functions is to interpolate in the frequency during the prediction (i.e. when calculating the ``MODEL_DATA``). WSClean can currently not do this, but will hopefully soon be implemented. Because of the already-available features mentioned above it has not been of very high priority, since the quality improvement on continuum imaging is rather small and the computational savings might only exist in some more exotic cases.

Relation to CASA's multi-term deconvolution
-------------------------------------------

To avoid confusion, WSClean's wide-band mode is not the same as the multi-frequency deconvolution algorithm in CASA, i.e., the algorithm described in `Sault & Wieringa (1994) <http://adsabs.harvard.edu/abs/1994A%26AS..108..585S>`_, further enhanced in `Rau and Cornwell (2011) <http://arxiv.org/abs/1106.2745>`_. The multi-frequency deconvolution algorithm in WSClean is called "joined-channel deconvolution" and is described in `Offringa and Smirnov (2017) <https://arxiv.org/abs/1706.06786>`_. As shown there, the multi-frequency implementation of WSClean is computationally *orders of magnitude* faster than MSMFS, which is the most important difference between the two algorithms. Accuracy wise, they produce similar results.

Although WSClean's mode and CASA's mode both mitigate the problem of spectral variation over the imaged bandwidth, the underlying algorithms are different and have different settings. Both CASA and WSClean can fit smooth functions over frequency. WSClean has additional functionality for the imaging of cubes that can have unsmooth spectral functions, such as for spectral lines or off-axis imaging (where the beam causes large/steep fluctuations over frequency). Smooth function fitting is implemented in :doc:`WSClean version 1.10 <changelogs/v1.10>`.
