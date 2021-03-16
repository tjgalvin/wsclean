A-Term correction
=================

A-term correction is available in WSClean when gridding with the Image Domain Gridder (IDG). Info on how to enable and effectively use IDG can be found on the :doc:`IDG page <image_domain_gridding>`.

Gridding with a-term screens
----------------------------

To use a-term correction or combine multiple corrections, a configuration file can be provided to WSClean to set up the a-term corrections to be applied by IDG. The configuration file can be specified with ``-aterm-config <filename>``, and has a structure like this:

.. code-block:: text

    # This is a test parset. Comments are made by starting a line with a hash symbol
    # The aterms option lists all the corrections that are made. For demonstration,
    # this parset makes all possible corrections:
    aterms = [ tec, dldm, diagonal, beam, paf ]
    
    # A tec correction has parameters 'images' and 'window' :
    # See the WSClean help for a description of the image format used.
    tec.images = [ aterms1-tec.fits aterms2-tec.fits ]
    # The window parameter is new since 2.7
    # It supports tukey, hann, raised_hann, rectangular or gaussian.
    # If not specified, raised_hann is used, which generally performs best.
    tec.window = raised_hann
    
    # The dldm correction (source shift) has parameters 'images', 'window' and 'update_interval'
    dldm.images = [ aterms1-dldm.fits aterms2-dldm.fits ]
    # How often (in seconds) to update the phases computed from the dl,dm values. 
    dldm.update_interval = 300
    dldm.window = raised_hann
    
    # The diagonal correction has parameters 'images' and 'window'.
    diagonal.images = [ aterms1-diag.fits aterms2-diag.fits ]
    diagonal.window = raised_hann
    
    # The beam correction has parameter 'update_interval'. It may also have
    # telescope-specific options, e.g. the lofar beam supports 'differential', and
    # 'usechannelfreq'.
    beam.differential = true
    beam.update_interval = 120
    beam.usechannelfreq = true
    beam.frequency_interpolation = true

    # The paf correction is an easy wrapper for multi-beaming arrays. In this mode, measured
    # beam gains are read from fits files. Each combination of antenna and beam  is stored in
    # a separate fits file with total power (scalar) values, and each fits file can have
    # multiple frequencies to support having a different beam at different frequencies.
    # Each measurement set is assumed to contain the data for a single beam, and
    # the order of the fits files specified here should match the order of the measurement set.
    paf.antenna_map = [ ant_0 ant_1 ant_2 ant_3 ant_4 ant_5 ant_6 ant_7 ant_8 ant_9 ant_10 ]
    paf.beam_map = [ 00 01 ]
    paf.beam_pointings = [ -09h25m00.0s 55d49m59.0s -09h35m33.0s 54d47m29.0s ]
    paf.file_template = beammodels/$ANT/CygA_191120_$BEAM_$ANT_I.fits
    paf.window = hann
 
The ``aterms`` keyword specifies a list of corrections to be applied, and each correction can have some parameters. The TEC and diagonal images are FITS image cubes with a special format; see the chapter below on TEC correction. Those parameters that are specified in the config file are no longer necessary on the command line (e.g. ``-grid-with-beam`` is no longer effective when a config file is specified with beam). The tec and diagonal image lists can take one or multiple fits files. The fits files are concatenated time-wise, so in an observation with 20 timesteps, the first one can specify the first 10 and the second one the last 10. This allows timegaps; e.g. if two observations are imaged that have some gap in time, different fits files can be made so that it is not required to "zero pad" the time in between.

The ``paf`` correction is described in the chapter on :doc:`combining pointings <combining_pointings>`. The other corrections are discussed below.

Kernel size
-----------

When gridding with a-terms, it is important to set the ``-aterm-kernel-size`` parameter to an appropriate value. A-term images are resampled onto a small grid, and they are made as smooth as necessary to fit within the specified aterm kernel size. If your aterm screens (beam / ionospheric / ...) are smooth, the value can be low, but if they need to have a higher resolution, some tuning might be required. In WSClean 2.10, the default size is 16 when specifying a aterm config file. Otherwise / before WSClean 2.10, the default size is 5. A size of five is very small, but enough for the LOFAR beam. It is probably too low for most TEC or gain correction screens. You can check the result of resampling the aterms to the specified size with the ``-save-aterms`` option that is described elsewhere on this page.

TEC correction
--------------

IDG can apply a spatially varying time-variable TEC term that can additionally be different for different antennas and output channels. To use this, the ``-aterm-config`` option described above should be used to supply a tec image. The provided TEC image should have 5 dimensions, ordered as follows:

- RA
- DEC
- Antenna
- Frequency
- Time

The number of antennas should either match with the imaged measurement set, or should have a dimension of one, in which case the same aterm is used for all antennas. The TIME dimension is optional: when not specified, the same corrections are applied to all times. The RA and DEC dimensions are interpolated on the IDG sub-grid by nearest neighbour interpolation. This is typically around 256 pixels, so providing image which are larger is not necessary.  The frequency and time axes are also interpolated. The RA and DEC dimensions should be in the standard radio imaging projection with appropriate CRPIX, CRVAL and CDELT settings. These parameters need also to be set appropriate for the TIME and FREQ axis, setting the time axis to have "aips time" values and the frequency axis to have values in Hz. The times in the FITS file have the same meaning as value in the TIME column in the measurement set; so they represent the time at the centre of the timestep. The screen is selected whose time is closest to that of the time in the TIME column.

Since TEC values are interpolated over frequency with its 1/nu relation, it is normally not required to have more than one channel in the image, unless higher order terms need to be corrected. The correction is constant per output channel, so the output channels have to be chosen such that they are fine enough to achieve the desired accuracy. The values in a TEC file are applied as "delta TEC terms", meaning that a value of zero implies no change to the gain of the antenna. The phase of the gain (in radians) is evaluated as:  ``phase = image[pixel] * -8.44797245e9 / frequency``, with frequency in Hz.

This is an example header of a aterm TEC fits file:

.. code-block:: text

    SIMPLE  =                    T / file does conform to FITS standard
    BITPIX  =                  -32 / number of bits per data pixel
    NAXIS   =                    5 / number of data axes
    NAXIS1  =                 1024 / length of RA axis
    NAXIS2  =                 1024 / length of DEC axis
    NAXIS3  =                   48 / length of ANTENNA axis
    NAXIS4  =                    1 / length of FREQ axis
    NAXIS5  =                   10 / length of TIME axis
    EXTEND  =                    T / FITS dataset may contain extensions
    [..]
    CTYPE1  = 'RA---SIN'           / Right ascension angle cosine
    CRPIX1  =                 513.
    CRVAL1  =          123.4002825
    CDELT1  =              -0.0125
    CUNIT1  = 'deg     '
    CTYPE2  = 'DEC--SIN'           / Declination angle cosine
    CRPIX2  =                 513.
    CRVAL2  =     48.2173836111111
    CDELT2  =               0.0125
    CUNIT2  = 'deg     '
    CTYPE3  = 'ANTENNA '
    CRPIX3  =                   1.
    CRVAL3  =                   0.
    CTYPE4  = 'FREQ    '           / Central frequency
    CRPIX4  =                   1.
    CRVAL4  =     138475036.621094
    CDELT4  =         183105.46875
    CUNIT4  = 'Hz      '
    CTYPE5  = 'TIME    '
    CRPIX5  =                   1.
    CRVAL5  =                   0. / Should be an AIPS time

dldm gain correction
--------------------

"Dl-dm" gain correction can apply a positionshift to correct the position of sources. This kind of correction works almost the same as TEC correction. It also requires a FITS file with 5 dimensions:

    RA, DEC, MATRIX, FREQ, TIME
    
Again, the TIME dimension is optional: when not specified, the same corrections are applied to all times. Like with TEC correction, the dimensions need to be given in this exact order. The dimension ``MATRIX`` should have 2 elements: one for the ``dl`` values, and one for the ``dm`` values. The other dimensions are as described for TEC correction. 

Diagonal gain correction
------------------------

Diagonal gain correction can correct the visibilities with a diagonal Jones matrix. Therefore, diagonal correction performs a correction with two complex values, one for XX and one for YY. Diagonal gain correction with IDG works almost the same as TEC correction. Instead of a FITS file with 5 dimensions, diagonal correction requires a FITS file with 6 dimensions:

    RA, DEC, MATRIX, ANTENNA, FREQ, TIME
    
Like with TEC correction, the dimensions need to be given in this exact order. Compared to the TEC aterms file, there's one extra dimension: ``MATRIX``. For diagonal gains, this matrix dimension has 4 elements: real XX, imaginary XX, real YY and imaginary YY. The other dimensions have their same use. The frequency axis is used to find the nearest image-frequency for each visibility (this works since :doc:`version 2.8 <changelogs/v2.8>`).

If you get images out with all NaNs, the gains might be all zero at some position. For TEC or dldm correction, this obviously is not a problem (zero phase=no correction), but for diagonal gains, a zero matrix leads to division by zero at some point. This can in particular happen because IDG pads the image -- so if one makes TEC aterm images that are exactly the size of the output image, they won't cover the border.

Analyzing / saving the a-terms
------------------------------

The ``-save-aterms`` can be useful for diagnostic output. It turns on saving of the TEC screen after resizing them to the IDG subgrid size and low-pass filtering them to the kernel size (see the kernel size section for more info). The output images are named "``aterm-ev0.fits``" and "``aterm-realxx0.fits``", with increasing numbers for the different aterms over time and counting further in subsequent cleaning iterations. Each image contains a mosaic of images, one image per antenna, starting counting in the bottom left. The images with "ev" in their name are the eigen value of the Jones matrix. These reflect e.g. the power of the beam when imaging with the beam. When imaging with only TEC aterm values, the values are all one, because a TEC change is just a phase change, and the eigenvalue of such a matrix is one: hence not very useful! The images with "realxx" in their names, are the real value of the first ("xx") element of the Jones matrix. These are more useful for assessing TEC aterm values.
