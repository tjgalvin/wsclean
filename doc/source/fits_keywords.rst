FITS keywords
=============

WSClean adds the following non-standard keywords to its output FITS files:

General
-------

.. code-block:: text

    WSCVERSI -- Version string for WSClean, e.g. "2.5.1"  (since WSClean 2.5).
    WSCVDATE -- Version date of WSClean, e.g. "2017-05-08"  (since WSClean 2.5).
    WSCNWLAY -- nwlayers
    WSCDATAC -- data column
    WSCWEIGH -- weighting (textual description, e.g. uniform or Briggs'(0.5))
    WSCGKRNL -- gridding (antialiasing) kernel size
    WSCCHANS / WSCCHANE -- channel range (only when specified)
    WSCTIMES / WSCTIMEE -- time range (only when specified)
    WSCFIELD -- imaged field number in MS
    WSCNVIS  -- gridded visibility count
    WSCENVIS -- effective nr of visibilities (only for natural weighting)
    WSCVWSUM -- sum of visibility weight (since WSClean 2.3)
    WSCIMGWG -- Weight to be used when averaging images together
    WSCNORMF -- Normalization factor that was applied to the image. The factor is
                useful to undo the normalization for e.g. conversion to Kelvins.
                WSCNORMF is available since WSClean 2.
    WSCTHRES -- Manually applied threshold

Cleaning
--------
   
.. code-block:: text

    WSCNITER -- max number of iterations specified
    WSCGAIN -- clean gain (minor iteration gain)
    WSCMGAIN -- major cleaning gain
    WSCNEGCM -- whether negative components were allowed to clean
    WSCNEGST -- whether stopping on first negative component during clean
    WSCMAJOR -- number of major iterations actually used before reaching stopping criteria
    WSCMINOR -- number of minor iterations actually used before reaching stopping criteria 

The actually used number of major/minor iterations are stored in ``WSCMAJOR`` and ``WSCMINOR``, whereas ``WSCMGAIN`` and ``WSCNITER`` describe the input parameters. 

No longer used
--------------

.. code-block:: text

    WSCSMPSF -- whether the PSF is made smaller prior to cleaning


    
