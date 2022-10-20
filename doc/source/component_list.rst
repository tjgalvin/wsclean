Component lists
===============

The benefit of a component list is that it allows representing the components of multi-scale cleaning without 'pixelizing' the larger kernels, thereby decreasing the number of required components in a sky model compared to representing all components with delta scales.

The option is enabled by adding ``-save-source-list`` on the command line

File format
-----------

The component list file format is compatible with the Blackboard Self-cal DP3 sky model formats, but be aware that versions after April 2017 of DP3 support the polynomial spectral function that is used.

When ``-save-source-list`` is added on the command line, a text file named ``<prefix>-sources.txt`` is saved. The file contains comma-separated values and looks as follows:

.. code-block:: text

    Format = Name, Type, Ra, Dec, I, SpectralIndex, LogarithmicSI, ReferenceFrequency='125584411.621094', MajorAxis, MinorAxis, Orientation
    s0c0,POINT,08:28:05.152,39.35.08.511,0.000748810650400475,[-0.00695379313004673,-0.0849693907803257],false,125584411.621094,,,
    s0c1,POINT,08:22:27.658,39.37.38.353,-0.000154968071120503,[-0.000898135869319762,0.0183710297781511],false,125584411.621094,,,
    s0c2,POINT,08:18:44.309,39.38.37.773,0.000233552686127518,[-0.000869089801859608,0.0828587947079702],false,125584411.621094,,,
    s0c3,POINT,08:03:07.538,39.37.02.717,0.000919058240247659,[0.001264109956439,0.0201438425344451],false,125584411.621094,,,
    [..]
    s1c0,GAUSSIAN,08:31:10.37,41.47.17.131,0.000723326710524984,[0.00344317919656096,-0.115990377833407],false,125584411.621094,83.6144111272856,83.6144111272856,0
    s1c1,GAUSSIAN,07:51:09.24,42.32.46.177,0.000660490865128381,[0.00404869217508666,-0.011844732049232],false,125584411.621094,83.6144111272856,83.6144111272856,0
    [..]
    
Each source (i.e., clean component) is one line in the file. The first line is a header starting with "``Format = ...``" that describes the columns. The header is allowed to specify a default value for this column, as is done above for the reference frequency. When a field is left empty, it should take the default value. See the BBS / DP3 documentation for more info on default values.

The columns
-----------

The **Name** contains a unique identifier for this component. In the current format, it consists of the letter 's' followed by the scale index of the component (0 meaning the smallest scale).

The **Type** column is either ``POINT`` or ``GAUSSIAN``. For a point component, the axes and orientation fields are empty. A Gaussian is a standard Gaussian as often used in sky models. When asking WSClean to output a source list, it will automatically enable Gaussian shaped components, instead of the default tapered quadratic. 

**RA** and **dec** are the central coordinates of the component, in notation of "hh:mm:ss.sss" and "dd.mm.ss.sss". 

The **I** column represents the flux density in Jy at the reference frequency.

The **SpectralIndex** column is an array of numbers surround by square brackets, that represent the coefficients of the logarithmic or ordinary polynomial (see **logarithmicSI**), when normalized to the reference frequency. The logarithmic polynomial function is given by

.. math::

    \log S(\nu) = \log (S_0) + c_0 \log \left(\frac{\nu}{\nu_0} \right) + c_1 \log \left( \frac{\nu}{\nu_0} \right)^2 + ...

The logarithms are "base 10" logarithms, such that it is equivalent to say :math:`S(\nu)=10^\left(\log (S_0) + ... \right)`.
Also note that :math:`c_0` represents the spectral index term.

An ordinary polynomial function is evaluated as

.. math::

    S(\nu) = S_0 + p_0 \left(\frac{\nu}{\nu_0} - 1\right) + p_1 \left(\frac{\nu}{\nu_0} - 1\right)^2 + ...

Note that the value 1 is subtracted in the base. This makes sure that :math:`S_0` (the "Stokes I" value) represents the flux density at the reference frequency.

**LogarithmicSI** is *true* or *false*, denoting that the spectral index column uses logarithmic or ordinary polynomials, respectively. Note that in absense of this column, logarithmic polynomials are used.

**ReferenceFrequency** gives the frequency in Hz at which the polynomial or logarithmic polynomial is normalized.

The **MajorAxis**, **MinorAxis** and **Orientation** columns define the shape of the Gaussian. The axes are given in units of arcseconds, and orientation is in degrees. Note that currently the major and minor axis of scales are always the same, and thus the orientation has no effect (and is therefore always zero).
 
History
-------

WSClean supports output of a component list since :doc:`version 2.4 <changelogs/v2.4>`. Note that :doc:`version 2.3 <changelogs/v2.3>` supported a different format which is now deprecated.
