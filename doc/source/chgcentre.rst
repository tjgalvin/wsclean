chgcentre
=========

The 'chgcentre' tool can be used to change the phase centre of a measurement set. It will recalculate the *uvw*-values (from the antenna locations, phase centre and time) and phase-rotate the visibilities. We found that the casa task 'fixvis' has several bugs (as of March 2014), and cannot be used to phase rotate a MWA measurement set properly.

Since :doc:`wsclean 2.10 <changelogs/v2.10>`, chgcentre is compiled and installed when wsclean is installed. See the :doc:`general installation instruction <installation>` for help.

Execute chgcentre without parameters to get help on the syntax:

.. code-block::

    A program to change the phase centre of a measurement set.
    Written by André Offringa (offringa@gmail.com).

    Syntax: chgcentre [options] <ms> <new ra> <new dec>

    The format of RA can either be 00h00m00.0s or 00:00:00.0
    The format of Dec can either be 00d00m00.0s or 00.00.00.0

    Example to rotate to HydA:
            chgcentre myset.ms 09h18m05.8s -12d05m44s

    Some options:
    -geozenith
            Will calculate the RA,dec of zenith for each timestep, and moves there. This make the set non-standard.
    -flipuvwsign
            Flips the UVW sign. Necessary for LOFAR, for unknown reasons.
    -minw
            Calculate the direction that gives the minimum w-values for the array.
    -zenith
            Shift to the average zenith value.
    -only-uvw
            Only update UVW values, do not apply the phase shift.
    -shiftback
            After changing the phase centre, project the visibilities back to the old phase centre. This is useful
            in WSClean for imaging with minimum w-values in a different projection.
    -f
            Force recalculation, even if destination is same as original phase direction.
    -datacolumn <name>
            Only phase-rotate the visibilities in the given column. Otherwise, the columns
            DATA, MODEL_DATA and CORRECTED_DATA will all be processed if they exist.
    -from-ms <ms>
            Rotate the measurement set to the same direction as specified
            in the provided measurement set.

When a measurement set contains multiple data columns (e.g., ``DATA``, ``MODEL_DATA``, ``CORRECTED_DATA``), each column will be updated (as long as they have a standard name).

If you do not provide a new RA and dec, chgcentre will give you some info about the measurement set:

.. code-block::

    $ chgcentre myobs.ms
    A program to change the phase centre of a measurement set.
    Written by André Offringa (offringa@gmail.com).

    Current phase direction:
      -00h59m31.7s -16d46m19s
    Zenith is at:
      -01h00m00.8s -26d46m44s
      (-01h00m49.0s -26d46m44s - -00h59m12.7s -26d46m44s)
    Min-w direction is at:
      -00h59m31.7s -26d46m19s

You can specify '``-zenith``' or '``-minw``' as option to rephase to the local array zenith or the direction orthogonal to the best-fit plane to the antennas. The latter is close to zenith, but provides slightly lower *w*-terms. This has not been tested on telescopes other than the MWA. The syntax for this is:

.. code-block::

    $ chgcentre -minw myobs.ms

This can be used in combination with WSClean's ``-shift`` parameter for *w*-snapshot imaging. In that case, the original phase centre should be specified with the ``-shift`` parameter. The net effect is that the measurement set is phase rotated to the sky direction with minimal *w*-terms, and shifted back along the tangent plane to the direction of interest. 

In :doc:`WSClean 2.11 <changelogs/v2.11>`, this approach replaced the ``-shiftback`` option of ``chgcentre`` for shifting the visibilities along the tangent plane. See the :doc:`w-snapshot algorithm page <w_snapshot_algorithm>` for more info. 

Legacy data with ``-shiftback`` applied
---------------------------------------

Before WSClean 2.11, it was possible to prepare w-snapshotting with the ``-shiftback`` option, e.g.:

.. code-block::

    $ chgcentre -minw -shiftback myobs.ms

A shifted measurement set uses wsclean-specific keywords. Support for this was removed in WSClean 2.11 (see the :doc:`changelog <changelogs/v2.11>` for details), and any observation for which these shifting keywords are detected will produce an error in WSClean 2.11. 

In case archival data to which ``-shiftback`` is applied needs to be imaged with WSClean 2.11, the option should be undone. A shifted measurement set can be restored by phase rotating it to its original RA/dec: ``chgcentre`` will detect the keywords in the measurement set, undo the shift and update the keywords.

A LOFAR bug
-----------

For unknown reasons, the uvw value needs to be flipped for LOFAR sets. As far as I know, this is not necessary for other telescopes, but LOFAR requires you to specify ``-flipuvwsign``.
