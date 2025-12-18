BDA Support in MeasurementSets
==============================

The information below was written during the implementation of baseline-dependent averaging in Dp3 and WSClean, and describes how BDA is stored in a measurement set.

Rationale
---------

The MeasurementSet (MS) format is generic enough to support Baseline
Dependent Averaging (BDA) in time. It is in fact generic enough to
support any kind of time averaging applied to every row in the
measurement set.

Averaging in frequency is also possible in the MS. However, it requires
creating a Spectral Window per frequency averaging (as the metadata
about frequency is stored in the ``SPECTRAL_WINDOW`` subtable).

No restrictions are imposed on averaging by the MeasurementSet format.
Therefore, time integrations can overlap, any ordering of rows is
allowed, and the integration times do not need to have any similarity
over the rows. In practice, however, no tools exist that allow this full
flexibility.

We suggest to add keywords to the MS that allow tools to determine
(without going through all the data) some of the characteristics of the
BD averaging method applied. In particular, tools that support BDA
should be quickly able to determine:

- If the data is BDA;
- How much memory is needed to hold all data for a certain interval of
  data (e.g. a 30 sec calibration interval); and
- Whether and how the data can be expanded back to a fully regular MS.
  That is, if the BDA is simple enough (as in the normal use-case) it
  should be possible to recover the original time and frequency axis.

A specific tool may work only on a subset of the MS, e.g. only on the
calibrator scans. Those may be BD averaged differently than e.g. the
target scans. We therefore suggest to record the regularity on a finer
level than on the entire MS.

Suggested implementation
------------------------

Time and frequency direction BDA are treated differently. We will
describe each.

We assume that the same averaging has been applied to all columns
(``DATA``, ``WEIGHT_SPECTRUM``, ``FLAGS`` but also ``CORRECTED_DATA`` and other
data columns). We think that this assumption is compatible with BDA
use-cases.

BDA_TIME_AXIS: Metadata for time table
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We suggest to add an optional structured set of meta information that
specify the regularity of the MS. These keywords are stored in a
subtable called ``BDA_TIME_AXIS`` of the main table, and specify for each
subset of the MS the characteristics of if and how BDA is applied.
Absence of the table implies no BDA is applied in the time direction.
The table has one non-optional keyword ``BDA_TIME_AXIS_VERSION``. The new
subtable specifies regularity only for the time axis, since for the
frequency axis metadata can be specified in the ``SPECTRAL_WINDOW``
subtable.

A measurement set selection might also be on ``FEED1``, ``FEED2``, ``SCAN_ID``
or ``PROCESSOR_ID``. Because this is very rare, we don't want to burden
all software with supporting such a selection, and therefore we
explicitly decided not to include these.

+--------------------------------+--------+-------+---------+---------------------+
| Name                           | Format | Units | Measure | Comments            |
+================================+========+=======+=========+=====================+
| **Columns**                    |        |       |         |                     |
+--------------------------------+--------+-------+---------+---------------------+
| *Keywords*                     |        |       |         |                     |
+--------------------------------+--------+-------+---------+---------------------+
| BDA_TIME_AXIS_VERSION          | String |       |         | Version tag         |
+--------------------------------+--------+-------+---------+---------------------+
| *Key*                          |        |       |         |                     |
+--------------------------------+--------+-------+---------+---------------------+
| FIELD_ID                       | Int    |       |         | Field id            |
+--------------------------------+--------+-------+---------+---------------------+
| BDA_FREQ_AXIS_ID               | Int    |       |         | Spectral window id. |
+--------------------------------+--------+-------+---------+---------------------+
| *(BDA_TIME_AXIS_ID)*           | Int    |       |         | Row id, Only in     |
|                                |        |       |         | MSv3                |
+--------------------------------+--------+-------+---------+---------------------+
| *Data*                         |        |       |         |                     |
+--------------------------------+--------+-------+---------+---------------------+
| IS_BDA_APPLIED                 | Bool   |       |         | BDA applied         |
+--------------------------------+--------+-------+---------+---------------------+
| *(SINGLE_FACTOR_PER_BASELINE)* | Bool   |       |         | Single factor used  |
+--------------------------------+--------+-------+---------+---------------------+
| *(MAX_TIME_INTERVAL)*          | Double | s     |         | Maximum time        |
|                                |        |       |         | interval            |
+--------------------------------+--------+-------+---------+---------------------+
| *(MIN_TIME_INTERVAL)*          | Double | s     |         | Minimum time        |
|                                |        |       |         | interval            |
+--------------------------------+--------+-------+---------+---------------------+
| *(UNIT_TIME_INTERVAL)*         | Double | s     |         | Time interval       |
+--------------------------------+--------+-------+---------+---------------------+
| *(INTEGER_INTERVAL_FACTORS)*   | Bool   |       |         | Interval factors    |
+--------------------------------+--------+-------+---------+---------------------+
| *(HAS_BDA_ORDERING)*           | Bool   |       |         | Ordered             |
+--------------------------------+--------+-------+---------+---------------------+

..

**FIELD_ID**

   Field identiﬁer (≥ 0).

**BDA_FREQ_AXIS_ID**

   Spectral window identifier (≥ 0). Together with the field identifier,
   this key uniquely identifies a row. The ``BDA_SET_ID``
   values in the ``SPECTRAL_WINDOW`` table refer to this key.

**BDA_TIME_AXIS_ID**

   Unique id.

**IS_BDA_APPLIED**

   BDA has been applied to the time axis.

**SINGLE_FACTOR_PER_BASELINE**

   For every baseline, the averaging factor is constant in time. If a
   specific baseline is averaged to 10 seconds for one timestep, it will
   be averaged to 10 seconds for every timestep over the selected data
   range. We expect this property to be true in the normal BDA cases.

**MAX_TIME_INTERVAL**

   Maximum ``TIME_INTERVAL`` over this subset. With BDA applied, this
   normally is the time interval of the smallest baselines. This value,
   together with ordering properties (discussed below), helps a software
   tool by telling it how long a baseline might have a contribution to
   an interval.

**MIN_TIME_INTERVAL**

   Minimum ``TIME_INTERVAL`` over this subset.

**UNIT_TIME_INTERVAL**

   This is basically the original ``TIME_INTERVAL``. If BDA is applied and
   the shortest baseline has not been averaged down, this value will be
   equal to ``MIN_TIME_INTERVAL``. If a high-time resolution measurement
   set is averaged down immediately with BDA, it might be that the
   shortest baseline is averaged and the longest baseline is not an
   integer factor of the shortest baseline averaging factor (and thus
   MIN_TIME_INTERVAL), whereas all baselines are still a multiple of
   some underlying time interval.

**INTEGER_INTERVAL_FACTORS**

   The TIME_INTERVAL – and therefore also the distance (difference
   between two TIMEs) between two consecutive timesteps – is an integer
   multiple of the UNIT_TIME_INTERVAL value. This implies that for all
   baselines, the first intervals starts at the same time, i.e. there is
   no offset.

**HAS_BDA_ORDERING**

   If a row *starts* at T_0 (where T_0 = ``TIME`` - 0.5 \*
   ``TIME_INTERVAL``) then all visibilities that *end* before T_0 are
   before this row. In other words, only overlapping intervals are
   allowed to not obey time ordering: non-overlapping intervals are
   strictly ordered.

BDA_FACTORS
~~~~~~~~~~~

This optional table describes fixed BDA time averaging factors for a
baseline. When this table is present, some values in the ``BDA_TIME_AXIS``
table are redundant and can be derived from this table.

================== ====== ===== ======= ==============================
\
Name               Format Units Measure Comments
**Columns**
*Key*
BDA_TIME_AXIS_ID   Int                  Reference to *BDA_TIME_AXIS*
SPECTRAL_WINDOW_ID Int                  Reference to *SPECTRAL_WINDOW*
ANTENNA1           Int                  Antenna 1
ANTENNA2           Int                  Antenna 2
*(ANTENNA3)*       Int                  Antenna 3
*Data*
TIME_FACTOR        Int                  Time averaging factor
================== ====== ===== ======= ==============================

**BDA_TIME_AXIS_ID**

   Refers to a row in the ``BDA_TIME_AXIS`` table, since a single baseline
   may have different fields. The ``BDA_TIME_AXIS`` row contains the
   ``FIELD_ID`` and common values for multiple baselines regarding time
   averaging.

**SPECTRAL_WINDOW_ID**

   Refers to a row in the ``SPECTRAL_WINDOW`` table, since a single
   baseline may have multiple spectral windows.

**ANTENNAn**

   Antenna identiﬁer, as indexed from ``ANTENNAn*\`` in ``MAIN``. Together,
   the antenna identifiers determine the baseline.

**TIME_FACTOR**

   Time averaging factor for the given baseline, field (via
   the ``BDA_TIME_AXIS_ID`` table) and spectral window. The effective time
   interval for the baseline is this ``TIME_FACTOR`` times the
   ``UNIT_TIME_INTERVAL`` value in the ``BDA_TIME_AXIS`` table.

SPECTRAL_WINDOW: continued (Metadata for frequency)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All other requirements (such as specified for the time metadata) can be
deduced from the metadata already in the ``SPECTRAL_WINDOW`` table.
Spectral BDA has been applied when: i) the BDA_SET_ID column
exists; and ii) multiple SPWs (rows in the table) have the same value in
the BDA_SET_ID column.

========================== ====== ===== ======= ==========
\
Name                       Format Units Measure Comments
*Data*
*(BDA_SET_ID)*             Int
========================== ====== ===== ======= ==========

**BDA_SET_ID**

   An id that links a set of spectral windows that cover the same
   (true/original) spectral window. It is the equal for all spectral
   windows where the only difference is the amount of frequency
   averaging. This id refers to the ``BDA_FREQ_AXIS_ID`` in the
   ``BDA_TIME_AXIS`` table, which contains information regarding the time
   averaging for the spectral window.

Ordering
--------

For BDA data, the keyword ``SORT_COLUMNS`` in the main table does not
suffice to determine the ordering (in time) of the rows. In the case of
BDA with non-strictly ordered overlapping intervals the ``SORT_COLUMNS``
should be "*None*" and more detail should be specified in the column
``HAS_BDA_ORDERING`` in the table ``BDA_TIME_AXIS``. In case the data is
still also strictly ordered in ``TIME``, this can of course still be
registered in the ``SORT_COLUMNS`` keyword.

Use cases
---------

1. Cut a BDA set half in time.
   If the set is cut half in time through a time averaged interval, the
   set becomes non-regular, and the integer factors do not make sense
   anymore. If the set is cut at a point where all intervals end, it can
   remain regular.
2. Cut a BDA set half in frequency.
   A frequency bin can be cut in half, then the spectral window mapping
   table needs updating. It may make more sense to cut at a point where
   all frequency bins have a matching end point.
3. Expand a BDA set back to what it was in a sensible way.
   If the table BDA_FACTOR is filled, this can be done. If
   INTEGER_INTERVAL_FACTORS is False, going back to original time
   resolution will probably not be possible. Frequency averaging has
   integer factors always, so expanding back to the original channels
   should always be possible.
4. Determine the maximum number of samples in a given interval.
   This can be done if INTEGER_INTERVAL_FACTORS is True, using
   MAX_TIME_INTERVAL and MIN_TIME_INTERVAL, and the metadata computed
   from the SPECTRAL_WINDOW table for each of the
   BDA_SET_IDs. If the BDA_FACTOR table is filled, the max
   and min factor in that table can be used.
5. Use different time averaging for different parts of the observation.
   E.g. if a specific part of the observation requires higher time
   accuracy (e.g. due to ionospheric behavior).
   For this, in the current design, a new FIELD_ID needs to be made for
   the different time windows. Each of the different FIELDs can then
   have regular BDA averaging. It may make sense to add back SCAN_ID to
   the selection criteria.

Example
-------

Original MS, before applying BDA:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SPECTRAL_WINDOW table:

+--------------------+---------------------------------+---------------------------+
| SPECTRAL_WINDOW_ID | CHAN_FREQ (MHz)                 | CHAN_WIDTH (MHz)          |
+====================+=================================+===========================+
| 0                  | [100,110,120,130]               | [10,10,10,10]             |
+--------------------+---------------------------------+---------------------------+
| 1                  | [1000,1100,1200,1300,1400,1500] | [100,100,100,100,100,100] |
+--------------------+---------------------------------+---------------------------+

Main table:

+-------------+--------------------+------+------+----------------+-------------+
| TIME        | SPECTRAL_WINDOW_ID | ANT1 | ANT2 | shape(DATA)[0] |             |
+=============+====================+======+======+================+=============+
| 0:00.0      | 0                  | 0    | 1    | 6              | (short      |
|             |                    |      |      |                | baseline)   |
+-------------+--------------------+------+------+----------------+-------------+
| 0:00.0      | 0                  | 0    | 2    | 6              | (long       |
|             |                    |      |      |                | baseline)   |
+-------------+--------------------+------+------+----------------+-------------+
| 0:00.0      | 1                  | 0    | 1    | 6              | (short      |
|             |                    |      |      |                | baseline)   |
+-------------+--------------------+------+------+----------------+-------------+
| 0:00.0      | 1                  | 0    | 2    | 6              | (long       |
|             |                    |      |      |                | baseline)   |
+-------------+--------------------+------+------+----------------+-------------+
| 0:01.0,     | Four rows, similar |      |      |                |             |
| 0:02.0, ... | the four rows      |      |      |                |             |
|             | above, for each    |      |      |                |             |
|             | time step.         |      |      |                |             |
+-------------+--------------------+------+------+----------------+-------------+

After applying BDA to the original MS:
--------------------------------------

In this example, the time averaging factor and the frequency averaging
factor are equal. In reality, they may differ.

BDA_FACTORS table:

+------------------+--------------------+----------+----------+-------------+-----------+
| BDA_TIME_AXIS_ID | BDA_FREQ_AXIS_ID   | ANTENNA1 | ANTENNA2 | TIME_FACTOR |           |
+==================+====================+==========+==========+=============+===========+
| 0                | 0                  | 0        | 1        | 4           | (short    |
|                  |                    |          |          |             | baseline) |
+------------------+--------------------+----------+----------+-------------+-----------+
| 0                | 1                  | 0        | 2        | 3           | (long     |
|                  |                    |          |          |             | baseline) |
+------------------+--------------------+----------+----------+-------------+-----------+
| 1                | 2                  | 0        | 1        | 3           | (short    |
|                  |                    |          |          |             | baseline) |
+------------------+--------------------+----------+----------+-------------+-----------+
| 1                | 3                  | 0        | 2        | 2           | (long     |
|                  |                    |          |          |             | baseline) |
+------------------+--------------------+----------+----------+-------------+-----------+

.. _spectral_window-table-1:

SPECTRAL_WINDOW table:

+--------------------+------------------+------------------+------------------------+
| SPECTRAL_WINDOW_ID | CHAN_FREQ (MHz)  | CHAN_WIDTH (MHz) | BDA_FREQ_AXIS_ID       |
+====================+==================+==================+========================+
| 0                  | [115]            | [40]             | 0                      |
+--------------------+------------------+------------------+------------------------+
| 1                  | [105,125]        | [20, 20]         | 0                      |
+--------------------+------------------+------------------+------------------------+
| 2                  | [1100,1400]      | [300,300]        | 1                      |
+--------------------+------------------+------------------+------------------------+
| 3                  | [1050,1250,1450] | [200,200,200]    | 1                      |
+--------------------+------------------+------------------+------------------------+

.. _main-table-1:

Main table:

This table shows the first two rows for each baseline + spectral window
combination:

+--------+-----------+--------------------+------+------+----------------+-----------+
| TIME   | (Original | SPECTRAL_WINDOW_ID | ANT1 | ANT2 | shape(DATA)[0] |           |
|        | TIMEs)    |                    |      |      |                |           |
+========+===========+====================+======+======+================+===========+
| 0:01.5 | 0,1,2,3   | 0                  | 0    | 1    | 1              | (short    |
|        |           |                    |      |      |                | baseline) |
+--------+-----------+--------------------+------+------+----------------+-----------+
| 0:05.5 | 4,5,6,7   | 0                  | 0    | 1    | 1              | (short    |
|        |           |                    |      |      |                | baseline) |
+--------+-----------+--------------------+------+------+----------------+-----------+
| 0:00.5 | 0,1       | 1                  | 0    | 2    | 2              | (long     |
|        |           |                    |      |      |                | baseline) |
+--------+-----------+--------------------+------+------+----------------+-----------+
| 0:02.5 | 2,3       | 1                  | 0    | 2    | 2              | (long     |
|        |           |                    |      |      |                | baseline) |
+--------+-----------+--------------------+------+------+----------------+-----------+
| 0:01.0 | 0,1,2     | 2                  | 0    | 1    | 2              | (short    |
|        |           |                    |      |      |                | baseline) |
+--------+-----------+--------------------+------+------+----------------+-----------+
| 0:04.0 | 3,4,5     | 2                  | 0    | 1    | 2              | (short    |
|        |           |                    |      |      |                | baseline) |
+--------+-----------+--------------------+------+------+----------------+-----------+
| 0:00.5 | 0,1       | 3                  | 0    | 2    | 3              | (long     |
|        |           |                    |      |      |                | baseline) |
+--------+-----------+--------------------+------+------+----------------+-----------+
| 0:02.5 | 2,3       | 3                  | 0    | 2    | 3              | (long     |
|        |           |                    |      |      |                | baseline) |
+--------+-----------+--------------------+------+------+----------------+-----------+

BDA_TIME_AXIS table (flipped):

========================== ====== ======
BDA_TIME_AXIS_ID           0      1
FIELD_ID                   0      0
BDA_FREQ_AXIS_ID           0      1
IS_BDA_APPLIED             True   True
SINGLE_FACTOR_PER_BASELINE True   True
MAX_TIME_INTERVAL          0:04.0 0:03.0
MIN_TIME_INTERVAL          0:02.0 0:02.0
UNIT_TIME_INTERVAL         0:01.0 0:01.0
INTEGER_INTERVAL_FACTORS   True   True
========================== ====== ======
