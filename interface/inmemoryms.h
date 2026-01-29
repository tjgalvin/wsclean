#ifndef WSCLEAN_INTERFACE_IN_MEMORY_MS_H_
#define WSCLEAN_INTERFACE_IN_MEMORY_MS_H_

#include <array>
#include <complex>
#include <cstdint>
#include <string>
#include <vector>

#include <aocommon/multibanddata.h>
#include <aocommon/polarization.h>

namespace wsclean {

/**
 * Single row of data, including meta data. Comparable to a row in a measurement
 * set. A row contains the data for one baseline, one timestep.
 */
struct InMemoryRow {
  /**
   * Flattened data array of size polarizations x frequencies, where
   * polarization is the fastest changing. The number of polarizations should
   * match with the length of InMemoryMs::polarizations. The number of channels
   * must match with the number of channels in InMemoryMs::bands[data_desc_id].
   * The number of channels may be variable, e.g. when baseline-dependent
   * averaging is used or when the observations covers multiple bands with
   * different channel counts.
   */
  std::vector<std::complex<float>> data;
  /**
   * Weights associated with the @c data values, with the same shape. For
   * flagged values, the weight should be set to zero.
   */
  std::vector<float> weights;
  std::array<double, 3> uvw{0.0, 0.0, 0.0};
  double time = 0.0;
  uint32_t data_desc_id = 0;
  uint32_t field_id = 0;
  uint32_t antenna1 = 0;
  uint32_t antenna2 = 0;
};

/**
 * An InMemoryMs is a memory buffer for the data for a measurement set. It
 * contains only the data that is relevant for WSClean, being the observed
 * visibility data and some of the meta data required for imaging.
 */
struct InMemoryMs {
  /**
   * Each row contain the data for one baseline, one timestep, with all channels
   * and all polarizations.
   */
  std::vector<InMemoryRow> rows;
  /**
   * Righ ascension in radians.
   */
  double phase_centre_ra = 0.0;
  /**
   * Declination in radians.
   */
  double phase_centre_dec = 0.0;
  /**
   * Bands that are used. The InMemoryRow::data_desc_id are indices into this
   * structure.
   */
  aocommon::MultiBandData bands;
  /**
   * Antenna names are currently only used when solutions are applied, to link
   * data in the h5 solution files to the right antennas. Note that the concept
   * of an "antenna" here is the same as an antenna in a measurement set. For
   * e.g. LOFAR, a station is an antenna. InMemoryRow::antenna1 and 2 are
   * indices into this structure.
   */
  std::vector<std::string> antenna_names;
  /**
   * List of polarizations in the data, e.g. {XX, XY, YX, YY}.
   */
  std::vector<aocommon::PolarizationEnum> polarizations;
  /**
   * Should be set to true when InMemoryRow::data_desc_id is used to distinguish
   * different averaging factors for frequency-direction baseline-dependent
   * averaging.
   */
  bool has_frequency_bda = false;
  /**
   * Could be in the future used as diagnostic; currently not used. Time BDA is
   * supported even if set to false.
   */
  bool has_time_bda = false;
  /**
   * Copied to output FITS files, e.g. "LOFAR" or "JVLA".
   */
  std::string telescope_name;
  /**
   * Copied to output FITS files, e.g. "M. Veldhuis".
   */
  std::string observer;
  /**
   * Copied to output FITS files, e.g. "3C196" or "LOTSS P113+24".
   */
  std::string field_name;
  /**
   * Copied to output FITS files, time in MJD seconds.
   */
  double start_time = 0.0;
  /**
   * Number of seconds between adjacent timesteps. Only relevant when
   * no time-direction baseline-dependent averaging is used. It is
   * used for determining time-smearing (only when explicitly enabled).
   */
  double interval = 0.0;
};

}  // namespace wsclean

#endif
