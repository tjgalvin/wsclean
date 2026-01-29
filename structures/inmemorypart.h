#ifndef WSCLEAN_INTERFACE_IN_MEMORY_PART_H_
#define WSCLEAN_INTERFACE_IN_MEMORY_PART_H_

#include <array>
#include <complex>
#include <cstdint>
#include <string>
#include <vector>

#include <aocommon/banddata.h>
#include <aocommon/multibanddata.h>
#include <aocommon/polarization.h>

#include "observationinfo.h"
#include "../msproviders/msprovider.h"  // for MSProvider::MetaData

namespace wsclean {

struct InMemoryPartRow {
  std::vector<std::complex<float>> data;
  std::vector<std::complex<float>> model_data;
  std::vector<float> weights;
};

struct InMemoryInfo {
  double start_time;
  double interval;
  std::vector<std::string> antenna_names;
  bool has_frequency_bda;
  ObservationInfo observation_info;
};

/**
 * An InMemoryPart is a memory buffer for the data from a measurement set,
 * after reordering and data selection. It contains one data part, i.e.
 * the data for a single output channel and single output polarization.
 * It might still hold multiple polarizations though.
 */
struct InMemoryPart {
  std::string description;
  std::vector<InMemoryPartRow> rows;
  std::shared_ptr<std::vector<MSProvider::MetaData>> meta_data;
  std::shared_ptr<InMemoryInfo> info;
  aocommon::MultiBandData selected_bands;
  bool is_regular;
  aocommon::PolarizationEnum polarization;
};

}  // namespace wsclean

#endif
