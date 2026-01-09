#include "msprovider.h"

#include "msreaders/msreader.h"

#include <aocommon/logger.h>

#include <schaapcommon/reordering/handledata.h>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include "../structures/msselection.h"

using aocommon::Logger;
using schaapcommon::reordering::MSSelection;
using schaapcommon::reordering::StorageManagerType;

namespace wsclean {

aocommon::MultiBandData MakeSelectedPartBands(
    const aocommon::MultiBandData& input,
    const aocommon::VectorMap<schaapcommon::reordering::ChannelRange>& ranges) {
  aocommon::MultiBandData result;
  for (const schaapcommon::reordering::ChannelRange& range : ranges) {
    if (!range.Empty()) {
      const aocommon::BandData& band = input[range.data_desc_id];
      result.SetBand(range.data_desc_id,
                     aocommon::BandData(band, range.start, range.end));
    }
  }
  return result;
}

std::vector<aocommon::MultiBandData> MakeSelectedBands(
    const aocommon::MultiBandData& input,
    const std::vector<
        aocommon::VectorMap<schaapcommon::reordering::ChannelRange>>&
        channels) {
  std::vector<aocommon::MultiBandData> result;
  for (const aocommon::VectorMap<schaapcommon::reordering::ChannelRange>&
           ranges : channels) {
    result.emplace_back(MakeSelectedPartBands(input, ranges));
  }
  return result;
}

bool HasFrequencyBda(const casacore::MeasurementSet& ms) {
  // This is a simple test that works in all current use-cases. Formally, we
  // should check if there are rows with the same BDA_SET_ID to conclude BDA is
  // enabled.
  return ms.spectralWindow().actualTableDesc().isColumn("BDA_SET_ID");
}

MSProvider::~MSProvider() = default;

void MSProvider::ResetModelColumn() {
  std::unique_ptr<MSReader> msReader = MakeReader();
  const std::vector<std::complex<float>> buffer(
      NMaxChannels() * NPolarizations(), {0.0f, 0.0f});
  while (msReader->CurrentRowAvailable()) {
    // Always overwrite
    const bool addToMS = false;
    WriteModel(buffer.data(), addToMS);
    NextOutputRow();
    msReader->NextInputRow();
  }
}

}  // namespace wsclean
