#include "msprovider.h"

#include "msreaders/msreader.h"

#include <aocommon/logger.h>

#include <schaapcommon/reordering/handledata.h>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include "../structures/msselection.h"

using aocommon::Logger;
using schaapcommon::reordering::MSSelection;

namespace wsclean {

MSProvider::~MSProvider() = default;

void MSProvider::ResetModelColumn() {
  std::unique_ptr<MSReader> msReader = MakeReader();
  const std::vector<std::complex<float>> buffer(NMaxChannels() * NPolarizations(), {0.0f, 0.0f});
  while (msReader->CurrentRowAvailable()) {
    // Always overwrite
    constexpr bool addToMS = false;
    WriteModel(buffer.data(), addToMS);
    NextOutputRow();
    msReader->NextInputRow();
  }
}

}  // namespace wsclean
