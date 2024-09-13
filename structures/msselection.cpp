#include "msselection.h"

#include <aocommon/logger.h>

#include <limits>

namespace wsclean {

bool SelectMsChannels(schaapcommon::reordering::MSSelection& selection,
                      const aocommon::MultiBandData& msBands, size_t dataDescId,
                      const ImagingTableEntry& entry) {
  const aocommon::BandData& band = msBands[dataDescId];
  double firstCh = band.ChannelFrequency(0);
  double lastCh = band.ChannelFrequency(band.ChannelCount() - 1);
  // Some mses have decreasing (i.e. reversed) channel frequencies in them
  bool isReversed = false;
  if (firstCh > lastCh) {
    std::swap(firstCh, lastCh);
    isReversed = true;
    aocommon::Logger::Debug
        << "Warning: MS has reversed channel frequencies.\n";
  }
  if (band.ChannelCount() != 0 && entry.lowestFrequency <= lastCh &&
      entry.highestFrequency >= firstCh) {
    size_t newStart, newEnd;
    if (isReversed) {
      aocommon::BandData::const_reverse_iterator lowPtr =
          std::lower_bound(band.rbegin(), band.rend(), entry.lowestFrequency);
      aocommon::BandData::const_reverse_iterator highPtr =
          std::lower_bound(lowPtr, band.rend(), entry.highestFrequency);

      if (highPtr == band.rend()) --highPtr;
      newStart = band.ChannelCount() - 1 - (highPtr - band.rbegin());
      newEnd = band.ChannelCount() - (lowPtr - band.rbegin());
    } else {
      const double *lowPtr, *highPtr;
      lowPtr =
          std::lower_bound(band.begin(), band.end(), entry.lowestFrequency);
      highPtr = std::lower_bound(lowPtr, band.end(), entry.highestFrequency);

      if (highPtr == band.end()) --highPtr;
      newStart = lowPtr - band.begin();
      newEnd = highPtr - band.begin() + 1;
    }

    selection.SetBandId(dataDescId);
    selection.SetChannelRange(newStart, newEnd);
    return true;
  } else {
    return false;
  }
}

}  // namespace wsclean
