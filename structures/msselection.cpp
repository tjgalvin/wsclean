#include "msselection.h"

#include <algorithm>

#include <aocommon/logger.h>
#include <aocommon/multibanddata.h>

using schaapcommon::reordering::ChannelRange;

namespace wsclean {

ChannelRange SelectMsChannels(const aocommon::MultiBandData& ms_bands,
                              size_t data_desc_id, double lowest_frequency,
                              double highest_frequency) {
  const aocommon::BandData& band = ms_bands[data_desc_id];
  double first_channel = band.ChannelFrequency(0);
  double last_channel = band.ChannelFrequency(band.ChannelCount() - 1);
  // Some mses have decreasing (i.e. reversed) channel frequencies in them
  bool is_reversed = false;
  if (first_channel > last_channel) {
    std::swap(first_channel, last_channel);
    is_reversed = true;
    aocommon::Logger::Debug
        << "Warning: MS has reversed channel frequencies.\n";
  }
  if (band.ChannelCount() != 0 && lowest_frequency <= last_channel &&
      highest_frequency >= first_channel) {
    size_t range_start, range_end;
    if (is_reversed) {
      aocommon::BandData::const_reverse_iterator low_ptr =
          std::lower_bound(band.rbegin(), band.rend(), lowest_frequency);
      aocommon::BandData::const_reverse_iterator high_ptr =
          std::upper_bound(low_ptr, band.rend(), highest_frequency);

      range_start = band.ChannelCount() - (high_ptr - band.rbegin());
      range_end = band.ChannelCount() - (low_ptr - band.rbegin());
    } else {
      const aocommon::BandData::const_iterator low_ptr =
          std::lower_bound(band.begin(), band.end(), lowest_frequency);
      aocommon::BandData::const_iterator high_ptr =
          std::upper_bound(low_ptr, band.end(), highest_frequency);

      range_start = low_ptr - band.begin();
      range_end = high_ptr - band.begin();
    }

    return ChannelRange{data_desc_id, range_start, range_end};
  } else {
    return ChannelRange(data_desc_id, 0, 0);
  }
}

}  // namespace wsclean
