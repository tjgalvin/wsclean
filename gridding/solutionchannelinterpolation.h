#ifndef WSCLEAN_GRIDDING_SOLUTION_CHANNEL_INTERPOLATION_H_
#define WSCLEAN_GRIDDING_SOLUTION_CHANNEL_INTERPOLATION_H_

#include <algorithm>

#include <aocommon/multibanddata.h>

namespace wsclean {

/**
 * Performs nearest neighbour interpolation in frequency direction.
 * The number of values in the from vector is expected to be
 * n_sets x from_band.ChannelCount() x n_values_per_channel. Typically (i.e., in
 * the visibility modifier), n_values_per_channel is set to the number of
 * solution parameter x n_antennas whereas n_sets is set to n_times.
 * @param [out] to is resized and filled with the interpolated values.
 */
template <typename ValueType>
inline void InterpolateChannels(const std::vector<ValueType>& from,
                                const aocommon::BandData& from_band,
                                std::vector<ValueType>& to,
                                const aocommon::BandData& to_band,
                                size_t n_values_per_channel, size_t n_sets) {
  const size_t n_from_channels = from_band.ChannelCount();
  const size_t n_to_channels = to_band.ChannelCount();
  // The 'from' band should not be empty, unless the to_band is also empty (in
  // which case the function does nothing).
  assert(n_from_channels > 0 || n_to_channels == 0);
  assert(n_sets * n_from_channels * n_values_per_channel == from.size());
  to.resize(n_sets * n_to_channels * n_values_per_channel);

  for (size_t set_index = 0; set_index != n_sets; ++set_index) {
    const ValueType* from_ptr =
        &from[set_index * n_from_channels * n_values_per_channel];
    ValueType* to_ptr = &to[set_index * n_to_channels * n_values_per_channel];
    size_t from_channel = 0;
    for (size_t to_channel = 0; to_channel != n_to_channels; ++to_channel) {
      // Continue increasing from_channel when the next channel is closer to the
      // destination channel.
      bool at_end = from_channel + 1 >= n_from_channels;
      const double to_freq = to_band.ChannelFrequency(to_channel);
      while (!at_end) {
        const double from_freq = from_band.ChannelFrequency(from_channel);
        const double next_from_freq =
            from_band.ChannelFrequency(from_channel + 1);
        if (std::abs(from_freq - to_freq) >
            std::abs(next_from_freq - to_freq)) {
          ++from_channel;
        } else {
          break;
        }
        at_end = from_channel + 1 >= n_from_channels;
      }
      const ValueType* from_location =
          from_ptr + from_channel * n_values_per_channel;
      to_ptr = std::copy_n(from_location, n_values_per_channel, to_ptr);
    }
  }
}

}  // namespace wsclean

#endif
