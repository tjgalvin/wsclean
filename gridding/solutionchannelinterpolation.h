#ifndef WSCLEAN_GRIDDING_SOLUTION_CHANNEL_INTERPOLATION_H_
#define WSCLEAN_GRIDDING_SOLUTION_CHANNEL_INTERPOLATION_H_

#include <algorithm>

#include <aocommon/jaggedarray.h>
#include <aocommon/multibanddata.h>

namespace wsclean {
namespace internal {
/**
 * Creates on row of the table created by @ref MakeChannelInterpolationTable(), which maps
 * channels from the @p to_band onto the @ref from_band.
 */
inline std::vector<uint32_t> MakeChannelInterpolationTableRow(const aocommon::BandData& from_band,
                                                              const aocommon::BandData& to_band) {
  const size_t n_from_channels = from_band.ChannelCount();
  const size_t n_to_channels = to_band.ChannelCount();
  assert(n_from_channels > 0 || n_to_channels == 0);

  size_t from_channel = 0;
  std::vector<uint32_t> table_row(n_to_channels);
  for (size_t to_channel = 0; to_channel != n_to_channels; ++to_channel) {
    // Continue increasing from_channel when the next channel is closer to the
    // destination channel.
    bool at_end = from_channel + 1 >= n_from_channels;
    const double to_freq = to_band.ChannelFrequency(to_channel);
    while (!at_end) {
      const double from_freq = from_band.ChannelFrequency(from_channel);
      const double next_from_freq = from_band.ChannelFrequency(from_channel + 1);
      if (std::abs(from_freq - to_freq) > std::abs(next_from_freq - to_freq)) {
        ++from_channel;
      } else {
        break;
      }
      at_end = from_channel + 1 >= n_from_channels;
    }
    table_row[to_channel] = from_channel;
  }
  return table_row;
}

}  // namespace internal

/**
 * Performs nearest neighbour interpolation for one row of data, with a pre-calculated
 * @p channel_interpolation_table_row. The number of channels in the @p from array is
 * determined from @p channel_interpolation_table_row. The number of parameters that
 * is interpolated is given by @p n_values_per_channel, and this is the fastest changing
 * index in @p to and @p from.
 * @param [in] from The input data, of size n_values_per_channel x n_highres_channels.
 * @param [out] to The interpolated output data, of minimum size
 * n_values_per_channel x n_to_channels.
 */
template <typename ValueType>
inline void InterpolateRowChannels(std::span<uint32_t> channel_interpolation_table_row,
                                   std::span<const ValueType> from, std::span<ValueType> to,
                                   size_t n_values_per_channel) {
  const size_t n_to_channels = channel_interpolation_table_row.size();
  assert(from.size() >= n_to_channels * n_values_per_channel);
  assert(to.size() >= n_to_channels * n_values_per_channel);

  const ValueType* from_ptr = from.data();
  ValueType* to_ptr = to.data();
  for (size_t to_channel = 0; to_channel != n_to_channels; ++to_channel) {
    const uint32_t from_channel = channel_interpolation_table_row[to_channel];
    const ValueType* from_location = from_ptr + from_channel * n_values_per_channel;
    to_ptr = std::copy_n(from_location, n_values_per_channel, to_ptr);
  }
}

/**
 * Specialization of @ref InterpolateRowChannels() for n_values_per_channel=1.
 */
template <typename ValueType>
inline void InterpolateRowChannels(std::span<uint32_t> channel_interpolation_table_row,
                                   std::span<const ValueType> from, std::span<ValueType> to) {
  const size_t n_to_channels = channel_interpolation_table_row.size();
  assert(from.size() >= n_to_channels);
  assert(to.size() >= n_to_channels);

  for (size_t to_channel = 0; to_channel != n_to_channels; ++to_channel) {
    const uint32_t from_channel = channel_interpolation_table_row[to_channel];
    to[to_channel] = from[from_channel];
  }
}

/**
 * Returns a jagged array that can map (data_desc_id, bda_channel_index) to a highres_channel_index.
 * This is used for nearest neighbour interpolation.
 */
inline aocommon::JaggedArray<uint32_t> MakeChannelInterpolationTable(
    const aocommon::MultiBandData& bands) {
  assert(bands.BandCount() > 0);
  const size_t from_data_desc_id = bands.DataDescIdWithMaxChannels();
  const aocommon::BandData from_band(bands[from_data_desc_id]);

  std::vector<std::vector<uint32_t>> table;
  for (size_t data_desc_id : bands.DataDescIds()) {
    if (table.size() <= data_desc_id) table.resize(data_desc_id + 1);
    table[data_desc_id] =
        internal::MakeChannelInterpolationTableRow(from_band, bands[data_desc_id]);
  }
  return aocommon::JaggedArray<uint32_t>(table);
}

}  // namespace wsclean

#endif
