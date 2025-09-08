#ifndef MS_SELECTION
#define MS_SELECTION

#include <schaapcommon/reordering/channelrange.h>

namespace aocommon {
class MultiBandData;
}  // namespace aocommon

namespace wsclean {

/**
 * Determines the channel index range that covers all channels of one
 * data_desc_id between @p lowest_frequency and @p highest_frequency, inclusive
 * on both sides. An empty range is returned when no channels match this
 * criterion.
 */
schaapcommon::reordering::ChannelRange SelectMsChannels(
    const aocommon::MultiBandData& ms_bands, size_t data_desc_id,
    double lowest_frequency, double highest_frequency);

}  // namespace wsclean

#endif
