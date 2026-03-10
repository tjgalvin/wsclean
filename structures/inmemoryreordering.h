#ifndef WSCLEAN_STRUCTURES_IN_MEMORY_REORDERING_H_
#define WSCLEAN_STRUCTURES_IN_MEMORY_REORDERING_H_

#include "inmemorypart.h"

#include "../interface/inmemoryms.h"

#include <cassert>
#include <map>
#include <string>
#include <vector>

#include <aocommon/vectormap.h>
#include <aocommon/polarization.h>

#include <schaapcommon/reordering/channelrange.h>

namespace wsclean {

struct InMemoryMs;
struct InMemoryPart;
class Settings;

class InMemoryHandle {
 public:
  InMemoryHandle() noexcept = default;
  InMemoryHandle(const InMemoryHandle&) = delete;
  InMemoryHandle(InMemoryHandle&&) noexcept = default;
  InMemoryHandle& operator=(const InMemoryHandle&) = delete;
  InMemoryHandle& operator=(InMemoryHandle&&) noexcept = default;

  /**
   * Get a part that does or does not exist. This function is typically
   * used during reordering. Use @ref GetExistingPart() after reordering.
   */
  InMemoryPart& GetPart(aocommon::PolarizationEnum polarization,
                        size_t out_channel) {
    return parts_[std::pair(polarization, out_channel)];
  }

  /**
   * Retrieve a part which which is known to exist. This function is typically
   * used by the MsProvider to obtain its already reordered data, and is
   * constant.
   */
  const InMemoryPart& GetExistingPart(aocommon::PolarizationEnum polarization,
                                      size_t out_channel) const {
    const auto iterator = parts_.find(std::pair(polarization, out_channel));
    assert(iterator != parts_.end());
    return iterator->second;
  }

  InMemoryPart& GetExistingPart(aocommon::PolarizationEnum polarization,
                                size_t out_channel) {
    return const_cast<InMemoryPart&>(
        static_cast<const InMemoryHandle*>(this)->GetExistingPart(polarization,
                                                                  out_channel));
  }

 private:
  std::map<std::pair<aocommon::PolarizationEnum, size_t>, InMemoryPart> parts_;
};

/**
 * @param channels has one element for every output channel, where each element
 * has all data desc ids for the corresponding output channel.
 */
InMemoryHandle ReorderInMemory(
    InMemoryMs&& data,
    const std::vector<
        aocommon::VectorMap<schaapcommon::reordering::ChannelRange>>& channels,
    const std::string& data_column_name, const std::string& model_column_name,
    bool set_model_to_zero, const Settings& settings);

}  // namespace wsclean

#endif
