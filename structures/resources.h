#ifndef STRUCTURES_RESOURCES_H_
#define STRUCTURES_RESOURCES_H_

#include <cstdint>
#include <cstddef>

class Resources {
 public:
  Resources() = default;
  Resources(std::size_t n_cpus, int64_t memory)
      : n_cpus_(n_cpus), memory_(memory) {}

  const size_t NCpus() const { return n_cpus_; }
  const int64_t Memory() const { return memory_; }

  const Resources GetPart(size_t part_size) const;

 private:
  std::size_t n_cpus_;
  int64_t memory_;
};

/**
 * Obtain available system memory, taking into account user-requested limits on
 * the usage. If both memory_fraction and abs_memory_limit reduce the amount of
 * available memory, they're combined by taking the strongest limit of the two.
 * They're not applied cummulative.
 * @param memory_fraction Value (0-1) to limit relative memory.
 * @param abs_memory_limit Absolute limit in gb, or 0 if no absolute limit is
 * set.
 * @returns available memory in bytes.
 */
int64_t GetAvailableMemory(double memory_fraction, double abs_memory_limit);

#endif
