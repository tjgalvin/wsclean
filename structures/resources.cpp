#include "resources.h"

#include <atomic>
#include <cassert>
#include <cmath>
#include <unistd.h>  // for sysconf

#include <aocommon/logger.h>

using aocommon::Logger;

namespace {
double RoundOneDecimal(double value) { return std::round(value * 10.0) / 10.0; }
}  // namespace

const Resources Resources::GetPart(size_t part_size) const {
  assert(part_size != 0);
  const size_t n_cpus =
      std::max<size_t>(1, (n_cpus_ + part_size - 1) / part_size);
  const int64_t memory = memory_ / part_size;
  return Resources(n_cpus, memory);
}

int64_t GetAvailableMemory(double memory_fraction, double abs_memory_limit) {
  assert(memory_fraction > 0.0 && memory_fraction <= 1.0);
  // During the first run of this function, some information is reported
  // This needs to be thread safe because gridders can call this function in
  // parallel.
  const int64_t gb = 1024.0 * 1024.0 * 1024.0;
  static std::atomic<bool> isFirst(true);
  const bool print_output = isFirst.exchange(false);
  const long page_count = sysconf(_SC_PHYS_PAGES);
  const long page_size = sysconf(_SC_PAGE_SIZE);
  int64_t memory = (int64_t)page_count * (int64_t)page_size;
  const double memory_size_in_gb = (double)memory / gb;
  if (memory_fraction == 1.0 && abs_memory_limit == 0.0) {
    if (print_output) {
      Logger::Info << "Detected " << RoundOneDecimal(memory_size_in_gb)
                   << " GB of system memory, usage not limited.\n";
    }
  } else {
    if (print_output) {
      double limit_in_gb = memory_size_in_gb * memory_fraction;
      if (abs_memory_limit != 0.0)
        limit_in_gb = std::min(limit_in_gb, abs_memory_limit);
      Logger::Info << "Detected " << RoundOneDecimal(memory_size_in_gb)
                   << " GB of system memory, usage limited to "
                   << RoundOneDecimal(limit_in_gb)
                   << " GB (frac=" << RoundOneDecimal(memory_fraction * 100.0)
                   << "%, ";
      if (abs_memory_limit == 0.0)
        Logger::Info << "no abs limit)\n";
      else
        Logger::Info << "abs limit=" << RoundOneDecimal(abs_memory_limit)
                     << "GB)\n";
    }

    memory = int64_t((double)page_count * (double)page_size * memory_fraction);
    if (abs_memory_limit != 0.0)
      memory = std::min<int64_t>(double(gb) * abs_memory_limit, memory);
  }
  return memory;
}
