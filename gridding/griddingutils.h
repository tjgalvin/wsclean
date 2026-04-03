#ifndef WSCLEAN_GRIDDING_GRIDDING_UTILS_H_
#define WSCLEAN_GRIDDING_GRIDDING_UTILS_H_

#include <cstddef>
#include <vector>

namespace wsclean {

/**
 * Calculate optimal chunk sizes for pipelined processing.
 *
 * @return A vector of row counts per chunk. The first chunk is sized as a smaller remainder
 * (controlled by memory_overlap_factor) to reduce startup latency. Subsequent chunks are maximized
 * to fit the available memory budget.
 */
std::vector<size_t> CalculateChunkSizes(size_t total_rows, size_t max_rows_in_memory,
                                        float memory_overlap_factor);

}  // namespace wsclean

#endif  // WSCLEAN_GRIDDING_GRIDDING_UTILS_H_
