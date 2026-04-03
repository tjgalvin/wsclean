#include "griddingutils.h"

#include <cassert>

namespace wsclean {

std::vector<size_t> CalculateChunkSizes(size_t total_rows, size_t max_rows_in_memory,
                                        float memory_overlap_factor) {
  assert(total_rows > 0);
  assert(max_rows_in_memory > 0);
  assert(memory_overlap_factor >= 1.0f);

  const size_t desired_rows = std::min(max_rows_in_memory, total_rows);

  const size_t ideal_rows_per_chunk = std::max<size_t>(desired_rows / memory_overlap_factor, 1);

  const size_t total_chunks = (total_rows + ideal_rows_per_chunk - 1) / ideal_rows_per_chunk;
  const size_t n_rows_in_first_chunk = total_rows - (ideal_rows_per_chunk * (total_chunks - 1));

  std::vector<size_t> chunk_sizes(total_chunks, ideal_rows_per_chunk);
  chunk_sizes[0] = n_rows_in_first_chunk;

  return chunk_sizes;
}

}  // namespace wsclean
