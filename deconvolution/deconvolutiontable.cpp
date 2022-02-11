#include "deconvolutiontable.h"

#include <cassert>

void DeconvolutionTable::AddEntry(
    std::unique_ptr<DeconvolutionTableEntry> entry) {
  const size_t original_channel_index = entry->original_channel_index;
  assert(original_channel_index < original_groups_.size());

  entry->index = entries_.size();
  entries_.push_back(std::move(entry));

  original_groups_[original_channel_index].push_back(entries_.back().get());
}
