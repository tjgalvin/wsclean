#include "deconvolutiontable.h"

#include <cassert>

void DeconvolutionTable::AddEntry(
    std::unique_ptr<DeconvolutionTableEntry> entry) {
  entry->index = entries_.size();
  entries_.push_back(std::move(entry));

  size_t group_id = entries_.back()->channel_group_id;
  assert(group_id >= entries_.front()->channel_group_id);
  group_id -= entries_.front()->channel_group_id;

  // As documented in deconvolutiontable.h, group_id may be at most 1 larger
  // than the largest existing group_id.
  assert(group_id <= channel_groups_.size());
  if (group_id == channel_groups_.size()) {
    // Create a new group for the entry.
    // The first entry of a group must have a psf image accessor.
    assert(entries_.back()->psf_accessor);
    channel_groups_.emplace_back(1, entries_.back().get());
  } else {
    // Add the entry to an existing group.
    channel_groups_[group_id].push_back(entries_.back().get());
  }
}
