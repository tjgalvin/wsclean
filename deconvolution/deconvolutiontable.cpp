#include "deconvolutiontable.h"

#include <cassert>

void DeconvolutionTable::AddEntry(
    std::unique_ptr<DeconvolutionTableEntry> entry) {
  entry->index = _entries.size();
  _entries.push_back(std::move(entry));

  size_t sqIndex = _entries.back()->squaredDeconvolutionIndex;
  assert(sqIndex >= _entries.front()->squaredDeconvolutionIndex);
  sqIndex -= _entries.front()->squaredDeconvolutionIndex;

  // As documented in deconvolutiontable.h, sqIndex may be at most 1 larger
  // than the largest existing sqIndex.
  assert(sqIndex <= _squaredGroups.size());
  if (sqIndex == _squaredGroups.size()) {
    // Create a new group for the entry.
    _squaredGroups.emplace_back(1, _entries.back().get());
  } else {
    // Add the entry to an existing group.
    _squaredGroups[sqIndex].push_back(_entries.back().get());
  }
}
