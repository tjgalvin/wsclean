#ifndef WSCLEAN_DECONVOLUTION_TABLE_H_
#define WSCLEAN_DECONVOLUTION_TABLE_H_

#include "deconvolutiontableentry.h"

#include <functional>
#include <memory>
#include <vector>

/**
 * The DeconvolutionTable contains DeconvolutionTableEntry's and groups entries
 * that have the same squaredDeconvolutionIndex.
 */
class DeconvolutionTable {
 public:
  using Entries = std::vector<std::unique_ptr<DeconvolutionTableEntry>>;
  using Group = std::vector<const DeconvolutionTableEntry*>;

  /**
   * Iterator class for looping over entries.
   *
   * Dereferencing this iterator yields a reference to the actual object instead
   * of a reference to the pointer for the object.
   */
  class EntryIterator {
    using BaseIterator = Entries::const_iterator;

   public:
    explicit EntryIterator(BaseIterator base_iterator)
        : base_iterator_(base_iterator) {}

    const DeconvolutionTableEntry& operator*() const {
      return **base_iterator_;
    }
    void operator++() { ++base_iterator_; }
    bool operator!=(const EntryIterator& other) const {
      return base_iterator_ != other.base_iterator_;
    }

   private:
    BaseIterator base_iterator_;
  };

  /**
   * @brief Constructs a new DeconvolutionTable object.
   *
   * @param n_original_groups The number of original channel groups. When adding
   * entries, their original channel index must be less than the number of
   * original groups.
   */
  explicit DeconvolutionTable(size_t n_original_groups)
      : entries_(), original_groups_(n_original_groups) {}

  /**
   * @return The table entries, grouped by their original channel index.
   * @see AddEntry()
   */
  const std::vector<Group>& OriginalGroups() const { return original_groups_; }

  EntryIterator begin() const { return EntryIterator(entries_.begin()); }
  EntryIterator end() const { return EntryIterator(entries_.end()); }

  /**
   * @brief Adds an entry to the table.
   *
   * The original channel index of the entry determines the original group for
   * the entry. It must be less than the number of original channel groups, as
   * given in the constructor.
   *
   * @param entry A new entry.
   */
  void AddEntry(std::unique_ptr<DeconvolutionTableEntry> entry);

  /**
   * @return A reference to the first entry.
   */
  const DeconvolutionTableEntry& Front() const { return *entries_.front(); }

  /**
   * @return The number of entries in the table.
   */
  size_t Size() const { return entries_.size(); }

 private:
  Entries entries_;

  /**
   * An original group has entries with equal original channel indices.
   */
  std::vector<Group> original_groups_;
};

#endif
