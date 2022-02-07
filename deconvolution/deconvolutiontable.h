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
    explicit EntryIterator(BaseIterator baseIt) : _baseIterator(baseIt) {}

    const DeconvolutionTableEntry& operator*() const { return **_baseIterator; }
    void operator++() { ++_baseIterator; }
    bool operator!=(const EntryIterator& other) const {
      return _baseIterator != other._baseIterator;
    }

   private:
    BaseIterator _baseIterator;
  };

  const std::vector<Group>& SquaredGroups() const { return _squaredGroups; }

  EntryIterator begin() const { return EntryIterator(_entries.begin()); }
  EntryIterator end() const { return EntryIterator(_entries.end()); }

  /**
   * @brief Adds an entry to the table.
   *
   * When adding multiple entries, these restrictions apply to the
   * the squaredDeconvolutionIndex ('sqIndex') of the entries:
   * - The sqIndex must be greater than or equal to the sqIndex of the first
   *   entry. (The first entry must have the smallest sqIndex of all entries.)
   * - The sqIndex must be at most 1 larger than the sqIndex of any of the
   *   previously added entries.
   * For example, valid indices are: 0-1-2-0-1-2-0-1-2 or 4-4-4-5-5-5-6-6-6.
   *
   * @param entry A new entry.
   */
  void AddEntry(std::unique_ptr<DeconvolutionTableEntry> entry);

  /**
   * @return A reference to the first entry.
   */
  const DeconvolutionTableEntry& Front() const { return *_entries.front(); }

  /**
   * @return The number of entries in the table.
   */
  size_t Size() const { return _entries.size(); }

 private:
  Entries _entries;
  std::vector<Group> _squaredGroups;
};

#endif
