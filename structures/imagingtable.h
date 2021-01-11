#ifndef WSCLEAN_IMAGING_TABLE_H
#define WSCLEAN_IMAGING_TABLE_H

#include "imagingtableentry.h"

#include <functional>
#include <memory>
#include <vector>

class ImagingTable {
 public:
  using EntryPtr = std::shared_ptr<ImagingTableEntry>;
  using Group = std::vector<EntryPtr>;
  using Groups = std::vector<Group>;

  /**
   * Iterator class for looping over entries.
   *
   * Dereferencing this iterator yields a reference to the actual object instead
   * of a reference to the shared pointer for the object.
   */
  class EntryIterator {
    using BaseIterator = Group::const_iterator;

   public:
    explicit EntryIterator(BaseIterator baseIt) : _baseIterator(baseIt) {}

    ImagingTableEntry& operator*() { return **_baseIterator; }
    const ImagingTableEntry& operator*() const { return **_baseIterator; }
    void operator++() { ++_baseIterator; }
    bool operator!=(const EntryIterator& other) const {
      return _baseIterator != other._baseIterator;
    }
    bool operator==(const EntryIterator& other) const {
      return _baseIterator == other._baseIterator;
    }

   private:
    BaseIterator _baseIterator;
  };

  ImagingTable() = default;

  size_t IndependentGroupCount() const { return _independentGroups.size(); }

  ImagingTable GetIndependentGroup(size_t index) const {
    return ImagingTable(_independentGroups[index]);
  }

  size_t SquaredGroupCount() const { return _squaredGroups.size(); }

  ImagingTable GetSquaredGroup(size_t index) const {
    return ImagingTable(_squaredGroups[index]);
  }

  const Groups& SquaredGroups() const { return _squaredGroups; }

  size_t FacetGroupCount() const { return _facetGroups.size(); }

  ImagingTable GetFacetGroup(size_t index) const {
    return ImagingTable(_facetGroups[index]);
  }

  const Groups& FacetGroups() const { return _facetGroups; }

  size_t EntryCount() const { return _entries.size(); }

  ImagingTableEntry& operator[](size_t index) { return *_entries[index]; }
  const ImagingTableEntry& operator[](size_t index) const {
    return *_entries[index];
  }

  const EntryIterator begin() const { return EntryIterator(_entries.begin()); }
  EntryIterator begin() { return EntryIterator(_entries.begin()); }

  const EntryIterator end() const { return EntryIterator(_entries.end()); }
  EntryIterator end() { return EntryIterator(_entries.end()); }

  void Clear() {
    _entries.clear();
    Update();
  }

  void AddEntry(std::unique_ptr<ImagingTableEntry> entry) {
    entry->index = _entries.size();
    _entries.push_back(std::move(entry));
  }

  void Update() {
    updateGroups(_independentGroups,
                 [](const ImagingTableEntry& e) { return e.joinedGroupIndex; });
    updateGroups(_facetGroups,
                 [](const ImagingTableEntry& e) { return e.facetIndex; });
    updateGroups(_squaredGroups, [](const ImagingTableEntry& e) {
      return e.squaredDeconvolutionIndex;
    });
  }

  void Print() const;

  ImagingTableEntry& Front() { return *_entries.front(); }
  const ImagingTableEntry& Front() const { return *_entries.front(); }

 private:
  explicit ImagingTable(const Group& entries);

  static void PrintEntry(const ImagingTableEntry& entry);
  void updateGroups(
      Groups& groups,
      std::function<size_t(const ImagingTableEntry&)> getIndex) const;

  Group _entries;

  Groups _independentGroups;
  Groups _facetGroups;
  Groups _squaredGroups;
};

#endif
