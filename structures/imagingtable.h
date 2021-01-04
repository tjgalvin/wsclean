#ifndef WSCLEAN_IMAGING_TABLE_H
#define WSCLEAN_IMAGING_TABLE_H

#include "imagingtableentry.h"

#include <functional>
#include <limits>
#include <memory>
#include <vector>

class ImagingTable {
 public:
  /**
   * Iterator class for looping over entries.
   *
   * Dereferencing this iterator yields a reference to the actual object instead
   * of a reference to the shared pointer for the object.
   */
  class EntryIterator {
    using BaseIterator =
        std::vector<std::shared_ptr<ImagingTableEntry>>::const_iterator;

   public:
    explicit EntryIterator(BaseIterator baseIt) : _baseIterator(baseIt) {}

    ImagingTableEntry& operator*() { return **_baseIterator; }
    const ImagingTableEntry& operator*() const { return **_baseIterator; }
    void operator++() { ++_baseIterator; }
    bool operator!=(const EntryIterator& other) const {
      return _baseIterator != other._baseIterator;
    }

   private:
    BaseIterator _baseIterator;
  };

  ImagingTable() = default;

  size_t IndependentGroupCount() const {
    return _independentGroupLookup.size();
  }

  ImagingTable GetIndependentGroup(size_t index) const {
    return ImagingTable(_independentGroupLookup[index]);
  }

  size_t FacetGroupCount() const { return _facetGroupLookup.size(); }

  ImagingTable GetFacetGroup(size_t index) const {
    return ImagingTable(_facetGroupLookup[index]);
  }

  size_t SquaredGroupCount() const { return _squaredGroupLookup.size(); }

  ImagingTable GetSquaredGroup(size_t index) const {
    return ImagingTable(_squaredGroupLookup[index]);
  }

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
    _entries.push_back(std::move(entry));
  }

  void Update() {
    updateGroupLookup(_independentGroupLookup, [](const ImagingTableEntry& e) {
      return e.joinedGroupIndex;
    });
    updateGroupLookup(_facetGroupLookup, [](const ImagingTableEntry& e) {
      return e.facetGroupIndex;
    });
    updateGroupLookup(_squaredGroupLookup, [](const ImagingTableEntry& e) {
      return e.squaredDeconvolutionIndex;
    });
  }

  void Print() const;

  ImagingTableEntry& Front() { return *_entries.front(); }
  const ImagingTableEntry& Front() const { return *_entries.front(); }

  const ImagingTableEntry* FirstWithHigherFrequency(double frequency) const {
    double currentDistance = std::numeric_limits<double>::max();
    const ImagingTableEntry* entry = nullptr;

    for (const auto& e : _entries) {
      if (e->CentralFrequency() > frequency &&
          e->CentralFrequency() - frequency < currentDistance) {
        currentDistance = e->CentralFrequency() - frequency;
        entry = &*e;
      }
    }
    return entry;
  }

  const ImagingTableEntry* FirstWithLowerFrequency(double frequency) const {
    double currentDistance = std::numeric_limits<double>::max();
    const ImagingTableEntry* entry = nullptr;

    for (const auto& e : _entries) {
      if (e->CentralFrequency() < frequency &&
          frequency - e->CentralFrequency() < currentDistance) {
        currentDistance = frequency - e->CentralFrequency();
        entry = &*e;
      }
    }
    return entry;
  }

 private:
  using ImagingTableEntryPtr = std::shared_ptr<ImagingTableEntry>;
  using GroupLookup = std::vector<std::vector<ImagingTableEntryPtr>>;

  explicit ImagingTable(const std::vector<ImagingTableEntryPtr>& entries);

  static void PrintEntry(const ImagingTableEntry& entry);
  void updateGroupLookup(
      GroupLookup& group,
      std::function<size_t(const ImagingTableEntry&)> getIndex) const;

  std::vector<ImagingTableEntryPtr> _entries;

  GroupLookup _independentGroupLookup;
  GroupLookup _facetGroupLookup;
  GroupLookup _squaredGroupLookup;
};

#endif
