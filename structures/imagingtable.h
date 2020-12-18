#ifndef WSCLEAN_IMAGING_TABLE_H
#define WSCLEAN_IMAGING_TABLE_H

#include "imagingtableentry.h"

#include <limits>
#include <memory>
#include <vector>

class ImagingTable {
 public:
  size_t IndependentGroupCount() const {
    return _independentGroupLookup.size();
  }

  ImagingTable GetIndependentGroup(size_t index) const;

  size_t SquaredGroupCount() const { return _squaredGroupLookup.size(); }

  ImagingTable GetSquaredGroup(size_t index) const;

  size_t EntryCount() const { return _entries.size(); }

  ImagingTableEntry& operator[](size_t index) { return *_entries[index]; }
  const ImagingTableEntry& operator[](size_t index) const {
    return *_entries[index];
  }

  void Clear() {
    _entries.clear();
    Update();
  }

  ImagingTableEntry& AddEntry() {
    _entries.emplace_back(new ImagingTableEntry());
    return *_entries.back();
  }

  void Update() {
    updateIndependentGroupLookup();
    updateSquaredGroupLookup();
  }

  void Print();

  ImagingTableEntry& Front() { return *_entries.front(); }
  const ImagingTableEntry& Front() const { return *_entries.front(); }

  const ImagingTableEntry* FirstWithHigherFrequency(double frequency) const {
    double currentDistance = std::numeric_limits<double>::max();
    ImagingTableEntry* entry = nullptr;

    for (auto& e : _entries) {
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
    ImagingTableEntry* entry = nullptr;

    for (auto& e : _entries) {
      if (e->CentralFrequency() < frequency &&
          frequency - e->CentralFrequency() < currentDistance) {
        currentDistance = frequency - e->CentralFrequency();
        entry = &*e;
      }
    }
    return entry;
  }

 private:
  void printIndependentGroup(bool isFinal);
  void updateIndependentGroupLookup();
  void updateSquaredGroupLookup();

  typedef std::shared_ptr<ImagingTableEntry> ImagingTableEntryPtr;
  std::vector<ImagingTableEntryPtr> _entries;

  std::vector<std::vector<ImagingTableEntryPtr>> _independentGroupLookup;
  std::vector<std::vector<ImagingTableEntryPtr>> _squaredGroupLookup;
};

#endif
