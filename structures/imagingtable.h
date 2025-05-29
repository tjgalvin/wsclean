#ifndef WSCLEAN_IMAGING_TABLE_H
#define WSCLEAN_IMAGING_TABLE_H

#include "imagingtableentry.h"

#include <functional>
#include <memory>
#include <vector>

#include <radler/work_table.h>

namespace wsclean {

/**
 * The ImagingTable contains ImagingTableEntry's and supports creating subtables
 * with entries that have the same attribute.
 * It supports the following grouping types:
 * - IndependentGroup : Entries have an equal joinedGroupIndex in each group.
 * - FacetGroup : Entries have an equal facetGroupIndex in each group.
 *   Each group thus contains all entries that form a single image.
 * - Facet : Entries have an equal facetIndex in each group.
 *   Each group thus contains all entries for a single facet.
 * - SquaredGroup : Entries have an equal squaredDeconvolutionIndex in each
 *   group.
 */
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
    ImagingTableEntry* operator->() { return _baseIterator->get(); }
    const ImagingTableEntry* operator->() const { return _baseIterator->get(); }
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

  ImagingTable(const ImagingTable& other,
               std::function<bool(const ImagingTableEntry&)> isSelected);

  explicit ImagingTable(const Group& entries);

  size_t IndependentGroupCount() const { return _independentGroups.size(); }

  ImagingTable GetIndependentGroup(size_t index) const {
    return ImagingTable(_independentGroups[index]);
  }

  size_t SquaredGroupCount() const { return _squaredGroups.size(); }

  ImagingTable GetSquaredGroup(size_t index) const {
    return ImagingTable(_squaredGroups[index]);
  }

  const Groups& SquaredGroups() const { return _squaredGroups; }

  // When an imagingtable is split into different tables, a facetGroupIndex may
  // be greater-equal than FacetGroupCount(). Since facetGroupIndices are used
  // for acquiring scheduler locks, always use (MaxFacetGroupIndex() + 1) for
  // determining the number of scheduler locks.
  size_t MaxFacetGroupIndex() const { return _entries.back()->facetGroupIndex; }

  /**
   * Groups the entries per facet group.
   * Each group holds all entries for one image.
   * @param is_selected Selection function. The default selection
   * function omits direction dependent PSF entries from the result.
   */
  Groups FacetGroups(std::function<bool(const ImagingTableEntry&)> is_selected =
                         [](const ImagingTableEntry& e) {
                           return !e.isDdPsf;
                         }) const {
    return CreateGroups(
        [](const ImagingTableEntry& e) { return e.facetGroupIndex; },
        std::move(is_selected));
  }

  size_t FacetCount() const { return _facets.size(); }

  ImagingTable GetFacet(size_t index) const {
    return ImagingTable(_facets[index]);
  }

  const Groups& Facets() const { return _facets; }

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
    _independentGroups.clear();
    _facets.clear();
    _squaredGroups.clear();
  }

  void AddEntry(std::unique_ptr<ImagingTableEntry> entry) {
    entry->index = _entries.size();
    _entries.push_back(std::move(entry));
  }

  void Update() {
    _independentGroups = CreateGroups(
        [](const ImagingTableEntry& e) { return e.joinedGroupIndex; });
    _facets =
        CreateGroups([](const ImagingTableEntry& e) { return e.facetIndex; },
                     [](const ImagingTableEntry& e) { return !e.isDdPsf; });
    _squaredGroups = CreateGroups(
        [](const ImagingTableEntry& e) { return e.squaredDeconvolutionIndex; });
  }

  void Print() const;

  ImagingTableEntry& Front() { return *_entries.front(); }
  const ImagingTableEntry& Front() const { return *_entries.front(); }

  /**
   * Copies all fields that are set by the gridding results from one
   * polarization to all other polarizations.
   */
  void AssignGridDataFromPolarization(aocommon::PolarizationEnum polarization);

  std::unique_ptr<radler::WorkTable> CreateDeconvolutionTable(
      int n_deconvolution_channels, CachedImageSet& psf_images,
      CachedImageSet& model_images, CachedImageSet& residual_images) const;

 private:
  static void PrintEntry(const ImagingTableEntry& entry);

  Groups CreateGroups(
      std::function<size_t(const ImagingTableEntry&)> getIndex,
      std::function<bool(const ImagingTableEntry&)> isSelected =
          [](const ImagingTableEntry&) { return true; }) const;

  Group _entries;

  Groups _independentGroups;
  Groups _facets;
  /**
   * Item i will contain a Group for which each entry has a
   * squaredDeconvolutionIndex equal to i.
   */
  Groups _squaredGroups;
};

}  // namespace wsclean

#endif
