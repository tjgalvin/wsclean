#include "imagingtable.h"
#include "../io/logger.h"

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <map>

#include <boost/make_unique.hpp>

ImagingTable::ImagingTable(const std::vector<EntryPtr>& entries)
    : _entries(entries),
      _independentGroups(),
      _facetGroups(),
      _facets(),
      _squaredGroups() {
  Update();
}

void ImagingTable::Print() const {
  Logger::Info << "=== IMAGING TABLE ===\n"
                  "       # Pol Ch JG ²G FG FI In Freq(MHz)\n";
  for (size_t i = 0; i != _independentGroups.size(); ++i) {
    Logger::Info << "| Independent group:\n";
    const ImagingTable independent(_independentGroups[i]);

    const ImagingTable::Groups& squaredGroups = independent.SquaredGroups();
    for (size_t s = 0; s != squaredGroups.size(); ++s) {
      const ImagingTable::Group& squared = squaredGroups[s];

      for (size_t e = 0; e != squared.size(); ++e) {
        if (s == 0 && e == 0)
          Logger::Info << "+-";
        else if ((i + 1) == _independentGroups.size())
          Logger::Info << "  ";
        else
          Logger::Info << "| ";

        if (e == 0)
          Logger::Info << "+-";
        else if ((s + 1) == squaredGroups.size())
          Logger::Info << "  ";
        else
          Logger::Info << "| ";

        PrintEntry(*squared[e]);
      }
    }
  }
}

void ImagingTable::PrintEntry(const ImagingTableEntry& entry) {
  std::ostringstream str;

  str << std::setw(2) << entry.index << "  ";
  str << aocommon::Polarization::TypeToShortString(entry.polarization) << "  ";
  str << std::setw(2) << entry.outputChannelIndex << " ";
  str << std::setw(2) << entry.joinedGroupIndex << " ";
  str << std::setw(2) << entry.squaredDeconvolutionIndex << " ";
  str << std::setw(2) << entry.facetGroupIndex << " ";
  str << std::setw(2) << entry.facetIndex << " ";
  str << std::setw(2) << entry.outputIntervalIndex << "  ";
  str << round(entry.bandStartFrequency * 1e-6) << "-"
      << round(entry.bandEndFrequency * 1e-6) << " (" << entry.inputChannelCount
      << ")";

  Logger::Info << "J-" << str.str() << '\n';
}

void ImagingTable::updateGroups(
    Groups& groups,
    std::function<size_t(const ImagingTableEntry&)> getIndex) const {
  std::map<size_t, Group> groupMap;

  for (const EntryPtr& e : _entries) {
    groupMap[getIndex(*e)].push_back(e);
  }

  groups.clear();
  for (auto& item : groupMap) {
    groups.emplace_back(std::move(item.second));
  }
}

void ImagingTable::AssignGridDataFromPolarization(
    aocommon::PolarizationEnum polarization) {
  for (Group& group : _squaredGroups) {
    const EntryPtr& sourceEntry = *std::find_if(
        group.begin(), group.end(), [polarization](const EntryPtr& e) {
          return e->polarization == polarization;
        });
    for (EntryPtr& entryPtr : group) {
      if (entryPtr != sourceEntry) {
        entryPtr->AssignGridData(*sourceEntry);
      }
    }
  }
}

std::unique_ptr<DeconvolutionTable> ImagingTable::CreateDeconvolutionTable(
    CachedImageSet& psf_images, CachedImageSet& model_images,
    CachedImageSet& residual_images) const {
  auto table = boost::make_unique<DeconvolutionTable>();

  int max_squared_index = -1;

  for (const EntryPtr& entry_ptr : _entries) {
    assert(entry_ptr);

    if (entry_ptr->imageCount >= 1) {
      CachedImageSet* psf_images_ptr = nullptr;

      // Only set psf_images_ptr for the first entry of each squared group.
      // This way, CreateDeconvolutionEntry() only creates psf accessors for
      // those entries.
      if (int(entry_ptr->squaredDeconvolutionIndex) > max_squared_index) {
        max_squared_index = entry_ptr->squaredDeconvolutionIndex;
        psf_images_ptr = &psf_images;
      }

      std::unique_ptr<DeconvolutionTableEntry> real_entry =
          entry_ptr->CreateDeconvolutionEntry(psf_images_ptr, model_images,
                                              residual_images, false);
      table->AddEntry(std::move(real_entry));
    }

    if (entry_ptr->imageCount == 2) {
      std::unique_ptr<DeconvolutionTableEntry> imaginary_entry =
          entry_ptr->CreateDeconvolutionEntry(nullptr, model_images,
                                              residual_images, true);
      table->AddEntry(std::move(imaginary_entry));
    }
  }
  return table;
}