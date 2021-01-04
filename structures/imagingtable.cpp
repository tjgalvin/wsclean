#include "imagingtable.h"
#include "../io/logger.h"

#include <iomanip>
#include <map>

ImagingTable::ImagingTable(const std::vector<ImagingTableEntryPtr>& entries)
    : _entries(entries),
      _independentGroupLookup(),
      _facetGroupLookup(),
      _squaredGroupLookup() {
  Update();
}

void ImagingTable::Print() const {
  Logger::Info << "=== IMAGING TABLE ===\n"
                  "         # Pol Ch JG FG ²G In Freq(MHz)\n";
  for (size_t i = 0; i != IndependentGroupCount(); ++i) {
    Logger::Info << "| Independent group:\n";
    const ImagingTable independent = GetIndependentGroup(i);

    for (size_t f = 0; f != independent.FacetGroupCount(); ++f) {
      const ImagingTable facet = independent.GetFacetGroup(f);

      for (size_t s = 0; s != facet.SquaredGroupCount(); ++s) {
        const ImagingTable squared = facet.GetSquaredGroup(s);

        for (size_t e = 0; e != squared._entries.size(); ++e) {
          if (f == 0 && s == 0 && e == 0)
            Logger::Info << "+-";
          else if ((i + 1) == IndependentGroupCount())
            Logger::Info << "  ";
          else
            Logger::Info << "| ";

          if (s == 0 && e == 0)
            Logger::Info << "+-";
          else if ((f + 1) == independent.FacetGroupCount())
            Logger::Info << "  ";
          else
            Logger::Info << "| ";

          if (e == 0)
            Logger::Info << "+-";
          else if ((s + 1) == facet.SquaredGroupCount())
            Logger::Info << "  ";
          else
            Logger::Info << "| ";

          PrintEntry(*squared._entries[e]);
        }
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
  str << std::setw(2) << entry.facetGroupIndex << " ";
  str << std::setw(2) << entry.squaredDeconvolutionIndex << " ";
  str << std::setw(2) << entry.outputIntervalIndex << "  ";
  str << round(entry.bandStartFrequency * 1e-6) << "-"
      << round(entry.bandEndFrequency * 1e-6) << " (" << entry.inputChannelCount
      << ")";

  Logger::Info << "J-" << str.str() << '\n';
}

void ImagingTable::updateGroupLookup(
    GroupLookup& group,
    std::function<size_t(const ImagingTableEntry&)> getIndex) const {
  std::map<size_t, size_t> indexToLookupSet;

  group.clear();
  for (const ImagingTableEntryPtr& e : _entries) {
    size_t groupIndex = getIndex(*e);
    std::map<size_t, size_t>::iterator found =
        indexToLookupSet.find(groupIndex);
    if (found == indexToLookupSet.end()) {
      indexToLookupSet.insert(std::make_pair(groupIndex, group.size()));
      group.emplace_back(1, e);
    } else {
      group[found->second].push_back(e);
    }
  }
}
