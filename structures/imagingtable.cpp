#include "imagingtable.h"
#include "../io/logger.h"

#include <iomanip>
#include <map>

ImagingTable::ImagingTable(const std::vector<EntryPtr>& entries)
    : _entries(entries),
      _independentGroups(),
      _facetGroups(),
      _squaredGroups() {
  Update();
}

void ImagingTable::Print() const {
  Logger::Info << "=== IMAGING TABLE ===\n"
                  "       # Pol Ch JG ²G FG In Freq(MHz)\n";
  for (size_t i = 0; i != IndependentGroupCount(); ++i) {
    Logger::Info << "| Independent group:\n";
    const ImagingTable independent = GetIndependentGroup(i);

    for (size_t s = 0; s != independent.SquaredGroupCount(); ++s) {
      const ImagingTable squared = independent.GetSquaredGroup(s);

      for (size_t e = 0; e != squared._entries.size(); ++e) {
        if (s == 0 && e == 0)
          Logger::Info << "+-";
        else if ((i + 1) == IndependentGroupCount())
          Logger::Info << "  ";
        else
          Logger::Info << "| ";

        if (e == 0)
          Logger::Info << "+-";
        else if ((s + 1) == independent.SquaredGroupCount())
          Logger::Info << "  ";
        else
          Logger::Info << "| ";

        PrintEntry(*squared._entries[e]);
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
