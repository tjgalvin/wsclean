#include "imagingtableentry.h"

ImagingTableEntry::ImagingTableEntry()
    : index(0),
      lowestFrequency(0.0),
      highestFrequency(0.0),
      bandStartFrequency(0.0),
      bandEndFrequency(0.0),
      inputChannelCount(0),
      polarization(aocommon::Polarization::StokesI),
      outputChannelIndex(0),
      outputIntervalIndex(0),
      msData(),
      squaredDeconvolutionIndex(0),
      joinedGroupIndex(0),
      imageCount(0),
      tmpFilePrefix(),
      imageWeight(0.0) {}

std::string ImagingTableEntry::ToString() {
  std::ostringstream str;
  if (index < 10) str << ' ';
  str << index << ' ';
  std::string polStr = aocommon::Polarization::TypeToShortString(polarization);
  if (polStr.size() < 2) str << ' ';
  str << polStr << "  ";
  if (outputChannelIndex < 10) str << ' ';
  str << outputChannelIndex << "  " << joinedGroupIndex << "  "
      << squaredDeconvolutionIndex << "  " << outputIntervalIndex << "  "
      << round(bandStartFrequency * 1e-6) << "-"
      << round(bandEndFrequency * 1e-6) << " (" << inputChannelCount << ")";

  return str.str();
}
