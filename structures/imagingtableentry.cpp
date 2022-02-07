#include "imagingtableentry.h"
#include "../deconvolution/deconvolutiontableentry.h"

#include <boost/make_unique.hpp>

ImagingTableEntry::ImagingTableEntry()
    : index(0),
      lowestFrequency(0.0),
      highestFrequency(0.0),
      bandStartFrequency(0.0),
      bandEndFrequency(0.0),
      inputChannelCount(0),
      polarization(aocommon::Polarization::StokesI),
      facetGroupIndex(0),
      facetIndex(0),
      facet(nullptr),
      centreShiftX(0),
      centreShiftY(0),
      outputChannelIndex(0),
      outputIntervalIndex(0),
      msData(),
      squaredDeconvolutionIndex(0),
      joinedGroupIndex(0),
      imageCount(0),
      tmpFilePrefix(),
      imageWeight(0.0) {}

std::unique_ptr<DeconvolutionTableEntry>
ImagingTableEntry::CreateDeconvolutionEntry(bool isImaginary) const {
  auto entry = boost::make_unique<DeconvolutionTableEntry>();

  entry->index = index;
  entry->bandStartFrequency = bandStartFrequency;
  entry->bandEndFrequency = bandEndFrequency;
  entry->polarization = polarization;
  entry->outputChannelIndex = outputChannelIndex;
  entry->outputIntervalIndex = outputIntervalIndex;
  entry->squaredDeconvolutionIndex = squaredDeconvolutionIndex;
  entry->isImaginary = isImaginary;
  entry->imageWeight = imageWeight;

  return entry;
}