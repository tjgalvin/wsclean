#include "imagingtableentry.h"
#include "../deconvolution/deconvolutiontableentry.h"
#include "../io/cachedimageaccessor.h"

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
ImagingTableEntry::CreateDeconvolutionEntry(CachedImageSet* psf_images,
                                            CachedImageSet& model_images,
                                            CachedImageSet& residual_images,
                                            bool is_imaginary) const {
  auto entry = boost::make_unique<DeconvolutionTableEntry>();

  entry->index = index;
  entry->band_start_frequency = bandStartFrequency;
  entry->band_end_frequency = bandEndFrequency;
  entry->polarization = polarization;
  entry->output_channel_index = outputChannelIndex;
  entry->output_interval_index = outputIntervalIndex;
  entry->channel_group_id = squaredDeconvolutionIndex;
  entry->image_weight = imageWeight;

  // A PSF accessor is only needed for the first entry of a squared group.
  if (psf_images) {
    entry->psf_accessor = boost::make_unique<CachedImageAccessor>(
        *psf_images, polarization, outputChannelIndex, is_imaginary);
  }
  entry->model_accessor = boost::make_unique<CachedImageAccessor>(
      model_images, polarization, outputChannelIndex, is_imaginary);
  entry->residual_accessor = boost::make_unique<CachedImageAccessor>(
      residual_images, polarization, outputChannelIndex, is_imaginary);

  return entry;
}