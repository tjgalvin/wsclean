#ifndef WSCLEAN_IMAGING_TABLE_ENTRY_H
#define WSCLEAN_IMAGING_TABLE_ENTRY_H

#include <aocommon/polarization.h>

#include <radler/work_table.h>

#include <memory>
#include <string>
#include <vector>

namespace radler {
struct DeconvolutionTableEntry;
}

namespace schaapcommon {
namespace facets {
class Facet;
}
}  // namespace schaapcommon

class CachedImageSet;

struct ImagingTableEntry {
  struct MSBandInfo {
    size_t bandIndex;
    size_t partIndex;
  };

  struct MSInfo {
    std::vector<MSBandInfo> bands;
  };

  ImagingTableEntry();

  /**
   * @brief Creates a DeconvolutionTableEntry for the ImagingTableEntry.
   *
   * Copies all necessary information from the ImagingTableEntry into a newly
   * created DeconvolutionTableEntry. Creates CachedImageAccessors for the new
   * entry, and initalizes them using the given CachedImageSets.
   * Creating the PSF image accessor is optional, since it is not always needed.
   * Image accessors for the model and residual images are always created.
   *
   * @param channel_index_offset Index of the first channel in the ImagingTable.
   * @param psf_images Pointer to a CachedImageSet for PSF images. If this
   * pointer is null, the created DeconvolutionTableEntry will have no
   * ImageAccessor for a PSF image.
   * @param model_images CachedImageSet for model images.
   * @param residual_images CachedImageSet for residual images.
   * @param is_imaginary False: Create a DeconvolutionTableEntry for an image
   * with real values. True: Create a DeconvolutionTableEntry for an image with
   * imaginary values.
   * @return A new DeconvolutionTableEntry.
   */
  std::unique_ptr<radler::WorkTableEntry> CreateDeconvolutionEntry(
      size_t channel_index_offset, CachedImageSet* psf_images,
      CachedImageSet& model_images, CachedImageSet& residual_images,
      bool is_imaginary) const;

  /**
   * Unique index of the entry within its ImagingTable.
   */
  size_t index;

  /**
   * Note that mses might have overlapping frequencies.
   */
  double lowestFrequency, highestFrequency;
  double bandStartFrequency, bandEndFrequency;
  double siCorrection;
  size_t inputChannelCount;

  aocommon::PolarizationEnum polarization;

  /**
   * All facets that belong to one image have the same facet group index.
   * ImagingTable uses this index for creating facet groups.
   */
  size_t facetGroupIndex;

  /**
   * Index of the entry within a facet group, which equals the index of 'facet'
   * in WSClean::_facets.
   */
  size_t facetIndex;

  /**
   * Pointer to a Facet. If it is null, faceting is not used.
   */
  std::shared_ptr<schaapcommon::facets::Facet> facet;

  /**
   * Difference from the centre to the facet centre to the facet centre.
   * Example, if the image centre was 100,100 and the facet centre 50,50, these
   * value will be -50,-50
   */
  int centreShiftX, centreShiftY;

  size_t outputChannelIndex;

  size_t outputIntervalIndex;

  /**
   * This vector links a filename index to MS data
   */
  std::vector<MSInfo> msData;

  /**
   * The group of entries with equal squaredDeconvolutionIndex should be
   * 'joinedly' deconvolved by adding their squared flux density values
   * together. Normally, all the polarizations from a single (output)channel /
   * timestep form such a group.
   */
  size_t squaredDeconvolutionIndex;

  /**
   * Entries with equal joinedGroupIndex are joinedly deconvolved.
   * Such a group of entries can be further split up in 'facet' and/or 'squared'
   * deconvolution groups.
   */
  size_t joinedGroupIndex;

  /**
   * A normal inversion results in '1' image. However, an XY
   * imaging run results in 2 (real and imaginary), while an
   * YX imaging run results in 0, as it is added to XY.
   */
  size_t imageCount;

  std::string tmpFilePrefix;

  double CentralFrequency() const {
    return 0.5 * (bandStartFrequency + bandEndFrequency);
  }

  /**
   * A number that scales with the estimated inverse-variance of the image. It
   * can be used when averaging images or fitting functions through the images
   * to get the optimal sensitivity. It is set after the first inversion.
   */
  double imageWeight;
  double normalizationFactor;

  void AssignGridData(const ImagingTableEntry& source) {
    imageWeight = source.imageWeight;
    normalizationFactor = source.normalizationFactor;
  }
};

#endif
