#ifndef WSCLEAN_IMAGING_TABLE_ENTRY_H
#define WSCLEAN_IMAGING_TABLE_ENTRY_H

#include <aocommon/polarization.h>

#include <string>
#include <vector>

namespace schaapcommon {
namespace facets {
class Facet;
}
}  // namespace schaapcommon

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
   * Index of the facet in WSClean::_facets.
   * ImagingTable uses this index for creating groups per facet.
   */
  size_t facetIndex;

  /**
   * Non-owning pointer to a Facet. If it is null, faceting is not used.
   */
  const schaapcommon::facets::Facet* facet;

  size_t outputChannelIndex;

  size_t outputIntervalIndex;

  /**
   * This vector links a filename index to MS data
   */
  std::vector<MSInfo> msData;

  /**
   * The group of entries with equal squaredDeconvolutionIndex
   * should be 'joinedly' deconvolved by adding their squared powers
   * together. Normally, all the polarizations from a single
   * (output)channel / timestep form such a group.
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
};

#endif
