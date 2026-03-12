#ifndef OUTPUT_CHANNEL_INFO_H
#define OUTPUT_CHANNEL_INFO_H

#include <cmath>
#include <cstring>
#include <vector>

#include "../gridding/averagecorrection.h"

namespace wsclean {

struct OutputChannelInfo {
  OutputChannelInfo(size_t n_facets = 0, size_t n_dd_psfs = 0)
      : averageFacetCorrection(n_facets),
        averageBeamFacetCorrection(n_facets),
        averageDdPsfCorrection(n_dd_psfs) {}
  double weight = 0.0;
  double normalizationFactor = 1.0;
  std::size_t wGridSize = 0;
  std::size_t visibilityCount = 0;
  double effectiveVisibilityCount = 0.0;
  double visibilityWeightSum = 0.0;
  double beamMaj = 0.0;
  double beamMin = 0.0;
  double beamPA = 0.0;
  // The beam size estimate is calculated from the longest baseline, and used
  // as initial value when fitting the (accurate) beam
  double beamSizeEstimate = 0.0;
  double theoreticBeamSize = 0.0;
  double psfNormalizationFactor = 1.0;
  // For dd psf mode, this is the facet index that holds the central
  // psf.
  std::size_t centralPsfIndex = 0;
  // See VisibilityModifier for an explanation
  std::vector<AverageCorrection> averageFacetCorrection;
  std::vector<AverageCorrection> averageBeamFacetCorrection;
  // Same as averageFacetCorrection, but then for the DD PSF facets
  std::vector<AverageCorrection> averageDdPsfCorrection;
};

/**
 * Calculate the smallest "theoretic" beam size in a vector of channels.
 * Non-finite/NaN values are skipped. If no finite values are present,
 * NaN is returned.
 */
inline double SmallestTheoreticBeamSize(const std::vector<OutputChannelInfo>& channels) {
  double smallest_theoretic_beam_size = std::numeric_limits<double>::max();
  for (const OutputChannelInfo& channel : channels) {
    const double value = channel.theoreticBeamSize;
    if (std::isfinite(value)) {
      smallest_theoretic_beam_size = std::min(smallest_theoretic_beam_size, value);
    }
  }
  return (smallest_theoretic_beam_size == std::numeric_limits<double>::max())
             ? std::numeric_limits<double>::quiet_NaN()
             : smallest_theoretic_beam_size;
}

}  // namespace wsclean

#endif
