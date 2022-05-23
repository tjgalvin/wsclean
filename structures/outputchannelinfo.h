#ifndef OUTPUT_CHANNEL_INFO_H
#define OUTPUT_CHANNEL_INFO_H

#include <cmath>
#include <cstring>
#include <vector>

struct OutputChannelInfo {
  OutputChannelInfo()
      : weight(0.0),
        normalizationFactor(1.0),
        wGridSize(0),
        visibilityCount(0),
        effectiveVisibilityCount(0.0),
        visibilityWeightSum(0.0),
        beamMaj(0.0),
        beamMin(0.0),
        beamPA(0.0),
        beamSizeEstimate(0.0),
        theoreticBeamSize(0.0),
        psfNormalizationFactor(1.0) {}
  double weight, normalizationFactor;
  std::size_t wGridSize, visibilityCount;
  double effectiveVisibilityCount, visibilityWeightSum;
  double beamMaj, beamMin, beamPA;
  // The beam size estimate is calculated from the longest baseline, and used
  // as initial value when fitting the (accurate) beam
  double beamSizeEstimate;
  double theoreticBeamSize, psfNormalizationFactor;
};

/**
 * Calculate the smallest "theoretic" beam size in a vector of channels.
 * Non-finite/NaN values are skipped. If no finite values are present,
 * NaN is returned.
 */
inline double SmallestTheoreticBeamSize(
    const std::vector<OutputChannelInfo>& channels) {
  double smallest_theoretic_beam_size = std::numeric_limits<double>::max();
  for (const OutputChannelInfo& channel : channels) {
    const double value = channel.theoreticBeamSize;
    if (std::isfinite(value)) {
      smallest_theoretic_beam_size =
          std::min(smallest_theoretic_beam_size, value);
    }
  }
  return (smallest_theoretic_beam_size == std::numeric_limits<double>::max())
             ? std::numeric_limits<double>::quiet_NaN()
             : smallest_theoretic_beam_size;
}

#endif
