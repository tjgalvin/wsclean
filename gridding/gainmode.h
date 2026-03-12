#ifndef GRIDDING_GAIN_MODE_H_
#define GRIDDING_GAIN_MODE_H_

#include <stdexcept>
#include <string>

#include <aocommon/polarization.h>

namespace wsclean {

/**
 * This enum summarizes the number of polarizations stored in the
 * measurement set provider for gridding, together with the type
 * of polarizations. It is mainly used for templating.
 */
enum class GainMode {
  /// Correct visibilities only with the X solution.
  kXX,
  /// Correct visibilities only with the Y solution.
  kYY,
  /// Correct Stokes I visibilities with the trace of the
  /// X and Y solutions.
  kTrace,
  /// Apply X and Y separately to the XX and YY visibilities.
  k2VisDiagonal,
  /// TODO to be added: Multiply the diagonal [X 0 ; 0 Y] matrix with the 2x2
  /// visibility matrix.
  // k4VisDiagonal,
  /// Correct visibilities with the full 2x2 complex matrix.
  kFull
};

constexpr std::size_t GetNVisibilities(GainMode mode) {
  switch (mode) {
    case GainMode::kXX:
    case GainMode::kYY:
    case GainMode::kTrace:
      return 1;
    case GainMode::k2VisDiagonal:
      return 2;
    case GainMode::kFull:
      return 4;
  }
  return 0;
}

constexpr bool AllowScalarCorrection(GainMode mode) { return mode == GainMode::kTrace; }

/**
 * @param polarization The polarization requested
 */
inline GainMode SelectGainMode(aocommon::PolarizationEnum polarization,
                               size_t n_visibility_polarizations) {
  switch (n_visibility_polarizations) {
    case 1:
      switch (polarization) {
        case aocommon::Polarization::XX:
          return GainMode::kXX;
        case aocommon::Polarization::YY:
          return GainMode::kYY;
        case aocommon::Polarization::StokesI:
          // polarization might also be RR or Stokes other than I. We still need
          // to provide a GainMode, so we also return trace in those cases.
        default:
          return GainMode::kTrace;
      }
      break;
    case 2:
      // When 2 polarizations are stored (XX, YY), it might be because Stokes I
      // is being imaged with diagonal solutions, but it might also be that both
      // polarizations are requested independently and imaged at the same time.
      // IDG could in theory do this, although it's currently not implemented.
      if (polarization == aocommon::Polarization::StokesI ||
          polarization == aocommon::Polarization::DiagonalInstrumental ||
          polarization == aocommon::Polarization::XX || polarization == aocommon::Polarization::YY)
        return GainMode::k2VisDiagonal;
      break;
    case 4:
      if (polarization == aocommon::Polarization::FullStokes ||
          polarization == aocommon::Polarization::Instrumental ||
          aocommon::Polarization::IsStokes(polarization) ||
          polarization == aocommon::Polarization::XX || polarization == aocommon::Polarization::YY)
        return GainMode::kFull;
      break;
  }
  throw std::runtime_error("Invalid combination of polarization (" +
                           aocommon::Polarization::TypeToFullString(polarization) +
                           ") and n_visibility_polarizations (" +
                           std::to_string(n_visibility_polarizations) + ") in GetGainMode()");
}

}  // namespace wsclean

#endif
