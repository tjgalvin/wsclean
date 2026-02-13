#include "msgridder.h"

#include "msprovidercollection.h"

#include "../math/calculatefftsize.h"

#include <aocommon/logger.h>
#include <aocommon/units/angle.h>

using aocommon::Logger;

namespace wsclean {

// Defined out of class to allow the class the be used in a std::unique_ptr.
MsGridder::~MsGridder() = default;

MsGridder::MsGridder(const Settings& settings,
                     MsProviderCollection& ms_provider_collection)
    : MsGridderData(settings),
      ms_data_vector_(ms_provider_collection.ms_data_vector_),
      w_grid_size_(settings.nWLayers),
      data_column_name_(settings.dataColumnName),
      small_inversion_(settings.minGridResolution),
      weighting_(settings.weightMode) {}

void MsGridder::CalculateOverallMetaData() {
  theoretical_beam_size_ = 1.0 / MaxBaseline();
  if (IsFirstTask()) {
    Logger::Info << "Theoretic beam = " +
                        aocommon::units::Angle::ToNiceString(
                            theoretical_beam_size_) +
                        "\n";
  }

  if (!HasTrimSize()) SetTrimSize(ImageWidth(), ImageHeight());

  actual_inversion_width_ = ImageWidth();
  actual_inversion_height_ = ImageHeight();
  actual_pixel_size_x_ = PixelSizeX();
  actual_pixel_size_y_ = PixelSizeY();

  if (SmallInversion()) {
    size_t optWidth, optHeight, minWidth, minHeight;
    CalculateFFTSize(actual_inversion_width_, actual_pixel_size_x_,
                     theoretical_beam_size_, minWidth, optWidth);
    CalculateFFTSize(actual_inversion_height_, actual_pixel_size_y_,
                     theoretical_beam_size_, minHeight, optHeight);
    if (optWidth < actual_inversion_width_ ||
        optHeight < actual_inversion_height_) {
      const size_t newWidth =
          std::max(std::min(optWidth, actual_inversion_width_), size_t(32));
      const size_t newHeight =
          std::max(std::min(optHeight, actual_inversion_height_), size_t(32));
      if (IsFirstTask()) {
        Logger::Info << "Minimal inversion size: " + std::to_string(minWidth) +
                            " x " + std::to_string(minHeight) +
                            ", using optimal: " + std::to_string(newWidth) +
                            " x " + std::to_string(newHeight) + "\n";
      }
      actual_pixel_size_x_ =
          (double(actual_inversion_width_) * actual_pixel_size_x_) /
          double(newWidth);
      actual_pixel_size_y_ =
          (double(actual_inversion_height_) * actual_pixel_size_y_) /
          double(newHeight);
      actual_inversion_width_ = newWidth;
      actual_inversion_height_ = newHeight;
    } else {
      if (IsFirstTask()) {
        Logger::Info
            << "Small inversion enabled, but inversion resolution already "
               "smaller than beam size: not using optimization.\n";
      }
    }
  }

  // Always call GetSuggestedWGridSize in the first iteration, since it then
  // logs the suggested wgrid size.
  const size_t suggestedGridSize =
      (IsFirstTask() || !hasWGridSize()) ? GetSuggestedWGridSize() : 0;
  actual_w_grid_size_ = hasWGridSize() ? w_grid_size_ : suggestedGridSize;

  // Allow derived gridders a chance to compute their own additional metadata.
  CalculateGridderMetaData();
}

}  // namespace wsclean
