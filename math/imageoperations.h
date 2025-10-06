#ifndef IMAGE_OPERATIONS_H
#define IMAGE_OPERATIONS_H

#include "../structures/outputchannelinfo.h"
#include "../io/imagefilename.h"

#include <aocommon/hmatrix4x4.h>
#include <aocommon/image.h>
#include <aocommon/polarization.h>

#include <vector>
#include <string>

namespace wsclean {

class Settings;

namespace math {

/**
 * Multiplies every pixel of a set of Stokes I,Q,U,V images with a Mueller
 * matrix. This function assumes that the Mueller matrix is a correction to
 * linear instrumental polarizations, i.e. it applies a correction on
 * XX/XY/YX/YY values. The Stokes IQUV image values are therefore first
 * converted to linear polarization, then the matrix is applied, and then the
 * values are converted and written back to Stokes images.
 */
void CorrectImagesForMuellerMatrix(const aocommon::HMC4x4& mueller_correction,
                                   std::array<aocommon::Image*, 4>& images);

/**
 * Like @ref CorrectImagesForMuellerMatrix(), but for XX/YY correction. It
 * is assumed that the XY/YX values are zero. This function does not do a
 * Stokes correction (unlike @ref CorrectImagesForMuellerMatrix()); instead
 * the input images are already assumed to have the same polarization system
 * as the Mueller matrix (i.e., linear for a linearly polarized telescope).
 */
void CorrectDualImagesForMuellerMatrix(
    const aocommon::HMC4x4& mueller_correction,
    std::array<aocommon::Image*, 2>& images);

void DetermineBeamSize(const Settings& settings, double& bMaj, double& bMin,
                       double& bPA, double& bTheoretical,
                       const aocommon::Image& image, double initialEstimate);

void MakeMFSImage(const Settings& settings,
                  const std::vector<OutputChannelInfo>& infoPerChannel,
                  OutputChannelInfo& mfsInfo, const std::string& suffix,
                  size_t intervalIndex, aocommon::PolarizationEnum pol,
                  ImageFilenameType image_type,
                  std::optional<size_t> directionIndex = std::nullopt);

void RenderMFSImage(const Settings& settings, const OutputChannelInfo& mfsInfo,
                    size_t intervalIndex, aocommon::PolarizationEnum pol,
                    bool isImaginary, bool isPBCorrected);

}  // namespace math
}  // namespace wsclean

#endif
