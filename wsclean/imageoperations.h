#ifndef IMAGE_OPERATIONS_H
#define IMAGE_OPERATIONS_H

#include "outputchannelinfo.h"

#include <aocommon/polarization.h>

#include <vector>
#include <string>

class ImageOperations {
 public:
  static void FitBeamSize(const class WSCleanSettings& settings, double& bMaj,
                          double& bMin, double& bPA, const double* image,
                          double beamEstimate);

  static void DetermineBeamSize(const class WSCleanSettings& settings,
                                double& bMaj, double& bMin, double& bPA,
                                double& bTheoretical, const double* image,
                                double initialEstimate);

  static void MakeMFSImage(const class WSCleanSettings& settings,
                           const std::vector<OutputChannelInfo>& infoPerChannel,
                           OutputChannelInfo& mfsInfo,
                           const std::string& suffix, size_t intervalIndex,
                           aocommon::PolarizationEnum pol, bool isImaginary,
                           bool isPSF = false);

  static void RenderMFSImage(const class WSCleanSettings& settings,
                             const OutputChannelInfo& mfsInfo,
                             size_t intervalIndex,
                             aocommon::PolarizationEnum pol, bool isImaginary,
                             bool isPBCorrected);
};

#endif
