#ifndef RMS_IMAGE_H
#define RMS_IMAGE_H

#include <aocommon/image.h>

class RMSImage {
 public:
  static void Make(aocommon::Image& rmsOutput,
                   const aocommon::Image& inputImage, double windowSize,
                   long double beamMaj, long double beamMin, long double beamPA,
                   long double pixelScaleL, long double pixelScaleM,
                   size_t threadCount);

  static void SlidingMinimum(aocommon::Image& output,
                             const aocommon::Image& input, size_t windowSize,
                             size_t threadCount);

  static void SlidingMaximum(aocommon::Image& output,
                             const aocommon::Image& input, size_t windowSize,
                             size_t threadCount);

  static void MakeWithNegativityLimit(aocommon::Image& rmsOutput,
                                      const aocommon::Image& inputImage,
                                      double windowSize, long double beamMaj,
                                      long double beamMin, long double beamPA,
                                      long double pixelScaleL,
                                      long double pixelScaleM,
                                      size_t threadCount);
};

#endif
