#ifndef RMS_IMAGE_H
#define RMS_IMAGE_H

#include "../structures/image.h"

class RMSImage {
 public:
  static void Make(Image& rmsOutput, const Image& inputImage, double windowSize,
                   long double beamMaj, long double beamMin, long double beamPA,
                   long double pixelScaleL, long double pixelScaleM);

  static void SlidingMinimum(Image& output, const Image& input,
                             size_t windowSize, size_t threadCount);

  static void SlidingMaximum(Image& output, const Image& input,
                             size_t windowSize, size_t threadCount) {
    Image flipped(input);
    flipped.Negate();
    SlidingMinimum(output, flipped, windowSize, threadCount);
    output.Negate();
  }

  static void MakeWithNegativityLimit(Image& rmsOutput, const Image& inputImage,
                                      double windowSize, long double beamMaj,
                                      long double beamMin, long double beamPA,
                                      long double pixelScaleL,
                                      long double pixelScaleM,
                                      size_t threadCount);
};

#endif
