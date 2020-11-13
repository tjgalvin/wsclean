#ifndef RMS_IMAGE_H
#define RMS_IMAGE_H

#include "image.h"

class RMSImage {
 public:
  static void Make(ImageF& rmsOutput, const ImageF& inputImage,
                   double windowSize, long double beamMaj, long double beamMin,
                   long double beamPA, long double pixelScaleL,
                   long double pixelScaleM);

  static void SlidingMinimum(ImageF& output, const ImageF& input,
                             size_t windowSize);

  static void SlidingMaximum(ImageF& output, const ImageF& input,
                             size_t windowSize) {
    ImageF flipped(input);
    flipped.Negate();
    SlidingMinimum(output, flipped, windowSize);
    output.Negate();
  }

  static void MakeWithNegativityLimit(ImageF& rmsOutput,
                                      const ImageF& inputImage,
                                      double windowSize, long double beamMaj,
                                      long double beamMin, long double beamPA,
                                      long double pixelScaleL,
                                      long double pixelScaleM);
};

#endif
