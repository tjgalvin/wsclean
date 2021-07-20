#include "rmsimage.h"

#include "modelrenderer.h"

#include <aocommon/staticfor.h>

void RMSImage::Make(Image& rmsOutput, const Image& inputImage,
                    double windowSize, long double beamMaj, long double beamMin,
                    long double beamPA, long double pixelScaleL,
                    long double pixelScaleM) {
  Image image(inputImage);
  rmsOutput = Image(image.Width(), image.Height(), 0.0);

  for (auto& val : image) val *= val;

  ModelRenderer::Restore(rmsOutput.data(), image.data(), image.Width(),
                         image.Height(), beamMaj * windowSize,
                         beamMin * windowSize, beamPA, pixelScaleL,
                         pixelScaleM);

  double s = std::sqrt(2.0 * M_PI);
  const long double sigmaMaj = beamMaj / (2.0L * sqrtl(2.0L * logl(2.0L)));
  const long double sigmaMin = beamMin / (2.0L * sqrtl(2.0L * logl(2.0L)));
  const double norm = 1.0 / (s * sigmaMaj / pixelScaleL * windowSize * s *
                             sigmaMin / pixelScaleL * windowSize);
  for (auto& val : rmsOutput) val = std::sqrt(val * norm);
}

void RMSImage::SlidingMinimum(Image& output, const Image& input,
                              size_t windowSize, size_t threadCount) {
  const size_t width = input.Width();
  output = Image(width, input.Height());
  Image temp(output);

  aocommon::StaticFor<size_t> loop(threadCount);

  loop.Run(0, input.Height(), [&](size_t yStart, size_t yEnd) {
    for (size_t y = yStart; y != yEnd; ++y) {
      float* outRowptr = &temp[y * width];
      const float* inRowptr = &input[y * width];
      for (size_t x = 0; x != width; ++x) {
        size_t left = std::max(x, windowSize / 2) - windowSize / 2;
        size_t right = std::min(x, width - windowSize / 2) + windowSize / 2;
        outRowptr[x] = *std::min_element(inRowptr + left, inRowptr + right);
      }
    }
  });

  loop.Run(0, width, [&](size_t xStart, size_t xEnd) {
    aocommon::UVector<float> vals;
    for (size_t x = xStart; x != xEnd; ++x) {
      for (size_t y = 0; y != input.Height(); ++y) {
        size_t top = std::max(y, windowSize / 2) - windowSize / 2;
        size_t bottom =
            std::min(y, input.Height() - windowSize / 2) + windowSize / 2;
        vals.clear();
        for (size_t winY = top; winY != bottom; ++winY)
          vals.push_back(temp[winY * width + x]);
        output[y * width + x] = *std::min_element(vals.begin(), vals.end());
      }
    }
  });
}

void RMSImage::MakeWithNegativityLimit(
    Image& rmsOutput, const Image& inputImage, double windowSize,
    long double beamMaj, long double beamMin, long double beamPA,
    long double pixelScaleL, long double pixelScaleM, size_t threadCount) {
  Make(rmsOutput, inputImage, windowSize, beamMaj, beamMin, beamPA, pixelScaleL,
       pixelScaleM);
  Image slidingMinimum(inputImage.Width(), inputImage.Height());
  double beamInPixels = std::max(beamMaj / pixelScaleL, 1.0L);
  SlidingMinimum(slidingMinimum, inputImage, windowSize * beamInPixels,
                 threadCount);
  for (size_t i = 0; i != rmsOutput.size(); ++i) {
    rmsOutput[i] = std::max<float>(rmsOutput[i],
                                   std::abs(slidingMinimum[i]) * (1.5 / 5.0));
  }
}
