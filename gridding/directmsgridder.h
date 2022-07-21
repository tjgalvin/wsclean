#ifndef DIRECT_MS_GRIDDER_H
#define DIRECT_MS_GRIDDER_H

#include <aocommon/image.h>
#include <aocommon/lane.h>

#include "msgridderbase.h"

#include "../structures/resources.h"

template <typename num_t>
class DirectMSGridder final : public MSGridderBase {
 public:
  const static size_t num_t_factor =
      (sizeof(num_t) + sizeof(double) - 1) / sizeof(double);

  DirectMSGridder(const Settings& settings, const Resources& resources);

  virtual void Invert() override;

  virtual void Predict(std::vector<aocommon::Image>&& images) override;

  virtual std::vector<aocommon::Image> ResultImages() override {
    return {std::move(_image)};
  }
  virtual size_t getSuggestedWGridSize() const override { return 1; }

 private:
  struct InversionSample {
    num_t uInLambda, vInLambda, wInLambda;
    std::complex<float> sample;
  };
  const Resources _resources;
  aocommon::Image _image;
  num_t* _sqrtLMTable;
  std::vector<num_t*> _layers;
  aocommon::Lane<InversionSample> _inversionLane;

  template <DDGainMatrix GainEntry>
  void invertMeasurementSet(const MSData& msData, class ProgressBar& progress,
                            size_t msIndex);
  void inversionWorker(size_t layer);
  void gridSample(const InversionSample& sample, size_t layer);
  void initializeSqrtLMLookupTable();

  num_t* allocate() { return new num_t[ImageWidth() * ImageHeight()]; }
  void freeImg(num_t* ptr) { delete[] ptr; }
};

#endif
