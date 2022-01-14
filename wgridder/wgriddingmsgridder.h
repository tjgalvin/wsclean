#ifndef WGRIDDING_MS_GRIDDER_H
#define WGRIDDING_MS_GRIDDER_H

#include "../gridding/msgridderbase.h"

#include <memory>

class WGriddingMSGridder final : public MSGridderBase {
 public:
  WGriddingMSGridder(const Settings& settings);

  virtual void Invert() override;

  virtual void Predict(std::vector<Image>&& images) override;

  virtual std::vector<Image> ResultImages() override {
    return {std::move(_image)};
  }

  virtual void FreeImagingData() override {}

  virtual size_t getSuggestedWGridSize() const override { return 1; }

 private:
  Image _image;

  template <DDGainMatrix GainEntry>
  void gridMeasurementSet(MSData& msData);

  template <DDGainMatrix GainEntry>
  void predictMeasurementSet(MSData& msData);

  size_t calculateMaxNRowsInMemory(size_t channelCount) const;

  void getActualTrimmedSize(size_t& trimmedWidth, size_t& trimmedHeight) const;

  size_t _cpuCount;
  int64_t _memSize;
  double _accuracy;
  std::unique_ptr<class WGriddingGridder_Simple> _gridder;
};

#endif
