#ifndef BUFFERED_MS_GRIDDER_H
#define BUFFERED_MS_GRIDDER_H

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

  virtual size_t ActualInversionWidth() const override {
    return _actualInversionWidth;
  }
  virtual size_t ActualInversionHeight() const override {
    return _actualInversionHeight;
  }

  virtual void FreeImagingData() override {}

  virtual size_t getSuggestedWGridSize() const override { return 1; }

 private:
  Image _image;

  void gridMeasurementSet(MSData& msData);

  void predictMeasurementSet(MSData& msData);

  size_t calculateMaxNRowsInMemory(size_t channelCount) const;

  void getTrimmedSize(size_t& trimmedWidth, size_t& trimmedHeight) const;

  size_t _cpuCount;
  int64_t _memSize;
  double _accuracy;
  std::unique_ptr<class WGriddingGridder_Simple> _gridder;
};

#endif
