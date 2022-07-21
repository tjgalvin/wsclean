#ifndef WGRIDDING_MS_GRIDDER_H
#define WGRIDDING_MS_GRIDDER_H

#include "../gridding/msgridderbase.h"
#include "../structures/resources.h"

#include <aocommon/image.h>

#include <memory>

class WGriddingGridder_Simple;

class WGriddingMSGridder final : public MSGridderBase {
 public:
  WGriddingMSGridder(const Settings& settings, const Resources& resources);
  ~WGriddingMSGridder();

  virtual void Invert() override;

  virtual void Predict(std::vector<aocommon::Image>&& images) override;

  virtual std::vector<aocommon::Image> ResultImages() override {
    return {std::move(_image)};
  }

  virtual void FreeImagingData() override {}

  virtual size_t getSuggestedWGridSize() const override { return 1; }

 private:
  aocommon::Image _image;

  template <DDGainMatrix GainEntry>
  void gridMeasurementSet(MSData& msData);

  template <DDGainMatrix GainEntry>
  void predictMeasurementSet(MSData& msData);

  size_t calculateMaxNRowsInMemory(size_t channelCount) const;

  void getActualTrimmedSize(size_t& trimmedWidth, size_t& trimmedHeight) const;

  const Resources _resources;
  double _accuracy;
  std::unique_ptr<WGriddingGridder_Simple> _gridder;
};

#endif
