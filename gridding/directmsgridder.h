#ifndef DIRECT_MS_GRIDDER_H
#define DIRECT_MS_GRIDDER_H

#include <aocommon/image.h>
#include <aocommon/lane.h>

#include "msgridder.h"

#include "../main/progressbar.h"
#include "../structures/resources.h"

namespace wsclean {

class ProgressBar;

template <typename num_t>
class DirectMSGridder final : public MsGridder {
 public:
  DirectMSGridder(const Settings& settings, const Resources& resources,
                  MsProviderCollection& ms_provider_collection);

  void StartInversion() final;
  size_t GridMeasurementSet(const MsProviderCollection::MsData& ms_data) final;
  void FinishInversion() final;

  void StartPredict(std::vector<aocommon::Image>&& images) final;
  size_t PredictMeasurementSet(const MsProviderCollection::MsData& ms_data) final;
  void FinishPredict() final;

  std::vector<aocommon::Image> ResultImages() final { return {std::move(_image)}; }
  size_t GetSuggestedWGridSize() const final { return 1; }

 private:
  struct InversionSample {
    num_t uInLambda, vInLambda, wInLambda;
    std::complex<float> sample;
  };
  const Resources _resources;
  aocommon::Image _image;
  aocommon::ImageBase<num_t> _sqrtLMTable;
  std::vector<aocommon::ImageBase<num_t>> _layers;
  aocommon::Lane<InversionSample> _inversionLane;

  void InvertMeasurementSet(const MsProviderCollection::MsData& ms_data, size_t ms_index);
  void gridSample(const InversionSample& sample, size_t layer);
  aocommon::ImageBase<num_t> GetSqrtLMLookupTable() const;
  std::unique_ptr<ProgressBar> progress_bar_;
};

}  // namespace wsclean

#endif
