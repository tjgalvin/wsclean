#ifndef WS_MS_GRIDDER_H
#define WS_MS_GRIDDER_H

#include "msgridder.h"
#include "wstackinggridder.h"

#include "../structures/resources.h"

#include <casacore/casa/Arrays/Array.h>
#include <casacore/tables/Tables/ArrayColumn.h>

#include <aocommon/image.h>
#include <aocommon/lane.h>
#include <aocommon/multibanddata.h>

#include <complex>
#include <memory>
#include <thread>

namespace wsclean {

class WSMSGridder final : public MsGridder {
 public:
  typedef WStackingGridder<float> GridderType;

  WSMSGridder(const Settings& settings, const Resources& resources,
              MsProviderCollection& ms_provider_collection);
  ~WSMSGridder() noexcept final;

  size_t GetNInversionPasses() const final { return _gridder->NPasses(); }
  void StartInversion() final;
  void StartInversionPass(size_t pass_index) final;
  size_t GridMeasurementSet(const MsProviderCollection::MsData& ms_data) final;
  void FinishInversionPass(size_t pass_index) final;
  void FinishInversion() final;

  size_t GetNPredictPasses() const final { return _gridder->NPasses(); }
  void StartPredict(std::vector<aocommon::Image>&& images) final;
  void StartPredictPass(size_t pass_index) final;
  size_t PredictMeasurementSet(
      const MsProviderCollection::MsData& ms_data) final;
  void FinishPredictPass(size_t pass_index) final;
  void FinishPredict() final;

  std::vector<aocommon::Image> ResultImages() final {
    if (IsComplex())
      return {std::move(_realImage), std::move(_imaginaryImage)};
    else
      return {std::move(_realImage)};
  }
  void FreeImagingData() final { _gridder.reset(); }

  size_t AntialiasingKernelSize() const { return _antialiasingKernelSize; }
  size_t OverSamplingFactor() const { return _overSamplingFactor; }

  bool HasNWSize() const { return NWWidth() != 0 || NWHeight() != 0; }
  size_t NWWidth() const { return GetSettings().widthForNWCalculation; }
  size_t NWHeight() const { return GetSettings().heightForNWCalculation; }
  double NWFactor() const { return GetSettings().nWLayersFactor; }

 private:
  struct InversionWorkSample {
    double uInLambda, vInLambda, wInLambda;
    std::complex<float> sample;
  };
  struct PredictionWorkItem {
    std::array<double, 3> uvw;
    std::unique_ptr<std::complex<float>[]> data;
    size_t data_desc_id;
    size_t rowId;
  };

  void countSamplesPerLayer(MsProviderCollection::MsData& msData);
  size_t GetSuggestedWGridSize() const final;

  void startInversionWorkThreads(size_t maxChannelCount);
  void finishInversionWorkThreads();
  void workThreadPerSample(aocommon::Lane<InversionWorkSample>* workLane);

  void predictCalcThread(aocommon::Lane<PredictionWorkItem>* input_lane,
                         aocommon::Lane<PredictionWorkItem>* output_lane,
                         const aocommon::MultiBandData* bands);

  void predictWriteThread(
      aocommon::Lane<PredictionWorkItem>* sampling_work_lane,
      const MsProviderCollection::MsData* ms_data, GainMode gain_mode);

  std::unique_ptr<GridderType> _gridder;
  std::vector<aocommon::Lane<InversionWorkSample>> _inversionCPULanes;
  std::vector<std::thread> _threadGroup;
  size_t _antialiasingKernelSize, _overSamplingFactor;
  const Resources _resources;
  size_t _laneBufferSize;
  aocommon::Image _realImage;
  aocommon::Image _imaginaryImage;
};

}  // namespace wsclean

#endif
