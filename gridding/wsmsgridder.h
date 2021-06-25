#ifndef WS_MS_GRIDDER_H
#define WS_MS_GRIDDER_H

#include "msgridderbase.h"
#include "wstackinggridder.h"

#include "../structures/multibanddata.h"
#include "../structures/image.h"

#include <casacore/casa/Arrays/Array.h>
#include <casacore/tables/Tables/ArrayColumn.h>

#include <aocommon/lane.h>

#include <complex>
#include <memory>
#include <thread>

class WSMSGridder final : public MSGridderBase {
 public:
  typedef WStackingGridder<float> GridderType;

  WSMSGridder(const Settings& settings);
  ~WSMSGridder() noexcept;

  virtual void Invert() override;

  virtual void Predict(std::vector<Image>&& images) override;

  virtual std::vector<Image> ResultImages() override {
    if (IsComplex())
      return {std::move(_realImage), std::move(_imaginaryImage)};
    else
      return {std::move(_realImage)};
  }
  virtual size_t ActualInversionWidth() const override {
    return _actualInversionWidth;
  }
  virtual size_t ActualInversionHeight() const override {
    return _actualInversionHeight;
  }

  virtual void FreeImagingData() override { _gridder.reset(); }

  size_t AntialiasingKernelSize() const { return _antialiasingKernelSize; }
  size_t OverSamplingFactor() const { return _overSamplingFactor; }

  bool HasNWSize() const { return _nwWidth != 0 || _nwHeight != 0; }
  size_t NWWidth() const { return _nwWidth; }
  size_t NWHeight() const { return _nwHeight; }
  double NWFactor() const { return _nwFactor; }

 private:
  struct InversionWorkSample {
    double uInLambda, vInLambda, wInLambda;
    std::complex<float> sample;
  };
  struct PredictionWorkItem {
    std::array<double, 3> uvw;
    std::unique_ptr<std::complex<float>[]> data;
    size_t rowId, dataDescId;
  };

  void gridMeasurementSet(MSData& msData);
  void countSamplesPerLayer(MSData& msData);
  virtual size_t getSuggestedWGridSize() const override;

  void predictMeasurementSet(MSData& msData);

  void workThread(aocommon::Lane<InversionRow>* workLane) {
    InversionRow workItem;
    while (workLane->read(workItem)) {
      _gridder->AddData(workItem.data, workItem.dataDescId, workItem.uvw[0],
                        workItem.uvw[1], workItem.uvw[2]);
      delete[] workItem.data;
    }
  }

  void startInversionWorkThreads(size_t maxChannelCount);
  void finishInversionWorkThreads();
  void workThreadPerSample(aocommon::Lane<InversionWorkSample>* workLane);

  void predictCalcThread(aocommon::Lane<PredictionWorkItem>* inputLane,
                         aocommon::Lane<PredictionWorkItem>* outputLane,
                         const MultiBandData* bandData);
  void predictWriteThread(aocommon::Lane<PredictionWorkItem>* samplingWorkLane,
                          const MSData* msData, const MultiBandData* bandData);

  std::unique_ptr<GridderType> _gridder;
  std::vector<aocommon::Lane<InversionWorkSample>> _inversionCPULanes;
  std::vector<std::thread> _threadGroup;
  size_t _nwWidth, _nwHeight;
  double _nwFactor;
  size_t _antialiasingKernelSize, _overSamplingFactor;
  size_t _cpuCount, _laneBufferSize;
  int64_t _memSize;
  Image _realImage, _imaginaryImage;
};

#endif
