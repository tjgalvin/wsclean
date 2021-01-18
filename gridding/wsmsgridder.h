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
  typedef WStackingGridderF GridderType;

  WSMSGridder(size_t threadCount, double memFraction, double absMemLimit);

  virtual void Invert() override;

  virtual void Predict(ImageF image) override {
    Predict(std::move(image), ImageF());
  }
  virtual void Predict(ImageF real, ImageF imaginary) override;

  virtual ImageF ImageRealResult() override { return std::move(_realImage); }
  virtual ImageF ImageImaginaryResult() override {
    if (!IsComplex())
      throw std::runtime_error(
          "No imaginary result available for non-complex inversion");
    return std::move(_imaginaryImage);
  }
  virtual size_t ActualInversionWidth() const override {
    return _actualInversionWidth;
  }
  virtual size_t ActualInversionHeight() const override {
    return _actualInversionHeight;
  }

  virtual void FreeImagingData() override { _gridder.reset(); }

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
                          const MSData* msData);

  std::unique_ptr<GridderType> _gridder;
  std::vector<aocommon::Lane<InversionWorkSample>> _inversionCPULanes;
  std::vector<std::thread> _threadGroup;
  size_t _cpuCount, _laneBufferSize;
  int64_t _memSize;
  ImageF _realImage, _imaginaryImage;
};

#endif
