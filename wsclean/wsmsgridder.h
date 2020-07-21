#ifndef WS_MS_GRIDDER_H
#define WS_MS_GRIDDER_H

#include "msgridderbase.h"
#include "wstackinggridder.h"

#include "../multibanddata.h"
#include "../image.h"

#include <casacore/casa/Arrays/Array.h>
#include <casacore/tables/Tables/ArrayColumn.h>

#include <aocommon/lane.h>

#include <complex>
#include <memory>
#include <thread>

namespace casacore {
class MeasurementSet;
}
class ImageBufferAllocator;

class WSMSGridder : public MSGridderBase {
 public:
  typedef WStackingGridderF GridderType;

  WSMSGridder(size_t threadCount, double memFraction, double absMemLimit);

  virtual void Invert() final override;

  virtual void Predict(Image image) final override {
    Predict(std::move(image), Image());
  }
  virtual void Predict(Image real, Image imaginary) final override;

  virtual Image ImageRealResult() final override {
    return std::move(_realImage);
  }
  virtual Image ImageImaginaryResult() final override {
    if (!IsComplex())
      throw std::runtime_error(
          "No imaginary result available for non-complex inversion");
    return std::move(_imaginaryImage);
  }
  virtual size_t ActualInversionWidth() const final override {
    return _actualInversionWidth;
  }
  virtual size_t ActualInversionHeight() const final override {
    return _actualInversionHeight;
  }

  virtual void FreeImagingData() final override { _gridder.reset(); }

 private:
  struct InversionWorkSample {
    double uInLambda, vInLambda, wInLambda;
    std::complex<float> sample;
  };
  struct PredictionWorkItem {
    double u, v, w;
    std::unique_ptr<std::complex<float>[]> data;
    size_t rowId, dataDescId;
  };

  void gridMeasurementSet(MSData& msData);
  void countSamplesPerLayer(MSData& msData);
  virtual size_t getSuggestedWGridSize() const final override;

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
                         aocommon::Lane<PredictionWorkItem>* outputLane);
  void predictWriteThread(aocommon::Lane<PredictionWorkItem>* samplingWorkLane,
                          const MSData* msData);

  std::unique_ptr<GridderType> _gridder;
  std::vector<aocommon::Lane<InversionWorkSample>> _inversionCPULanes;
  std::vector<std::thread> _threadGroup;
  size_t _cpuCount, _laneBufferSize;
  int64_t _memSize;
  Image _realImage, _imaginaryImage;
};

#endif
