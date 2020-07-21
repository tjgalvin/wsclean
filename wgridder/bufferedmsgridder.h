#ifndef BUFFERED_MS_GRIDDER_H
#define BUFFERED_MS_GRIDDER_H

#include "../wsclean/msgridderbase.h"

#include "../multibanddata.h"

#include <aocommon/lane.h>

#include <casacore/casa/Arrays/Array.h>
#include <casacore/tables/Tables/ArrayColumn.h>

#include <complex>
#include <memory>
#include <thread>

namespace casacore {
class MeasurementSet;
}
class ImageBufferAllocator;

class BufferedMSGridder : public MSGridderBase {
 public:
  BufferedMSGridder(size_t threadCount, double memFraction, double absMemLimit);

  virtual void Invert() final override;

  virtual void Predict(Image image) final override;
  virtual void Predict(Image, Image) final override {
    throw std::runtime_error("Can not do imaginary imaging in this mode");
  }

  virtual Image ImageRealResult() final override { return std::move(_image); }
  virtual Image ImageImaginaryResult() final override {
    throw std::runtime_error("Can not do imaginary imaging in this mode");
  }

  virtual size_t ActualInversionWidth() const final override {
    return _actualInversionWidth;
  }
  virtual size_t ActualInversionHeight() const final override {
    return _actualInversionHeight;
  }

  virtual void FreeImagingData() final override {}

  virtual size_t getSuggestedWGridSize() const final override { return 1; }

 private:
  Image _image;

  void gridMeasurementSet(MSData& msData);

  void predictMeasurementSet(MSData& msData);

  size_t calculateMaxNRowsInMemory(size_t channelCount) const;

  void getTrimmedSize(size_t& trimmedWidth, size_t& trimmedHeight) const;

  size_t _cpuCount;
  int64_t _memSize;
  std::unique_ptr<class WGriddingGridder_Simple> _gridder;
};

#endif
