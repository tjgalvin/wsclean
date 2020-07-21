#ifndef NOT_IMPLEMENTED_GRIDDER_H
#define NOT_IMPLEMENTED_GRIDDER_H

#include "msgridderbase.h"

#include <stdexcept>
#include <string>

class UnavailableGridder : public MSGridderBase {
 public:
  UnavailableGridder(const class WSCleanSettings&) { doThrow(); }

  virtual ~UnavailableGridder() final override { doThrow(); }

  virtual void Invert() final override { doThrow(); }

  virtual void Predict(Image) final override { doThrow(); }

  virtual void Predict(Image, Image) final override { doThrow(); }

  virtual Image ImageRealResult() final override {
    doThrow();
    return Image();
  }

  virtual Image ImageImaginaryResult() final override {
    doThrow();
    return Image();
  }

  static void SavePBCorrectedImages(class FitsWriter& /*writer*/,
                                    class ImageFilename& /*filename*/,
                                    const std::string& /*filenameKind*/,
                                    const WSCleanSettings&) {}

  static void SaveBeamImage(const class ImagingTableEntry& /*entry*/,
                            class ImageFilename& /*filename*/,
                            const WSCleanSettings&, double, double, double,
                            double, const MetaDataCache&) {}

 private:
  virtual size_t getSuggestedWGridSize() const final override {
    doThrow();
    return 0;
  }

  void doThrow() const {
    throw std::runtime_error(
        "This gridder is not available, because WSClean was not compiled to "
        "have this gridder. Use a different gridder or recompile WSClean and "
        "make sure the necessary prerequisites are satisfied.");
  }
};

#endif
