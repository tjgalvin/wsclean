#ifndef NOT_IMPLEMENTED_GRIDDER_H
#define NOT_IMPLEMENTED_GRIDDER_H

#include "msgridderbase.h"

#include <stdexcept>
#include <string>

class UnavailableGridder final : public MSGridderBase {
 public:
  UnavailableGridder(const class Settings&) { doThrow(); }

  virtual ~UnavailableGridder() override { doThrow(); }

  virtual void Invert() override { doThrow(); }

  virtual void Predict(ImageF) override { doThrow(); }

  virtual void Predict(ImageF, ImageF) override { doThrow(); }

  virtual ImageF ImageRealResult() override {
    doThrow();
    return ImageF();
  }

  virtual ImageF ImageImaginaryResult() override {
    doThrow();
    return ImageF();
  }

  static void SavePBCorrectedImages(class FitsWriter& /*writer*/,
                                    class ImageFilename& /*filename*/,
                                    const std::string& /*filenameKind*/,
                                    const Settings&) {}

  static void SaveBeamImage(const struct ImagingTableEntry& /*entry*/,
                            class ImageFilename& /*filename*/, const Settings&,
                            double, double, double, double,
                            const MetaDataCache&) {}

 private:
  virtual size_t getSuggestedWGridSize() const override {
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
