#ifndef NOT_IMPLEMENTED_GRIDDER_H
#define NOT_IMPLEMENTED_GRIDDER_H

#include "msgridder.h"
#include "../structures/resources.h"

#include <aocommon/image.h>
#include <aocommon/fits/fitswriter.h>

#include <stdexcept>
#include <string>

namespace wsclean {

class UnavailableGridder final : public MsGridder {
 public:
  UnavailableGridder(const class Settings& settings, const Resources& resources,
                     MsProviderCollection& ms_provider_collection)
      : MsGridder(settings, ms_provider_collection) {
    doThrow();
  }

  ~UnavailableGridder() = default;

  void StartInversion() final { doThrow(); }
  size_t GridMeasurementSet(const MsProviderCollection::MsData& ms_data) final {
    doThrow();
    return 0;
  }
  void FinishInversion() final { doThrow(); }

  void StartPredict(std::vector<aocommon::Image>&& images) final { doThrow(); }
  size_t PredictMeasurementSet(
      const MsProviderCollection::MsData& /*ms_data*/) final {
    doThrow();
    return 0;
  }
  void FinishPredict() final { doThrow(); }

  std::vector<aocommon::Image> ResultImages() final {
    doThrow();
    return {};
  }

  static void SavePBCorrectedImages(aocommon::FitsWriter& /*writer*/,
                                    class ImageFilename& /*filename*/,
                                    const std::string& /*filenameKind*/,
                                    const Settings&) {}

  static void SaveBeamImage(const struct ImagingTableEntry& /*entry*/,
                            class ImageFilename& /*filename*/, const Settings&,
                            double, double, double, double,
                            const AverageBeam&) {}

  void SetAverageBeam(std::unique_ptr<AverageBeam>) { doThrow(); }

  std::unique_ptr<AverageBeam> ReleaseAverageBeam() {
    doThrow();
    return {};
  }

 private:
  size_t GetSuggestedWGridSize() const final {
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

}  // namespace wsclean

#endif
