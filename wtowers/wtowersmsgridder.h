#ifndef WSCLEAN_WTOWERS_MS_GRIDDER_H
#define WSCLEAN_WTOWERS_MS_GRIDDER_H

#include <aocommon/image.h>

#include "../gridding/msgridder.h"
#include "../structures/resources.h"

// Main known outstanding issues:
// TODO: Currently we apply padding twice; once externally from the generic
// WSClean code and then a second time in the W-towers code
// TODO: Currently we assume square images, further changes required to support
// faceting
// TODO: Currently W-towers uses all threads for a single instance, further
// changes required to support parallel imaging
// TODO: Currently W-towers cannot do on the fly corrections, so doesn't support
// -shared-facet-reads
// TODO: Currently we reuse wgridder memory estimates, we need to implement
// W-towers specific estimates
// TODO: Other general cleanup, profiling, testing etc.
// TODO: We need at least one unit test for this gridder

namespace wsclean {

class WTowersGridderBase;

class WTowersMsGridder final : public MsGridder {
 public:
  WTowersMsGridder(const Settings& settings, const Resources& resources,
                   MsProviderCollection& ms_provider_collection);
  ~WTowersMsGridder() final;

  void CalculateGridderMetaData() final;

  void StartInversion() final;
  size_t GridMeasurementSet(const MsProviderCollection::MsData& ms_data) final;
  void FinishInversion() final;

  void StartPredict(std::vector<aocommon::Image>&& images) final;
  size_t PredictMeasurementSet(
      const MsProviderCollection::MsData& ms_data) final;
  void FinishPredict() final;

  std::vector<aocommon::Image> ResultImages() final {
    return {std::move(image_)};
  }

  void FreeImagingData() final {}

  size_t GetSuggestedWGridSize() const final { return 1; }

 private:
  std::unique_ptr<WTowersGridderBase> MakeGridder(size_t width,
                                                  size_t height) const;

  size_t CalculateConstantMemory() const final;
  size_t CalculateMaxRowsInMemory(int64_t available_memory,
                                  size_t constant_memory,
                                  double additional_per_row_consumption,
                                  size_t per_row_uvw_consumption,
                                  size_t channel_count,
                                  size_t num_polarizations_stored) const final;

  void GetActualTrimmedSize(size_t& trimmed_width,
                            size_t& trimmed_height) const;

  aocommon::Image image_;
  const Resources resources_;
  int wtowers_subgrid_size_;
  int wtowers_support_;
  int wtowers_w_support_;
  double wtowers_padding_;
  double wtowers_w_padding_;
  double accuracy_;
  double max_abs_w_;
  std::unique_ptr<WTowersGridderBase> gridder_;
};

}  // namespace wsclean

#endif  // WSCLEAN_WTOWERS_MS_GRIDDER_H
