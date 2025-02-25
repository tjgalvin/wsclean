#ifndef WGRIDDING_MS_GRIDDER_H
#define WGRIDDING_MS_GRIDDER_H

#include "../gridding/msgridder.h"
#include "../structures/resources.h"

#include <aocommon/image.h>

#include <memory>

namespace wsclean {

class WGridderBase;

class WGriddingMSGridder final : public MsGridder {
 public:
  WGriddingMSGridder(const Settings& settings, const Resources& resources,
                     MsProviderCollection& ms_provider_collection,
                     bool use_tuned_wgridder);
  ~WGriddingMSGridder() final;

  void StartInversion() final;
  size_t GridMeasurementSet(const MsProviderCollection::MsData& ms_data) final;
  void GridSharedMeasurementSetChunk(
      bool apply_corrections, size_t n_polarizations, size_t n_rows,
      const double* uvws, const double* frequencies,
      const aocommon::BandData& selected_band,
      const std::pair<size_t, size_t>* antennas,
      const std::complex<float>* visibilities, const size_t* time_offsets,
      size_t n_antennas,
      const std::vector<std::complex<float>>& parm_response) final;
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
  aocommon::Image image_;

  std::unique_ptr<WGridderBase> MakeGridder(size_t width, size_t height) const;

  size_t CalculateConstantMemory() const final;
  size_t CalculateMaxRowsInMemory(int64_t available_memory,
                                  size_t constant_memory,
                                  size_t additional_per_row_consumption,
                                  size_t channel_count,
                                  size_t num_polarizations_stored) const final;

  void GetActualTrimmedSize(size_t& trimmedWidth, size_t& trimmedHeight) const;

  const Resources resources_;
  double accuracy_;
  bool use_tuned_wgridder_;
  std::unique_ptr<WGridderBase> gridder_;
};

}  // namespace wsclean

#endif
