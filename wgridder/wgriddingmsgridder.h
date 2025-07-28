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
      const aocommon::BandData& selected_band, size_t data_desc_id,
      const std::pair<size_t, size_t>* antennas,
      const std::complex<float>* visibilities, const size_t* time_offsets,
      size_t n_antennas, const std::vector<std::complex<float>>& parm_response,
      const BeamResponseCacheChunk& beam_response) final;
  void FinishInversion() final;

  void StartPredict(std::vector<aocommon::Image>&& images) final;
  size_t PredictMeasurementSet(
      const MsProviderCollection::MsData& ms_data) final;
  void PredictChunk(size_t n_rows, size_t n_channels, const double* frequencies,
                    const double* uvws,
                    std::complex<float>* visibilities) const final;
  void FinishPredict() final;

  std::vector<aocommon::Image> ResultImages() final {
    return {std::move(image_)};
  }

  void FreeImagingData() final {}

  size_t GetSuggestedWGridSize() const final { return 1; }

 private:
  std::unique_ptr<WGridderBase> MakeGridder(size_t width, size_t height) const;

  size_t CalculateConstantMemory() const final;
  /**
   * Returns a suggested maximum nr rows that should be stored in memory. This
   * function is for regular data.
   * @param available_memory System memory available for this gridder.
   * @param constant_memory Constant memory size used by the gridder, e.g. to
   * store the image. Typically the value returned by @ref
   * CalculateConstantMemory().
   * @param additional_per_row_consumption Any external over heads to storing a
   * row of data (not being the data or their uvws).
   * @param per_row_uvw_consumption Data size of storing the uvws (typically 3 x
   * sizeof(double)).
   * @param channel_count Channels per row.
   * @param num_polarizations_stored Polarizations per row (typically 1).
   */
  size_t CalculateMaxRowsInMemory(int64_t available_memory,
                                  size_t constant_memory,
                                  double additional_per_row_consumption,
                                  size_t per_row_uvw_consumption,
                                  size_t channel_count,
                                  size_t num_polarizations_stored) const final;
  /**
   * Returns a suggested maximum nr rows that should be stored in memory. This
   * function is for irregular data that needs to be flattened. As a
   * consequence, each value will have its own uvw values. The parameters are
   * otherwise similar to @ref CalculateMaxRowsInMemory().
   */
  size_t CalculateMaxVisibilitiesInMemory(
      int64_t available_memory, size_t constant_memory,
      double additional_per_visibility_consumption,
      size_t per_visibility_uvw_consumption,
      size_t num_polarizations_stored) const;

  void GetActualTrimmedSize(size_t& trimmedWidth, size_t& trimmedHeight) const;

  void WritePredictChunk(MSProvider* ms_provider, size_t n_rows,
                         size_t n_antennas,
                         const aocommon::BandData& selected_band,
                         const std::vector<MSProvider::MetaData>& metadatas,
                         std::complex<float>* visibilities);

  size_t GridRegularMeasurementSet(const MsProviderCollection::MsData& ms_data);
  size_t GridBdaMeasurementSet(const MsProviderCollection::MsData& ms_data);
  size_t PredictRegularMeasurementSet(
      const MsProviderCollection::MsData& ms_data);
  size_t PredictBdaMeasurementSet(const MsProviderCollection::MsData& ms_data);

  aocommon::Image image_;
  const Resources resources_;
  double accuracy_;
  bool use_tuned_wgridder_;
  std::unique_ptr<WGridderBase> gridder_;
};

}  // namespace wsclean

#endif
