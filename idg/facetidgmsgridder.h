#ifndef FACET_IDG_MS_GRIDDER_H_
#define FACET_IDG_MS_GRIDDER_H_

#ifdef HAVE_IDG

#include "../gridding/msgridder.h"
#include "../structures/resources.h"

#include <idg-api.h>

#include <aocommon/image.h>
#include <aocommon/uvector.h>

#include "../main/stopwatch.h"

#include "idgcommon.h"

namespace wsclean {

struct ImagingTableEntry;
class ImageFilename;
class Settings;

class FacetIdgMsGridder final : public MsGridder {
 public:
  FacetIdgMsGridder(const Settings& settings, const Resources& resources,
                    MsProviderCollection& measurement_sets);

  ~FacetIdgMsGridder() final = default;

  // Unlike the IdgMsGridder, the FacetIdgMsGridder only requires one inversion
  // pass since we forgo of the average beam.
  size_t GetNInversionPasses() const final { return 1; };
  void StartInversion() final;
  void StartInversionPass(size_t pass_index) final {
    ResetVisibilityCounters();
  }
  size_t GridMeasurementSet(const MsProviderCollection::MsData& ms_data) final;
  void FinishInversionPass(size_t pass_index) final {}
  void FinishInversion() final;

  void StartPredict(std::vector<aocommon::Image>&& images) final;
  size_t PredictMeasurementSet(
      const MsProviderCollection::MsData& ms_data) final;
  void FinishPredict() final {}

  std::vector<aocommon::Image> ResultImages() final;

 private:
  size_t GetSuggestedWGridSize() const final { return 1; }

  bool PrepareForMeasurementSet(const MsProviderCollection::MsData& ms_data,
                                idg::api::BufferSetType);

  void PredictRow(IDGPredictionRow& row,
                  const std::vector<std::string>& antennaNames);
  void ComputePredictionBuffer(const std::vector<std::string>& antennaNames);

  std::unique_ptr<idg::api::BufferSet> buffer_set_;
  aocommon::UVector<double> image_;
  MSProvider* output_provider_;
  aocommon::MultiBandData selected_bands_;
  idg::api::Type proxy_type_;
  int buffer_size_;
  idg::api::options_type options_;
  const Resources resources_;
};

void init_optimal_taper_1D(int subgridsize, int gridsize, float kernelsize,
                           float padding, float* taper_subgrid,
                           float* taper_grid);
void init_optimal_gridding_taper_1D(int subgridsize, int gridsize,
                                    float kernelsize, float* taper_subgrid,
                                    float* taper_grid);

}  // namespace wsclean

#else

#include "../gridding/unavailablegridder.h"

#define FacetIdgMsGridder UnavailableGridder

#endif  // HAVE IDG

#endif  // FACET_IDG_MS_GRIDDER_H
