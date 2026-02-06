#ifndef IDG_MS_GRIDDER_H
#define IDG_MS_GRIDDER_H

#ifdef HAVE_IDG

#include "../gridding/msgridder.h"
#include "../structures/resources.h"

#include <idg-api.h>

#include <aocommon/image.h>
#include <aocommon/lane.h>
#include <aocommon/uvector.h>
#include <aocommon/fits/fitswriter.h>

#ifdef HAVE_EVERYBEAM
#include <EveryBeam/aterms/atermbase.h>
#endif  // HAVE_EVERYBEAM

#include "../main/stopwatch.h"

#include "idgcommon.h"

namespace wsclean {

struct ImagingTableEntry;
class ImageFilename;
class Settings;

class IdgMsGridder final : public MsGridder {
 public:
  IdgMsGridder(const Settings& settings, const Resources& resources,
               MsProviderCollection& measurement_sets);

  ~IdgMsGridder() final = default;

  // If we are computing a PSF or have a cached average beam then we only do one
  // pass. Otherwise we do an additional first pass to compute the average beam.
  size_t GetNInversionPasses() const final;
  void StartInversion() final;
  void StartInversionPass(size_t pass_index) final;
  size_t GridMeasurementSet(const MsProviderCollection::MsData& ms_data) final;
  void FinishInversionPass(size_t pass_index) final;
  void FinishInversion() final;

  void StartPredict(std::vector<aocommon::Image>&& images) final;
  size_t PredictMeasurementSet(
      const MsProviderCollection::MsData& ms_data) final;
  void FinishPredict() final;

  std::vector<aocommon::Image> ResultImages() final;

  static void SavePBCorrectedImages(aocommon::FitsWriter& writer,
                                    const ImageFilename& filename,
                                    const std::string& filenameKind,
                                    const Settings& settings);

  static void SaveBeamImage(const ImagingTableEntry& entry,
                            ImageFilename& filename, const Settings& settings,
                            double ra, double dec, double pdl, double pdm,
                            const AverageBeam& average_beam);

  void SetAverageBeam(std::unique_ptr<AverageBeam> average_beam);

  std::unique_ptr<AverageBeam> ReleaseAverageBeam();

 private:
  std::unique_ptr<AverageBeam> _averageBeam;

  size_t GetSuggestedWGridSize() const final {
    return 1;  // TODO
  }

  void readConfiguration();

#ifdef HAVE_EVERYBEAM
  std::unique_ptr<class everybeam::aterms::ATermBase> getATermMaker(
      const MsProviderCollection::MsData& msData);
  bool prepareForMeasurementSet(
      const MsProviderCollection::MsData& ms_data,
      std::unique_ptr<everybeam::aterms::ATermBase>& aTermMaker,
      aocommon::UVector<std::complex<float>>& aTermBuffer,
      idg::api::BufferSetType);
#else
  bool prepareForMeasurementSet(
      const MsProviderCollection::MsData& ms_data,
      aocommon::UVector<std::complex<float>>& aTermBuffer,
      idg::api::BufferSetType);
#endif  // HAVE_EVERYBEAM

  void predictRow(IDGPredictionRow& row,
                  const std::vector<std::string>& antennaNames);
  void computePredictionBuffer(const std::vector<std::string>& antennaNames);

  std::unique_ptr<idg::api::BufferSet> _bufferset;
  size_t _subgridSize;
  aocommon::UVector<double> _image;
  aocommon::UVector<float> _taper_subgrid;
  aocommon::UVector<float> _taper_grid;
  MSProvider* _outputProvider;
  aocommon::MultiBandData _selectedBands;
  idg::api::Type _proxyType;
  int _buffersize;
  idg::api::options_type _options;
  const Resources _resources;
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

#define IdgMsGridder UnavailableGridder

#endif  // HAVE IDG

#endif
