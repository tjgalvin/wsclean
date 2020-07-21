#ifndef IDG_MS_GRIDDER_H
#define IDG_MS_GRIDDER_H

#ifdef HAVE_IDG

#include "../wsclean/msgridderbase.h"

//#include "interface.h"
#include <idg-api.h>

#include <aocommon/lane.h>
#include <aocommon/uvector.h>

#include <boost/thread/mutex.hpp>

class IdgMsGridder : public MSGridderBase {
 public:
  IdgMsGridder(const class WSCleanSettings& settings);

  virtual ~IdgMsGridder() final override {}

  virtual void Invert() final override;

  virtual void Predict(Image real) final override;

  virtual void Predict(Image real, Image imaginary) final override;

  virtual Image ImageRealResult() final override;

  virtual Image ImageImaginaryResult() final override;

  static void SavePBCorrectedImages(class FitsWriter& writer,
                                    const class ImageFilename& filename,
                                    const std::string& filenameKind,
                                    const WSCleanSettings& settings);

  static void SaveBeamImage(const class ImagingTableEntry& entry,
                            class ImageFilename& filename,
                            const WSCleanSettings& settings, double ra,
                            double dec, double pdl, double pdm,
                            const MetaDataCache& cache);

 private:
  AverageBeam* _averageBeam;

  virtual size_t getSuggestedWGridSize() const override final {
    return 1;  // TODO
  }

  void gridMeasurementSet(MSGridderBase::MSData& msData);
  void gridThreadFunction();

  void predictMeasurementSet(MSGridderBase::MSData& msData);
  void readConfiguration();

  void setIdgType();

  std::unique_ptr<class ATermBase> getATermMaker(MSGridderBase::MSData& msData);
  bool prepareForMeasurementSet(
      MSGridderBase::MSData& msData, std::unique_ptr<ATermBase>& aTermMaker,
      aocommon::UVector<std::complex<float>>& aTermBuffer,
      idg::api::BufferSetType);

  struct IDGInversionRow : public MSGridderBase::InversionRow {
    size_t antenna1, antenna2, timeIndex;
  };
  struct IDGPredictionRow {
    double uvw[3];
    size_t dataDescId, antenna1, antenna2, timeIndex, rowId;
  };
  void predictRow(IDGPredictionRow& row);
  void computePredictionBuffer(size_t dataDescId);

  std::unique_ptr<idg::api::BufferSet> _bufferset;
  size_t _subgridSize;
  aocommon::UVector<double> _image;
  aocommon::UVector<float> _taper_subgrid;
  aocommon::UVector<float> _taper_grid;
  MSProvider* _outputProvider;
  MultiBandData _selectedBands;
  const WSCleanSettings& _settings;
  idg::api::Type _proxyType;
  int _buffersize;
  idg::api::options_type _options;
};

void init_optimal_taper_1D(int subgridsize, int gridsize, float kernelsize,
                           float padding, float* taper_subgrid,
                           float* taper_grid);
void init_optimal_gridding_taper_1D(int subgridsize, int gridsize,
                                    float kernelsize, float* taper_subgrid,
                                    float* taper_grid);

#else

#include "../wsclean/unavailablegridder.h"

#define IdgMsGridder UnavailableGridder

#endif  // HAVE IDG

#endif
