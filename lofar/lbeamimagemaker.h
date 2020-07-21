#ifndef BEAM_IMAGE_MAKER
#define BEAM_IMAGE_MAKER

#include "../hmatrix4x4.h"
#include <aocommon/polarization.h>
#include "../wsclean/imagingtable.h"
#include "../wsclean/primarybeamimageset.h"

#include <aocommon/uvector.h>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/measures/Measures/MDirection.h>

#ifdef HAVE_LOFAR_BEAM
#include <StationResponse/Station.h>
#endif

#include <set>

class LBeamImageMaker {
 public:
  LBeamImageMaker(const ImagingTableEntry* tableEntry)
      : _tableEntry(tableEntry),
        _undersample(8),
        _secondsBeforeBeamUpdate(1800),
        _useDifferentialBeam(false) {}

  void AddMS(class MSProvider* msProvider, const MSSelection* selection,
             size_t msIndex) {
    _msProviders.push_back(MSProviderInfo(msProvider, selection, msIndex));
  }

  void SetImageWeight(std::shared_ptr<class ImageWeights> imageWeights) {
    _imageWeights = imageWeights;
  }

  void SetImageDetails(size_t width, size_t height, double pixelSizeX,
                       double pixelSizeY, double phaseCentreRA,
                       double phaseCentreDec, double phaseCentreDL,
                       double phaseCentreDM) {
    _width = width;
    _height = height;
    _pixelSizeX = pixelSizeX;
    _pixelSizeY = pixelSizeY;
    _phaseCentreRA = phaseCentreRA;
    _phaseCentreDec = phaseCentreDec;
    _phaseCentreDL = phaseCentreDL;
    _phaseCentreDM = phaseCentreDM;
  }

  PrimaryBeamImageSet Make();

  void SetUseDifferentialBeam(bool useDifferentialBeam) {
    _useDifferentialBeam = useDifferentialBeam;
  }

  void SetUndersampling(size_t undersamplingFactor) {
    _undersample = undersamplingFactor;
  }

  void SetSaveIntermediateImages(bool saveIntermediateImages) {
    _saveIntermediateImages = saveIntermediateImages;
  }

  void SetSecondsBeforeBeamUpdate(size_t seconds) {
    _secondsBeforeBeamUpdate = seconds;
  }

 private:
#ifdef HAVE_LOFAR_BEAM
  class WeightMatrix {
   public:
    explicit WeightMatrix(size_t nAntenna)
        : _nAntenna(nAntenna), _weights(nAntenna * nAntenna, 0) {}
    double& Value(size_t a1, size_t a2) {
      if (a1 < a2)
        return _weights[a1 * _nAntenna + a2];
      else
        return _weights[a1 + a2 * _nAntenna];
    }
    const double& Value(size_t a1, size_t a2) const {
      if (a1 < a2)
        return _weights[a1 * _nAntenna + a2];
      else
        return _weights[a1 + a2 * _nAntenna];
    }

   private:
    size_t _nAntenna;
    aocommon::UVector<double> _weights;
  };

  void makeBeamForMS(std::vector<HMC4x4>& _matrices, MSProvider& msProvider,
                     const MSSelection& selection, double centralFrequency);

  void makeBeamSnapshot(
      const std::vector<LOFAR::StationResponse::Station::Ptr>& stations,
      const WeightMatrix& weights, HMC4x4* matrices, double time,
      double frequency, double subbandFrequency,
      const casacore::MeasFrame& frame);

  void calculateStationWeights(const class ImageWeights& imageWeights,
                               double& totalWeight,
                               WeightMatrix& baselineWeights,
                               SynchronizedMS& ms, MSProvider& msProvider,
                               const MSSelection& selection, double endTime);

  void logWeights(casacore::MeasurementSet& ms,
                  const aocommon::UVector<double>& weights);
#endif

  struct MSProviderInfo {
    MSProviderInfo(MSProvider* _provider, const MSSelection* _selection,
                   size_t _msIndex)
        : provider(_provider), selection(_selection), msIndex(_msIndex) {}
    MSProvider* provider;
    const MSSelection* selection;
    size_t msIndex;
  };

  const ImagingTableEntry* _tableEntry;
  std::vector<MSProviderInfo> _msProviders;

  std::shared_ptr<class ImageWeights> _imageWeights;
  class ImageBufferAllocator* _allocator;

  size_t _width, _height, _sampledWidth, _sampledHeight;
  size_t _undersample, _secondsBeforeBeamUpdate;
  double _pixelSizeX, _pixelSizeY, _phaseCentreRA, _phaseCentreDec,
      _phaseCentreDL, _phaseCentreDM;
  double _sPixelSizeX, _sPixelSizeY, _totalWeightSum;
  bool _useDifferentialBeam, _saveIntermediateImages;
  casacore::MDirection _delayDir, _preappliedDir, _tileBeamDir;
};

#endif
