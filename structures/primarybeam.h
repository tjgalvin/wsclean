#ifndef PRIMARY_BEAM_H
#define PRIMARY_BEAM_H

#include <string>

#include "primarybeamimageset.h"

#include "../io/imagefilename.h"

#include "../structures/imagingtable.h"

#include "../msproviders/msprovider.h"

#include <aocommon/fits/fitswriter.h>
#include <aocommon/fits/fitsreader.h>
#include <aocommon/polarization.h>
#include <aocommon/uvector.h>

namespace everybeam {
namespace coords {
struct CoordinateSystem;
}  // namespace coords
}  // namespace everybeam

class PrimaryBeam {
 public:
  PrimaryBeam(const class Settings& settings);
  ~PrimaryBeam();

  void SetPhaseCentre(double ra, double dec, double dl, double dm) {
    _phaseCentreRA = ra;
    _phaseCentreDec = dec;
    _phaseCentreDL = dl;
    _phaseCentreDM = dm;
  }

  void CorrectImages(const ImageFilename& imageName,
                     std::vector<float*>& images) {
    PrimaryBeamImageSet beamImages = load(imageName, _settings);
    if (_settings.polarizations.size() == 1 &&
        *_settings.polarizations.begin() == aocommon::Polarization::StokesI) {
      beamImages.ApplyStokesI(images[0]);
    } else if (_settings.polarizations.size() == 4 &&
               aocommon::Polarization::HasFullStokesPolarization(
                   _settings.polarizations)) {
      beamImages.ApplyFullStokes(images.data());
    }
  }

  PrimaryBeamImageSet Load(const ImageFilename& imageName) {
    return load(imageName, _settings);
  }

  void AddMS(std::unique_ptr<class MSDataDescription> description);

  void MakeBeamImages(const ImageFilename& imageName,
                      const ImagingTableEntry& entry,
                      std::shared_ptr<class ImageWeights> imageWeights);

  void CorrectImages(class aocommon::FitsWriter& writer,
                     const ImageFilename& imageName,
                     const std::string& filenameKind);

 private:
  const Settings& _settings;
  double _phaseCentreRA, _phaseCentreDec, _phaseCentreDL, _phaseCentreDM;
  size_t _undersample, _secondsBeforeBeamUpdate;
  std::vector<std::unique_ptr<class MSDataDescription>> _msList;
  struct MSProviderInfo {
    MSProviderInfo(MSProvider* _provider, const MSSelection* _selection,
                   size_t _msIndex)
        : provider(_provider), selection(_selection), msIndex(_msIndex) {}
    MSProvider* provider;
    const MSSelection* selection;
    size_t msIndex;
  };
  std::vector<MSProviderInfo> _msProviders;

  static PrimaryBeamImageSet load(const ImageFilename& imageName,
                                  const Settings& settings);
#ifdef HAVE_EVERYBEAM
  /**
   * @brief Lower triangular matrix representation of baseline weights.
   *
   */
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

    /**
     * @brief Get the weights per baseline
     *
     * @return aocommon::UVector<double>
     */
    aocommon::UVector<double> GetBaselineWeights() const {
      int nbaselines = _nAntenna * (_nAntenna + 1) / 2;
      aocommon::UVector<double> baseline_weights(nbaselines, 0);

      int index = 0;
      for (size_t a1 = 0; a1 != _nAntenna; ++a1) {
        for (size_t a2 = a1; a2 != _nAntenna; ++a2) {
          baseline_weights[index] = _weights[a1 * _nAntenna + a2];
          ++index;
        }
      }
      return baseline_weights;
    }

   private:
    size_t _nAntenna;
    aocommon::UVector<double> _weights;
  };

  PrimaryBeamImageSet MakeImage(const ImagingTableEntry& entry,
                                std::shared_ptr<ImageWeights> imageWeights);

  double MakeBeamForMS(
      aocommon::UVector<double>& buffer, MSProvider& msProvider,
      const MSSelection& selection, const ImageWeights& imageWeights,
      const everybeam::coords::CoordinateSystem& coordinateSystem,
      double centralFrequency);

  std::tuple<double, double, size_t> GetTimeInfo(MSProvider& msProvider);

  void CalculateStationWeights(const class ImageWeights& imageWeights,
                               WeightMatrix& baselineWeights,
                               SynchronizedMS& ms, MSProvider& msProvider,
                               const MSSelection& selection, double endTime);
#endif  // HAVE_EVERYBEAM
};

#endif  // PRIMARY_BEAM_H
