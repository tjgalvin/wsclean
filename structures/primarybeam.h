#ifndef PRIMARY_BEAM_H
#define PRIMARY_BEAM_H

#include <string>

#include "primarybeamimageset.h"

#include "../io/imagefilename.h"

#include "../structures/imagingtable.h"

#include "../msproviders/msprovider.h"

#include "../scheduling/metadatacache.h"

#include <aocommon/fits/fitswriter.h>
#include <aocommon/fits/fitsreader.h>
#include <aocommon/polarization.h>
#include <aocommon/uvector.h>

#ifdef HAVE_EVERYBEAM
#include <EveryBeam/beammode.h>
#include <EveryBeam/beamnormalisationmode.h>
#endif

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
      beamImages.ApplyStokesI(images[0], _settings.primaryBeamLimit);
    } else if (_settings.polarizations.size() == 4 &&
               aocommon::Polarization::HasFullStokesPolarization(
                   _settings.polarizations)) {
      beamImages.ApplyFullStokes(images.data(), _settings.primaryBeamLimit);
    }
  }

  PrimaryBeamImageSet Load(const ImageFilename& imageName) {
    return load(imageName, _settings);
  }

  void AddMS(std::unique_ptr<class MSDataDescription> description);

  void MakeBeamImages(const ImageFilename& imageName,
                      const ImagingTableEntry& entry,
                      std::shared_ptr<class ImageWeights> imageWeights);

  /**
   * @brief Correct images for the primary beam by multiplying the input image
   * by the (simplified) inverse of the beam. Before the beam is applied, the
   * beam is corrected by solutions obtained from an H5 solution file if @param
   * requiresH5Correction is true. In that case, the beam images are overwritten
   * by their corrected counterparts.
   *
   * @param writer FitsWriter
   * @param imageName Image name object from which prefixes or polarization can
   * be derived.
   * @param filenameKind string specifying which image will be corrected
   * @param table Imaging table of a single FacetGroup
   * @param metaCache MSGridder meta data cache, containing image weights an
   * (summed) H5 facet solutions.
   * @param requiresH5Correction Correct beam images for piecewise constant h5
   * solution?
   */
  void CorrectImages(
      class aocommon::FitsWriter& writer, const ImageFilename& imageName,
      const std::string& filenameKind, const ImagingTable& table,
      const std::map<size_t, std::unique_ptr<MetaDataCache>>& metaCache,
      bool requiresH5Correction);

  size_t GetUndersamplingFactor() const { return _undersample; };
  size_t GetBeamUpdateTime() const { return _secondsBeforeBeamUpdate; };

 private:
  // Compute undersampling factor from the primaryBeamGridSize.
  // In case of rectangular images, the undersampling factor is derived
  // from the shortest dimension.
  static size_t computeUndersamplingFactor(const Settings& settings);

  const Settings& _settings;
  double _phaseCentreRA, _phaseCentreDec, _phaseCentreDL, _phaseCentreDM;
  const size_t _undersample;
  const size_t _secondsBeforeBeamUpdate;
#ifdef HAVE_EVERYBEAM
  const everybeam::BeamMode _beamMode;
  const everybeam::BeamNormalisationMode _beamNormalisationMode;
#endif
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
