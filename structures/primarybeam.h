#ifndef PRIMARY_BEAM_H
#define PRIMARY_BEAM_H

#include <string>

#include "primarybeamimageset.h"

#include "../io/imagefilename.h"

#include "../structures/imagingtable.h"

#include "../msproviders/msprovider.h"

#include "../scheduling/metadatacache.h"

#include <aocommon/coordinatesystem.h>
#include <aocommon/fits/fitswriter.h>
#include <aocommon/fits/fitsreader.h>
#include <aocommon/polarization.h>
#include <aocommon/uvector.h>

#ifdef HAVE_EVERYBEAM
#include <EveryBeam/beammode.h>
#include <EveryBeam/beamnormalisationmode.h>
#include <EveryBeam/load.h>
#endif

class PrimaryBeam {
 public:
  PrimaryBeam(const class Settings& settings);
  ~PrimaryBeam();

  void SetPhaseCentre(double ra, double dec, double l_shift, double m_shift) {
    _phaseCentreRA = ra;
    _phaseCentreDec = dec;
    _lShift = l_shift;
    _mShift = m_shift;
  }

  void CorrectImages(const ImageFilename& imageName,
                     std::vector<float*>& images) {
    if (_settings.polarizations.size() == 1 &&
        *_settings.polarizations.begin() == aocommon::Polarization::StokesI) {
      const PrimaryBeamImageSet beamImages = LoadStokesI(imageName);
      beamImages.ApplyStokesI(images[0], _settings.primaryBeamLimit);
    } else if (_settings.polarizations.size() == 4 &&
               aocommon::Polarization::HasFullStokesPolarization(
                   _settings.polarizations)) {
      const PrimaryBeamImageSet beamImages = LoadFull(imageName);
      beamImages.ApplyFullStokes(images.data(), _settings.primaryBeamLimit);
    } else
      throw std::runtime_error("Unsupported primary beam correction");
  }

  PrimaryBeamImageSet LoadFull(const ImageFilename& imageName) {
    const std::set<size_t> kFullIndices = {0, 1, 2,  3,  4,  5,  6,  7,
                                           8, 9, 10, 11, 12, 13, 14, 15};
    return Load(imageName, kFullIndices);
  }
  PrimaryBeamImageSet LoadDiagonal(const ImageFilename& imageName) {
    const std::set<size_t> kDiagonalIndices = {0, 3, 8, 15};
    return Load(imageName, kDiagonalIndices);
  }
  PrimaryBeamImageSet LoadStokesI(const ImageFilename& imageName) {
    const std::set<size_t> kStokesIIndices = {0, 15};
    return Load(imageName, kStokesIIndices);
  }

  void AddMS(std::unique_ptr<class MSDataDescription> description);

  void MakeOrReuse(const ImageFilename& imageName,
                   const ImagingTableEntry& entry,
                   std::shared_ptr<class ImageWeights> imageWeights,
                   size_t field_id);

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
      aocommon::FitsWriter& writer, const ImageFilename& imageName,
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
  double _phaseCentreRA, _phaseCentreDec, _lShift, _mShift;
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

  PrimaryBeamImageSet Load(const ImageFilename& imageName,
                           const std::set<size_t>& elements);
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

  void MakeImage(const ImageFilename& imageName, const ImagingTableEntry& entry,
                 std::shared_ptr<ImageWeights> imageWeights, size_t field_id);

  /**
   * Calculate the average beam for one measurement set.
   * @param result The average beam values are assigned to this vector.
   */
  double MakeBeamForMS(std::vector<aocommon::HMC4x4>& result,
                       MSProvider& msProvider, const MSSelection& selection,
                       const ImageWeights& imageWeights,
                       const aocommon::CoordinateSystem& coordinateSystem,
                       double centralFrequency, size_t fieldId);

  std::tuple<double, double, size_t> GetTimeInfo(MSProvider& msProvider);

  void CalculateStationWeights(const ImageWeights& imageWeights,
                               WeightMatrix& baselineWeights,
                               SynchronizedMS& ms, MSProvider& msProvider,
                               const MSSelection& selection, double endTime);
#endif  // HAVE_EVERYBEAM
};

#endif  // PRIMARY_BEAM_H
