#ifndef STRUCTURES_PRIMARY_BEAM_H_
#define STRUCTURES_PRIMARY_BEAM_H_

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

class ImageWeights;
class MSDataDescription;
class Settings;

class PrimaryBeam {
 public:
  PrimaryBeam(const Settings& settings);
  ~PrimaryBeam();

  void SetPhaseCentre(double ra, double dec, double l_shift, double m_shift) {
    phase_centre_ra_ = ra;
    phase_centre_dec_ = dec;
    l_shift_ = l_shift;
    m_shift_ = m_shift;
  }

  void CorrectImages(const ImageFilename& image_name,
                     std::vector<float*>& images) {
    if (settings_.polarizations.size() == 1 &&
        *settings_.polarizations.begin() == aocommon::Polarization::StokesI) {
      const PrimaryBeamImageSet beam_images = LoadStokesI(image_name);
      beam_images.ApplyStokesI(images[0], settings_.primaryBeamLimit);
    } else if (settings_.polarizations.size() == 4 &&
               aocommon::Polarization::HasFullStokesPolarization(
                   settings_.polarizations)) {
      const PrimaryBeamImageSet beam_images = LoadFull(image_name);
      beam_images.ApplyFullStokes(images.data(), settings_.primaryBeamLimit);
    } else
      throw std::runtime_error("Unsupported primary beam correction");
  }

  PrimaryBeamImageSet LoadFull(const ImageFilename& image_name) {
    const std::set<size_t> kFullIndices = {0, 1, 2,  3,  4,  5,  6,  7,
                                           8, 9, 10, 11, 12, 13, 14, 15};
    return Load(image_name, kFullIndices);
  }
  PrimaryBeamImageSet LoadDiagonal(const ImageFilename& image_name) {
    const std::set<size_t> kDiagonalIndices = {0, 3, 8, 15};
    return Load(image_name, kDiagonalIndices);
  }
  PrimaryBeamImageSet LoadStokesI(const ImageFilename& image_name) {
    const std::set<size_t> kStokesIIndices = {0, 15};
    return Load(image_name, kStokesIIndices);
  }

  void AddMS(std::unique_ptr<MSDataDescription> description);

  void MakeOrReuse(const ImageFilename& image_name,
                   const ImagingTableEntry& entry,
                   std::shared_ptr<ImageWeights> image_weights,
                   size_t field_id);

  /**
   * @brief Correct images for the primary beam by multiplying the input image
   * by the (simplified) inverse of the beam. Before the beam is applied, the
   * beam is corrected by solutions obtained from an H5 solution file if @param
   * requiresH5Correction is true. In that case, the beam images are overwritten
   * by their corrected counterparts.
   *
   * @param writer FitsWriter
   * @param image_name Image name object from which prefixes or polarization can
   * be derived.
   * @param filename_kind string specifying which image will be corrected
   * @param table Imaging table of a single FacetGroup
   * @param meta_cache MSGridder meta data cache, containing image weights an
   * (summed) H5 facet solutions.
   * @param requires_gain_correction Correct beam images for piecewise constant
   * h5 solution?
   */
  void CorrectImages(
      aocommon::FitsWriter& writer, const ImageFilename& image_name,
      const std::string& filename_kind, const ImagingTable& table,
      const std::map<size_t, std::unique_ptr<MetaDataCache>>& meta_cache,
      bool requires_gain_correction);

  size_t GetUndersamplingFactor() const { return undersample_; };
  size_t GetBeamUpdateTime() const { return seconds_before_beam_update_; };

 private:
  // Compute undersampling factor from the primaryBeamGridSize.
  // In case of rectangular images, the undersampling factor is derived
  // from the shortest dimension.
  static size_t computeUndersamplingFactor(const Settings& settings);

  const Settings& settings_;
  double phase_centre_ra_, phase_centre_dec_, l_shift_, m_shift_;
  const size_t undersample_;
  const size_t seconds_before_beam_update_;
#ifdef HAVE_EVERYBEAM
  const everybeam::BeamMode beam_mode_;
  const everybeam::BeamNormalisationMode beam_normalisation_mode_;
#endif
  std::vector<std::unique_ptr<MSDataDescription>> ms_list_;
  struct MSProviderInfo {
    MSProviderInfo(MSProvider* _provider, const MSSelection* _selection,
                   size_t _ms_index)
        : provider(_provider), selection(_selection), ms_index(_ms_index) {}
    MSProvider* provider;
    const MSSelection* selection;
    size_t ms_index;
  };
  std::vector<MSProviderInfo> ms_providers_;

  PrimaryBeamImageSet Load(const ImageFilename& image_name,
                           const std::set<size_t>& elements);
#ifdef HAVE_EVERYBEAM
  /**
   * @brief Lower triangular matrix representation of baseline weights.
   *
   */
  class WeightMatrix {
   public:
    explicit WeightMatrix(size_t n_antenna)
        : n_antenna_(n_antenna), weights_(n_antenna * n_antenna, 0) {}
    double& Value(size_t a1, size_t a2) {
      if (a1 < a2)
        return weights_[a1 * n_antenna_ + a2];
      else
        return weights_[a1 + a2 * n_antenna_];
    }
    const double& Value(size_t a1, size_t a2) const {
      if (a1 < a2)
        return weights_[a1 * n_antenna_ + a2];
      else
        return weights_[a1 + a2 * n_antenna_];
    }

    /**
     * @brief Get the weights per baseline
     *
     * @return aocommon::UVector<double>
     */
    aocommon::UVector<double> GetBaselineWeights() const {
      int nbaselines = n_antenna_ * (n_antenna_ + 1) / 2;
      aocommon::UVector<double> baseline_weights(nbaselines, 0);

      int index = 0;
      for (size_t a1 = 0; a1 != n_antenna_; ++a1) {
        for (size_t a2 = a1; a2 != n_antenna_; ++a2) {
          baseline_weights[index] = weights_[a1 * n_antenna_ + a2];
          ++index;
        }
      }
      return baseline_weights;
    }

   private:
    size_t n_antenna_;
    aocommon::UVector<double> weights_;
  };

  void MakeImage(const ImageFilename& image_name,
                 const ImagingTableEntry& entry,
                 std::shared_ptr<ImageWeights> image_weights, size_t field_id);

  /**
   * Calculate the average beam for one measurement set.
   * @param result The average beam values are assigned to this vector.
   */
  double MakeBeamForMS(std::vector<aocommon::HMC4x4>& result,
                       MSProvider& ms_provider, const MSSelection& selection,
                       const ImageWeights& image_weights,
                       const aocommon::CoordinateSystem& coordinate_system,
                       double central_frequency, size_t field_id);

  std::tuple<double, double, size_t> GetTimeInfo(MSProvider& ms_provider);

  void CalculateStationWeights(const ImageWeights& image_weights,
                               WeightMatrix& baseline_weights,
                               SynchronizedMS& ms, MSProvider& ms_provider,
                               const MSSelection& selection, double end_time);
#endif  // HAVE_EVERYBEAM
};

#endif  // STRUCTURES_PRIMARY_BEAM_H_
