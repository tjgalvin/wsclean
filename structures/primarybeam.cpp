#include "primarybeam.h"
#include "../msproviders/msreaders/msreader.h"

#include "../main/settings.h"

#include "../structures/imageweights.h"

#include "../msproviders/msdatadescription.h"

#include "../io/findmwacoefffile.h"

#include <aocommon/image.h>
#include <aocommon/logger.h>
#include <aocommon/multibanddata.h>

#include <schaapcommon/facets/facetimage.h>

#include <boost/filesystem/operations.hpp>
#include <boost/algorithm/string/case_conv.hpp>

#include <optional>
#include <stdexcept>

#ifdef HAVE_EVERYBEAM
#include <EveryBeam/aterms/atermconfig.h>
#include <EveryBeam/griddedresponse/griddedresponse.h>

using everybeam::ATermSettings;
using everybeam::aterms::ATermConfig;
#endif

using aocommon::CoordinateSystem;
using aocommon::Image;
using aocommon::Logger;
using aocommon::Polarization;
using aocommon::PolarizationEnum;

using schaapcommon::reordering::MSSelection;

namespace wsclean {

namespace {

/// Returns a fitswriter initialized for the given coordinates and entry.
aocommon::FitsWriter MakeWriter(const CoordinateSystem& coordinates,
                                const ImagingTableEntry& entry) {
  aocommon::FitsWriter writer;
  writer.SetImageDimensions(coordinates.width, coordinates.height,
                            coordinates.ra, coordinates.dec, coordinates.dl,
                            coordinates.dm);
  writer.SetPhaseCentreShift(coordinates.l_shift, coordinates.m_shift);
  writer.SetFrequency(entry.CentralFrequency(),
                      entry.bandEndFrequency - entry.bandStartFrequency);
  return writer;
}

std::string BeamFilename(const Settings& settings,
                         const ImageFilename& image_name,
                         std::optional<size_t> element_index) {
  const std::string prefix = image_name.GetBeamPrefix(settings);
  return element_index ? prefix + "-" + std::to_string(*element_index) + ".fits"
                       : prefix + ".fits";
}

/**
 * @param element_index when given, the index is used in the filename (e.g.
 * "wsclean-beam-3.fits"). Otherwise, a filename without index is used (e.g.
 * "wsclean-beam.fits").
 */
void WriteBeamElement(const ImageFilename& image_name, const Image& beam_image,
                      const Settings& settings,
                      std::optional<size_t> element_index,
                      const aocommon::FitsWriter& writer) {
  const std::string filename =
      BeamFilename(settings, image_name, element_index);
  writer.Write(filename, beam_image.Data());
}

aocommon::Image Load(const Settings& settings, const ImageFilename& image_name,
                     std::optional<size_t> element_index) {
  const std::string filename =
      BeamFilename(settings, image_name, element_index);
  if (boost::filesystem::exists(filename)) {
    aocommon::FitsReader reader(filename);
    Image image(reader.ImageWidth(), reader.ImageHeight());
    reader.Read(image.Data());
    return image;
  } else {
    return aocommon::Image();
  }
}

#ifdef HAVE_EVERYBEAM
void WriteBeamImages(const ImageFilename& image_name,
                     std::vector<aocommon::HMC4x4>&& beam,
                     const Settings& settings, const ImagingTableEntry& entry,
                     const CoordinateSystem& coordinates,
                     size_t undersampling_factor) {
  using everybeam::griddedresponse::GriddedResponse;
  const bool use_squared_beam = settings.UseFacetCorrections();
  const aocommon::FitsWriter writer = MakeWriter(coordinates, entry);
  Image upsampled(coordinates.width, coordinates.height);
  if (use_squared_beam) {
    for (aocommon::HMC4x4& pixel : beam) {
      if (pixel.Invert()) {
        const double stokes_i_sqrt = std::sqrt(
            0.5 * (pixel.Data(0) + 2.0 * pixel.Data(9) + pixel.Data(15)));
        pixel.Data(0) = stokes_i_sqrt == 0.0 ? 0.0 : 1.0 / stokes_i_sqrt;
      } else {
        pixel.Data(0) = 0.0;
      }
    }
    GriddedResponse::UpsampleCorrection(upsampled.Data(), 0, coordinates.width,
                                        coordinates.height, beam,
                                        undersampling_factor);
    WriteBeamElement(image_name, upsampled, settings, {}, writer);
  } else {
    std::vector<size_t> required_elements;
    const bool pseudo_correction = settings.polarizations.size() == 1 &&
                                   (entry.polarization == Polarization::RR ||
                                    entry.polarization == Polarization::LL);
    const bool stokes_i_correction =
        settings.polarizations.size() == 1 &&
        entry.polarization == Polarization::StokesI;
    const bool diagonal_correction =
        settings.polarizations ==
        std::set{aocommon::Polarization::XX, aocommon::Polarization::YY};
    if (pseudo_correction || stokes_i_correction) {
      // Require: m_00, m_03, m_30 and m_33 ; see aocommon::HMC4x4::Data()
      // m_03 is complex conjugate of m_30. The imaginary value is not
      // necessary for Stokes I correction as it cancels out.
      required_elements = {0, 9, 15};
    } else if (diagonal_correction) {
      // In diagonal (=xx,yy) correction, a 2x2 matrix with the xx-to-xx,
      // xx-to-yy, yy-to-xx and yy-to-yy values needs to be inverted. This
      // matrix inversion is affected (I think) by the imaginary value of
      // the off-diagonal, so we need image 10 too.
      required_elements = {0, 9, 10, 15};
    } else {
      required_elements = {0, 1, 2,  3,  4,  5,  6,  7,
                           8, 9, 10, 11, 12, 13, 14, 15};
    }
    for (size_t element : required_elements) {
      Logger::Debug << "Upsampling beam element " << element << "...\n";
      GriddedResponse::UpsampleCorrection(upsampled.Data(), element,
                                          coordinates.width, coordinates.height,
                                          beam, undersampling_factor);
      WriteBeamElement(image_name, upsampled, settings, element, writer);
    }
  }
}
#endif

void ApplyFacetCorrections(const ImageFilename& image_name,
                           const Settings& settings,
                           const CoordinateSystem& coordinates,
                           const ImagingTable::Group& group,
                           const OutputChannelInfo& channel_info) {
  schaapcommon::facets::FacetImage facet_image(coordinates.width,
                                               coordinates.height, 1);

  // group.front() can be used, because the central frequency and start/end
  // frequency are equal inside a FacetGroup
  const aocommon::FitsWriter writer = MakeWriter(coordinates, *group.front());

  // When facet corrections are applied, only a scalar correction is left to
  // be applied on the images.
  Image beam_image = Load(settings, image_name, {});
  if (!beam_image.Empty()) {
    std::vector<float*> image_pointer = {beam_image.Data()};
    for (const std::shared_ptr<ImagingTableEntry>& entry : group) {
      // The visibilities are weighted by the beam and h5 facet solutions
      // when gridding, and the visibilities themselves are "apparent"
      // causing the images to be weighted by the square of those gains.
      // The images are then per facet fully Mueller corrected to make them
      // flat gain (=true instrinsic flux), after which they are again divided
      // by the sqrt of the Stokes I gain (a scalar correction) to make them
      // approximately flat noise for deconvolution. What's left to do here is
      // therefore to take out that sqrt of Stokes I gain. This is done in
      // two steps: an estimate of the average squared beam is corrected as a
      // smooth image correction, whereas the residual facet solution gains
      // are facet-based. In case the facet gains are near unity, this would
      // result in a smooth image. Otherwise, the beam and solution
      // contributions aren't easily separable (as they are averages of
      // squares), so we just use (full_correction / beam_correction) as
      // factor. The final image needs to be divided by this factor. However,
      // since we are scaling the beam images, we need to apply the inverse of
      // that.
      const double beam_factor =
          channel_info.averageBeamFacetCorrection[entry->facetIndex]
              .GetStokesIValue();
      const double full_factor =
          channel_info.averageFacetCorrection[entry->facetIndex]
              .GetStokesIValue();
      const double correction =
          beam_factor == 0.0 ? full_factor : full_factor / beam_factor;
      facet_image.SetFacet(*entry->facet, true);
      facet_image.MultiplyImageInsideFacet(image_pointer,
                                           std::sqrt(correction));
    }

    WriteBeamElement(image_name, beam_image, settings, {}, writer);
  }
}

#ifdef HAVE_EVERYBEAM
std::unique_ptr<everybeam::telescope::Telescope> PrepareEveryBeam(
    SynchronizedMS& ms, const Settings& settings,
    everybeam::TelescopeType telescope_type) {
  // Pass the settings to EveryBeam::Options struct
  const bool frequency_interpolation = true;
  const bool use_channel_frequency = true;
  const std::string element_response_model = settings.beamModel;

  const std::string coefficients_path =
      (telescope_type == everybeam::TelescopeType::kMWATelescope)
          ? wsclean::mwa::FindCoeffFile(settings.mwaPath)
          : "";

  ATermSettings aterm_settings;
  aterm_settings.coeff_path = coefficients_path;
  aterm_settings.data_column_name = settings.dataColumnName;
  const everybeam::Options options = ATermConfig::ConvertToEBOptions(
      *ms, aterm_settings, frequency_interpolation,
      settings.beamNormalisationMode, use_channel_frequency,
      element_response_model, settings.beamMode);

  // Make telescope
  return everybeam::Load(ms.MS(), options);
}

#endif  // HAVE_EVERYBEAM

}  // namespace

/**
 * Returns the intervals to which a beam is calculated over.
 * The returned vector holds for each interval the start and
 * end row indices, as well as the central time of the interval.
 */
std::vector<BeamInterval> GetBeamIntervals(MSProvider& ms_provider,
                                           double seconds_before_beam_update) {
  double start_time = 0.0;
  double previous_time = 0.0;
  std::unique_ptr<MSReader> ms_reader = ms_provider.MakeReader();
  size_t row = 0;
  size_t start_row = 0;
  std::vector<BeamInterval> result;
  if (ms_reader->CurrentRowAvailable()) {
    MSProvider::MetaData meta;
    ms_reader->ReadMeta(meta);
    start_time = meta.time;
    previous_time = meta.time;
    ms_reader->NextInputRow();
    ++row;
    while (ms_reader->CurrentRowAvailable()) {
      ms_reader->ReadMeta(meta);
      if (std::abs(meta.time - start_time) > seconds_before_beam_update) {
        result.emplace_back(BeamInterval{start_row, row - 1,
                                         (start_time + previous_time) * 0.5});
        start_row = row;
        start_time = meta.time;
      }
      previous_time = meta.time;
      ms_reader->NextInputRow();
      ++row;
    }
    if (row - 1 != start_row)
      result.emplace_back(
          BeamInterval{start_row, row - 1, (start_time + previous_time) * 0.5});
  }

  return result;
}

PrimaryBeam::PrimaryBeam(const Settings& settings)
    : settings_(settings),
      phase_centre_ra_(0.0),
      phase_centre_dec_(0.0),
      l_shift_(0.0),
      m_shift_(0.0),
      undersample_(computeUndersamplingFactor(settings)),
      seconds_before_beam_update_(settings.primaryBeamUpdateTime)
#ifdef HAVE_EVERYBEAM
      ,
      beam_mode_(everybeam::ParseBeamMode(settings.beamMode)),
      beam_normalisation_mode_(
          everybeam::ParseBeamNormalisationMode(settings.beamNormalisationMode))
#endif
{
}

PrimaryBeam::~PrimaryBeam() = default;

void PrimaryBeam::AddMS(std::unique_ptr<MSDataDescription> description) {
  ms_list_.emplace_back(std::move(description));
}

void PrimaryBeam::CorrectBeamForFacetGain(
    const ImageFilename& image_name, const ImagingTable::Group& group,
    const OutputChannelInfo& channel_info) {
  const CoordinateSystem coordinates{settings_.trimmedImageWidth,
                                     settings_.trimmedImageHeight,
                                     phase_centre_ra_,
                                     phase_centre_dec_,
                                     settings_.pixelScaleX,
                                     settings_.pixelScaleY,
                                     l_shift_,
                                     m_shift_};
  ApplyFacetCorrections(image_name, settings_, coordinates, group,
                        channel_info);
}

void PrimaryBeam::CorrectImages(aocommon::FitsWriter& writer,
                                const ImageFilename& image_name,
                                const std::string& filename_kind) {
  if (settings_.polarizations.size() == 1 || filename_kind == "psf") {
    const PrimaryBeamImageSet beam_images = LoadStokesI(image_name);
    PolarizationEnum pol = *settings_.polarizations.begin();

    const bool pseudo_correction =
        settings_.polarizations.size() == 1 &&
        (pol == Polarization::RR || pol == Polarization::LL);
    if (pseudo_correction)
      Logger::Warn
          << "Warning: not all polarizations are available for full beam "
             "correction, performing pseudo-Stokes I beam correction.\n";
    if (pol == Polarization::StokesI || pseudo_correction) {
      ImageFilename stokes_i_name(image_name);
      stokes_i_name.SetPolarization(pol);
      std::string prefix;
      if (filename_kind == "psf")
        prefix = stokes_i_name.GetPSFPrefix(settings_);
      else
        prefix = stokes_i_name.GetPrefix(settings_);
      aocommon::FitsReader reader(prefix + "-" + filename_kind + ".fits");
      Image image(reader.ImageWidth(), reader.ImageHeight());
      reader.Read(image.Data());
      beam_images.ApplyStokesI(image.Data(), settings_.primaryBeamLimit);
      writer.Write(prefix + "-" + filename_kind + "-pb.fits", image.Data());
    } else {
      throw std::runtime_error(
          "Primary beam correction is requested, but this is not supported "
          "when imaging a single polarization that is not Stokes I. Either "
          "image all four polarizations or turn off beam correction.");
    }
  } else if (settings_.polarizations ==
             std::set<aocommon::PolarizationEnum>{aocommon::Polarization::XX,
                                                  aocommon::Polarization::YY}) {
    const PrimaryBeamImageSet beam_images = LoadDiagonal(image_name);
    Image images[2];
    std::unique_ptr<aocommon::FitsReader> reader;
    for (size_t pol_index = 0; pol_index != 2; ++pol_index) {
      const aocommon::PolarizationEnum pol = (pol_index == 0)
                                                 ? aocommon::Polarization::XX
                                                 : aocommon::Polarization::YY;
      ImageFilename name(image_name);
      name.SetPolarization(pol);
      reader = std::make_unique<aocommon::FitsReader>(
          name.GetPrefix(settings_) + "-" + filename_kind + ".fits");
      images[pol_index] = Image(reader->ImageWidth(), reader->ImageHeight());
      reader->Read(images[pol_index].Data());
    }

    float* image_ptrs[2] = {images[0].Data(), images[1].Data()};
    beam_images.ApplyDiagonal(image_ptrs, settings_.primaryBeamLimit);

    for (size_t pol_index = 0; pol_index != 2; ++pol_index) {
      const aocommon::PolarizationEnum pol = (pol_index == 0)
                                                 ? aocommon::Polarization::XX
                                                 : aocommon::Polarization::YY;
      ImageFilename name(image_name);
      name.SetPolarization(pol);
      writer.SetPolarization(pol);
      writer.Write(name.GetPrefix(settings_) + "-" + filename_kind + "-pb.fits",
                   images[pol_index].Data());
    }
  } else if (aocommon::Polarization::HasFullStokesPolarization(
                 settings_.polarizations)) {
    const PrimaryBeamImageSet beam_images = LoadFull(image_name);
    Image images[4];
    std::unique_ptr<aocommon::FitsReader> reader;
    for (size_t pol_index = 0; pol_index != 4; ++pol_index) {
      aocommon::PolarizationEnum pol =
          aocommon::Polarization::IndexToStokes(pol_index);
      ImageFilename name(image_name);
      name.SetPolarization(pol);
      reader = std::make_unique<aocommon::FitsReader>(
          name.GetPrefix(settings_) + "-" + filename_kind + ".fits");
      images[pol_index] = Image(reader->ImageWidth(), reader->ImageHeight());
      reader->Read(images[pol_index].Data());
    }

    float* image_ptrs[4] = {images[0].Data(), images[1].Data(),
                            images[2].Data(), images[3].Data()};
    beam_images.ApplyFullStokes(image_ptrs, settings_.primaryBeamLimit);
    for (size_t pol_index = 0; pol_index != 4; ++pol_index) {
      aocommon::PolarizationEnum pol =
          aocommon::Polarization::IndexToStokes(pol_index);
      ImageFilename name(image_name);
      name.SetPolarization(pol);
      writer.SetPolarization(pol);
      writer.Write(name.GetPrefix(settings_) + "-" + filename_kind + "-pb.fits",
                   images[pol_index].Data());
    }
  } else {
    throw std::runtime_error(
        "Primary beam correction can only be performed on Stokes I, "
        "polarizations (XX,YY) or when "
        "imaging all four polarizations.");
  }
}

PrimaryBeamImageSet PrimaryBeam::Load(const ImageFilename& image_name,
                                      const std::set<size_t>& elements) {
  assert(!elements.empty());
  // This function will be called for Stokes I, diagonal or Full Jones
  // correction, so we can assume that the first element is always required:
  assert(*elements.begin() == 0);
  const bool use_squared_beam = settings_.UseFacetCorrections();
  if (use_squared_beam || settings_.gridderType == GridderType::IDG) {
    PrimaryBeamImageSet beam_images;
    // IDG and facet-based imaging produce only a Stokes I beam, and images
    // have already been corrected for the rest. Currently we just load that
    // beam into the 4 diagonal entries. This is a bit wasteful so might
    // require a better strategy for big images.
    ImageFilename pol_name(image_name);
    pol_name.SetPolarization(aocommon::Polarization::StokesI);
    aocommon::FitsReader reader(pol_name.GetBeamPrefix(settings_) + ".fits");
    beam_images[0] =
        Image(settings_.trimmedImageWidth, settings_.trimmedImageHeight);
    reader.Read(beam_images[0].Data());

    // Copy image zero to images on the diagonal (see aocommon::HMC4x4)
    const std::array<size_t, 3> diagonal_entries = {3, 8, 15};
    for (size_t element : diagonal_entries) {
      if (elements.count(element)) beam_images[element] = beam_images[0];
    }
    return beam_images;
  } else {
    PrimaryBeamImageSet beam_images;
    for (size_t element = 0; element != beam_images.NImages(); ++element) {
      if (elements.count(element)) {
        beam_images[element] = wsclean::Load(settings_, image_name, element);
      }
    }
    return beam_images;
  }
}

void PrimaryBeam::MakeUnitary(const ImagingTableEntry& entry,
                              const ImageFilename& image_name,
                              const Settings& settings) {
  const size_t width = settings.trimmedImageWidth;
  const size_t height = settings.trimmedImageHeight;
  const aocommon::CoordinateSystem coordinates{width,
                                               height,
                                               phase_centre_ra_,
                                               phase_centre_dec_,
                                               settings_.pixelScaleX,
                                               settings_.pixelScaleY,
                                               l_shift_,
                                               m_shift_};
  const aocommon::FitsWriter writer = MakeWriter(coordinates, entry);
  const Image image(width, height, 1.0f);
  Logger::Debug << "Writing unitary beam...\n";
  WriteBeamElement(image_name, image, settings, {}, writer);
}

#ifndef HAVE_EVERYBEAM
void PrimaryBeam::MakeOrReuse(const ImageFilename& image_name,
                              const ImagingTableEntry& entry,
                              std::shared_ptr<ImageWeights> image_weights,
                              size_t field_id) {
  throw std::runtime_error(
      "PrimaryBeam correction requested, but the software has been compiled "
      "without EveryBeam. Recompile your software and make sure that "
      "cmake finds the EveryBeam library.");
}
#else
void PrimaryBeam::MakeOrReuse(const ImageFilename& image_name,
                              const ImagingTableEntry& entry,
                              std::shared_ptr<ImageWeights> image_weights,
                              size_t field_id) {
  bool use_existing_beam = false;
  if (settings_.reusePrimaryBeam) {
    ImageFilename first_pol_name(image_name);
    first_pol_name.SetPolarization(image_name.GetPolarization());
    first_pol_name.SetIsImaginary(false);
    std::string f(first_pol_name.GetBeamPrefix(settings_) + "-0.fits");
    if (boost::filesystem::exists(f)) {
      aocommon::FitsReader reader(f);
      if (reader.ImageWidth() == settings_.trimmedImageWidth &&
          reader.ImageHeight() == settings_.trimmedImageHeight) {
        use_existing_beam = true;
        Logger::Info << "File '" << f
                     << "' exists on disk -- reusing files for primary beam.\n";
      } else {
        Logger::Info << "File '" << f
                     << "' exists on disk but has different dimensions. Beam "
                        "will be recreated.\n";
      }
    } else {
      Logger::Info << "Primary beam not yet available (file '" << f
                   << "' does not exist). Beam will be created.\n";
    }
  }
  if (!use_existing_beam) {
    Logger::Info << " == Constructing primary beam ==\n";
    MakeImage(image_name, entry, image_weights, field_id);
  }
}

void PrimaryBeam::MakeImage(const ImageFilename& image_name,
                            const ImagingTableEntry& entry,
                            std::shared_ptr<ImageWeights> image_weights,
                            size_t field_id) {
  const size_t width(settings_.trimmedImageWidth);
  const size_t height(settings_.trimmedImageHeight);

  std::vector<std::unique_ptr<MSProvider>> providers;
  for (size_t i = 0; i != ms_list_.size(); ++i) {
    providers.emplace_back(ms_list_[i]->GetProvider());
    ms_providers_.push_back(
        MSProviderInfo(providers.back().get(), &ms_list_[i]->Selection(), i));
  }

  aocommon::CoordinateSystem coordinates{width,
                                         height,
                                         phase_centre_ra_,
                                         phase_centre_dec_,
                                         settings_.pixelScaleX,
                                         settings_.pixelScaleY,
                                         l_shift_,
                                         m_shift_};

  std::vector<aocommon::HMC4x4> result;
  double ms_weight_sum = 0;
  for (const MSProviderInfo& ms_provider_info : ms_providers_) {
    const MSSelection& selection = *ms_provider_info.selection;
    double central_frequency;
    {
      aocommon::BandData band(ms_provider_info.provider->Band(),
                              selection.ChannelRangeStart(),
                              selection.ChannelRangeEnd());
      central_frequency = band.CentreFrequency();
    }

    std::vector<aocommon::HMC4x4> ms_beam;
    const double ms_weight =
        MakeBeamForMS(ms_beam, *ms_provider_info.provider, selection,
                      *image_weights, coordinates, central_frequency, field_id);
    if (result.empty()) {
      result = std::move(ms_beam);
      if (ms_weight > 0.0) {
        for (aocommon::HMC4x4& m : result) {
          m *= ms_weight;
        }
        ms_weight_sum += ms_weight;
      } else {
        // In case of a zero weight beam, the beam might contain NaNs or be
        // otherwise bad, so explicitly set it to zero.
        result.assign(result.size(), aocommon::HMC4x4::Zero());
      }
    } else if (ms_weight > 0.0) {
      assert(ms_beam.size() == result.size());
      for (size_t i = 0; i != result.size(); ++i) {
        result[i] += ms_beam[i] * ms_weight;
      }
      ms_weight_sum += ms_weight;
    }
  }

  // Apply MS weights
  for (size_t i = 0; i != result.size(); ++i) {
    result[i] /= ms_weight_sum;
  }

  WriteBeamImages(image_name, std::move(result), settings_, entry, coordinates,
                  undersample_);
}

double PrimaryBeam::MakeBeamForMS(
    std::vector<aocommon::HMC4x4>& result, MSProvider& ms_provider,
    const MSSelection& selection, const ImageWeights& image_weights,
    const aocommon::CoordinateSystem& coordinateSystem,
    double central_frequency, size_t field_id) {
  Logger::Debug << "Counting timesteps...\n";
  const std::vector<BeamInterval> intervals =
      GetBeamIntervals(ms_provider, seconds_before_beam_update_);
  Logger::Debug << "Dividing MS in " << intervals.size() << " intervals.\n";

  SynchronizedMS ms = ms_provider.MS();
  const everybeam::TelescopeType telescope_type =
      everybeam::GetTelescopeType(*ms);
  std::unique_ptr<everybeam::telescope::Telescope> telescope =
      PrepareEveryBeam(ms, settings_, telescope_type);

  casacore::MEpoch::ScalarColumn time_column(
      *ms, ms->columnName(casacore::MSMainEnums::TIME));
  std::size_t n_baselines =
      telescope->GetNrStations() * (telescope->GetNrStations() + 1) / 2;

  std::unique_ptr<everybeam::griddedresponse::GriddedResponse> grid_response =
      telescope->GetGriddedResponse(coordinateSystem);

  // Time array and baseline weights only relevant for LOFAR, MWA and SKA.
  // MWA beam needs scrutiny, this telescope might be amenable to a
  // more efficient implementation
  double ms_weight = 0;
  switch (telescope_type) {
    // These are the telescopes that require time information OR are
    // heterogenous (not all antennas have the same response)
    case everybeam::TelescopeType::kLofarTelescope:
    case everybeam::TelescopeType::kAARTFAAC:
    case everybeam::TelescopeType::kMWATelescope:
    case everybeam::TelescopeType::kOSKARTelescope:
    case everybeam::TelescopeType::kSkaMidTelescope:
    case everybeam::TelescopeType::kOvroLwaTelescope: {
      std::vector<double> baseline_weights(n_baselines * intervals.size(), 0.0);
      std::vector<double> time_array(intervals.size(), 0.0);
      // Loop over the intervalCounts
      std::unique_ptr<MSReader> reader = ms_provider.MakeReader();
      size_t current_row = 0;

      for (size_t interval_index = 0; interval_index != intervals.size();
           ++interval_index) {
        const BeamInterval& interval = intervals[interval_index];

        // Skip to the start row
        while (current_row < interval.start_row) {
          reader->NextInputRow();
          ++current_row;
        }

        // Set value in time array
        time_array[interval_index] = interval.central_time;

        WeightMatrix weights(telescope->GetNrStations());
        CalculateStationWeights(image_weights, weights, ms, *reader, selection,
                                current_row, interval.end_row);

        // Get the baseline weights from the baseline_weight matrix
        aocommon::UVector<double> interval_weights =
            weights.GetBaselineWeights();
        std::copy(interval_weights.begin(), interval_weights.end(),
                  baseline_weights.begin() + n_baselines * interval_index);
      }
      // Compute MS weight
      ms_weight = std::accumulate(baseline_weights.begin(),
                                  baseline_weights.end(), 0.0);
      const bool use_squared_beam = settings_.UseFacetCorrections();
      result = grid_response->UndersampledIntegratedCorrection(
          beam_mode_, time_array, central_frequency, field_id, undersample_,
          baseline_weights, use_squared_beam);
    } break;
    // Using 'default:' gives compatibility with different EveryBeam versions
    default: {
      if (telescope_type == everybeam::TelescopeType::kATCATelescope ||
          telescope_type == everybeam::TelescopeType::kGMRTTelescope) {
        Logger::Warn << "Warning: ATCA and GMRT primary beam corrections have "
                        "not yet been tested extensively!\n";
      }
      // The dish response is time independent, so leaving zero is fine:
      std::vector<double> time_array(1, 0);
      // A weight of 1 is used for these time independent telescopes
      ms_weight = 1.0;
      // baseline weights have no effect on homogeneous arrays, so leave at 1
      std::vector<double> baseline_weights(n_baselines, 1);
      if (settings_.fieldIds[0] == MSSelection::kAllFields) {
        Logger::Warn
            << "Warning: primary beam correction together with '-fields "
               "ALL' is not properly supported\n";
        Logger::Warn << "       : The beam will be calculated only for the "
                        "first field!\n";
      }
      result = grid_response->UndersampledIntegratedCorrection(
          beam_mode_, time_array, central_frequency, field_id, undersample_,
          baseline_weights, false);
    } break;
    case everybeam::TelescopeType::kUnknownTelescope:
      throw std::runtime_error(
          "Unknown telescope type! If your telescope is supposed to be "
          "supported, try upgrading EveryBeam.");
  }

  return ms_weight;
}

void PrimaryBeam::CalculateStationWeights(const ImageWeights& imageWeights,
                                          WeightMatrix& baselineWeights,
                                          SynchronizedMS& ms,
                                          MSReader& ms_reader,
                                          const MSSelection& selection,
                                          size_t& current_row, size_t end_row) {
  casacore::MSAntenna antenna_table(ms->antenna());
  aocommon::UVector<double> per_antenna_weights(antenna_table.nrow(), 0.0);

  aocommon::MultiBandData multi_band(ms->spectralWindow(),
                                     ms->dataDescription());
  const size_t n_channels =
      selection.ChannelRangeEnd() - selection.ChannelRangeStart();
  const size_t n_polarizations = ms_reader.NPolarizations();
  aocommon::UVector<float> weight_array(n_channels * n_polarizations);
  const aocommon::BandData band = multi_band[ms_reader.DataDescId()];
  while (ms_reader.CurrentRowAvailable() && current_row <= end_row) {
    MSProvider::MetaData meta_data;
    ms_reader.ReadMeta(meta_data);
    ms_reader.ReadWeights(weight_array.data());

    for (size_t ch = 0; ch != n_channels; ++ch) {
      const double u = meta_data.uInM / band.ChannelWavelength(ch);
      const double v = meta_data.vInM / band.ChannelWavelength(ch);
      const double iw = imageWeights.GetWeight(u, v);
      const double w = weight_array[ch * n_polarizations] * iw;
      baselineWeights.Value(meta_data.antenna1, meta_data.antenna2) += w;
    }
    ms_reader.NextInputRow();
    ++current_row;
  }
}

#endif  // HAVE_EVERYBEAM

size_t PrimaryBeam::computeUndersamplingFactor(const Settings& settings) {
  return std::max(
      std::min(settings.trimmedImageWidth / settings.primaryBeamGridSize,
               settings.trimmedImageHeight / settings.primaryBeamGridSize),
      (size_t)1);
}

}  // namespace wsclean
