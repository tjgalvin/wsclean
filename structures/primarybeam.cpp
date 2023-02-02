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

void WriteBeamElement(const ImageFilename& image_name, const Image& beam_image,
                      const Settings& settings, size_t element_index,
                      const aocommon::FitsWriter& writer) {
  writer.Write(image_name.GetBeamPrefix(settings) + "-" +
                   std::to_string(element_index) + ".fits",
               beam_image.Data());
}

aocommon::Image Load(const Settings& settings, const ImageFilename& image_name,
                     size_t element) {
  const std::string filename = image_name.GetBeamPrefix(settings) + "-" +
                               std::to_string(element) + ".fits";
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
                     const std::vector<aocommon::HMC4x4>& beam,
                     const Settings& settings, const ImagingTableEntry& entry,
                     const CoordinateSystem& coordinates,
                     size_t undersampling_factor) {
  std::vector<size_t> required_elements;
  const bool pseudo_correction = settings.polarizations.size() == 1 &&
                                 (entry.polarization == Polarization::RR ||
                                  entry.polarization == Polarization::LL);
  const bool stokes_i_correction = settings.polarizations.size() == 1 &&
                                   entry.polarization == Polarization::StokesI;
  if (pseudo_correction || stokes_i_correction) {
    // m_00 and m_33 ; see aocommon::HMC4x4::Data()
    required_elements = {0, 15};
  } else if (settings.polarizations ==
             std::set<aocommon::PolarizationEnum>{aocommon::Polarization::XX,
                                                  aocommon::Polarization::YY}) {
    // m_00, m_11, m_22 and m_33
    required_elements = {0, 3, 8, 15};
  } else {
    required_elements = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  }
  const aocommon::FitsWriter writer = MakeWriter(coordinates, entry);
  Image upsampled(coordinates.width, coordinates.height);
  for (size_t element : required_elements) {
    Logger::Debug << "Upsampling beam element " << element << "...\n";
    using everybeam::griddedresponse::GriddedResponse;
    GriddedResponse::UpsampleResponse(upsampled.Data(), element,
                                      coordinates.width, coordinates.height,
                                      beam, undersampling_factor);
    WriteBeamElement(image_name, upsampled, settings, element, writer);
  }
}
#endif

void ApplyFacetCorrections(
    const ImageFilename& imageName, const Settings& settings,
    const CoordinateSystem& coordinates, const ImagingTable& table,
    const std::map<size_t, std::unique_ptr<MetaDataCache>>& metaCache) {
  if (settings.polarizations ==
      std::set<PolarizationEnum>{Polarization::XX, Polarization::YY}) {
    // FIXME: to be implemented
    // This should multiply the 16 images (representing a Hermitian 4x4
    // matrix) with the diagonal 4x4 matrix with diagonal entries [1/sqrt(mx*
    // mx) ; 0 ; 0 ; 1/sqrt(my* my)] where mx the weighted h5 sum for the
    // XX-polarization and my the weighted h5sum for the YY-polarization the
    // result is, however, not Hermitian anymore.
    throw std::runtime_error(
        "Correcting the restored image both for H5Parm solutions and beam "
        "effects is not yet implemented for XX/YY imaging.");
  } else {
    schaapcommon::facets::FacetImage facet_image(coordinates.width,
                                                 coordinates.height, 1);

    // table.Front() can be used, because the central frequency and start/end
    // frequency are equal inside a FacetGroup
    const aocommon::FitsWriter writer = MakeWriter(coordinates, table.Front());

    // Process the images one by one to avoid loading all of them in memory
    // at the same time.
    for (size_t i = 0; i != 16; ++i) {
      Image beam_image = Load(settings, imageName, i);
      if (!beam_image.Empty()) {
        std::vector<float*> imagePtr{beam_image.Data()};
        for (const ImagingTableEntry& entry : table) {
          const float m = metaCache.at(entry.index)->h5Sum / entry.imageWeight;
          const float factor = 1.0 / std::sqrt(m);
          facet_image.SetFacet(*entry.facet, true);
          facet_image.MultiplyImageInsideFacet(imagePtr, factor);
        }

        WriteBeamElement(imageName, beam_image, settings, i, writer);
      }
    }
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

PrimaryBeam::~PrimaryBeam() {}

void PrimaryBeam::AddMS(std::unique_ptr<MSDataDescription> description) {
  ms_list_.emplace_back(std::move(description));
}

void PrimaryBeam::CorrectImages(
    aocommon::FitsWriter& writer, const ImageFilename& image_name,
    const std::string& filename_kind, const ImagingTable& table,
    const std::map<size_t, std::unique_ptr<MetaDataCache>>& meta_cache,
    bool requires_gain_correction) {
  if (requires_gain_correction) {
    const CoordinateSystem coordinates{settings_.trimmedImageWidth,
                                       settings_.trimmedImageHeight,
                                       phase_centre_ra_,
                                       phase_centre_dec_,
                                       settings_.pixelScaleX,
                                       settings_.pixelScaleY,
                                       l_shift_,
                                       m_shift_};
    ApplyFacetCorrections(image_name, settings_, coordinates, table,
                          meta_cache);
  }

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
      reader.reset(new aocommon::FitsReader(name.GetPrefix(settings_) + "-" +
                                            filename_kind + ".fits"));
      images[pol_index] = Image(reader->ImageWidth(), reader->ImageHeight());
      reader->Read(images[pol_index].Data());
    }

    float* image_ptrs[4] = {images[0].Data(), images[1].Data(),
                            images[2].Data(), images[3].Data()};
    beam_images.ApplyFullStokes(image_ptrs, settings_.primaryBeamLimit);
    for (size_t pol_index = 0; pol_index != 2; ++pol_index) {
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
  if (settings_.gridderType == GridderType::IDG) {
    PrimaryBeamImageSet beam_images;
    // IDG produces only a Stokes I beam, and has already corrected for the
    // rest. Currently we just load that beam into the diagonal entries of the
    // real component of XX and YY. This is
    // a bit wasteful so might require a better strategy for big images.
    ImageFilename pol_name(image_name);
    pol_name.SetPolarization(aocommon::Polarization::StokesI);
    aocommon::FitsReader reader(pol_name.GetBeamPrefix(settings_) + ".fits");
    beam_images[0] =
        Image(settings_.trimmedImageWidth, settings_.trimmedImageHeight);
    reader.Read(beam_images[0].Data());
    for (size_t i = 0;
         i != settings_.trimmedImageWidth * settings_.trimmedImageHeight; ++i)
      beam_images[0][i] = std::sqrt(beam_images[0][i]);

    // Copy zero entry to images on the diagonal (see aocommon::HMC4x4)
    std::array<size_t, 3> diagonal_entries = {3, 8, 15};
    for (size_t element : diagonal_entries) {
      if (elements.count(element)) beam_images[element] = beam_images[0];
    }
    return beam_images;
  } else {
    PrimaryBeamImageSet beam_images(settings_.trimmedImageWidth,
                                    settings_.trimmedImageHeight);
    for (size_t element = 0; element != beam_images.NImages(); ++element) {
      if (elements.count(element)) {
        beam_images[element] = ::Load(settings_, image_name, element);
      }
    }
    return beam_images;
  }
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
    // TODO: channelFrequency calculation might be telescope specific?
    const ImagingTableEntry::MSInfo& ms_info =
        entry.msData[ms_provider_info.ms_index];
    const MSSelection& selection = *ms_provider_info.selection;
    aocommon::MultiBandData band;
    {
      SynchronizedMS ms = ms_provider_info.provider->MS();
      band =
          aocommon::MultiBandData(ms->spectralWindow(), ms->dataDescription());
    }
    double central_frequency = 0.0;
    for (size_t data_desc_id = 0; data_desc_id != band.DataDescCount();
         ++data_desc_id) {
      const aocommon::BandData sub_band(band[data_desc_id],
                                        selection.ChannelRangeStart(),
                                        selection.ChannelRangeEnd());
      central_frequency += sub_band.CentreFrequency();
    }
    central_frequency /= ms_info.bands.size();

    std::vector<aocommon::HMC4x4> ms_beam;
    const double ms_weight =
        MakeBeamForMS(ms_beam, *ms_provider_info.provider, selection,
                      *image_weights, coordinates, central_frequency, field_id);
    if (result.empty()) {
      result = std::move(ms_beam);
      for (aocommon::HMC4x4& m : result) {
        m *= ms_weight;
      }
    } else {
      assert(ms_beam.size() == result.size());
      for (size_t i = 0; i != result.size(); ++i) {
        result[i] += ms_beam[i] * ms_weight;
      }
    }
    ms_weight_sum += ms_weight;
  }

  // Apply MS weights
  for (size_t i = 0; i != result.size(); ++i) {
    result[i] /= ms_weight_sum;
  }

  WriteBeamImages(image_name, result, settings_, entry, coordinates,
                  undersample_);
}

double PrimaryBeam::MakeBeamForMS(
    std::vector<aocommon::HMC4x4>& result, MSProvider& ms_provider,
    const MSSelection& selection, const ImageWeights& image_weights,
    const aocommon::CoordinateSystem& coordinateSystem,
    double central_frequency, size_t field_id) {
  // Get time info
  double start_time;
  double end_time;
  size_t interval_count;
  std::tie(start_time, end_time, interval_count) = GetTimeInfo(ms_provider);

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

  // Time array and baseline weights only relevant for LOFAR, MWA (and probably
  // SKA-LOW). MWA beam needs scrutiny, this telescope might be amenable to a
  // more efficient implementation
  double ms_weight = 0;
  switch (telescope_type) {
    case everybeam::TelescopeType::kLofarTelescope:
    case everybeam::TelescopeType::kAARTFAAC:
    case everybeam::TelescopeType::kMWATelescope:
    case everybeam::TelescopeType::kOSKARTelescope:
    case everybeam::TelescopeType::kSkaMidTelescope: {
      std::vector<double> baseline_weights(n_baselines * interval_count, 0);
      std::vector<double> time_array(interval_count, 0);
      // Loop over the intervalCounts
      ms_provider.ResetWritePosition();
      for (size_t interval_index = 0; interval_index != interval_count;
           ++interval_index) {
        // Find the mid time step
        double first_time = start_time + (end_time - start_time) *
                                             interval_index / interval_count;
        double last_time = start_time + (end_time - start_time) *
                                            (interval_index + 1) /
                                            interval_count;
        casacore::MEpoch time_epoch = casacore::MEpoch(
            casacore::MVEpoch((0.5 / 86400.0) * (first_time + last_time)),
            time_column(0).getRef());

        // Set value in time array
        time_array[interval_index] = time_epoch.getValue().get() * 86400.0;

        WeightMatrix weights(telescope->GetNrStations());
        CalculateStationWeights(image_weights, weights, ms, ms_provider,
                                selection, last_time);

        // Get the baseline weights from the baseline_weight matrix
        aocommon::UVector<double> interval_weights =
            weights.GetBaselineWeights();
        std::copy(interval_weights.begin(), interval_weights.end(),
                  baseline_weights.begin() + n_baselines * interval_index);
      }
      // Compute MS weight
      ms_weight = std::accumulate(baseline_weights.begin(),
                                  baseline_weights.end(), 0.0);
      result = grid_response->UndersampledIntegratedResponse(
          beam_mode_, time_array, central_frequency, field_id, undersample_,
          baseline_weights);
    } break;
    case everybeam::TelescopeType::kVLATelescope:
    case everybeam::TelescopeType::kATCATelescope:
    case everybeam::TelescopeType::kGMRTTelescope: {
      if (telescope_type == everybeam::TelescopeType::kATCATelescope ||
          telescope_type == everybeam::TelescopeType::kGMRTTelescope) {
        Logger::Warn << "Warning: ATCA and GMRT primary beam corrections have "
                        "not yet been tested!\n";
      }
      // The dish response is time independent, so leaving zero is fine:
      std::vector<double> time_array(1, 0);
      // A weight of 1 is used for these time independent telescopes
      ms_weight = 1.0;
      // baseline weights have no effect on homogeneous arrays, so leave at 1
      std::vector<double> baseline_weights(n_baselines, 1);
      if (settings_.fieldIds[0] == MSSelection::ALL_FIELDS) {
        Logger::Warn
            << "Warning: primary beam correction together with '-fields "
               "ALL' is not properly supported\n";
        Logger::Warn << "       : The beam will be calculated only for the "
                        "first field!\n";
      }
      result = grid_response->UndersampledIntegratedResponse(
          beam_mode_, time_array, central_frequency, field_id, undersample_,
          baseline_weights);
    } break;
    case everybeam::TelescopeType::kUnknownTelescope:
      throw std::runtime_error("Warning: Unknown telescope type!");
  }

  return ms_weight;
}

void PrimaryBeam::CalculateStationWeights(const ImageWeights& imageWeights,
                                          WeightMatrix& baselineWeights,
                                          SynchronizedMS& ms,
                                          MSProvider& msProvider,
                                          const MSSelection& selection,
                                          double endTime) {
  casacore::MSAntenna antenna_table(ms->antenna());
  aocommon::UVector<double> per_antenna_weights(antenna_table.nrow(), 0.0);

  aocommon::MultiBandData multi_band(ms->spectralWindow(),
                                     ms->dataDescription());
  size_t n_channels =
      selection.ChannelRangeEnd() - selection.ChannelRangeStart();
  size_t n_polarizations =
      (msProvider.Polarization() == aocommon::Polarization::Instrumental) ? 4
                                                                          : 1;
  aocommon::UVector<float> weight_array(n_channels * n_polarizations);
  std::unique_ptr<MSReader> ms_reader = msProvider.MakeReader();
  const aocommon::BandData band = multi_band[msProvider.DataDescId()];
  while (ms_reader->CurrentRowAvailable()) {
    MSProvider::MetaData meta_data;
    ms_reader->ReadMeta(meta_data);
    if (meta_data.time >= endTime) break;
    ms_reader->ReadWeights(weight_array.data());

    for (size_t ch = 0; ch != n_channels; ++ch) {
      const double u = meta_data.uInM / band.ChannelWavelength(ch);
      const double v = meta_data.vInM / band.ChannelWavelength(ch);
      const double iw = imageWeights.GetWeight(u, v);
      const double w = weight_array[ch * n_polarizations] * iw;
      baselineWeights.Value(meta_data.antenna1, meta_data.antenna2) += w;
    }
    ms_reader->NextInputRow();
  }
}

std::tuple<double, double, size_t> PrimaryBeam::GetTimeInfo(
    MSProvider& ms_provider) {
  Logger::Debug << "Counting timesteps...\n";
  ms_provider.ResetWritePosition();
  size_t n_timesteps = 0;
  double start_time = 0.0;
  double end_time = 0.0;
  std::unique_ptr<MSReader> ms_reader = ms_provider.MakeReader();
  if (ms_reader->CurrentRowAvailable()) {
    MSProvider::MetaData meta;
    ms_reader->ReadMeta(meta);
    start_time = meta.time;
    end_time = meta.time;
    ++n_timesteps;
    ms_reader->NextInputRow();
    while (ms_reader->CurrentRowAvailable()) {
      ms_reader->ReadMeta(meta);
      if (end_time != meta.time) {
        ++n_timesteps;
        end_time = meta.time;
      }
      ms_reader->NextInputRow();
    }
  }
  if (start_time == end_time) {
    ++end_time;
    --start_time;
  }
  const double total_seconds = end_time - start_time;
  size_t n_intervals =
      std::max<size_t>(1, (total_seconds + seconds_before_beam_update_ - 1) /
                              seconds_before_beam_update_);
  if (n_intervals > n_timesteps) n_intervals = n_timesteps;

  Logger::Debug << "MS spans " << total_seconds << " seconds, dividing in "
                << n_intervals << " intervals.\n";
  return std::make_tuple(start_time, end_time, n_intervals);
}
#endif  // HAVE_EVERYBEAM

size_t PrimaryBeam::computeUndersamplingFactor(const Settings& settings) {
  return std::max(
      std::min(settings.trimmedImageWidth / settings.primaryBeamGridSize,
               settings.trimmedImageHeight / settings.primaryBeamGridSize),
      (size_t)1);
}
