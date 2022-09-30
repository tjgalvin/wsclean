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

void WriteBeamElement(const ImageFilename& imageName, const Image& beam_image,
                      const Settings& settings, size_t element_index,
                      const aocommon::FitsWriter& writer) {
  writer.Write(imageName.GetBeamPrefix(settings) + "-" +
                   std::to_string(element_index) + ".fits",
               beam_image.Data());
}

aocommon::Image Load(const Settings& settings, const ImageFilename& imageName,
                     size_t element) {
  const std::string filename = imageName.GetBeamPrefix(settings) + "-" +
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
void WriteBeamImages(const ImageFilename& imageName,
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
    WriteBeamElement(imageName, upsampled, settings, element, writer);
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
    schaapcommon::facets::FacetImage facetImage(coordinates.width,
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
          facetImage.SetFacet(*entry.facet, true);
          facetImage.MultiplyImageInsideFacet(imagePtr, factor);
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
  const bool frequencyInterpolation = true;
  const bool useChannelFrequency = true;
  const std::string elementResponseModel = settings.beamModel;

  const std::string coeff_path =
      (telescope_type == everybeam::TelescopeType::kMWATelescope)
          ? wsclean::mwa::FindCoeffFile(settings.mwaPath)
          : "";

  ATermSettings aterm_settings;
  aterm_settings.coeff_path = coeff_path;
  aterm_settings.data_column_name = settings.dataColumnName;
  everybeam::Options options = ATermConfig::ConvertToEBOptions(
      *ms, aterm_settings, frequencyInterpolation,
      settings.beamNormalisationMode, useChannelFrequency, elementResponseModel,
      settings.beamMode);

  // Make telescope
  return everybeam::Load(ms.MS(), options);
}
#endif  // HAVE_EVERYBEAM

}  // namespace

PrimaryBeam::PrimaryBeam(const Settings& settings)
    : _settings(settings),
      _phaseCentreRA(0.0),
      _phaseCentreDec(0.0),
      _lShift(0.0),
      _mShift(0.0),
      _undersample(computeUndersamplingFactor(settings)),
      _secondsBeforeBeamUpdate(settings.primaryBeamUpdateTime)
#ifdef HAVE_EVERYBEAM
      ,
      _beamMode(everybeam::ParseBeamMode(settings.beamMode)),
      _beamNormalisationMode(
          everybeam::ParseBeamNormalisationMode(settings.beamNormalisationMode))
#endif
{
}

PrimaryBeam::~PrimaryBeam() {}

void PrimaryBeam::AddMS(std::unique_ptr<MSDataDescription> description) {
  _msList.emplace_back(std::move(description));
}

void PrimaryBeam::CorrectImages(
    aocommon::FitsWriter& writer, const ImageFilename& imageName,
    const std::string& filenameKind, const ImagingTable& table,
    const std::map<size_t, std::unique_ptr<MetaDataCache>>& metaCache,
    bool requiresH5Correction) {
  if (requiresH5Correction) {
    const CoordinateSystem coordinates{_settings.trimmedImageWidth,
                                       _settings.trimmedImageHeight,
                                       _phaseCentreRA,
                                       _phaseCentreDec,
                                       _settings.pixelScaleX,
                                       _settings.pixelScaleY,
                                       _lShift,
                                       _mShift};
    ApplyFacetCorrections(imageName, _settings, coordinates, table, metaCache);
  }

  if (_settings.polarizations.size() == 1 || filenameKind == "psf") {
    const PrimaryBeamImageSet beamImages = LoadStokesI(imageName);
    PolarizationEnum pol = *_settings.polarizations.begin();

    const bool pseudo_correction =
        _settings.polarizations.size() == 1 &&
        (pol == Polarization::RR || pol == Polarization::LL);
    if (pseudo_correction)
      Logger::Warn
          << "Warning: not all polarizations are available for full beam "
             "correction, performing pseudo-Stokes I beam correction.\n";
    if (pol == Polarization::StokesI || pseudo_correction) {
      ImageFilename stokesIName(imageName);
      stokesIName.SetPolarization(pol);
      std::string prefix;
      if (filenameKind == "psf")
        prefix = stokesIName.GetPSFPrefix(_settings);
      else
        prefix = stokesIName.GetPrefix(_settings);
      aocommon::FitsReader reader(prefix + "-" + filenameKind + ".fits");
      Image image(reader.ImageWidth(), reader.ImageHeight());
      reader.Read(image.Data());
      beamImages.ApplyStokesI(image.Data(), _settings.primaryBeamLimit);
      writer.Write(prefix + "-" + filenameKind + "-pb.fits", image.Data());
    } else {
      throw std::runtime_error(
          "Primary beam correction is requested, but this is not supported "
          "when imaging a single polarization that is not Stokes I. Either "
          "image all four polarizations or turn off beam correction.");
    }
  } else if (_settings.polarizations ==
             std::set<aocommon::PolarizationEnum>{aocommon::Polarization::XX,
                                                  aocommon::Polarization::YY}) {
    const PrimaryBeamImageSet beamImages = LoadDiagonal(imageName);
    Image images[2];
    std::unique_ptr<aocommon::FitsReader> reader;
    for (size_t polIndex = 0; polIndex != 2; ++polIndex) {
      const aocommon::PolarizationEnum pol = (polIndex == 0)
                                                 ? aocommon::Polarization::XX
                                                 : aocommon::Polarization::YY;
      ImageFilename name(imageName);
      name.SetPolarization(pol);
      reader.reset(new aocommon::FitsReader(name.GetPrefix(_settings) + "-" +
                                            filenameKind + ".fits"));
      images[polIndex] = Image(reader->ImageWidth(), reader->ImageHeight());
      reader->Read(images[polIndex].Data());
    }

    float* imagePtrs[2] = {images[0].Data(), images[1].Data()};
    beamImages.ApplyDiagonal(imagePtrs, _settings.primaryBeamLimit);

    for (size_t polIndex = 0; polIndex != 2; ++polIndex) {
      const aocommon::PolarizationEnum pol = (polIndex == 0)
                                                 ? aocommon::Polarization::XX
                                                 : aocommon::Polarization::YY;
      ImageFilename name(imageName);
      name.SetPolarization(pol);
      writer.SetPolarization(pol);
      writer.Write(name.GetPrefix(_settings) + "-" + filenameKind + "-pb.fits",
                   images[polIndex].Data());
    }
  } else if (aocommon::Polarization::HasFullStokesPolarization(
                 _settings.polarizations)) {
    const PrimaryBeamImageSet beamImages = LoadFull(imageName);
    Image images[4];
    std::unique_ptr<aocommon::FitsReader> reader;
    for (size_t polIndex = 0; polIndex != 4; ++polIndex) {
      aocommon::PolarizationEnum pol =
          aocommon::Polarization::IndexToStokes(polIndex);
      ImageFilename name(imageName);
      name.SetPolarization(pol);
      reader.reset(new aocommon::FitsReader(name.GetPrefix(_settings) + "-" +
                                            filenameKind + ".fits"));
      images[polIndex] = Image(reader->ImageWidth(), reader->ImageHeight());
      reader->Read(images[polIndex].Data());
    }

    float* imagePtrs[4] = {images[0].Data(), images[1].Data(), images[2].Data(),
                           images[3].Data()};
    beamImages.ApplyFullStokes(imagePtrs, _settings.primaryBeamLimit);
    for (size_t polIndex = 0; polIndex != 4; ++polIndex) {
      aocommon::PolarizationEnum pol =
          aocommon::Polarization::IndexToStokes(polIndex);
      ImageFilename name(imageName);
      name.SetPolarization(pol);
      writer.SetPolarization(pol);
      writer.Write(name.GetPrefix(_settings) + "-" + filenameKind + "-pb.fits",
                   images[polIndex].Data());
    }
  } else {
    throw std::runtime_error(
        "Primary beam correction can only be performed on Stokes I, "
        "polarizations (XX,YY) or when "
        "imaging all four polarizations.");
  }
}

PrimaryBeamImageSet PrimaryBeam::Load(const ImageFilename& imageName,
                                      const std::set<size_t>& elements) {
  assert(!elements.empty());
  // This function will be called for Stokes I, diagonal or Full Jones
  // correction, so we can assume that the first element is always required:
  assert(*elements.begin() == 0);
  if (_settings.gridderType == GridderType::IDG) {
    PrimaryBeamImageSet beamImages;
    // IDG produces only a Stokes I beam, and has already corrected for the
    // rest. Currently we just load that beam into the diagonal entries of the
    // real component of XX and YY. This is
    // a bit wasteful so might require a better strategy for big images.
    ImageFilename polName(imageName);
    polName.SetPolarization(aocommon::Polarization::StokesI);
    aocommon::FitsReader reader(polName.GetBeamPrefix(_settings) + ".fits");
    reader.Read(beamImages[0].Data());
    for (size_t i = 0;
         i != _settings.trimmedImageWidth * _settings.trimmedImageHeight; ++i)
      beamImages[0][i] = std::sqrt(beamImages[0][i]);

    // Copy zero entry to images on the diagonal (see aocommon::HMC4x4)
    std::array<size_t, 3> diagonal_entries = {3, 8, 15};
    for (size_t element : diagonal_entries) {
      if (elements.count(element)) beamImages[element] = beamImages[0];
    }
    return beamImages;
  } else {
    PrimaryBeamImageSet beamImages(_settings.trimmedImageWidth,
                                   _settings.trimmedImageHeight);
    for (size_t element = 0; element != beamImages.NImages(); ++element) {
      if (elements.count(element)) {
        beamImages[element] = ::Load(_settings, imageName, element);
      }
    }
    return beamImages;
  }
}

#ifndef HAVE_EVERYBEAM
void PrimaryBeam::MakeOrReuse(const ImageFilename& imageName,
                              const ImagingTableEntry& entry,
                              std::shared_ptr<ImageWeights> imageWeights,
                              size_t field_id) {
  throw std::runtime_error(
      "PrimaryBeam correction requested, but the software has been compiled "
      "without EveryBeam. Recompile your software and make sure that "
      "cmake finds the EveryBeam library.");
}
#else
void PrimaryBeam::MakeOrReuse(const ImageFilename& imageName,
                              const ImagingTableEntry& entry,
                              std::shared_ptr<ImageWeights> imageWeights,
                              size_t field_id) {
  bool useExistingBeam = false;
  if (_settings.reusePrimaryBeam) {
    ImageFilename firstPolName(imageName);
    firstPolName.SetPolarization(imageName.GetPolarization());
    firstPolName.SetIsImaginary(false);
    std::string f(firstPolName.GetBeamPrefix(_settings) + "-0.fits");
    if (boost::filesystem::exists(f)) {
      aocommon::FitsReader reader(f);
      if (reader.ImageWidth() == _settings.trimmedImageWidth &&
          reader.ImageHeight() == _settings.trimmedImageHeight) {
        useExistingBeam = true;
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
  if (!useExistingBeam) {
    Logger::Info << " == Constructing primary beam ==\n";
    MakeImage(imageName, entry, imageWeights, field_id);
  }
}

void PrimaryBeam::MakeImage(const ImageFilename& imageName,
                            const ImagingTableEntry& entry,
                            std::shared_ptr<ImageWeights> imageWeights,
                            size_t field_id) {
  const size_t width(_settings.trimmedImageWidth);
  const size_t height(_settings.trimmedImageHeight);

  std::vector<std::unique_ptr<MSProvider>> providers;
  for (size_t i = 0; i != _msList.size(); ++i) {
    providers.emplace_back(_msList[i]->GetProvider());
    _msProviders.push_back(
        MSProviderInfo(providers.back().get(), &_msList[i]->Selection(), i));
  }

  aocommon::CoordinateSystem coordinates{width,
                                         height,
                                         _phaseCentreRA,
                                         _phaseCentreDec,
                                         _settings.pixelScaleX,
                                         _settings.pixelScaleY,
                                         _lShift,
                                         _mShift};

  std::vector<aocommon::HMC4x4> result;
  double ms_weight_sum = 0;
  for (const MSProviderInfo& msProviderInfo : _msProviders) {
    // TODO: channelFrequency calculation might be telescope specific?
    const ImagingTableEntry::MSInfo& msInfo =
        entry.msData[msProviderInfo.msIndex];
    const MSSelection& selection = *msProviderInfo.selection;
    aocommon::MultiBandData band;
    {
      SynchronizedMS ms = msProviderInfo.provider->MS();
      band =
          aocommon::MultiBandData(ms->spectralWindow(), ms->dataDescription());
    }
    double centralFrequency = 0.0;
    for (size_t dataDescId = 0; dataDescId != band.DataDescCount();
         ++dataDescId) {
      aocommon::BandData subBand(band[dataDescId],
                                 selection.ChannelRangeStart(),
                                 selection.ChannelRangeEnd());
      centralFrequency += subBand.CentreFrequency();
    }
    centralFrequency /= msInfo.bands.size();

    std::vector<aocommon::HMC4x4> ms_beam;
    const double ms_weight =
        MakeBeamForMS(ms_beam, *msProviderInfo.provider, selection,
                      *imageWeights, coordinates, centralFrequency, field_id);
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

  WriteBeamImages(imageName, result, _settings, entry, coordinates,
                  _undersample);
}

double PrimaryBeam::MakeBeamForMS(
    std::vector<aocommon::HMC4x4>& result, MSProvider& msProvider,
    const MSSelection& selection, const ImageWeights& imageWeights,
    const aocommon::CoordinateSystem& coordinateSystem, double centralFrequency,
    size_t fieldId) {
  // Get time info
  double startTime, endTime;
  size_t intervalCount;
  std::tie(startTime, endTime, intervalCount) = GetTimeInfo(msProvider);

  SynchronizedMS ms = msProvider.MS();
  const everybeam::TelescopeType telescope_type =
      everybeam::GetTelescopeType(*ms);
  std::unique_ptr<everybeam::telescope::Telescope> telescope =
      PrepareEveryBeam(ms, _settings, telescope_type);

  casacore::MEpoch::ScalarColumn timeColumn(
      *ms, ms->columnName(casacore::MSMainEnums::TIME));
  std::size_t nbaselines =
      telescope->GetNrStations() * (telescope->GetNrStations() + 1) / 2;
  std::vector<double> baseline_weights(nbaselines * intervalCount, 0);
  std::vector<double> time_array(intervalCount, 0);

  // Time array and baseline weights only relevant for LOFAR, MWA (and probably
  // SKA-LOW). MWA beam needs scrutiny, this telescope might be amenable to a
  // more efficient implementation
  double ms_weight = 0;
  switch (telescope_type) {
    case everybeam::TelescopeType::kLofarTelescope:
    case everybeam::TelescopeType::kAARTFAAC:
    case everybeam::TelescopeType::kMWATelescope:
    case everybeam::TelescopeType::kOSKARTelescope:
    case everybeam::TelescopeType::kSkaMidTelescope:
      // Loop over the intervalCounts
      msProvider.ResetWritePosition();
      for (size_t intervalIndex = 0; intervalIndex != intervalCount;
           ++intervalIndex) {
        // Find the mid time step
        double firstTime =
            startTime + (endTime - startTime) * intervalIndex / intervalCount;
        double lastTime = startTime + (endTime - startTime) *
                                          (intervalIndex + 1) / intervalCount;
        casacore::MEpoch timeEpoch = casacore::MEpoch(
            casacore::MVEpoch((0.5 / 86400.0) * (firstTime + lastTime)),
            timeColumn(0).getRef());

        // Set value in time array
        time_array[intervalIndex] = timeEpoch.getValue().get() * 86400.0;

        WeightMatrix baselineWeights(telescope->GetNrStations());
        CalculateStationWeights(imageWeights, baselineWeights, ms, msProvider,
                                selection, lastTime);

        // Get the baseline weights from the baselineWeight matrix
        aocommon::UVector<double> interval_weights =
            baselineWeights.GetBaselineWeights();
        std::copy(interval_weights.begin(), interval_weights.end(),
                  baseline_weights.begin() + nbaselines * intervalIndex);
      }
      // Compute MS weight
      ms_weight = std::accumulate(baseline_weights.begin(),
                                  baseline_weights.end(), 0.0);
      break;
    case everybeam::TelescopeType::kVLATelescope:
    case everybeam::TelescopeType::kATCATelescope:
    case everybeam::TelescopeType::kGMRTTelescope:
      if (telescope_type == everybeam::TelescopeType::kATCATelescope ||
          telescope_type == everybeam::TelescopeType::kGMRTTelescope) {
        Logger::Warn << "Warning: ATCA and GMRT primary beam corrections have "
                        "not yet been tested!\n";
      }
      // Assign weight of 1 for these "time independent" telescopes
      ms_weight = 1.0;
      if (_settings.fieldIds[0] == MSSelection::ALL_FIELDS) {
        Logger::Warn
            << "Warning: primary beam correction together with '-fields "
               "ALL' is not properly supported\n";
        Logger::Warn << "       : The beam will be calculated only for the "
                        "first field!\n";
      }
      break;
    case everybeam::TelescopeType::kUnknownTelescope:
      Logger::Warn << "Warning: Unknown telescope type!\n";
      break;
  }

  std::unique_ptr<everybeam::griddedresponse::GriddedResponse> grid_response =
      telescope->GetGriddedResponse(coordinateSystem);

  result = grid_response->UndersampledIntegratedResponse(
      _beamMode, time_array, centralFrequency, fieldId, _undersample,
      baseline_weights);
  return ms_weight;
}

void PrimaryBeam::CalculateStationWeights(const ImageWeights& imageWeights,
                                          WeightMatrix& baselineWeights,
                                          SynchronizedMS& ms,
                                          MSProvider& msProvider,
                                          const MSSelection& selection,
                                          double endTime) {
  casacore::MSAntenna antTable(ms->antenna());
  aocommon::UVector<double> perAntennaWeights(antTable.nrow(), 0.0);

  aocommon::MultiBandData multiband(ms->spectralWindow(),
                                    ms->dataDescription());
  size_t channelCount =
      selection.ChannelRangeEnd() - selection.ChannelRangeStart();
  size_t polarizationCount =
      (msProvider.Polarization() == aocommon::Polarization::Instrumental) ? 4
                                                                          : 1;
  aocommon::UVector<float> weightArr(channelCount * polarizationCount);
  std::unique_ptr<MSReader> msReader = msProvider.MakeReader();
  const aocommon::BandData band = multiband[msProvider.DataDescId()];
  while (msReader->CurrentRowAvailable()) {
    MSProvider::MetaData metaData;
    msReader->ReadMeta(metaData);
    if (metaData.time >= endTime) break;
    msReader->ReadWeights(weightArr.data());

    for (size_t ch = 0; ch != channelCount; ++ch) {
      double u = metaData.uInM / band.ChannelWavelength(ch),
             v = metaData.vInM / band.ChannelWavelength(ch);
      double iw = imageWeights.GetWeight(u, v);
      double w = weightArr[ch * polarizationCount] * iw;
      baselineWeights.Value(metaData.antenna1, metaData.antenna2) += w;
    }
    msReader->NextInputRow();
  }
}

std::tuple<double, double, size_t> PrimaryBeam::GetTimeInfo(
    MSProvider& msProvider) {
  Logger::Debug << "Counting timesteps...\n";
  msProvider.ResetWritePosition();
  size_t timestepCount = 0;
  double startTime = 0.0, endTime = 0.0;
  std::unique_ptr<MSReader> msReader = msProvider.MakeReader();
  if (msReader->CurrentRowAvailable()) {
    MSProvider::MetaData meta;
    msReader->ReadMeta(meta);
    startTime = meta.time;
    endTime = meta.time;
    ++timestepCount;
    msReader->NextInputRow();
    while (msReader->CurrentRowAvailable()) {
      msReader->ReadMeta(meta);
      if (endTime != meta.time) {
        ++timestepCount;
        endTime = meta.time;
      }
      msReader->NextInputRow();
    }
  }
  if (startTime == endTime) {
    ++endTime;
    --startTime;
  }
  const double totalSeconds = endTime - startTime;
  size_t intervalCount =
      std::max<size_t>(1, (totalSeconds + _secondsBeforeBeamUpdate - 1) /
                              _secondsBeforeBeamUpdate);
  if (intervalCount > timestepCount) intervalCount = timestepCount;

  Logger::Debug << "MS spans " << totalSeconds << " seconds, dividing in "
                << intervalCount << " intervals.\n";
  return std::make_tuple(startTime, endTime, intervalCount);
}
#endif  // HAVE_EVERYBEAM

size_t PrimaryBeam::computeUndersamplingFactor(const Settings& settings) {
  return std::max(
      std::min(settings.trimmedImageWidth / settings.primaryBeamGridSize,
               settings.trimmedImageHeight / settings.primaryBeamGridSize),
      (size_t)1);
}
