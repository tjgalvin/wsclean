#include "primarybeam.h"
#include "../msproviders/msreaders/msreader.h"

#include "../main/settings.h"

#include "../io/logger.h"

#include "../structures/imageweights.h"

#include "../msproviders/msdatadescription.h"

#include "../io/findmwacoefffile.h"

#include <boost/filesystem/operations.hpp>
#include <boost/algorithm/string/case_conv.hpp>

#include <stdexcept>

#ifdef HAVE_EVERYBEAM
#include <EveryBeam/load.h>
#include <EveryBeam/aterms/atermconfig.h>
#include <EveryBeam/griddedresponse/griddedresponse.h>
#include <EveryBeam/coords/coordutils.h>

using everybeam::ATermSettings;
using everybeam::aterms::ATermConfig;
#endif

PrimaryBeam::PrimaryBeam(const Settings& settings)
    : _settings(settings),
      _phaseCentreRA(0.0),
      _phaseCentreDec(0.0),
      _phaseCentreDL(0.0),
      _phaseCentreDM(0.0) {
  // Undersampling factor (default 8) and beam update time (default 1800s), see
  // settings.h
  _undersample = _settings.primaryBeamUndersampling;
  _secondsBeforeBeamUpdate = _settings.primaryBeamUpdateTime;
}

PrimaryBeam::~PrimaryBeam() {}

void PrimaryBeam::AddMS(std::unique_ptr<class MSDataDescription> description) {
  _msList.emplace_back(std::move(description));
}

void PrimaryBeam::CorrectImages(aocommon::FitsWriter& writer,
                                const ImageFilename& imageName,
                                const std::string& filenameKind) {
  PrimaryBeamImageSet beamImages = load(imageName, _settings);
  if (_settings.polarizations.size() == 1 || filenameKind == "psf") {
    aocommon::PolarizationEnum pol = *_settings.polarizations.begin();

    if (pol == aocommon::Polarization::StokesI) {
      ImageFilename stokesIName(imageName);
      stokesIName.SetPolarization(pol);
      std::string prefix;
      if (filenameKind == "psf")
        prefix = stokesIName.GetPSFPrefix(_settings);
      else
        prefix = stokesIName.GetPrefix(_settings);
      aocommon::FitsReader reader(prefix + "-" + filenameKind + ".fits");
      Image image(reader.ImageWidth(), reader.ImageHeight());
      reader.Read(image.data());

      beamImages.ApplyStokesI(image.data(), _settings.primaryBeamLimit);
      writer.Write(prefix + "-" + filenameKind + "-pb.fits", image.data());
    } else {
      throw std::runtime_error(
          "Primary beam correction is requested, but this is not supported "
          "when imaging a single polarization that is not Stokes I. Either "
          "image all four polarizations or turn off beam correction.");
    }
  } else if (aocommon::Polarization::HasFullStokesPolarization(
                 _settings.polarizations)) {
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
      reader->Read(images[polIndex].data());
    }

    float* imagePtrs[4] = {images[0].data(), images[1].data(), images[2].data(),
                           images[3].data()};
    beamImages.ApplyFullStokes(imagePtrs, _settings.primaryBeamLimit);
    for (size_t polIndex = 0; polIndex != 4; ++polIndex) {
      aocommon::PolarizationEnum pol =
          aocommon::Polarization::IndexToStokes(polIndex);
      ImageFilename name(imageName);
      name.SetPolarization(pol);
      writer.SetPolarization(pol);
      writer.Write(name.GetPrefix(_settings) + "-" + filenameKind + "-pb.fits",
                   images[polIndex].data());
    }
  } else {
    throw std::runtime_error(
        "Primary beam correction can only be performed on Stokes I or when "
        "imaging all four polarizations.");
  }
}

PrimaryBeamImageSet PrimaryBeam::load(const ImageFilename& imageName,
                                      const Settings& settings) {
  if (settings.useIDG) {
    PrimaryBeamImageSet beamImages(settings.trimmedImageWidth,
                                   settings.trimmedImageHeight, 8);
    // IDG produces only a Stokes I beam, and has already corrected for the
    // rest. Currently we just load that beam into real component of XX and YY,
    // and set the other 6 images to zero. This is a bit wasteful so might
    // require a better strategy for big images.
    ImageFilename polName(imageName);
    polName.SetPolarization(aocommon::Polarization::StokesI);
    aocommon::FitsReader reader(polName.GetBeamPrefix(settings) + ".fits");
    reader.Read(beamImages[0].data());
    for (size_t i = 0;
         i != settings.trimmedImageWidth * settings.trimmedImageHeight; ++i)
      beamImages[0][i] = std::sqrt(beamImages[0][i]);
    std::copy_n(beamImages[0].data(),
                settings.trimmedImageWidth * settings.trimmedImageHeight,
                beamImages[6].data());
    for (size_t i = 1; i != 8; ++i) {
      if (i != 6)
        std::fill_n(beamImages[i].data(),
                    settings.trimmedImageWidth * settings.trimmedImageHeight,
                    0.0);
    }
    return beamImages;
  } else {
    try {
      PrimaryBeamImageSet beamImages(settings.trimmedImageWidth,
                                     settings.trimmedImageHeight, 8);
      aocommon::PolarizationEnum linPols[4] = {
          aocommon::Polarization::XX, aocommon::Polarization::XY,
          aocommon::Polarization::YX, aocommon::Polarization::YY};
      for (size_t i = 0; i != 8; ++i) {
        aocommon::PolarizationEnum p = linPols[i / 2];
        ImageFilename polName(imageName);
        polName.SetPolarization(p);
        polName.SetIsImaginary(i % 2 != 0);
        aocommon::FitsReader reader(polName.GetBeamPrefix(settings) + ".fits");
        reader.Read(beamImages[i].data());
      }
      return beamImages;
    } catch (std::exception&) {
      PrimaryBeamImageSet beamImages(settings.trimmedImageWidth,
                                     settings.trimmedImageHeight, 16);
      for (size_t i = 0; i != 16; ++i) {
        aocommon::FitsReader reader(imageName.GetBeamPrefix(settings) + "-" +
                                    std::to_string(i) + ".fits");
        reader.Read(beamImages[i].data());
      }
      return beamImages;
    }
  }
}

#ifndef HAVE_EVERYBEAM
void PrimaryBeam::MakeBeamImages(const ImageFilename& imageName,
                                 const ImagingTableEntry& entry,
                                 std::shared_ptr<ImageWeights> imageWeights) {
  throw std::runtime_error(
      "PrimaryBeam correction requested, but the software has been compiled "
      "without EveryBeam. Recompile your software and make sure that "
      "cmake finds the EveryBeam library.");
}
#else
void PrimaryBeam::MakeBeamImages(const ImageFilename& imageName,
                                 const ImagingTableEntry& entry,
                                 std::shared_ptr<ImageWeights> imageWeights) {
  if (entry.facet)
    throw std::runtime_error(
        "Making beam images for facets is not implemented");

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

    PrimaryBeamImageSet beamImages;
    beamImages = MakeImage(entry, imageWeights);

    // Save the 16 beam images as fits files
    aocommon::FitsWriter writer;
    writer.SetImageDimensions(_settings.trimmedImageWidth,
                              _settings.trimmedImageHeight, _phaseCentreRA,
                              _phaseCentreDec, _settings.pixelScaleX,
                              _settings.pixelScaleY);
    writer.SetPhaseCentreShift(_phaseCentreDL, _phaseCentreDM);
    for (size_t i = 0; i != 16; ++i) {
      writer.SetFrequency(entry.CentralFrequency(),
                          entry.bandEndFrequency - entry.bandStartFrequency);
      writer.Write(imageName.GetBeamPrefix(_settings) + "-" +
                       std::to_string(i) + ".fits",
                   beamImages[i].data());
    }
  }
}

PrimaryBeamImageSet PrimaryBeam::MakeImage(
    const ImagingTableEntry& entry,
    std::shared_ptr<ImageWeights> imageWeights) {
  std::vector<std::unique_ptr<MSProvider> > providers;
  for (size_t i = 0; i != _msList.size(); ++i) {
    providers.emplace_back(_msList[i]->GetProvider());
    _msProviders.push_back(
        MSProviderInfo(providers.back().get(), &_msList[i]->Selection(), i));
  }

  size_t width(_settings.trimmedImageWidth),
      height(_settings.trimmedImageHeight);
  everybeam::coords::CoordinateSystem coordinateSystem{width,
                                                       height,
                                                       _phaseCentreRA,
                                                       _phaseCentreDec,
                                                       _settings.pixelScaleX,
                                                       _settings.pixelScaleY,
                                                       _phaseCentreDL,
                                                       _phaseCentreDM};

  aocommon::UVector<double> buffer_total(width * height * 16, 0);
  double ms_weight_sum = 0;
  for (const MSProviderInfo& msProviderInfo : _msProviders) {
    // TODO: channelFrequency calculation might be telescope specific?
    const ImagingTableEntry::MSInfo& msInfo =
        entry.msData[msProviderInfo.msIndex];
    const MSSelection& selection = *msProviderInfo.selection;

    SynchronizedMS ms = msProviderInfo.provider->MS();
    MultiBandData band(ms->spectralWindow(), ms->dataDescription());
    ms.Reset();
    double centralFrequency = 0.0;
    for (size_t dataDescId = 0; dataDescId != band.DataDescCount();
         ++dataDescId) {
      BandData subBand(band[dataDescId], selection.ChannelRangeStart(),
                       selection.ChannelRangeEnd());
      centralFrequency += subBand.CentreFrequency();
    }
    centralFrequency /= msInfo.bands.size();

    aocommon::UVector<double> buffer(width * height * 16, 0);
    double ms_weight =
        MakeBeamForMS(buffer, *msProviderInfo.provider, selection,
                      *imageWeights, coordinateSystem, centralFrequency);
    for (size_t i = 0; i != buffer_total.size(); ++i) {
      buffer_total[i] += ms_weight * buffer[i];
    }
    ms_weight_sum += ms_weight;
  }

  // Apply MS weights
  for (size_t i = 0; i != buffer_total.size(); ++i) {
    buffer_total[i] /= ms_weight_sum;
  }

  PrimaryBeamImageSet beamImages(width, height, 16);
  beamImages.SetToZero();

  // Copy buffer_total data into beam_images
  for (size_t p = 0; p != 16; ++p) {
    std::copy_n(buffer_total.data() + p * width * height, width * height,
                &beamImages[p][0]);
  }
  return beamImages;
}

double PrimaryBeam::MakeBeamForMS(
    aocommon::UVector<double>& buffer, MSProvider& msProvider,
    const MSSelection& selection, const ImageWeights& imageWeights,
    const everybeam::coords::CoordinateSystem& coordinateSystem,
    double centralFrequency) {
  SynchronizedMS ms = msProvider.MS();

  // Get time info
  double startTime, endTime;
  size_t intervalCount;
  std::tie(startTime, endTime, intervalCount) = GetTimeInfo(msProvider);

  casacore::MEpoch::ScalarColumn timeColumn(
      *ms, ms->columnName(casacore::MSMainEnums::TIME));

  // Pass the settings to EveryBeam::Options struct
  const bool frequencyInterpolation = true;
  const bool useChannelFrequency = true;
  const std::string elementResponseModel = _settings.beamModel;

  everybeam::TelescopeType telescope_type = everybeam::GetTelescopeType(*ms);
  const std::string coeff_path =
      (telescope_type == everybeam::TelescopeType::kMWATelescope)
          ? wsclean::mwa::FindCoeffFile(_settings.mwaPath)
          : "";

  ATermSettings aterm_settings;
  aterm_settings.coeff_path = coeff_path;
  aterm_settings.data_column_name = _settings.dataColumnName;
  everybeam::Options options = ATermConfig::ConvertToEBOptions(
      *ms, aterm_settings, frequencyInterpolation,
      _settings.useDifferentialLofarBeam, useChannelFrequency,
      elementResponseModel);

  // Make telescope
  std::unique_ptr<everybeam::telescope::Telescope> telescope =
      everybeam::Load(ms.MS(), options);

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
      if (telescope_type == everybeam::TelescopeType::kATCATelescope) {
        Logger::Warn << "ATCA primary beam correction not yet tested!\n";
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

  // Note: field id is hard coded to 0
  grid_response->CalculateIntegratedResponse(buffer.data(), time_array,
                                             centralFrequency, 0, _undersample,
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

  MultiBandData multiband(ms->spectralWindow(), ms->dataDescription());
  size_t channelCount =
      selection.ChannelRangeEnd() - selection.ChannelRangeStart();
  size_t polarizationCount =
      (msProvider.Polarization() == aocommon::Polarization::Instrumental) ? 4
                                                                          : 1;
  aocommon::UVector<float> weightArr(channelCount * polarizationCount);
  std::unique_ptr<MSReader> msReader = msProvider.MakeReader();
  while (msReader->CurrentRowAvailable()) {
    MSProvider::MetaData metaData;
    msReader->ReadMeta(metaData);
    if (metaData.time >= endTime) break;
    const BandData& band(multiband[metaData.dataDescId]);
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
