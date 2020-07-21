#include "primarybeam.h"
#include "wscleansettings.h"
#include "logger.h"

#include "../fitswriter.h"
#include <aocommon/matrix2x2.h>

#include "../fitsreader.h"

#include "../lofar/lbeamimagemaker.h"

#include "../msproviders/msdatadescription.h"
#include "../mwa/mwabeam.h"

#include "../primarybeam/atcabeam.h"
#include "../primarybeam/vlabeam.h"

#include <boost/filesystem/operations.hpp>

#include <boost/algorithm/string/case_conv.hpp>

PrimaryBeam::PrimaryBeam(const WSCleanSettings& settings)
    : _settings(settings),
      _phaseCentreRA(0.0),
      _phaseCentreDec(0.0),
      _phaseCentreDL(0.0),
      _phaseCentreDM(0.0) {}

PrimaryBeam::~PrimaryBeam() {}

void PrimaryBeam::AddMS(std::unique_ptr<class MSDataDescription> description) {
  _msList.emplace_back(std::move(description));
}

void PrimaryBeam::MakeBeamImages(const ImageFilename& imageName,
                                 const ImagingTableEntry& entry,
                                 std::shared_ptr<ImageWeights> imageWeights) {
  bool useExistingBeam = false;
  if (_settings.reusePrimaryBeam) {
    ImageFilename firstPolName(imageName);
    firstPolName.SetPolarization(aocommon::Polarization::XX);
    firstPolName.SetIsImaginary(false);
    std::string f(firstPolName.GetBeamPrefix(_settings) + ".fits");
    if (boost::filesystem::exists(f)) {
      FitsReader reader(f);
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

    {
      // A single outputimage/entry might consist of multiple measurement sets.
      // In that case we only use the first measurement set:
      std::unique_ptr<MSProvider> provider = _msList.front()->GetProvider();
      SynchronizedMS ms(provider->MS());
      Telescope::TelescopeType type = Telescope::GetType(*ms);
      ms.Reset();
      provider.reset();
      switch (type) {
        case Telescope::LOFAR:
        case Telescope::AARTFAAC:
          beamImages = makeLOFARImage(entry, imageWeights);
          break;
        case Telescope::MWA:
          beamImages = PrimaryBeamImageSet(_settings.trimmedImageWidth,
                                           _settings.trimmedImageHeight, 8);
          beamImages.SetToZero();
          makeMWAImage(beamImages, entry);
          break;
        case Telescope::ATCA:
          beamImages = PrimaryBeamImageSet(_settings.trimmedImageWidth,
                                           _settings.trimmedImageHeight, 8);
          beamImages.SetToZero();
          makeATCAImage(beamImages, entry);
          break;
        case Telescope::VLA:
          beamImages = PrimaryBeamImageSet(_settings.trimmedImageWidth,
                                           _settings.trimmedImageHeight, 8);
          beamImages.SetToZero();
          makeVLAImage(beamImages, entry);
          break;
        default:
          throw std::runtime_error("Can't make beam for this telescope");
      }
    }

    if (beamImages.NImages() == 8) {
      // Save the 8 beam images as fits files
      aocommon::PolarizationEnum linPols[4] = {
          aocommon::Polarization::XX, aocommon::Polarization::XY,
          aocommon::Polarization::YX, aocommon::Polarization::YY};
      FitsWriter writer;
      writer.SetImageDimensions(_settings.trimmedImageWidth,
                                _settings.trimmedImageHeight, _phaseCentreRA,
                                _phaseCentreDec, _settings.pixelScaleX,
                                _settings.pixelScaleY);
      writer.SetPhaseCentreShift(_phaseCentreDL, _phaseCentreDM);
      for (size_t i = 0; i != 8; ++i) {
        aocommon::PolarizationEnum p = linPols[i / 2];
        ImageFilename polName(imageName);
        polName.SetPolarization(p);
        polName.SetIsImaginary(i % 2 != 0);
        writer.SetPolarization(p);
        writer.SetFrequency(entry.CentralFrequency(),
                            entry.bandEndFrequency - entry.bandStartFrequency);
        writer.Write<double>(polName.GetBeamPrefix(_settings) + ".fits",
                             beamImages[i].data());
      }
    } else {
      // Save the 16 beam images as fits files
      FitsWriter writer;
      writer.SetImageDimensions(_settings.trimmedImageWidth,
                                _settings.trimmedImageHeight, _phaseCentreRA,
                                _phaseCentreDec, _settings.pixelScaleX,
                                _settings.pixelScaleY);
      writer.SetPhaseCentreShift(_phaseCentreDL, _phaseCentreDM);
      for (size_t i = 0; i != 16; ++i) {
        writer.SetFrequency(entry.CentralFrequency(),
                            entry.bandEndFrequency - entry.bandStartFrequency);
        writer.Write<double>(imageName.GetBeamPrefix(_settings) + "-" +
                                 std::to_string(i) + ".fits",
                             beamImages[i].data());
      }
    }
  }
}

void PrimaryBeam::CorrectImages(FitsWriter& writer,
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
      FitsReader reader(prefix + "-" + filenameKind + ".fits");
      Image image(reader.ImageWidth(), reader.ImageHeight());
      reader.Read(image.data());

      beamImages.ApplyStokesI(image.data());
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
    std::unique_ptr<FitsReader> reader;
    for (size_t polIndex = 0; polIndex != 4; ++polIndex) {
      aocommon::PolarizationEnum pol =
          aocommon::Polarization::IndexToStokes(polIndex);
      ImageFilename name(imageName);
      name.SetPolarization(pol);
      reader.reset(new FitsReader(name.GetPrefix(_settings) + "-" +
                                  filenameKind + ".fits"));
      images[polIndex] = Image(reader->ImageWidth(), reader->ImageHeight());
      reader->Read(images[polIndex].data());
    }

    double* imagePtrs[4] = {images[0].data(), images[1].data(),
                            images[2].data(), images[3].data()};
    beamImages.ApplyFullStokes(imagePtrs);
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
                                      const WSCleanSettings& settings) {
  if (settings.useIDG) {
    PrimaryBeamImageSet beamImages(settings.trimmedImageWidth,
                                   settings.trimmedImageHeight, 8);
    // IDG produces only a Stokes I beam, and has already corrected for the
    // rest. Currently we just load that beam into real component of XX and YY,
    // and set the other 6 images to zero. This is a bit wasteful so might
    // require a better strategy for big images.
    ImageFilename polName(imageName);
    polName.SetPolarization(aocommon::Polarization::StokesI);
    FitsReader reader(polName.GetBeamPrefix(settings) + ".fits");
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
        FitsReader reader(polName.GetBeamPrefix(settings) + ".fits");
        reader.Read(beamImages[i].data());
      }
      return beamImages;
    } catch (std::exception&) {
      PrimaryBeamImageSet beamImages(settings.trimmedImageWidth,
                                     settings.trimmedImageHeight, 16);
      for (size_t i = 0; i != 16; ++i) {
        FitsReader reader(imageName.GetBeamPrefix(settings) + "-" +
                          std::to_string(i) + ".fits");
        reader.Read(beamImages[i].data());
      }
      return beamImages;
    }
  }
}

PrimaryBeamImageSet PrimaryBeam::makeLOFARImage(
    const ImagingTableEntry& entry,
    std::shared_ptr<ImageWeights> imageWeights) {
  LBeamImageMaker lbeam(&entry);
  std::vector<std::unique_ptr<MSProvider>> providers;
  for (size_t i = 0; i != _msList.size(); ++i) {
    providers.emplace_back(_msList[i]->GetProvider());
    lbeam.AddMS(providers.back().get(), &_msList[i]->Selection(), i);
  }
  lbeam.SetUseDifferentialBeam(_settings.useDifferentialLofarBeam);
  lbeam.SetSaveIntermediateImages(_settings.saveATerms);
  lbeam.SetImageDetails(_settings.trimmedImageWidth,
                        _settings.trimmedImageHeight, _settings.pixelScaleX,
                        _settings.pixelScaleY, _phaseCentreRA, _phaseCentreDec,
                        _phaseCentreDL, _phaseCentreDM);
  lbeam.SetImageWeight(std::move(imageWeights));
  lbeam.SetUndersampling(_settings.primaryBeamUndersampling);
  lbeam.SetSecondsBeforeBeamUpdate(_settings.primaryBeamUpdateTime);
  return lbeam.Make();
}

void PrimaryBeam::makeMWAImage(PrimaryBeamImageSet& beamImages,
                               const ImagingTableEntry& entry) {
  MWABeam mwaBeam(&entry);
  std::vector<std::unique_ptr<MSProvider>> providers;
  for (size_t i = 0; i != _msList.size(); ++i) {
    providers.emplace_back(_msList[i]->GetProvider());
    mwaBeam.AddMS(providers.back().get(), &_msList[i]->Selection(), i);
  }
  mwaBeam.SetImageDetails(_settings.trimmedImageWidth,
                          _settings.trimmedImageHeight, _settings.pixelScaleX,
                          _settings.pixelScaleY, _phaseCentreRA,
                          _phaseCentreDec, _phaseCentreDL, _phaseCentreDM);
  mwaBeam.SetUndersampling(_settings.primaryBeamUndersampling);
  if (!_settings.mwaPath.empty()) mwaBeam.SetSearchPath(_settings.mwaPath);
  mwaBeam.Make(beamImages);
}

void PrimaryBeam::makeFromVoltagePattern(PrimaryBeamImageSet& beamImages,
                                         const ImagingTableEntry& entry,
                                         VoltagePattern& vp) {
  std::unique_ptr<MSProvider> msProvider = _msList.front()->GetProvider();
  SynchronizedMS ms(msProvider->MS());
  casacore::ArrayColumn<double> pointingDirCol(
      ms->field(), casacore::MSField::columnName(casacore::MSField::DELAY_DIR));
  size_t fieldRow = _settings.fieldIds[0];
  if (fieldRow == MSSelection::ALL_FIELDS) {
    Logger::Warn << "Warning: primary beam correction together with '-fields "
                    "ALL' is not properly supported\n";
    Logger::Warn
        << "       : The beam will be calculated only for the first field!\n";
    fieldRow = 0;
  }
  casacore::Array<double> pDir = pointingDirCol(fieldRow);
  double pDirRA = *pDir.cbegin();
  double pDirDec = *(pDir.cbegin() + 1);
  vp.Render(beamImages, _settings.pixelScaleX, _settings.pixelScaleY,
            _phaseCentreRA, _phaseCentreDec, pDirRA, pDirDec, _phaseCentreDL,
            _phaseCentreDM, entry.CentralFrequency());
}

void PrimaryBeam::makeATCAImage(PrimaryBeamImageSet& beamImages,
                                const ImagingTableEntry& entry) {
  Logger::Info << "Calculating ATCA primary beam...\n";
  ATCABeam::Band band = ATCABeam::GetBand(entry.CentralFrequency() * 1e-9);
  VoltagePattern vp = ATCABeam::CalculateVoltagePattern(band);
  makeFromVoltagePattern(beamImages, entry, vp);
}

void PrimaryBeam::makeVLAImage(PrimaryBeamImageSet& beamImages,
                               const ImagingTableEntry& entry) {
  Logger::Info << "Calculating VLA primary beam...\n";
  std::string bandName;
  {
    std::unique_ptr<MSProvider> msProvider = _msList.front()->GetProvider();
    SynchronizedMS ms(msProvider->MS());
    casacore::ScalarColumn<casacore::String> bandNameCol(
        ms->spectralWindow(), casacore::MSSpectralWindow::columnName(
                                  casacore::MSSpectralWindow::NAME));
    // TODO for now, just pick the first band. Of course, this should be an
    // integral over all bands!
    size_t bandIndex = entry.msData[0].bands[0].bandIndex;
    bandName = bandNameCol(bandIndex);
  }
  std::array<double, 5> coefs =
      VLABeam::GetCoefficients(bandName, entry.CentralFrequency());
  VoltagePattern vp;
  vp.frequencies.assign(1, entry.CentralFrequency());
  vp.maximumRadiusArcMin = 53.0;
  aocommon::UVector<double> coefsVec(coefs.begin(), coefs.end());
  vp.EvaluatePolynomial(coefsVec, false);
  makeFromVoltagePattern(beamImages, entry, vp);
}
