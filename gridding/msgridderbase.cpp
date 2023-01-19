#include "msgridderbase.h"

#include "../math/calculatefftsize.h"

#include "../msproviders/msprovider.h"
#include "../msproviders/msreaders/msreader.h"

#include "../structures/imageweights.h"

#include <aocommon/logger.h>
#include <aocommon/units/angle.h>

#include <limits>
#include <aocommon/multibanddata.h>

#include <schaapcommon/h5parm/h5parm.h>
#include <schaapcommon/h5parm/jonesparameters.h>
#include <schaapcommon/h5parm/soltab.h>

#include <casacore/casa/Arrays/Cube.h>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MCPosition.h>
#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>
#include <casacore/ms/MeasurementSets/MSAntennaColumns.h>

#include <casacore/tables/Tables/ArrColDesc.h>
#include <casacore/tables/Tables/TableRecord.h>

using aocommon::Logger;
using schaapcommon::h5parm::JonesParameters;

namespace {

/**
 * @brief Select unique times from a given MSProvider
 */
std::vector<double> SelectUniqueTimes(MSProvider& msProvider) {
  std::unique_ptr<MSReader> msReader = msProvider.MakeReader();
  std::vector<double> msTimes;
  while (msReader->CurrentRowAvailable()) {
    MSProvider::MetaData metaData;
    msReader->ReadMeta(metaData);
    // Assumes that the time instants in the MS are in ascending order.
    // In case this is violated, the returned vector will contain redundant
    // entries.
    if (msTimes.empty() || metaData.time != msTimes.back()) {
      msTimes.push_back(metaData.time);
    }
    msReader->NextInputRow();
  }
  return msTimes;
}
}  // namespace

// Defined out of class to allow the class the be used in a std::unique_ptr.
MSGridderBase::~MSGridderBase() = default;

MSGridderBase::MSData::MSData()
    : msIndex(0),
      dataDescId(0),
      matchingRows(0),
      totalRowsProcessed(0),
      antennaNames() {}

MSGridderBase::MSGridderBase(const Settings& settings)
    : _metaDataCache(nullptr),
      _settings(settings),
      _actualInversionWidth(0),
      _actualInversionHeight(0),
      _actualPixelSizeX(0),
      _actualPixelSizeY(0),
      _phaseCentreRA(0.0),
      _phaseCentreDec(0.0),
      _l_shift(0.0),
      _m_shift(0.0),
      _mainImageDL(0.0),
      _mainImageDM(0.0),
      _facetIndex(0),
      _facetGroupIndex(0),
      _msIndex(0),
      _isFacet(false),
      _imagePadding(1.0),
      _imageWidth(0),
      _imageHeight(0),
      _trimWidth(0),
      _trimHeight(0),
      _pixelSizeX(settings.pixelScaleX),
      _pixelSizeY(settings.pixelScaleY),
      _wGridSize(settings.nWLayers),
      _actualWGridSize(0),
      _measurementSets(),
      _dataColumnName(settings.dataColumnName),
      _psfMode(PsfMode::kNone),
      _doSubtractModel(false),
      _smallInversion(settings.smallInversion),
      _wLimit(settings.wLimit / 100.0),
      _precalculatedWeightInfo(nullptr),
      _polarization(aocommon::Polarization::StokesI),
      _isComplex(false),
      _weighting(settings.weightMode),
      _isFirstIteration(false),
      _visibilityWeightingMode(settings.visibilityWeightingMode),
      _gridMode(GriddingKernelMode::KaiserBessel),
      _storeImagingWeights(false),
      _theoreticalBeamSize(0.0),
      _hasFrequencies(false),
      _freqHigh(0.0),
      _freqLow(0.0),
      _bandStart(0.0),
      _bandEnd(0.0),
      _startTime(0.0),
      _griddedVisibilityCount(0),
      _totalWeight(0.0),
      _maxGriddedWeight(0.0),
      _visibilityWeightSum(0.0),
      _predictReader(nullptr) {
#ifdef HAVE_EVERYBEAM
  _visibilityModifier.SetBeamInfo(settings.beamMode,
                                  settings.beamNormalisationMode);
#endif
}

std::vector<std::string> MSGridderBase::getAntennaNames(
    const casacore::MSAntenna& msAntenna) {
  const casacore::ScalarColumn<casacore::String> antennaNameColumn(
      msAntenna, msAntenna.columnName(casacore::MSAntenna::NAME));

  std::vector<std::string> antennaNames;
  antennaNames.reserve(antennaNameColumn.nrow());
  for (size_t i = 0; i < antennaNameColumn.nrow(); ++i) {
    antennaNames.push_back(antennaNameColumn(i));
  }
  return antennaNames;
}

void MSGridderBase::initializePointResponse(
    const MSGridderBase::MSData& msData) {
#ifdef HAVE_EVERYBEAM
  if (_settings.applyFacetBeam || _settings.gridWithBeam) {
    const std::string element_response_string =
        !_settings.beamModel.empty() ? _settings.beamModel : "DEFAULT";
    _visibilityModifier.InitializePointResponse(
        msData.msProvider->MS(), _settings.facetBeamUpdateTime,
        element_response_string, msData.bandData.ChannelCount(),
        _settings.dataColumnName, _settings.mwaPath);
  } else {
    _visibilityModifier.SetNoPointResponse();
  }
#else
  if (_settings.applyFacetBeam || _settings.gridWithBeam) {
    throw std::runtime_error(
        "-apply-facet-beam or -grid-with-beam was set, but wsclean was not "
        "compiled "
        "with EveryBeam. Please compile wsclean with EveryBeam to "
        "use the Facet Beam functionality");
  }
#endif
}

void MSGridderBase::StartMeasurementSet(const MSGridderBase::MSData& msData,
                                        bool isPredict) {
  initializePointResponse(msData);
  _msIndex = msData.msIndex;
  if (isPredict) initializePredictReader(*msData.msProvider);
}

void MSGridderBase::initializePredictReader(MSProvider& msProvider) {
  _predictReader = msProvider.MakeReader();
}

void MSGridderBase::initializeBandData(const casacore::MeasurementSet& ms,
                                       MSGridderBase::MSData& msData) {
  msData.bandData = aocommon::MultiBandData(ms)[msData.dataDescId];
  if (Selection(msData.msIndex).HasChannelRange()) {
    msData.startChannel = Selection(msData.msIndex).ChannelRangeStart();
    msData.endChannel = Selection(msData.msIndex).ChannelRangeEnd();
    Logger::Debug << "Selected channels: " << msData.startChannel << '-'
                  << msData.endChannel << '\n';
    if (msData.startChannel >= msData.bandData.ChannelCount() ||
        msData.endChannel > msData.bandData.ChannelCount() ||
        msData.startChannel == msData.endChannel) {
      std::ostringstream str;
      str << "An invalid channel range was specified! Measurement set only has "
          << msData.bandData.ChannelCount()
          << " channels, requested imaging range is " << msData.startChannel
          << " -- " << msData.endChannel << '.';
      throw std::runtime_error(str.str());
    }
  } else {
    msData.startChannel = 0;
    msData.endChannel = msData.bandData.ChannelCount();
  }
}

template <size_t NPolInMSProvider>
void MSGridderBase::calculateWLimits(MSGridderBase::MSData& msData) {
  Logger::Info << "Determining min and max w & theoretical beam size... ";
  Logger::Info.Flush();
  msData.maxW = 0.0;
  msData.maxWWithFlags = 0.0;
  msData.minW = 1e100;
  msData.maxBaselineUVW = 0.0;
  msData.maxBaselineInM = 0.0;
  const aocommon::BandData selectedBand = msData.SelectedBand();
  std::vector<float> weightArray(selectedBand.ChannelCount() *
                                 NPolInMSProvider);
  double curTimestep = -1, firstTime = -1, lastTime = -1;
  size_t nTimesteps = 0;
  std::unique_ptr<MSReader> msReader = msData.msProvider->MakeReader();
  const double smallestWavelength = selectedBand.SmallestWavelength();
  const double longestWavelength = selectedBand.LongestWavelength();
  while (msReader->CurrentRowAvailable()) {
    MSProvider::MetaData metaData;
    msReader->ReadMeta(metaData);

    if (curTimestep != metaData.time) {
      curTimestep = metaData.time;
      ++nTimesteps;
      if (firstTime == -1) firstTime = curTimestep;
      lastTime = curTimestep;
    }

    const double wHi = std::fabs(metaData.wInM / smallestWavelength);
    const double wLo = std::fabs(metaData.wInM / longestWavelength);
    const double baselineInM = std::sqrt(metaData.uInM * metaData.uInM +
                                         metaData.vInM * metaData.vInM +
                                         metaData.wInM * metaData.wInM);
    const double halfWidth = 0.5 * ImageWidth();
    const double halfHeight = 0.5 * ImageHeight();
    if (wHi > msData.maxW || wLo < msData.minW ||
        baselineInM / selectedBand.SmallestWavelength() >
            msData.maxBaselineUVW) {
      msReader->ReadWeights(weightArray.data());
      const float* weightPtr = weightArray.data();
      for (size_t ch = 0; ch != selectedBand.ChannelCount(); ++ch) {
        const double wavelength = selectedBand.ChannelWavelength(ch);
        double wInL = metaData.wInM / wavelength;
        msData.maxWWithFlags = std::max(msData.maxWWithFlags, fabs(wInL));
        if (*weightPtr != 0.0) {
          double uInL = metaData.uInM / wavelength,
                 vInL = metaData.vInM / wavelength,
                 x = uInL * PixelSizeX() * ImageWidth(),
                 y = vInL * PixelSizeY() * ImageHeight(),
                 imagingWeight = GetImageWeights()->GetWeight(uInL, vInL);
          if (imagingWeight != 0.0) {
            if (std::floor(x) > -halfWidth && std::ceil(x) < halfWidth &&
                std::floor(y) > -halfHeight && std::ceil(y) < halfHeight) {
              msData.maxW = std::max(msData.maxW, std::fabs(wInL));
              msData.minW = std::min(msData.minW, std::fabs(wInL));
              msData.maxBaselineUVW =
                  std::max(msData.maxBaselineUVW, baselineInM / wavelength);
              msData.maxBaselineInM =
                  std::max(msData.maxBaselineInM, baselineInM);
            }
          }
        }
        weightPtr += NPolInMSProvider;
      }
    }

    msReader->NextInputRow();
  }

  if (msData.minW == 1e100) {
    msData.minW = 0.0;
    msData.maxWWithFlags = 0.0;
    msData.maxW = 0.0;
  }

  Logger::Info << "DONE (w=[" << msData.minW << ":" << msData.maxW
               << "] lambdas, maxuvw=" << msData.maxBaselineUVW << " lambda)\n";
  if (msData.maxWWithFlags != msData.maxW) {
    Logger::Debug << "Discarded data has higher w value of "
                  << msData.maxWWithFlags << " lambda.\n";
  }

  if (lastTime == firstTime || nTimesteps < 2)
    msData.integrationTime = 1;
  else
    msData.integrationTime = (lastTime - firstTime) / (nTimesteps - 1);
}

template void MSGridderBase::calculateWLimits<1>(MSGridderBase::MSData& msData);
template void MSGridderBase::calculateWLimits<2>(MSGridderBase::MSData& msData);
template void MSGridderBase::calculateWLimits<4>(MSGridderBase::MSData& msData);

void MSGridderBase::initializeMSDataVector(
    std::vector<MSGridderBase::MSData>& msDataVector) {
  if (MeasurementSetCount() == 0)
    throw std::runtime_error(
        "Something is wrong during inversion: no measurement sets given to "
        "inversion algorithm");
  msDataVector = std::vector<MSGridderBase::MSData>(MeasurementSetCount());

  resetMetaData();
  // FIXME: migrate data members to GriddingResult
  _metaDataCache->h5Sum = 0.0;
  _metaDataCache->correctionSum = 0.0;

  bool hasCache = !_metaDataCache->msDataVector.empty();
  if (!hasCache) _metaDataCache->msDataVector.resize(MeasurementSetCount());

  if (!_settings.facetSolutionFiles.empty()) {
    _visibilityModifier.ResetCache(MeasurementSetCount(),
                                   _settings.facetSolutionFiles,
                                   _settings.facetSolutionTables);
  }

  for (size_t i = 0; i != MeasurementSetCount(); ++i) {
    msDataVector[i].msIndex = i;
    initializeMeasurementSet(msDataVector[i], _metaDataCache->msDataVector[i],
                             hasCache);
  }
  calculateOverallMetaData(msDataVector.data());
}

void MSGridderBase::initializeMeasurementSet(MSGridderBase::MSData& msData,
                                             MetaDataCache::Entry& cacheEntry,
                                             bool isCacheInitialized) {
  MSProvider& msProvider = MeasurementSet(msData.msIndex);
  msData.msProvider = &msProvider;
  SynchronizedMS ms(msProvider.MS());
  if (ms->nrow() == 0) throw std::runtime_error("Table has no rows (no data)");

  msData.antennaNames = getAntennaNames(ms->antenna());
  msData.dataDescId = msProvider.DataDescId();

  initializeBandData(*ms, msData);

  if (HasDenormalPhaseCentre())
    Logger::Debug << "Set has denormal phase centre: dl=" << _l_shift
                  << ", dm=" << _m_shift << '\n';

  calculateMSLimits(msData.SelectedBand(), msProvider.StartTime());

  if (isCacheInitialized) {
    msData.maxW = cacheEntry.maxW;
    msData.maxWWithFlags = cacheEntry.maxWWithFlags;
    msData.minW = cacheEntry.minW;
    msData.maxBaselineUVW = cacheEntry.maxBaselineUVW;
    msData.maxBaselineInM = cacheEntry.maxBaselineInM;
    msData.integrationTime = cacheEntry.integrationTime;
  } else {
    if (msProvider.Polarization() == aocommon::Polarization::Instrumental)
      calculateWLimits<4>(msData);
    else if (msProvider.Polarization() ==
             aocommon::Polarization::DiagonalInstrumental)
      calculateWLimits<2>(msData);
    else
      calculateWLimits<1>(msData);
    cacheEntry.maxW = msData.maxW;
    cacheEntry.maxWWithFlags = msData.maxWWithFlags;
    cacheEntry.minW = msData.minW;
    cacheEntry.maxBaselineUVW = msData.maxBaselineUVW;
    cacheEntry.maxBaselineInM = msData.maxBaselineInM;
    cacheEntry.integrationTime = msData.integrationTime;
  }

  if (!_settings.facetSolutionFiles.empty()) {
    _visibilityModifier.SetMSTimes(msData.msIndex,
                                   SelectUniqueTimes(*msData.msProvider));
  }
}

void MSGridderBase::calculateOverallMetaData(const MSData* msDataVector) {
  _maxW = 0.0;
  _minW = std::numeric_limits<double>::max();
  double maxBaseline = 0.0;

  for (size_t i = 0; i != MeasurementSetCount(); ++i) {
    const MSData& msData = msDataVector[i];

    maxBaseline = std::max(maxBaseline, msData.maxBaselineUVW);
    _maxW = std::max(_maxW, msData.maxW);
    _minW = std::min(_minW, msData.minW);
  }
  if (_minW > _maxW) {
    _minW = _maxW;
    Logger::Error
        << "*** Error! ***\n"
           "*** Calculating maximum and minimum w values failed! Make sure the "
           "data selection and scale settings are correct!\n"
           "***\n";
  }

  _theoreticalBeamSize = 1.0 / maxBaseline;
  if (IsFirstIteration()) {
    Logger::Info << "Theoretic beam = "
                 << aocommon::units::Angle::ToNiceString(_theoreticalBeamSize)
                 << "\n";
  }
  if (_wLimit != 0.0) {
    _maxW *= (1.0 - _wLimit);
    if (_maxW < _minW) _maxW = _minW;
  }

  if (!HasTrimSize()) SetTrimSize(ImageWidth(), ImageHeight());

  _actualInversionWidth = ImageWidth();
  _actualInversionHeight = ImageHeight();
  _actualPixelSizeX = PixelSizeX();
  _actualPixelSizeY = PixelSizeY();

  if (SmallInversion()) {
    size_t optWidth, optHeight, minWidth, minHeight;
    CalculateFFTSize(_actualInversionWidth, _actualPixelSizeX,
                     _theoreticalBeamSize, minWidth, optWidth);
    CalculateFFTSize(_actualInversionHeight, _actualPixelSizeY,
                     _theoreticalBeamSize, minHeight, optHeight);
    if (optWidth < _actualInversionWidth ||
        optHeight < _actualInversionHeight) {
      const size_t newWidth =
          std::max(std::min(optWidth, _actualInversionWidth), size_t(32));
      const size_t newHeight =
          std::max(std::min(optHeight, _actualInversionHeight), size_t(32));
      if (IsFirstIteration()) {
        Logger::Info << "Minimal inversion size: " << minWidth << " x "
                     << minHeight << ", using optimal: " << newWidth << " x "
                     << newHeight << "\n";
      }
      _actualPixelSizeX = (double(_actualInversionWidth) * _actualPixelSizeX) /
                          double(newWidth);
      _actualPixelSizeY = (double(_actualInversionHeight) * _actualPixelSizeY) /
                          double(newHeight);
      _actualInversionWidth = newWidth;
      _actualInversionHeight = newHeight;
    } else {
      if (IsFirstIteration()) {
        Logger::Info
            << "Small inversion enabled, but inversion resolution already "
               "smaller than beam size: not using optimization.\n";
      }
    }
  }

  // Always call getSuggestedWGridSize in the first iteration, since it then
  // logs the suggested wgrid size.
  const size_t suggestedGridSize =
      (IsFirstIteration() || !hasWGridSize()) ? getSuggestedWGridSize() : 0;
  _actualWGridSize = hasWGridSize() ? _wGridSize : suggestedGridSize;
}

template <size_t PolarizationCount, DDGainMatrix GainEntry>
void MSGridderBase::writeVisibilities(
    MSProvider& msProvider, const std::vector<std::string>& antennaNames,
    const aocommon::BandData& curBand, std::complex<float>* buffer) {
  assert(GetPsfMode() == PsfMode::kNone);  // The PSF is never predicted.

  if (_visibilityModifier.HasH5Parm()) {
    assert(!_settings.facetRegionFilename.empty());
    MSProvider::MetaData metaData;
    _predictReader->ReadMeta(metaData);
    // When the facet beam is applied, the row will be incremented later in this
    // function
    if (!_settings.applyFacetBeam) {
      _predictReader->NextInputRow();
    }

    _visibilityModifier.CacheParmResponse(metaData.time, antennaNames, curBand,
                                          _msIndex);

    _visibilityModifier.CorrectParmResponse<PolarizationCount, GainEntry>(
        buffer, _msIndex, curBand, antennaNames.size(), metaData.antenna1,
        metaData.antenna2);
  }

#ifdef HAVE_EVERYBEAM
  if (_settings.applyFacetBeam) {
    assert(!_settings.facetRegionFilename.empty());
    MSProvider::MetaData metaData;
    _predictReader->ReadMeta(metaData);
    _predictReader->NextInputRow();

    _visibilityModifier.CacheBeamResponse(metaData.time, metaData.fieldId,
                                          curBand);

    _visibilityModifier.ApplyBeamResponse<PolarizationCount, GainEntry>(
        buffer, curBand, metaData.antenna1, metaData.antenna2);
  }
#endif

  {
    WriterLockManager::LockGuard guard = _writerLockManager->GetLock(
        _facetGroupIndex * MeasurementSetCount() + _msIndex);
    msProvider.WriteModel(buffer, IsFacet());
  }
  msProvider.NextOutputRow();
}

template void MSGridderBase::writeVisibilities<1, DDGainMatrix::kXX>(
    MSProvider& msProvider, const std::vector<std::string>& antennaNames,
    const aocommon::BandData& curBand, std::complex<float>* buffer);

template void MSGridderBase::writeVisibilities<1, DDGainMatrix::kYY>(
    MSProvider& msProvider, const std::vector<std::string>& antennaNames,
    const aocommon::BandData& curBand, std::complex<float>* buffer);

template void MSGridderBase::writeVisibilities<1, DDGainMatrix::kTrace>(
    MSProvider& msProvider, const std::vector<std::string>& antennaNames,
    const aocommon::BandData& curBand, std::complex<float>* buffer);

template void MSGridderBase::writeVisibilities<2, DDGainMatrix::kFull>(
    MSProvider& msProvider, const std::vector<std::string>& antennaNames,
    const aocommon::BandData& curBand, std::complex<float>* buffer);

template void MSGridderBase::writeVisibilities<4, DDGainMatrix::kFull>(
    MSProvider& msProvider, const std::vector<std::string>& antennaNames,
    const aocommon::BandData& curBand, std::complex<float>* buffer);

template <size_t PolarizationCount, DDGainMatrix GainEntry>
void MSGridderBase::readAndWeightVisibilities(
    MSReader& msReader, const std::vector<std::string>& antennaNames,
    InversionRow& rowData, const aocommon::BandData& curBand,
    float* weightBuffer, std::complex<float>* modelBuffer,
    const bool* isSelected) {
  const std::size_t dataSize = curBand.ChannelCount() * PolarizationCount;
  if (GetPsfMode() != PsfMode::kNone) {
    // Visibilities for a point source at the phase centre are all ones
    std::fill_n(rowData.data, dataSize, 1.0);
    double dl = 0.0;
    double dm = 0.0;
    if (GetPsfMode() == PsfMode::kSingle) {
      // The point source is shifted to the centre of the main image
      dl = MainImageDL();
      dm = MainImageDM();
    } else {  // GetPsfMode() == PsfMode::kDirectionDependent
      // The point source is shifted to the centre of the current DdPsf position
      dl = LShift();
      dm = MShift();
    }
    if (dl != 0.0 || dm != 0.0) {
      const double dn = std::sqrt(1.0 - dl * dl - dm * dm) - 1.0;
      const double shiftFactor =
          2.0 * M_PI *
          (rowData.uvw[0] * dl + rowData.uvw[1] * dm + rowData.uvw[2] * dn);
      rotateVisibilities<PolarizationCount>(curBand, shiftFactor, rowData.data);
    }
  } else {
    msReader.ReadData(rowData.data);
  }
  rowData.rowId = msReader.RowId();

  if (DoSubtractModel()) {
    msReader.ReadModel(modelBuffer);
    std::complex<float>* modelIter = modelBuffer;
    for (std::complex<float>* iter = rowData.data;
         iter != rowData.data + dataSize; ++iter) {
      *iter -= *modelIter;
      modelIter++;
    }
  }

  msReader.ReadWeights(weightBuffer);

  // Any visibilities that are not gridded in this pass
  // should not contribute to the weight sum, so set these
  // to have zero weight.
  for (size_t i = 0; i != dataSize; ++i) {
    if (!isSelected[i]) weightBuffer[i] = 0.0;
  }

  switch (GetVisibilityWeightingMode()) {
    case VisibilityWeightingMode::NormalVisibilityWeighting:
      // The weight buffer already contains the visibility weights: do nothing
      break;
    case VisibilityWeightingMode::SquaredVisibilityWeighting:
      // Square the visibility weights
      for (size_t i = 0; i != dataSize; ++i) weightBuffer[i] *= weightBuffer[i];
      break;
    case VisibilityWeightingMode::UnitVisibilityWeighting:
      // Set the visibility weights to one
      for (size_t i = 0; i != dataSize; ++i) {
        if (weightBuffer[i] != 0.0) weightBuffer[i] = 1.0f;
      }
      break;
  }

  // Precompute imaging weights
  _scratchImageWeights.resize(curBand.ChannelCount());
  for (size_t ch = 0; ch != curBand.ChannelCount(); ++ch) {
    const double u = rowData.uvw[0] / curBand.ChannelWavelength(ch);
    const double v = rowData.uvw[1] / curBand.ChannelWavelength(ch);
    _scratchImageWeights[ch] = GetImageWeights()->GetWeight(u, v);
  }
  if (StoreImagingWeights())
    msReader.WriteImagingWeights(_scratchImageWeights.data());

  if (IsFacet() && (GetPsfMode() != PsfMode::kSingle)) {
    const bool apply_forward = GetPsfMode() == PsfMode::kDirectionDependent;
    if ((_settings.applyFacetBeam || _settings.gridWithBeam) &&
        _visibilityModifier.HasH5Parm()) {
#ifdef HAVE_EVERYBEAM
      // Load and apply (in conjugate) both the beam and the h5parm solutions
      MSProvider::MetaData metaData;
      msReader.ReadMeta(metaData);
      _visibilityModifier.CacheBeamResponse(metaData.time, metaData.fieldId,
                                            curBand);
      _visibilityModifier.CacheParmResponse(metaData.time, antennaNames,
                                            curBand, _msIndex);
      const VisibilityModifier::DualResult result =
          _visibilityModifier.ApplyConjugatedDual<PolarizationCount, GainEntry>(
              rowData.data, weightBuffer, _scratchImageWeights.data(), curBand,
              antennaNames, metaData.antenna1, metaData.antenna2, _msIndex,
              apply_forward);
      _metaDataCache->h5Sum += result.h5Sum;
      _metaDataCache->correctionSum += result.correctionSum;
    } else if (_settings.applyFacetBeam || _settings.gridWithBeam) {
      // Load and apply only the conjugate beam
      MSProvider::MetaData metaData;
      msReader.ReadMeta(metaData);
      _visibilityModifier.CacheBeamResponse(metaData.time, metaData.fieldId,
                                            curBand);
      _metaDataCache->correctionSum +=
          _visibilityModifier
              .ApplyConjugatedBeamResponse<PolarizationCount, GainEntry>(
                  rowData.data, weightBuffer, _scratchImageWeights.data(),
                  curBand, metaData.antenna1, metaData.antenna2, apply_forward);
#endif  // HAVE_EVERYBEAM
    } else if (_visibilityModifier.HasH5Parm()) {
      MSProvider::MetaData metaData;
      msReader.ReadMeta(metaData);

      _visibilityModifier.CacheParmResponse(metaData.time, antennaNames,
                                            curBand, _msIndex);

      _metaDataCache->correctionSum +=
          _visibilityModifier
              .CorrectConjugatedParmResponse<PolarizationCount, GainEntry>(
                  rowData.data, weightBuffer, _scratchImageWeights.data(),
                  _msIndex, curBand, antennaNames.size(), metaData.antenna1,
                  metaData.antenna2, apply_forward);
    }
  }

  // Apply visibility and imaging weights
  std::complex<float>* dataIter = rowData.data;
  float* weightIter = weightBuffer;
  for (size_t ch = 0; ch != curBand.ChannelCount(); ++ch) {
    for (size_t p = 0; p != PolarizationCount; ++p) {
      const double cumWeight = *weightIter * _scratchImageWeights[ch];
      if (p == 0 && cumWeight != 0.0) {
        // Visibility weight sum is the sum of weights excluding imaging weights
        _visibilityWeightSum += *weightIter;
        _maxGriddedWeight = std::max(cumWeight, _maxGriddedWeight);
        ++_griddedVisibilityCount;
        // Total weight includes imaging weights
        _totalWeight += cumWeight;
      }
      *weightIter = cumWeight;
      *dataIter *= cumWeight;
      ++dataIter;
      ++weightIter;
    }
  }
}

template void MSGridderBase::readAndWeightVisibilities<1, DDGainMatrix::kXX>(
    MSReader& msReader, const std::vector<std::string>& antennaNames,
    InversionRow& newItem, const aocommon::BandData& curBand,
    float* weightBuffer, std::complex<float>* modelBuffer,
    const bool* isSelected);

template void MSGridderBase::readAndWeightVisibilities<1, DDGainMatrix::kYY>(
    MSReader& msReader, const std::vector<std::string>& antennaNames,
    InversionRow& newItem, const aocommon::BandData& curBand,
    float* weightBuffer, std::complex<float>* modelBuffer,
    const bool* isSelected);

template void MSGridderBase::readAndWeightVisibilities<1, DDGainMatrix::kTrace>(
    MSReader& msReader, const std::vector<std::string>& antennaNames,
    InversionRow& newItem, const aocommon::BandData& curBand,
    float* weightBuffer, std::complex<float>* modelBuffer,
    const bool* isSelected);

template void MSGridderBase::readAndWeightVisibilities<2, DDGainMatrix::kFull>(
    MSReader& msReader, const std::vector<std::string>& antennaNames,
    InversionRow& newItem, const aocommon::BandData& curBand,
    float* weightBuffer, std::complex<float>* modelBuffer,
    const bool* isSelected);

template void MSGridderBase::readAndWeightVisibilities<4, DDGainMatrix::kFull>(
    MSReader& msReader, const std::vector<std::string>& antennaNames,
    InversionRow& newItem, const aocommon::BandData& curBand,
    float* weightBuffer, std::complex<float>* modelBuffer,
    const bool* isSelected);

template <size_t PolarizationCount>
void MSGridderBase::rotateVisibilities(const aocommon::BandData& bandData,
                                       double shiftFactor,
                                       std::complex<float>* dataIter) {
  for (size_t ch = 0; ch != bandData.ChannelCount(); ++ch) {
    const double wShiftRad = shiftFactor / bandData.ChannelWavelength(ch);
    const std::complex<float> phasor(std::cos(wShiftRad), std::sin(wShiftRad));
    for (size_t p = 0; p != PolarizationCount; ++p) {
      *dataIter *= phasor;
      ++dataIter;
    }
  }
}

template void MSGridderBase::rotateVisibilities<1>(
    const aocommon::BandData& bandData, double shiftFactor,
    std::complex<float>* dataIter);
template void MSGridderBase::rotateVisibilities<2>(
    const aocommon::BandData& bandData, double shiftFactor,
    std::complex<float>* dataIter);
template void MSGridderBase::rotateVisibilities<4>(
    const aocommon::BandData& bandData, double shiftFactor,
    std::complex<float>* dataIter);
