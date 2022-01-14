#include "msgridderbase.h"

#include "../io/logger.h"

#include "../math/calculatefftsize.h"

#include "../msproviders/msprovider.h"
#include "../msproviders/msreaders/msreader.h"

#include "../structures/imageweights.h"

#include <aocommon/units/angle.h>

#ifdef HAVE_EVERYBEAM
#include <EveryBeam/load.h>
#include <EveryBeam/aterms/atermconfig.h>
#include <EveryBeam/pointresponse/phasedarraypoint.h>
#endif

// Only needed for EB/H5Parm related options
#include "../io/findmwacoefffile.h"
#include <limits>
#include <aocommon/matrix2x2.h>
#include <aocommon/multibanddata.h>

#include <schaapcommon/h5parm/h5parm.h>
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

#include <atomic>

using schaapcommon::h5parm::JonesParameters;

namespace {
/**
 * @brief Apply conjugated gains to the visibilities.
 *
 * @tparam PolarizationCount polarization count, 4 for IDG, 1 for all other
 * gridders.
 * @tparam GainEntry Which entry or entries from the gain matrices should be
 * taken into account when correcting the visibilities? See also the
 * documentation of DDGainMatrix.
 */
template <size_t PolarizationCount, DDGainMatrix GainEntry>
void ApplyConjugatedGain(std::complex<float>* visibilities,
                         const aocommon::MC2x2F& gain1,
                         const aocommon::MC2x2F& gain2);

template <>
void ApplyConjugatedGain<1, DDGainMatrix::kXX>(
    std::complex<float>* visibilities, const aocommon::MC2x2F& gain1,
    const aocommon::MC2x2F& gain2) {
  *visibilities = conj(gain1[0]) * (*visibilities) * gain2[0];
}

template <>
void ApplyConjugatedGain<1, DDGainMatrix::kYY>(
    std::complex<float>* visibilities, const aocommon::MC2x2F& gain1,
    const aocommon::MC2x2F& gain2) {
  *visibilities = conj(gain1[3]) * (*visibilities) * gain2[3];
}

template <>
void ApplyConjugatedGain<1, DDGainMatrix::kTrace>(
    std::complex<float>* visibilities, const aocommon::MC2x2F& gain1,
    const aocommon::MC2x2F& gain2) {
  // Stokes-I
  *visibilities *= 0.5f * gain2.DoubleDot(gain1.Conjugate());
}

template <>
void ApplyConjugatedGain<4, DDGainMatrix::kFull>(
    std::complex<float>* visibilities, const aocommon::MC2x2F& gain1,
    const aocommon::MC2x2F& gain2) {
  // All polarizations
  const aocommon::MC2x2F visibilities_mc2x2(visibilities);
  const aocommon::MC2x2F result =
      gain1.HermThenMultiply(visibilities_mc2x2).Multiply(gain2);
  result.AssignTo(visibilities);
}

/**
 * @brief Apply gains to the visibilities.
 *
 * @tparam PolarizationCount polarization count, 4 for IDG, 1 for all other
 * gridders.
 * @tparam GainEntry Which entry or entries from the gain matrices should be
 * taken into account when correcting the visibilities? See also the
 * documentation of DDGainMatrix.
 */
template <size_t PolarizationCount, DDGainMatrix GainEntry>
void ApplyGain(std::complex<float>* visibilities, const aocommon::MC2x2F& gain1,
               const aocommon::MC2x2F& gain2);

template <>
void ApplyGain<1, DDGainMatrix::kXX>(std::complex<float>* visibilities,
                                     const aocommon::MC2x2F& gain1,
                                     const aocommon::MC2x2F& gain2) {
  *visibilities = gain1[0] * (*visibilities) * conj(gain2[0]);
}

template <>
void ApplyGain<1, DDGainMatrix::kYY>(std::complex<float>* visibilities,
                                     const aocommon::MC2x2F& gain1,
                                     const aocommon::MC2x2F& gain2) {
  *visibilities = gain1[3] * (*visibilities) * conj(gain2[3]);
}

template <>
void ApplyGain<1, DDGainMatrix::kTrace>(std::complex<float>* visibilities,
                                        const aocommon::MC2x2F& gain1,
                                        const aocommon::MC2x2F& gain2) {
  // Stokes-I.
  *visibilities *= 0.5f * gain1.DoubleDot(gain2.Conjugate());
}

template <>
void ApplyGain<4, DDGainMatrix::kFull>(std::complex<float>* visibilities,
                                       const aocommon::MC2x2F& gain1,
                                       const aocommon::MC2x2F& gain2) {
  // All polarizations
  const aocommon::MC2x2F visibilities_mc2x2(visibilities);
  const aocommon::MC2x2F result =
      gain1.Multiply(visibilities_mc2x2).MultiplyHerm(gain2);
  result.AssignTo(visibilities);
}

/**
 * @brief Compute the gain based from the given gain matrices.
 *
 * @tparam GainEntry Which entry or entries from the gain matrices should be
 * taken into account? See DDGainMatrix for further documentation.
 */
template <DDGainMatrix GainEntry>
std::complex<float> ComputeGain(const aocommon::MC2x2F& gain1,
                                const aocommon::MC2x2F& gain2);

template <>
std::complex<float> ComputeGain<DDGainMatrix::kXX>(
    const aocommon::MC2x2F& gain1, const aocommon::MC2x2F& gain2) {
  return gain2[0] * conj(gain1[0]);
}

template <>
std::complex<float> ComputeGain<DDGainMatrix::kYY>(
    const aocommon::MC2x2F& gain1, const aocommon::MC2x2F& gain2) {
  return gain2[3] * conj(gain1[3]);
}

template <>
std::complex<float> ComputeGain<DDGainMatrix::kTrace>(
    const aocommon::MC2x2F& gain1, const aocommon::MC2x2F& gain2) {
  return 0.5f * gain2.DoubleDot(gain1.Conjugate());
}

template <>
std::complex<float> ComputeGain<DDGainMatrix::kFull>(
    const aocommon::MC2x2F& gain1, const aocommon::MC2x2F& gain2) {
  throw std::runtime_error("Not implemented!");
}

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
      _phaseCentreDL(0.0),
      _phaseCentreDM(0.0),
      _facetDirectionRA(0.0),
      _facetDirectionDec(0.0),
      _facetIndex(0),
      _facetGroupIndex(0),
      _msIndex(0),
      _additivePredict(false),
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
      _doImagePSF(false),
      _doSubtractModel(false),
      _smallInversion(settings.smallInversion),
      _wLimit(settings.wLimit / 100.0),
      _precalculatedWeightInfo(nullptr),
      _polarization(aocommon::Polarization::StokesI),
      _isComplex(false),
      _weighting(settings.weightMode),
      _isFirstIteration(false),
      _visibilityWeightingMode(settings.visibilityWeightingMode),
      _gridMode(GridMode::KaiserBesselKernel),
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
      _predictReader(nullptr),
#ifdef HAVE_EVERYBEAM
      _beamMode(everybeam::ParseBeamMode(settings.beamMode)),
      _beamNormalisationMode(everybeam::ParseBeamNormalisationMode(
          settings.beamNormalisationMode)),
#endif
      _cachedParmResponse(),
      _h5parms(),
      _h5SolTabs(),
      _correctType(),
      _cachedMSTimes(),
      _timeOffset() {
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

int64_t MSGridderBase::getAvailableMemory(double memFraction,
                                          double absMemLimit) {
  static std::atomic<bool> isFirst(true);
  bool printOutput = isFirst.exchange(false);
  long int pageCount = sysconf(_SC_PHYS_PAGES),
           pageSize = sysconf(_SC_PAGE_SIZE);
  int64_t memory = (int64_t)pageCount * (int64_t)pageSize;
  double memSizeInGB = (double)memory / (1024.0 * 1024.0 * 1024.0);
  if (memFraction == 1.0 && absMemLimit == 0.0 && printOutput) {
    Logger::Info << "Detected " << round(memSizeInGB * 10.0) / 10.0
                 << " GB of system memory, usage not limited.\n";
  } else {
    double limitInGB = memSizeInGB * memFraction;
    if (absMemLimit != 0.0 && limitInGB > absMemLimit) limitInGB = absMemLimit;
    if (printOutput) {
      Logger::Info << "Detected " << round(memSizeInGB * 10.0) / 10.0
                   << " GB of system memory, usage limited to "
                   << round(limitInGB * 10.0) / 10.0
                   << " GB (frac=" << round(memFraction * 1000.0) / 10.0
                   << "%, ";
      if (absMemLimit == 0.0)
        Logger::Info << "no limit)\n";
      else
        Logger::Info << "limit=" << round(absMemLimit * 10.0) / 10.0 << "GB)\n";
    }

    memory = int64_t((double)pageCount * (double)pageSize * memFraction);
    if (absMemLimit != 0.0 &&
        double(memory) > double(1024.0 * 1024.0 * 1024.0) * absMemLimit)
      memory = int64_t(double(absMemLimit) * double(1024.0 * 1024.0 * 1024.0));
  }
  return memory;
}

void MSGridderBase::initializePointResponse(
    const MSGridderBase::MSData& msData) {
  SynchronizedMS ms(msData.msProvider->MS());
#ifdef HAVE_EVERYBEAM
  if (_settings.applyFacetBeam && !_settings.facetRegionFilename.empty()) {
    // Hard-coded for now
    const bool frequency_interpolation = true;
    const bool use_channel_frequency = true;
    const std::string element_response_string =
        !_settings.beamModel.empty() ? _settings.beamModel : "DEFAULT";

    // Get path to coefficient file for MWA telescope
    everybeam::TelescopeType telescope_type = everybeam::GetTelescopeType(*ms);
    const std::string coeff_path =
        (telescope_type == everybeam::TelescopeType::kMWATelescope)
            ? wsclean::mwa::FindCoeffFile(_settings.mwaPath)
            : "";

    everybeam::ATermSettings aterm_settings;
    aterm_settings.coeff_path = coeff_path;
    aterm_settings.data_column_name = _settings.dataColumnName;

    everybeam::Options options =
        everybeam::aterms::ATermConfig::ConvertToEBOptions(
            *ms, aterm_settings, frequency_interpolation,
            _settings.beamNormalisationMode, use_channel_frequency,
            element_response_string, _settings.beamMode);

    _telescope = everybeam::Load(*ms, options);
    _pointResponse =
        _telescope->GetPointResponse(msData.msProvider->StartTime());
    _pointResponse->SetUpdateInterval(_settings.facetBeamUpdateTime);
    _cachedBeamResponse.resize(msData.bandData.ChannelCount() *
                               _pointResponse->GetAllStationsBufferSize());
  } else {
    if (_settings.applyFacetBeam) {
      throw std::runtime_error(
          "-apply-facet-beam was set, but no corresponding facet "
          "regions file was specified.");
    }
    _pointResponse = nullptr;
    _cachedBeamResponse.resize(0);
  }
#else
  if (_settings.applyFacetBeam && !_settings.facetRegionFilename.empty()) {
    throw std::runtime_error(
        "-apply-facet-beam was set, but wsclean was not compiled "
        "with EveryBeam. Please compile wsclean with EveryBeam to "
        "use the Facet Beam functionality");
  }
#endif
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
  _metaDataCache->beamSum = 0.0;
  _metaDataCache->h5Sum = 0.0;

  bool hasCache = !_metaDataCache->msDataVector.empty();
  if (!hasCache) _metaDataCache->msDataVector.resize(MeasurementSetCount());

  if (!DoImagePSF() && !_settings.facetSolutionFiles.empty()) {
    // Assign, rather than a resize here to make sure that
    // caches are re-initialized - even in the case an MSGridderBase
    // object would be re-used for a multiple gridding tasks.
    _cachedParmResponse.assign(MeasurementSetCount(), {});
    _cachedMSTimes.assign(MeasurementSetCount(), {});
    _timeOffset.assign(MeasurementSetCount(), 0u);
  }

  for (size_t i = 0; i != MeasurementSetCount(); ++i) {
    msDataVector[i].msIndex = i;
    initializeMeasurementSet(msDataVector[i], _metaDataCache->msDataVector[i],
                             hasCache);
    if (!DoImagePSF() && !_settings.facetSolutionFiles.empty()) {
      _cachedMSTimes[i] = SelectUniqueTimes(*msDataVector[i].msProvider);
    }
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
    Logger::Debug << "Set has denormal phase centre: dl=" << _phaseCentreDL
                  << ", dm=" << _phaseCentreDM << '\n';

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
    _h5parms.resize(MeasurementSetCount());
    _h5SolTabs.resize(MeasurementSetCount());
    _correctType.resize(MeasurementSetCount());

    for (size_t i = 0; i < _h5parms.size(); ++i) {
      if (_settings.facetSolutionFiles.size() > 1) {
        _h5parms[i].reset(
            new schaapcommon::h5parm::H5Parm(_settings.facetSolutionFiles[i]));
      } else {
        _h5parms[i].reset(
            new schaapcommon::h5parm::H5Parm(_settings.facetSolutionFiles[0]));
      }

      if (_settings.facetSolutionTables.size() == 1) {
        _h5SolTabs[i] = std::make_pair(
            &_h5parms[i]->GetSolTab(_settings.facetSolutionTables[0]), nullptr);
        _correctType[i] = JonesParameters::StringToCorrectType(
            _h5SolTabs[i].first->GetType());
      } else if (_settings.facetSolutionTables.size() == 2) {
        const std::array<std::string, 2> solTabTypes{
            _h5parms[i]->GetSolTab(_settings.facetSolutionTables[0]).GetType(),
            _h5parms[i]->GetSolTab(_settings.facetSolutionTables[1]).GetType()};

        auto itrA =
            std::find(solTabTypes.begin(), solTabTypes.end(), "amplitude");
        auto itrP = std::find(solTabTypes.begin(), solTabTypes.end(), "phase");

        if (itrA == solTabTypes.end() || itrP == solTabTypes.end()) {
          throw std::runtime_error(
              "WSClean expected solution tables with name 'amplitude' and "
              "'phase', but received " +
              solTabTypes[0] + " and " + solTabTypes[1]);
        } else {
          const size_t idxA = std::distance(solTabTypes.begin(), itrA);
          const size_t idxP = std::distance(solTabTypes.begin(), itrP);
          _h5SolTabs[i] = std::make_pair(
              &_h5parms[i]->GetSolTab(_settings.facetSolutionTables[idxA]),
              &_h5parms[i]->GetSolTab(_settings.facetSolutionTables[idxP]));
        }

        const size_t npol1 = _h5SolTabs[i].first->HasAxis("pol")
                                 ? _h5SolTabs[i].first->GetAxis("pol").size
                                 : 1;
        const size_t npol2 = _h5SolTabs[i].second->HasAxis("pol")
                                 ? _h5SolTabs[i].second->GetAxis("pol").size
                                 : 1;
        if (npol1 == 1 && npol2 == 1) {
          _correctType[i] = JonesParameters::CorrectType::SCALARGAIN;
        } else if (npol1 == 2 && npol2 == 2) {
          _correctType[i] = JonesParameters::CorrectType::GAIN;
        } else if (npol1 == 4 && npol2 == 4) {
          _correctType[i] = JonesParameters::CorrectType::FULLJONES;
        } else {
          throw std::runtime_error(
              "Incorrect or mismatching number of polarizations in the "
              "provided soltabs. Number of polarizations should be either "
              "all 1, 2 or 4, but received " +
              std::to_string(npol1) + " and " + std::to_string(npol2));
        }
      } else {
        throw std::runtime_error(
            "Specify the solution table name(s) with "
            "-soltab-names=soltabname1[OPTIONAL,soltabname2]");
      }
    }
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
  assert(!DoImagePSF());  // The PSF is never predicted.

  if (!_h5parms.empty()) {
    MSProvider::MetaData metaData;
    _predictReader->ReadMeta(metaData);
    // When the facet beam is applied, the row will be incremented later in this
    // function
    if (!_settings.applyFacetBeam) {
      _predictReader->NextInputRow();
    }

    const size_t nparms =
        (_correctType[_msIndex] == JonesParameters::CorrectType::FULLJONES) ? 4
                                                                            : 2;

    // Only extract DD solutions if the corresponding cache entry is empty.
    if (_cachedParmResponse[_msIndex].empty()) {
      const std::vector<double> freqs(curBand.begin(), curBand.end());
      const size_t responseSize = _cachedMSTimes[_msIndex].size() *
                                  freqs.size() * antennaNames.size() * nparms;
      const std::string dirName = _h5parms[_msIndex]->GetNearestSource(
          _facetDirectionRA, _facetDirectionDec);
      const size_t dirIndex = _h5SolTabs[_msIndex].first->GetDirIndex(dirName);
      JonesParameters jonesParameters(
          freqs, _cachedMSTimes[_msIndex], antennaNames, _correctType[_msIndex],
          JonesParameters::InterpolationType::NEAREST, dirIndex,
          _h5SolTabs[_msIndex].first, _h5SolTabs[_msIndex].second, false, 0.0f,
          0u, JonesParameters::MissingAntennaBehavior::kUnit);
      // parms (Casacore::Cube) is column major
      const auto parms = jonesParameters.GetParms();
      _cachedParmResponse[_msIndex].assign(&parms(0, 0, 0),
                                           &parms(0, 0, 0) + responseSize);
    }

    const size_t nchannels = curBand.ChannelCount();
    auto it =
        std::find(_cachedMSTimes[_msIndex].begin() + _timeOffset[_msIndex],
                  _cachedMSTimes[_msIndex].end(), metaData.time);
    if (it != _cachedMSTimes[_msIndex].end()) {
      _timeOffset[_msIndex] =
          std::distance(_cachedMSTimes[_msIndex].begin(), it);
    } else {
      throw std::runtime_error(
          "Time not found in cached times. A potential reason could be that "
          "the "
          "time values in the provided MS are not in ascending order.");
    }

    std::complex<float>* iter = buffer;
    if (nparms == 2) {
      for (size_t ch = 0; ch < nchannels; ++ch) {
        // Column major indexing
        const size_t offset = (_timeOffset[_msIndex] * nchannels + ch) *
                              antennaNames.size() * nparms;
        const size_t offset1 = offset + metaData.antenna1 * nparms;
        const size_t offset2 = offset + metaData.antenna2 * nparms;
        const aocommon::MC2x2F gain1(
            _cachedParmResponse[_msIndex][offset1], 0, 0,
            _cachedParmResponse[_msIndex][offset1 + 1]);
        const aocommon::MC2x2F gain2(
            _cachedParmResponse[_msIndex][offset2], 0, 0,
            _cachedParmResponse[_msIndex][offset2 + 1]);
        ApplyGain<PolarizationCount, GainEntry>(iter, gain1, gain2);
        iter += PolarizationCount;
      }
    } else {
      for (size_t ch = 0; ch < nchannels; ++ch) {
        // Column major indexing
        const size_t offset = (_timeOffset[_msIndex] * nchannels + ch) *
                              antennaNames.size() * nparms;
        const size_t offset1 = offset + metaData.antenna1 * nparms;
        const size_t offset2 = offset + metaData.antenna2 * nparms;
        const aocommon::MC2x2F gain1(&_cachedParmResponse[_msIndex][offset1]);
        const aocommon::MC2x2F gain2(&_cachedParmResponse[_msIndex][offset2]);
        ApplyGain<PolarizationCount, GainEntry>(iter, gain1, gain2);
        iter += PolarizationCount;
      }
    }
  }

#ifdef HAVE_EVERYBEAM
  if (_settings.applyFacetBeam && !_settings.facetRegionFilename.empty()) {
    MSProvider::MetaData metaData;
    _predictReader->ReadMeta(metaData);
    _predictReader->NextInputRow();
    _pointResponse->UpdateTime(metaData.time);
    if (_pointResponse->HasTimeUpdate()) {
      for (size_t ch = 0; ch < curBand.ChannelCount(); ++ch) {
        _pointResponse->ResponseAllStations(
            _beamMode,
            &_cachedBeamResponse[ch *
                                 _pointResponse->GetAllStationsBufferSize()],
            _facetDirectionRA, _facetDirectionDec, curBand.ChannelFrequency(ch),
            metaData.fieldId);
      }
    }

    std::complex<float>* iter = buffer;
    for (size_t ch = 0; ch < curBand.ChannelCount(); ++ch) {
      const size_t offset = ch * _pointResponse->GetAllStationsBufferSize();
      const size_t offset1 = offset + metaData.antenna1 * 4u;
      const size_t offset2 = offset + metaData.antenna2 * 4u;

      const aocommon::MC2x2F gain1(&_cachedBeamResponse[offset1]);
      const aocommon::MC2x2F gain2(&_cachedBeamResponse[offset2]);
      ApplyGain<PolarizationCount, GainEntry>(iter, gain1, gain2);
      iter += PolarizationCount;
    }
  }
#endif

  {
    WriterLockManager::LockGuard guard = _writerLockManager->GetLock(
        _facetGroupIndex * MeasurementSetCount() + _msIndex);
    msProvider.WriteModel(buffer, _additivePredict);
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

template void MSGridderBase::writeVisibilities<4, DDGainMatrix::kFull>(
    MSProvider& msProvider, const std::vector<std::string>& antennaNames,
    const aocommon::BandData& curBand, std::complex<float>* buffer);

#ifdef HAVE_EVERYBEAM
template <size_t PolarizationCount, DDGainMatrix GainEntry>
void MSGridderBase::ApplyConjugatedFacetBeam(MSReader& msReader,
                                             InversionRow& rowData,
                                             const aocommon::BandData& curBand,
                                             float* weightBuffer) {
  MSProvider::MetaData metaData;
  msReader.ReadMeta(metaData);

  _pointResponse->UpdateTime(metaData.time);
  if (_pointResponse->HasTimeUpdate()) {
    for (size_t ch = 0; ch < curBand.ChannelCount(); ++ch) {
      _pointResponse->ResponseAllStations(
          _beamMode,
          &_cachedBeamResponse[ch * _pointResponse->GetAllStationsBufferSize()],
          _facetDirectionRA, _facetDirectionDec, curBand.ChannelFrequency(ch),
          metaData.fieldId);
    }
  }

  // rowData.data contains the visibilities
  std::complex<float>* iter = rowData.data;
  float* weightIter = weightBuffer;
  for (size_t ch = 0; ch < curBand.ChannelCount(); ++ch) {
    const size_t offset = ch * _pointResponse->GetAllStationsBufferSize();
    const size_t offset1 = offset + metaData.antenna1 * 4u;
    const size_t offset2 = offset + metaData.antenna2 * 4u;

    const aocommon::MC2x2F gain1(&_cachedBeamResponse[offset1]);
    const aocommon::MC2x2F gain2(&_cachedBeamResponse[offset2]);
    ApplyConjugatedGain<PolarizationCount, GainEntry>(iter, gain1, gain2);
    const std::complex<float> g = ComputeGain<GainEntry>(gain1, gain2);

    const float weight = *weightIter * _scratchWeights[ch];
    _metaDataCache->beamSum += (conj(g) * weight * g).real();

    // Only admissible PolarizationCount for applying the facet beam is 1.
    iter += PolarizationCount;
    weightIter += PolarizationCount;
  }
}
#endif

template <size_t PolarizationCount, DDGainMatrix GainEntry>
void MSGridderBase::ApplyConjugatedH5Parm(
    MSReader& msReader, const std::vector<std::string>& antennaNames,
    InversionRow& rowData, const aocommon::BandData& curBand,
    float* weightBuffer) {
  MSProvider::MetaData metaData;
  msReader.ReadMeta(metaData);

  const size_t nparms =
      (_correctType[_msIndex] == JonesParameters::CorrectType::FULLJONES) ? 4
                                                                          : 2;

  // Only extract DD solutions if the corresponding cache entry is empty.
  if (_cachedParmResponse[_msIndex].empty()) {
    const std::vector<double> freqs(curBand.begin(), curBand.end());
    const size_t responseSize = _cachedMSTimes[_msIndex].size() * freqs.size() *
                                antennaNames.size() * nparms;
    const std::string dirName = _h5parms[_msIndex]->GetNearestSource(
        _facetDirectionRA, _facetDirectionDec);
    const size_t dirIndex = _h5SolTabs[_msIndex].first->GetDirIndex(dirName);
    JonesParameters jonesParameters(
        freqs, _cachedMSTimes[_msIndex], antennaNames, _correctType[_msIndex],
        JonesParameters::InterpolationType::NEAREST, dirIndex,
        _h5SolTabs[_msIndex].first, _h5SolTabs[_msIndex].second, false, 0.0f,
        0u, JonesParameters::MissingAntennaBehavior::kUnit);
    // parms (Casacore::Cube) is column major
    const auto parms = jonesParameters.GetParms();
    _cachedParmResponse[_msIndex].assign(&parms(0, 0, 0),
                                         &parms(0, 0, 0) + responseSize);
  }

  const size_t nchannels = curBand.ChannelCount();
  auto it = std::find(_cachedMSTimes[_msIndex].begin() + _timeOffset[_msIndex],
                      _cachedMSTimes[_msIndex].end(), metaData.time);
  if (it != _cachedMSTimes[_msIndex].end()) {
    // Update _timeOffset value with index
    _timeOffset[_msIndex] = std::distance(_cachedMSTimes[_msIndex].begin(), it);
  } else {
    throw std::runtime_error(
        "Time not found in cached times. A potential reason could be that the "
        "time values in the provided MS are not in ascending order.");
  }

  // Conditional could be templated once C++ supports partial function
  // specialization
  std::complex<float>* iter = rowData.data;
  float* weightIter = weightBuffer;
  if (nparms == 2) {
    for (size_t ch = 0; ch < nchannels; ++ch) {
      // Column major indexing
      const size_t offset = (_timeOffset[_msIndex] * nchannels + ch) *
                            antennaNames.size() * nparms;
      const size_t offset1 = offset + metaData.antenna1 * nparms;
      const size_t offset2 = offset + metaData.antenna2 * nparms;
      const aocommon::MC2x2F gain1(_cachedParmResponse[_msIndex][offset1], 0, 0,
                                   _cachedParmResponse[_msIndex][offset1 + 1]);
      const aocommon::MC2x2F gain2(_cachedParmResponse[_msIndex][offset2], 0, 0,
                                   _cachedParmResponse[_msIndex][offset2 + 1]);
      ApplyConjugatedGain<PolarizationCount, GainEntry>(iter, gain1, gain2);
      const std::complex<float> g = ComputeGain<GainEntry>(gain1, gain2);

      const float weight = *weightIter * _scratchWeights[ch];
      _metaDataCache->h5Sum += (conj(g) * weight * g).real();

      // Only admissible PolarizationCount for applying gains from solution
      // file is 1.
      iter += PolarizationCount;
      weightIter += PolarizationCount;
    }
  } else {
    for (size_t ch = 0; ch < nchannels; ++ch) {
      // Column major indexing
      const size_t offset = (_timeOffset[_msIndex] * nchannels + ch) *
                            antennaNames.size() * nparms;
      const size_t offset1 = offset + metaData.antenna1 * nparms;
      const size_t offset2 = offset + metaData.antenna2 * nparms;
      const aocommon::MC2x2F gain1(&_cachedParmResponse[_msIndex][offset1]);
      const aocommon::MC2x2F gain2(&_cachedParmResponse[_msIndex][offset2]);
      ApplyConjugatedGain<PolarizationCount, GainEntry>(iter, gain1, gain2);
      const std::complex<float> g = ComputeGain<GainEntry>(gain1, gain2);

      const float weight = *weightIter * _scratchWeights[ch];
      _metaDataCache->h5Sum += (conj(g) * weight * g).real();

      // Only admissible PolarizationCount for applying gains from solution
      // file is 1.
      iter += PolarizationCount;
      weightIter += PolarizationCount;
    }
  }
}

template <size_t PolarizationCount, DDGainMatrix GainEntry>
void MSGridderBase::readAndWeightVisibilities(
    MSReader& msReader, const std::vector<std::string>& antennaNames,
    InversionRow& rowData, const aocommon::BandData& curBand,
    float* weightBuffer, std::complex<float>* modelBuffer,
    const bool* isSelected) {
  const std::size_t dataSize = curBand.ChannelCount() * PolarizationCount;
  if (DoImagePSF()) {
    std::fill_n(rowData.data, dataSize, 1.0);
    if (HasDenormalPhaseCentre() && _settings.facetRegionFilename.empty()) {
      const double lmsqrt = std::sqrt(1.0 - PhaseCentreDL() * PhaseCentreDL() -
                                      PhaseCentreDM() * PhaseCentreDM());
      const double shiftFactor = 2.0 * M_PI *
                                 ((rowData.uvw[0] * PhaseCentreDL() +
                                   rowData.uvw[1] * PhaseCentreDM()) +
                                  rowData.uvw[2] * (lmsqrt - 1.0));
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
  for (size_t ch = 0; ch != dataSize; ++ch) {
    if (!isSelected[ch]) weightBuffer[ch] = 0.0;
  }

  switch (GetVisibilityWeightingMode()) {
    case VisibilityWeightingMode::NormalVisibilityWeighting:
      // The weight buffer already contains the visibility weights: do nothing
      break;
    case VisibilityWeightingMode::SquaredVisibilityWeighting:
      // Square the visibility weights
      for (size_t chp = 0; chp != dataSize; ++chp)
        weightBuffer[chp] *= weightBuffer[chp];
      break;
    case VisibilityWeightingMode::UnitVisibilityWeighting:
      // Set the visibility weights to one
      for (size_t chp = 0; chp != dataSize; ++chp) {
        if (weightBuffer[chp] != 0.0) weightBuffer[chp] = 1.0f;
      }
      break;
  }

  // Precompute imaging weights
  _scratchWeights.resize(curBand.ChannelCount());
  for (size_t ch = 0; ch != curBand.ChannelCount(); ++ch) {
    const double u = rowData.uvw[0] / curBand.ChannelWavelength(ch);
    const double v = rowData.uvw[1] / curBand.ChannelWavelength(ch);
    _scratchWeights[ch] = GetImageWeights()->GetWeight(u, v);
  }
  if (StoreImagingWeights())
    msReader.WriteImagingWeights(_scratchWeights.data());

  if (!DoImagePSF()) {
#ifdef HAVE_EVERYBEAM
    if (_settings.applyFacetBeam && !_settings.facetRegionFilename.empty()) {
      ApplyConjugatedFacetBeam<PolarizationCount, GainEntry>(
          msReader, rowData, curBand, weightBuffer);
    }
#endif

    if (!_h5parms.empty()) {
      ApplyConjugatedH5Parm<PolarizationCount, GainEntry>(
          msReader, antennaNames, rowData, curBand, weightBuffer);
    }
  }

  // Calculate imaging weights
  std::complex<float>* dataIter = rowData.data;
  float* weightIter = weightBuffer;
  for (size_t ch = 0; ch != curBand.ChannelCount(); ++ch) {
    for (size_t p = 0; p != PolarizationCount; ++p) {
      const double cumWeight = *weightIter * _scratchWeights[ch];
      if (p == 0 && cumWeight != 0.0) {
        // Visibility weight sum is the sum of weights excluding imaging weights
        _visibilityWeightSum += *weightIter;
        _maxGriddedWeight = std::max(cumWeight, _maxGriddedWeight);
        ++_griddedVisibilityCount;
        // Total weight includes imaging weights
        _totalWeight += cumWeight;
      }
      *weightIter = cumWeight;
      *dataIter *= *weightIter;
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
template void MSGridderBase::rotateVisibilities<4>(
    const aocommon::BandData& bandData, double shiftFactor,
    std::complex<float>* dataIter);
