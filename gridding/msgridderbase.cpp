#include "msgridderbase.h"

#include "../io/logger.h"

#include "../math/calculatefftsize.h"

#include "../msproviders/msprovider.h"

#include "../structures/imageweights.h"

#include "../units/angle.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MCPosition.h>
#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>

#include <casacore/tables/Tables/ArrColDesc.h>
#include <casacore/tables/Tables/TableRecord.h>

#include <atomic>

MSGridderBase::MSData::MSData()
    : msIndex(0), matchingRows(0), totalRowsProcessed(0) {}

MSGridderBase::MSData::~MSData() {}

MSGridderBase::MSGridderBase()
    : MeasurementSetGridder(),
      _theoreticalBeamSize(0.0),
      _actualInversionWidth(0),
      _actualInversionHeight(0),
      _actualPixelSizeX(0),
      _actualPixelSizeY(0),
      _metaDataCache(nullptr),
      _hasFrequencies(false),
      _freqHigh(0.0),
      _freqLow(0.0),
      _bandStart(0.0),
      _bandEnd(0.0),
      _startTime(0.0),
      _griddedVisibilityCount(0),
      _totalWeight(0.0),
      _maxGriddedWeight(0.0),
      _visibilityWeightSum(0.0) {}

MSGridderBase::~MSGridderBase() {}

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

void MSGridderBase::initializeBandData(casacore::MeasurementSet& ms,
                                       MSGridderBase::MSData& msData) {
  msData.bandData = MultiBandData(ms.spectralWindow(), ms.dataDescription());
  if (Selection(msData.msIndex).HasChannelRange()) {
    msData.startChannel = Selection(msData.msIndex).ChannelRangeStart();
    msData.endChannel = Selection(msData.msIndex).ChannelRangeEnd();
    Logger::Debug << "Selected channels: " << msData.startChannel << '-'
                  << msData.endChannel << '\n';
    if (msData.startChannel >= msData.bandData.MaxChannels() ||
        msData.endChannel > msData.bandData.MaxChannels() ||
        msData.startChannel == msData.endChannel) {
      std::ostringstream str;
      str << "An invalid channel range was specified! Measurement set only has "
          << msData.bandData.MaxChannels()
          << " channels, requested imaging range is " << msData.startChannel
          << " -- " << msData.endChannel << '.';
      throw std::runtime_error(str.str());
    }
  } else {
    msData.startChannel = 0;
    msData.endChannel = msData.bandData.MaxChannels();
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
  MultiBandData selectedBand = msData.SelectedBand();
  std::vector<float> weightArray(selectedBand.MaxChannels() * NPolInMSProvider);
  msData.msProvider->Reset();
  double curTimestep = -1, firstTime = -1, lastTime = -1;
  size_t nTimesteps = 0;
  while (msData.msProvider->CurrentRowAvailable()) {
    MSProvider::MetaData metaData;
    msData.msProvider->ReadMeta(metaData);

    if (curTimestep != metaData.time) {
      curTimestep = metaData.time;
      ++nTimesteps;
      if (firstTime == -1) firstTime = curTimestep;
      lastTime = curTimestep;
    }

    const BandData& curBand = selectedBand[metaData.dataDescId];
    double wHi = fabs(metaData.wInM / curBand.SmallestWavelength());
    double wLo = fabs(metaData.wInM / curBand.LongestWavelength());
    double baselineInM =
        sqrt(metaData.uInM * metaData.uInM + metaData.vInM * metaData.vInM +
             metaData.wInM * metaData.wInM);
    double halfWidth = 0.5 * ImageWidth(), halfHeight = 0.5 * ImageHeight();
    if (wHi > msData.maxW || wLo < msData.minW ||
        baselineInM / curBand.SmallestWavelength() > msData.maxBaselineUVW) {
      msData.msProvider->ReadWeights(weightArray.data());
      const float* weightPtr = weightArray.data();
      for (size_t ch = 0; ch != curBand.ChannelCount(); ++ch) {
        const double wavelength = curBand.ChannelWavelength(ch);
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

    msData.msProvider->NextRow();
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

  bool hasCache = !_metaDataCache->msDataVector.empty();
  if (!hasCache) _metaDataCache->msDataVector.resize(MeasurementSetCount());
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
                 << Angle::ToNiceString(_theoreticalBeamSize) << "\n";
  }
  if (HasWLimit()) {
    _maxW *= (1.0 - WLimit());
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
      size_t newWidth =
          std::max(std::min(optWidth, _actualInversionWidth), size_t(32));
      size_t newHeight =
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

  if (IsFirstIteration() || !HasWGridSize()) {
    size_t suggestedGridSize = getSuggestedWGridSize();
    if (!HasWGridSize())
      SetActualWGridSize(suggestedGridSize);
    else
      SetActualWGridSize(WGridSize());
  } else
    SetActualWGridSize(WGridSize());
}

template <size_t PolarizationCount>
void MSGridderBase::readAndWeightVisibilities(MSProvider& msProvider,
                                              InversionRow& rowData,
                                              const BandData& curBand,
                                              float* weightBuffer,
                                              std::complex<float>* modelBuffer,
                                              const bool* isSelected) {
  const std::size_t dataSize = curBand.ChannelCount() * PolarizationCount;
  if (DoImagePSF()) {
    std::fill_n(rowData.data, dataSize, 1.0);
    if (HasDenormalPhaseCentre()) {
      const double lmsqrt = std::sqrt(1.0 - PhaseCentreDL() * PhaseCentreDL() -
                                      PhaseCentreDM() * PhaseCentreDM());
      const double shiftFactor = 2.0 * M_PI *
                                 ((rowData.uvw[0] * PhaseCentreDL() +
                                   rowData.uvw[1] * PhaseCentreDM()) +
                                  rowData.uvw[2] * (lmsqrt - 1.0));
      rotateVisibilities<PolarizationCount>(curBand, shiftFactor, rowData.data);
    }
  } else {
    msProvider.ReadData(rowData.data);
  }
  rowData.rowId = msProvider.RowId();

  if (DoSubtractModel()) {
    msProvider.ReadModel(modelBuffer);
    std::complex<float>* modelIter = modelBuffer;
    for (std::complex<float>* iter = rowData.data;
         iter != rowData.data + dataSize; ++iter) {
      *iter -= *modelIter;
      modelIter++;
    }
  }

  msProvider.ReadWeights(weightBuffer);

  // Any visibilities that are not gridded in this pass
  // should not contribute to the weight sum, so set these
  // to have zero weight.
  for (size_t ch = 0; ch != dataSize; ++ch) {
    if (!isSelected[ch]) weightBuffer[ch] = 0.0;
  }

  switch (VisibilityWeightingMode()) {
    case NormalVisibilityWeighting:
      // The weight buffer already contains the visibility weights: do nothing
      break;
    case SquaredVisibilityWeighting:
      // Square the visibility weights
      for (size_t chp = 0; chp != dataSize; ++chp)
        weightBuffer[chp] *= weightBuffer[chp];
      break;
    case UnitVisibilityWeighting:
      // Set the visibility weights to one
      for (size_t chp = 0; chp != dataSize; ++chp) {
        if (weightBuffer[chp] != 0.0) weightBuffer[chp] = 1.0f;
      }
      break;
  }

  // Calculate imaging weights
  std::complex<float>* dataIter = rowData.data;
  float* weightIter = weightBuffer;
  _scratchWeights.resize(curBand.ChannelCount());
  for (size_t ch = 0; ch != curBand.ChannelCount(); ++ch) {
    double u = rowData.uvw[0] / curBand.ChannelWavelength(ch),
           v = rowData.uvw[1] / curBand.ChannelWavelength(ch),
           imageWeight = GetImageWeights()->GetWeight(u, v);
    _scratchWeights[ch] = imageWeight;

    for (size_t p = 0; p != PolarizationCount; ++p) {
      double cumWeight = *weightIter * imageWeight;
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
  if (StoreImagingWeights())
    msProvider.WriteImagingWeights(rowData.rowId, _scratchWeights.data());
}

template void MSGridderBase::readAndWeightVisibilities<1>(
    MSProvider& msProvider, InversionRow& newItem, const BandData& curBand,
    float* weightBuffer, std::complex<float>* modelBuffer,
    const bool* isSelected);

template void MSGridderBase::readAndWeightVisibilities<4>(
    MSProvider& msProvider, InversionRow& newItem, const BandData& curBand,
    float* weightBuffer, std::complex<float>* modelBuffer,
    const bool* isSelected);

template <size_t PolarizationCount>
void MSGridderBase::rotateVisibilities(const BandData& bandData,
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
    const BandData& bandData, double shiftFactor,
    std::complex<float>* dataIter);
template void MSGridderBase::rotateVisibilities<4>(
    const BandData& bandData, double shiftFactor,
    std::complex<float>* dataIter);
