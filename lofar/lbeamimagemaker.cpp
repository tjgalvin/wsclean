#include "lbeamimagemaker.h"

#ifndef HAVE_LOFAR_BEAM
#include <stdexcept>
PrimaryBeamImageSet LBeamImageMaker::Make() {
  throw std::runtime_error(
      "LOFAR beam imager maker called, but the software has been compiled "
      "without the LOFAR beam. Recompile your software and make sure that "
      "cmake finds the LOFAR station response library.");
}
#else

#include "../units/angle.h"
#include "../units/radeccoord.h"

#include "../fftresampler.h"
#include "../fitsreader.h"
#include "../fitswriter.h"
#include "../imageweights.h"

#include "../progressbar.h"

#include "../wsclean/imageweightcache.h"
#include "../wsclean/logger.h"
#include "../msproviders/msprovider.h"
#include "../multibanddata.h"

#include "lofarbeamkeywords.h"

#include <StationResponse/LofarMetaDataUtil.h>
#include <StationResponse/ITRFConverter.h>

#include <casacore/ms/MeasurementSets/MSField.h>

#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/tables/Tables/ScalarColumn.h>

#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MCPosition.h>

#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>
#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>

#include <aocommon/banddata.h>
#include <aocommon/imagecoordinates.h>
#include <aocommon/uvector.h>
#include <aocommon/matrix2x2.h>

#include <stdexcept>

using namespace LOFAR::StationResponse;
using aocommon::MC2x2;
using aocommon::MC4x4;

static void dirToITRFVector(const casacore::MDirection& dir,
                            ITRFConverter& convert, vector3r_t& itrf);

PrimaryBeamImageSet LBeamImageMaker::Make() {
  _sampledWidth = _width / _undersample;
  _sampledHeight = _width / _undersample;
  _sPixelSizeX = _pixelSizeX * _undersample;
  _sPixelSizeY = _pixelSizeX * _undersample;

  _totalWeightSum = 0.0;

  std::vector<HMC4x4> matrices(_sampledWidth * _sampledHeight, HMC4x4::Zero());

  Logger::Debug << "Making beam for " << _msProviders.size() << " parts ("
                << _tableEntry->msData.size() << " ms)\n";
  for (const MSProviderInfo& msProviderInfo : _msProviders) {
    const ImagingTableEntry::MSInfo& msInfo =
        _tableEntry->msData[msProviderInfo.msIndex];
    const MSSelection& selection = *msProviderInfo.selection;
    MultiBandData band;
    {
      SynchronizedMS ms = msProviderInfo.provider->MS();
      band = MultiBandData(ms->spectralWindow(), ms->dataDescription());
    }
    double centralFrequency = 0.0;
    for (size_t dataDescId = 0; dataDescId != band.DataDescCount();
         ++dataDescId) {
      BandData subBand(band[dataDescId], selection.ChannelRangeStart(),
                       selection.ChannelRangeEnd());
      centralFrequency += subBand.CentreFrequency();
    }
    centralFrequency /= msInfo.bands.size();
    makeBeamForMS(matrices, *msProviderInfo.provider, selection,
                  centralFrequency);
  }

  for (size_t j = 0; j != _sampledWidth * _sampledHeight; ++j) {
    matrices[j] *= 1.0 / _totalWeightSum;
  }

  PrimaryBeamImageSet beamImages(_width, _height, 16);
  beamImages.SetToZero();

  if (_width != _sampledWidth || _height != _sampledHeight) {
    FFTResampler resampler(_sampledWidth, _sampledHeight, _width, _height, 1);
    Image scratchA(_sampledWidth, _sampledHeight), scratchB(_width, _height);
    for (size_t p = 0; p != 16; ++p) {
      for (size_t i = 0; i != _sampledWidth * _sampledHeight; ++i)
        scratchA[i] = matrices[i].Data(p);
      resampler.Resample(scratchA.data(), scratchB.data());
      std::copy_n(scratchB.data(), _width * _height, &beamImages[p][0]);
    }
  } else {
    for (size_t p = 0; p != 16; ++p) {
      for (size_t i = 0; i != _sampledWidth * _sampledHeight; ++i)
        beamImages[p][i] = matrices[i].Data(p);
    }
  }

  return beamImages;
}

void LBeamImageMaker::makeBeamForMS(std::vector<HMC4x4>& matrices,
                                    MSProvider& msProvider,
                                    const MSSelection& selection,
                                    double centralFrequency) {
  /**
   * Read some meta data from the measurement set
   */
  SynchronizedMS ms = msProvider.MS();

  casacore::MSAntenna aTable = ms->antenna();
  if (aTable.nrow() == 0) throw std::runtime_error("No antennae in set");
  casacore::MPosition::ScalarColumn antPosColumn(
      aTable, aTable.columnName(casacore::MSAntennaEnums::POSITION));
  casacore::MPosition arrayPos = antPosColumn(0);

  Logger::Debug << "Making beam for frequency " << centralFrequency * 1e-6
                << " MHz.\n";

  casacore::MSField fieldTable(ms->field());
  casacore::ScalarMeasColumn<casacore::MDirection> delayDirColumn(
      fieldTable,
      casacore::MSField::columnName(casacore::MSFieldEnums::DELAY_DIR));
  if (fieldTable.nrow() != 1)
    throw std::runtime_error("Set has multiple fields");
  _delayDir = delayDirColumn(0);

  _useDifferentialBeam = LOFARBeamKeywords::GetPreappliedBeamDirection(
      *ms, msProvider.DataColumnName(), _useDifferentialBeam, _preappliedDir);

  if (fieldTable.tableDesc().isColumn("LOFAR_TILE_BEAM_DIR")) {
    casacore::ArrayMeasColumn<casacore::MDirection> tileBeamDirColumn(
        fieldTable, "LOFAR_TILE_BEAM_DIR");
    _tileBeamDir = *(tileBeamDirColumn(0).data());
  } else {
    _tileBeamDir = _delayDir;
  }

  Logger::Debug << "Using delay direction: "
                << RaDecCoord::RaDecToString(_delayDir.getAngle().getValue()[0],
                                             _delayDir.getAngle().getValue()[1])
                << '\n';
  Logger::Debug << "Using tile direction: "
                << RaDecCoord::RaDecToString(
                       _tileBeamDir.getAngle().getValue()[0],
                       _tileBeamDir.getAngle().getValue()[1])
                << '\n';

  std::vector<Station::Ptr> stations(aTable.nrow());
  readStations(*ms, stations.begin());

  Logger::Debug << "Counting timesteps...\n";
  msProvider.Reset();
  size_t timestepCount = 0;
  double startTime = 0.0, endTime = 0.0;
  if (msProvider.CurrentRowAvailable()) {
    MSProvider::MetaData meta;
    msProvider.ReadMeta(meta);
    startTime = meta.time;
    endTime = meta.time;
    ++timestepCount;
    msProvider.NextRow();
    while (msProvider.CurrentRowAvailable()) {
      msProvider.ReadMeta(meta);
      if (endTime != meta.time) {
        ++timestepCount;
        endTime = meta.time;
      }
      msProvider.NextRow();
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

  casacore::MEpoch::ScalarColumn timeColumn(
      *ms, ms->columnName(casacore::MSMainEnums::TIME));
  casacore::MEpoch midTime(
      casacore::MVEpoch((0.5 / 86400.0) * (startTime + endTime)),
      timeColumn(0).getRef());
  Logger::Debug << "Mid time for full selection: " << midTime << '\n';
  casacore::MeasFrame midFrame(arrayPos, midTime);
  const casacore::MDirection::Ref hadecRef(casacore::MDirection::HADEC,
                                           midFrame);
  const casacore::MDirection::Ref azelgeoRef(casacore::MDirection::AZELGEO,
                                             midFrame);
  const casacore::MDirection::Ref midJ2000Ref(casacore::MDirection::J2000,
                                              midFrame);

  double refIntervalWeight = 0.0;
  msProvider.Reset();
  for (size_t intervalIndex = 0; intervalIndex != intervalCount;
       ++intervalIndex) {
    // Find the mid time step
    double firstTime =
        startTime + (endTime - startTime) * intervalIndex / intervalCount;
    double lastTime =
        startTime + (endTime - startTime) * (intervalIndex + 1) / intervalCount;
    casacore::MEpoch timeEpoch = casacore::MEpoch(
        casacore::MVEpoch((0.5 / 86400.0) * (firstTime + lastTime)),
        timeColumn(0).getRef());
    Logger::Debug << "Mid time for this interval: " << timeEpoch << '\n';

    casacore::MeasFrame frame(arrayPos, timeEpoch);
    const casacore::MDirection::Ref j2000Ref(casacore::MDirection::J2000,
                                             frame);

    if (_useDifferentialBeam)
      Logger::Debug << "Making differential snapshot beam for " << timeEpoch
                    << "\n";
    else
      Logger::Debug << "Making snapshot beam for " << timeEpoch << "\n";

    double intervalWeight = 0.0;
    WeightMatrix baselineWeights(stations.size());
    calculateStationWeights(*_imageWeights, intervalWeight, baselineWeights, ms,
                            msProvider, selection, lastTime);

    if (refIntervalWeight == 0.0) refIntervalWeight = intervalWeight;
    Logger::Debug << "Relative weight of this interval: "
                  << ((intervalWeight == 0.0)
                          ? 0.0
                          : (intervalWeight * 100.0 / refIntervalWeight))
                  << " %\n";

    if (intervalWeight != 0.0) {
      std::vector<HMC4x4> singleMatrices(_sampledWidth * _sampledHeight);

      makeBeamSnapshot(stations, baselineWeights, singleMatrices.data(),
                       timeEpoch.getValue().get() * 86400.0, centralFrequency,
                       centralFrequency, frame);

      _totalWeightSum += intervalWeight;
      for (size_t j = 0; j != _sampledWidth * _sampledHeight; ++j) {
        matrices[j] += singleMatrices[j] * intervalWeight;
      }
      if (_saveIntermediateImages) {
        static size_t index = 0;
        std::ostringstream str;
        str << "beam-intermediate-" << index << ".fits";
        aocommon::UVector<double> img(_sampledWidth * _sampledHeight);
        for (size_t px = 0; px != _sampledWidth * _sampledHeight; ++px) {
          double norm = singleMatrices[px].Norm();
          img[px] = norm;
          if (px == _sampledWidth / 2 + (_sampledHeight / 2) * _sampledWidth)
            Logger::Debug << "Central pixel: " << img[px] << '\n';
        }
        FitsWriter writer;
        writer.SetImageDimensions(_sampledWidth, _sampledHeight);
        writer.Write(str.str(), img.data());
        ++index;
      }
    }
  }
}

static void dirToITRFVector(const casacore::MDirection& dir,
                            ITRFConverter& convert, vector3r_t& itrf) {
  casacore::MDirection itrfDir = convert.toDirection(dir);
  casacore::Vector<double> itrfVal = itrfDir.getValue().getValue();
  itrf[0] = itrfVal[0];
  itrf[1] = itrfVal[1];
  itrf[2] = itrfVal[2];
}

void LBeamImageMaker::makeBeamSnapshot(
    const std::vector<Station::Ptr>& stations, const WeightMatrix& weights,
    HMC4x4* matrices, double time, double frequency, double subbandFrequency,
    const casacore::MeasFrame& frame) {
  static const casacore::Unit radUnit("rad");
  const casacore::MDirection::Ref j2000Ref(casacore::MDirection::J2000, frame);
  ITRFConverter converter(time);

  vector3r_t station0, tile0, diffBeamCentre;
  dirToITRFVector(_delayDir, converter, station0);
  dirToITRFVector(_tileBeamDir, converter, tile0);

  Logger::Debug << "Time=" << time << '\n';

  std::vector<MC2x2> inverseCentralGain;
  if (_useDifferentialBeam) {
    dirToITRFVector(_preappliedDir, converter, diffBeamCentre);
    inverseCentralGain.resize(stations.size());
    for (size_t a = 0; a != stations.size(); ++a) {
      matrix22c_t gainMatrix = stations[a]->response(
          time, frequency, diffBeamCentre, subbandFrequency, station0, tile0);
      inverseCentralGain[a][0] = gainMatrix[0][0];
      inverseCentralGain[a][1] = gainMatrix[0][1];
      inverseCentralGain[a][2] = gainMatrix[1][0];
      inverseCentralGain[a][3] = gainMatrix[1][1];
      if (!inverseCentralGain[a].Invert()) {
        inverseCentralGain[a] = MC2x2::NaN();
      }
    }
  }

  ProgressBar progressBar("Constructing beam");
  for (size_t y = 0; y != _sampledHeight; ++y) {
    for (size_t x = 0; x != _sampledWidth; ++x) {
      double l, m, ra, dec;
      aocommon::ImageCoordinates::XYToLM(x, y, _sPixelSizeX, _sPixelSizeY,
                                         _sampledWidth, _sampledHeight, l, m);
      l += _phaseCentreDL;
      m += _phaseCentreDM;
      aocommon::ImageCoordinates::LMToRaDec(l, m, _phaseCentreRA,
                                            _phaseCentreDec, ra, dec);

      casacore::MDirection imageDir(
          casacore::MVDirection(casacore::Quantity(ra, radUnit),
                                casacore::Quantity(dec, radUnit)),
          j2000Ref);

      vector3r_t itrfDirection;
      dirToITRFVector(imageDir, converter, itrfDirection);

      std::vector<MC2x2> stationGains(stations.size());
      for (size_t a = 0; a != stations.size(); ++a) {
        matrix22c_t gainMatrix = stations[a]->response(
            time, frequency, itrfDirection, subbandFrequency, station0, tile0);
        stationGains[a][0] = gainMatrix[0][0];
        stationGains[a][1] = gainMatrix[0][1];
        stationGains[a][2] = gainMatrix[1][0];
        stationGains[a][3] = gainMatrix[1][1];

        if (_useDifferentialBeam) {
          // The data have been premultiplied with the central beam C, and we
          // want to return a matrix that corrects the data for the full beam B.
          // Given our data:
          //    data = Ci^-1 vis Cj^-*
          // we want to multiple data with a differential beam matrix D such
          // that
          //    Di^-1 data Dj^-* = Bi^-1 vis Bj^-*
          // With B the full beam matrix. We can solve for D:
          //    Di^-1 Ci^-1 = Bi^-1  (and Cj^-* Dj^-* = Bj^-* )
          //          Di^-1 = Bi^-1 Ci
          //          Di    = Ci^-1 Bi
          // which means: diffBeam = inverseCentralGain * stationGain
          MC2x2 diffBeam;
          MC2x2::ATimesB(diffBeam, inverseCentralGain[a], stationGains[a]);
          stationGains[a] = diffBeam;
        }
      }

      HMC4x4 gain = HMC4x4::Zero();
      double totalWeight = 0.0;
      for (size_t a1 = 0; a1 != stations.size(); ++a1) {
        for (size_t a2 = a1; a2 != stations.size(); ++a2) {
          MC4x4 baselineBeam = MC4x4::KroneckerProduct(
                                   stationGains[a2].HermTranspose().Transpose(),
                                   stationGains[a1]) +
                               MC4x4::KroneckerProduct(
                                   stationGains[a1].HermTranspose().Transpose(),
                                   stationGains[a2]);
          double w = weights.Value(a1, a2);
          gain += HMC4x4(baselineBeam) * (0.5 * w);
          totalWeight += w;
        }
        matrices[y * _sampledWidth + x] = gain * (1.0 / totalWeight);
      }
    }
    progressBar.SetProgress(y, _sampledHeight);
  }
  progressBar.SetProgress(_sampledHeight, _sampledHeight);
}

void LBeamImageMaker::calculateStationWeights(
    const ImageWeights& imageWeights, double& totalWeight,
    WeightMatrix& baselineWeights, SynchronizedMS& ms, MSProvider& msProvider,
    const MSSelection& selection, double endTime) {
  casacore::MSAntenna antTable(ms->antenna());
  totalWeight = 0.0;
  aocommon::UVector<double> perAntennaWeights(antTable.nrow(), 0.0);

  MultiBandData multiband(ms->spectralWindow(), ms->dataDescription());
  size_t channelCount =
      selection.ChannelRangeEnd() - selection.ChannelRangeStart();
  size_t polarizationCount =
      (msProvider.Polarization() == aocommon::Polarization::Instrumental) ? 4
                                                                          : 1;
  aocommon::UVector<float> weightArr(channelCount * polarizationCount);

  while (msProvider.CurrentRowAvailable()) {
    MSProvider::MetaData metaData;
    msProvider.ReadMeta(metaData);
    if (metaData.time >= endTime) break;
    const BandData& band(multiband[metaData.dataDescId]);
    msProvider.ReadWeights(weightArr.data());

    for (size_t ch = 0; ch != channelCount; ++ch) {
      double u = metaData.uInM / band.ChannelWavelength(ch),
             v = metaData.vInM / band.ChannelWavelength(ch);
      double iw = imageWeights.GetWeight(u, v);
      double w = weightArr[ch * polarizationCount] * iw;
      totalWeight += w;
      perAntennaWeights[metaData.antenna1] += w;
      perAntennaWeights[metaData.antenna2] += w;
      baselineWeights.Value(metaData.antenna1, metaData.antenna2) += w;
    }
    msProvider.NextRow();
  }
  if (Logger::IsVerbose()) logWeights(*ms, perAntennaWeights);
}

void LBeamImageMaker::logWeights(casacore::MeasurementSet& ms,
                                 const aocommon::UVector<double>& weights) {
  casacore::MSAntenna antTable(ms.antenna());
  Logger::Debug << "Weights:";
  casacore::ScalarColumn<casacore::String> antNameCol(
      antTable, casacore::MSAntenna::columnName(casacore::MSAntenna::NAME));
  double maxWeight = 0.0;
  for (size_t a = 0; a != antTable.nrow(); ++a)
    maxWeight = std::max(weights[a], maxWeight);

  for (size_t a = 0; a != antTable.nrow(); ++a) {
    std::string name = antNameCol(a);
    double w =
        (maxWeight == 0.0) ? 0.0 : round(weights[a] * 1000.0 / maxWeight) * 0.1;
    Logger::Debug << ' ' << name << '=' << w;
  }
  Logger::Debug << '\n';
}

#endif
