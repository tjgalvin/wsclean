#ifndef ATERM_CONFIG_H
#define ATERM_CONFIG_H

#include <cassert>
#include <string>

#include "atermbase.h"
#include "atermbeam.h"
#include "atermstub.h"
#include "dishaterm.h"
#include "dldmaterm.h"
#include "fitsaterm.h"
#include "lofarbeamterm.h"
#include "mwabeamterm.h"
#include "pafbeamterm.h"
#include "telescope.h"

#include <aocommon/matrix2x2.h>

#include "../parsetreader.h"

#include "../wsclean/wscleansettings.h"
#include "../units/radeccoord.h"

using aocommon::Matrix2x2;

class ATermConfig : public ATermBase {
 public:
  ATermConfig(size_t nAntenna, const CoordinateSystem& coordinateSystem,
              const WSCleanSettings& settings)
      : _nAntenna(nAntenna),
        _coordinateSystem(coordinateSystem),
        _settings(settings) {}

  void Read(const std::string& msFilename, casacore::MeasurementSet ms,
            const std::string& parset) {
    auto iter = std::find(_settings.filenames.begin(),
                          _settings.filenames.end(), msFilename);
    assert(iter != _settings.filenames.end());
    size_t filenameIndex = iter - _settings.filenames.begin();

    ParsetReader reader(parset);
    std::vector<std::string> aterms = reader.GetStringList("aterms");
    if (aterms.empty())
      throw std::runtime_error(
          "No a-term correction given in parset (aterms key is an empty list)");

    for (const std::string atermName : aterms) {
      std::string atermType =
          reader.GetStringOr(atermName + ".type", atermName);
      if (atermType == "tec") {
        std::vector<std::string> tecFiles =
            reader.GetStringList(atermName + ".images");
        std::unique_ptr<FitsATerm> f(
            new FitsATerm(_nAntenna, _coordinateSystem));
        f->OpenTECFiles(tecFiles);
        std::string windowStr =
            reader.GetStringOr(atermName + ".window", "raised-hann");
        WindowFunction::Type window = WindowFunction::GetType(windowStr);
        if (window == WindowFunction::Tukey)
          f->SetTukeyWindow(double(_settings.paddedImageWidth) /
                            _settings.trimmedImageWidth);
        else
          f->SetWindow(window);
        f->SetDownSample(reader.GetBoolOr(atermName + ".downsample", true));
        _aterms.emplace_back(std::move(f));
      } else if (atermType == "diagonal") {
        std::vector<std::string> diagFiles =
            reader.GetStringList(atermName + ".images");
        std::unique_ptr<FitsATerm> f(
            new FitsATerm(_nAntenna, _coordinateSystem));
        f->OpenDiagGainFiles(diagFiles);
        std::string windowStr =
            reader.GetStringOr(atermName + ".window", "raised-hann");
        WindowFunction::Type window = WindowFunction::GetType(windowStr);
        if (window == WindowFunction::Tukey)
          f->SetTukeyWindow(double(_settings.paddedImageWidth) /
                            _settings.trimmedImageWidth);
        else
          f->SetWindow(window);
        f->SetDownSample(reader.GetBoolOr(atermName + ".downsample", true));
        _aterms.emplace_back(std::move(f));
      } else if (atermType == "dldm") {
        std::vector<std::string> dldmFiles =
            reader.GetStringList(atermName + ".images");
        std::unique_ptr<DLDMATerm> f(
            new DLDMATerm(_nAntenna, _coordinateSystem));
        f->Open(dldmFiles);
        f->SetUpdateInterval(
            reader.GetDoubleOr("dldm.update_interval", 5.0 * 60.0));
        std::string windowStr =
            reader.GetStringOr(atermName + ".window", "raised-hann");
        WindowFunction::Type window = WindowFunction::GetType(windowStr);
        if (window == WindowFunction::Tukey)
          f->SetTukeyWindow(double(_settings.paddedImageWidth) /
                            _settings.trimmedImageWidth);
        else
          f->SetWindow(window);
        f->SetDownSample(reader.GetBoolOr(atermName + ".downsample", true));
        _aterms.emplace_back(std::move(f));
      } else if (atermType == "beam") {
        std::unique_ptr<ATermBeam> beam;
        switch (Telescope::GetType(ms)) {
          case Telescope::AARTFAAC:
          case Telescope::LOFAR: {
            bool differential = reader.GetBoolOr("beam.differential", false);
            bool useChannelFrequency =
                reader.GetBoolOr("beam.usechannelfreq", true);
            std::unique_ptr<LofarBeamTerm> lofarBeam(new LofarBeamTerm(
                ms, _coordinateSystem, _settings.dataColumnName));
            lofarBeam->SetUseDifferentialBeam(differential);
            lofarBeam->SetUseChannelFrequency(useChannelFrequency);
            beam = std::move(lofarBeam);
            break;
          }
          case Telescope::MWA: {
            std::unique_ptr<MWABeamTerm> mwaTerm(
                new MWABeamTerm(ms, _coordinateSystem));
            mwaTerm->SetSearchPath(_settings.mwaPath);
            beam = std::move(mwaTerm);
            break;
          }
          case Telescope::VLA: {
            beam.reset(new DishATerm(ms, _coordinateSystem));
            break;
          }
          default: {
            // This is here to make sure ATermStub compiles. This call should be
            // the same as the call for LofarBeamTerm(..)
            beam.reset(
                new ATermStub(ms, _coordinateSystem, _settings.dataColumnName));
            throw std::runtime_error("Can't make beam for this telescope");
          }
        }
        double updateInterval = reader.GetDoubleOr(
            "beam.update_interval", _settings.beamAtermUpdateTime);
        beam->SetUpdateInterval(updateInterval);
        _aterms.emplace_back(std::move(beam));
      } else if (atermType == "paf") {
        std::vector<std::string> antennaMap =
            reader.GetStringList(atermName + ".antenna_map");
        if (antennaMap.size() != _nAntenna)
          throw std::runtime_error(
              "Antenna map in paf term of aterm config contains " +
              std::to_string(antennaMap.size()) +
              " antennas, whereas the measurement set consists of " +
              std::to_string(_nAntenna));
        std::vector<std::string> beamMap =
            reader.GetStringList(atermName + ".beam_map");
        if (beamMap.size() != _settings.filenames.size())
          throw std::runtime_error(
              "Number of beams specified in aterm config (" +
              std::to_string(beamMap.size()) +
              ") should match the number of measurement sets specified on the "
              "command line (" +
              std::to_string(_settings.filenames.size()) + ")");
        std::vector<std::string> beamPointings =
            reader.GetStringList(atermName + ".beam_pointings");
        if (beamPointings.size() != _settings.filenames.size() * 2)
          throw std::runtime_error("Size of beam pointings is invalid");
        double beamRA = RaDecCoord::ParseRA(beamPointings[filenameIndex * 2]);
        double beamDec =
            RaDecCoord::ParseDec(beamPointings[filenameIndex * 2 + 1]);
        std::string fileTemplate =
            reader.GetString(atermName + ".file_template");
        std::unique_ptr<PAFBeamTerm> f(new PAFBeamTerm(_coordinateSystem));
        f->Open(fileTemplate, antennaMap, beamMap[filenameIndex], beamRA,
                beamDec);
        std::string windowStr =
            reader.GetStringOr(atermName + ".window", "raised-hann");
        WindowFunction::Type window = WindowFunction::GetType(windowStr);
        if (window == WindowFunction::Tukey)
          f->SetTukeyWindow(double(_settings.paddedImageWidth) /
                            _settings.trimmedImageWidth);
        else
          f->SetWindow(window);
        f->SetDownSample(reader.GetBoolOr(atermName + ".downsample", true));
        f->SetReferenceFrequency(
            reader.GetDoubleOr(atermName + ".reference_frequency", 0.0));
        _aterms.emplace_back(std::move(f));
      }
      _aterms.back()->SetSaveATerms(
          false, _settings.prefixName);  // done by config after combining
    }
    Logger::Debug << "Constructed an a-term configuration with "
                  << _aterms.size() << " terms.\n";
    if (_aterms.empty()) {
      throw std::runtime_error(
          "The specified a-term configuration does not define any terms to "
          "apply");
    }
    if (_aterms.size() > 1) {
      _previousAterms.resize(_aterms.size());
      for (aocommon::UVector<std::complex<float>>& buf : _previousAterms)
        buf.resize(_coordinateSystem.width * _coordinateSystem.height *
                   _nAntenna * 4);
    }
  }

  virtual bool Calculate(std::complex<float>* buffer, double time,
                         double frequency, size_t fieldId,
                         const double* uvwInM) final override {
    if (_aterms.size() == 1) {
      bool result =
          _aterms.front()->Calculate(buffer, time, frequency, fieldId, uvwInM);
      if (result)
        saveATermsIfNecessary(buffer, _nAntenna, _coordinateSystem.width,
                              _coordinateSystem.height);
      return result;
    } else {
      bool isUpdated = false;
      for (size_t i = 0; i != _aterms.size(); ++i) {
        bool atermUpdated = _aterms[i]->Calculate(
            _previousAterms[i].data(), time, frequency, fieldId, uvwInM);
        isUpdated = isUpdated || atermUpdated;
      }

      if (isUpdated) {
        std::copy(_previousAterms[0].begin(), _previousAterms[0].end(), buffer);
        for (size_t i = 1; i != _aterms.size(); ++i) {
          for (size_t j = 0; j != _coordinateSystem.width *
                                      _coordinateSystem.height * _nAntenna * 4;
               j += 4) {
            std::complex<float> scratch[4];
            Matrix2x2::ATimesB(scratch, &_previousAterms[i][j], &buffer[j]);
            Matrix2x2::Assign(&buffer[j], scratch);
          }
        }
        saveATermsIfNecessary(buffer, _nAntenna, _coordinateSystem.width,
                              _coordinateSystem.height);
      }

      return isUpdated;
    }
  }

  virtual double AverageUpdateTime() const override {
    double avgTime = _aterms.front()->AverageUpdateTime();
    for (size_t i = 1; i < _aterms.size(); ++i)
      avgTime = std::min(avgTime, _aterms[i]->AverageUpdateTime());
    return avgTime;
  }

 private:
  size_t _nAntenna;
  CoordinateSystem _coordinateSystem;
  std::vector<std::unique_ptr<ATermBase>> _aterms;
  std::vector<aocommon::UVector<std::complex<float>>> _previousAterms;
  const WSCleanSettings& _settings;
};

#endif
