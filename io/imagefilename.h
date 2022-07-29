#ifndef IMAGE_FILENAME_H
#define IMAGE_FILENAME_H

#include "../main/settings.h"

#include <aocommon/polarization.h>

class ImageFilename {
 public:
  ImageFilename()
      : _polarization(aocommon::Polarization::StokesI),
        _channelIndex(0),
        _intervalIndex(0),
        _isImaginary(false) {}

  ImageFilename(size_t channelIndex, size_t intervalIndex)
      : _polarization(aocommon::Polarization::StokesI),
        _channelIndex(channelIndex),
        _intervalIndex(intervalIndex),
        _isImaginary(false) {}

  std::string GetPrefix(const Settings& settings) const {
    return GetPrefix(settings, _polarization, _channelIndex, _intervalIndex,
                     _isImaginary);
  }

  std::string GetPSFPrefix(const Settings& settings) const {
    return GetPSFPrefix(settings, _channelIndex, _intervalIndex);
  }

  std::string GetBeamPrefix(const Settings& settings) const {
    return GetBeamPrefix(settings, _channelIndex, _intervalIndex);
  }

  aocommon::PolarizationEnum GetPolarization() const { return _polarization; }
  void SetPolarization(aocommon::PolarizationEnum p) { _polarization = p; }
  void SetIsImaginary(bool isImaginary) { _isImaginary = isImaginary; }

  static std::string GetPSFPrefix(const Settings& settings, size_t channelIndex,
                                  size_t intervalIndex, size_t ddPsfIndex) {
    std::ostringstream partPrefixNameStr;
    partPrefixNameStr << settings.prefixName;
    if (settings.intervalsOut != 1)
      partPrefixNameStr << "-t" << fourDigitStr(intervalIndex);
    if (settings.ddPsfGridHeight > 1 || settings.ddPsfGridWidth > 1)
      partPrefixNameStr << "-d" << fourDigitStr(ddPsfIndex);
    if (settings.channelsOut != 1)
      partPrefixNameStr << '-' << fourDigitStr(channelIndex);
    return partPrefixNameStr.str();
  }

  static std::string GetPSFPrefix(const Settings& settings, size_t channelIndex,
                                  size_t intervalIndex) {
    std::ostringstream partPrefixNameStr;
    partPrefixNameStr << settings.prefixName;
    if (settings.intervalsOut != 1)
      partPrefixNameStr << "-t" << fourDigitStr(intervalIndex);
    if (settings.channelsOut != 1)
      partPrefixNameStr << '-' << fourDigitStr(channelIndex);
    return partPrefixNameStr.str();
  }

  static std::string GetPrefix(const Settings& settings,
                               aocommon::PolarizationEnum polarization,
                               size_t channelIndex, size_t intervalIndex,
                               bool isImaginary) {
    std::ostringstream partPrefixNameStr;
    partPrefixNameStr << settings.prefixName;
    if (settings.intervalsOut != 1)
      partPrefixNameStr << "-t" << fourDigitStr(intervalIndex);
    if (settings.channelsOut != 1)
      partPrefixNameStr << '-' << fourDigitStr(channelIndex);
    if (settings.polarizations.size() != 1) {
      partPrefixNameStr << '-'
                        << aocommon::Polarization::TypeToShortString(
                               polarization);
      if (isImaginary) partPrefixNameStr << 'i';
    }
    return partPrefixNameStr.str();
  }

  static std::string GetBeamPrefix(const Settings& settings,
                                   size_t channelIndex, size_t intervalIndex) {
    std::ostringstream partPrefixNameStr;
    partPrefixNameStr << settings.prefixName;
    if (settings.intervalsOut != 1)
      partPrefixNameStr << "-t" << fourDigitStr(intervalIndex);
    if (settings.channelsOut != 1)
      partPrefixNameStr << '-' << fourDigitStr(channelIndex);
    partPrefixNameStr << "-beam";
    return partPrefixNameStr.str();
  }

  static std::string GetMFSPrefix(const Settings& settings,
                                  aocommon::PolarizationEnum polarization,
                                  size_t intervalIndex, bool isImaginary,
                                  bool isPSF) {
    std::ostringstream partPrefixNameStr;
    partPrefixNameStr << settings.prefixName;
    if (settings.intervalsOut != 1)
      partPrefixNameStr << "-t" << fourDigitStr(intervalIndex);
    if (settings.channelsOut != 1) partPrefixNameStr << "-MFS";
    if (settings.polarizations.size() != 1 && !isPSF) {
      partPrefixNameStr << '-'
                        << aocommon::Polarization::TypeToShortString(
                               polarization);
      if (isImaginary) partPrefixNameStr << 'i';
    }
    return partPrefixNameStr.str();
  }

 private:
  aocommon::PolarizationEnum _polarization;
  size_t _channelIndex;
  size_t _intervalIndex;
  bool _isImaginary;

  static std::string fourDigitStr(size_t val) {
    std::ostringstream str;
    if (val < 1000) str << '0';
    if (val < 100) str << '0';
    if (val < 10) str << '0';
    str << val;
    return str.str();
  }
};

#endif
