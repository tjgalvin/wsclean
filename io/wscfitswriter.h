#ifndef WSC_FITS_WRITER_H
#define WSC_FITS_WRITER_H

#include <optional>
#include <string>

#include <aocommon/fits/fitsreader.h>
#include <aocommon/fits/fitswriter.h>
#include <aocommon/polarization.h>

#include <radler/radler.h>

#include "imagefilename.h"
#include "../main/settings.h"
#include "../structures/observationinfo.h"
#include "../structures/outputchannelinfo.h"
#include "../structures/imagingtable.h"

class WSCFitsWriter {
 public:
  WSCFitsWriter(const ImagingTableEntry& entry, bool isImaginary,
                const Settings& settings,
                const std::optional<radler::Radler>& deconvolution,
                const ObservationInfo& observationInfo, size_t majorIterationNr,
                const std::string& commandLine,
                const OutputChannelInfo& channelInfo, bool isModel,
                double startTime);

  WSCFitsWriter(const ImagingTableEntry& entry,
                aocommon::PolarizationEnum polarization, bool isImaginary,
                const Settings& settings,
                const std::optional<radler::Radler>& deconvolution,
                const ObservationInfo& observationInfo, size_t majorIterationNr,
                const std::string& commandLine,
                const OutputChannelInfo& channelInfo, bool isModel,
                double startTime);

  explicit WSCFitsWriter(aocommon::FitsReader& templateReader);

  aocommon::FitsWriter& Writer() { return _writer; }

  template <typename NumT>
  void WriteImage(const std::string& suffix, const NumT* image);

  template <typename NumT>
  void WriteUV(const std::string& suffix, const NumT* image);

  template <typename NumT>
  void WritePSF(const std::string& fullname, const NumT* image);

  /**
   * Restore an elliptical beam using a FFT deconvolution directly from images.
   */
  static void Restore(const class Settings& settings);

  static void RestoreList(const class Settings& settings);

  static ObservationInfo ReadObservationInfo(
      const aocommon::FitsReader& reader);

 private:
  void setSettingsKeywords(const Settings& settings,
                           const std::string& commandLine);

  void setGridderConfiguration(const Settings& settings,
                               const ObservationInfo& observationInfo,
                               double startTime);

  void setDeconvolutionKeywords(const Settings& settings);

  void setDeconvolutionResultKeywords(size_t minorIterationNr,
                                      size_t majorIterationNr);

  void setChannelKeywords(const ImagingTableEntry& entry,
                          aocommon::PolarizationEnum polarization,
                          const OutputChannelInfo& channelInfo);

  void copyWSCleanKeywords(aocommon::FitsReader& reader);

 private:
  aocommon::FitsWriter _writer;
  std::string _filenamePrefix;
};

#endif
