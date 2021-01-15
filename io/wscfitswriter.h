#ifndef WSC_FITS_WRITER_H
#define WSC_FITS_WRITER_H

#include "imagefilename.h"

#include "../structures/observationinfo.h"
#include "../structures/outputchannelinfo.h"
#include "../structures/imagingtable.h"

#include "../main/settings.h"

#include <aocommon/fits/fitswriter.h>
#include <aocommon/polarization.h>

#include <string>

using aocommon::FitsReader;
using aocommon::FitsWriter;

class WSCFitsWriter {
 public:
  WSCFitsWriter(const ImagingTableEntry& entry, bool isImaginary,
                const Settings& settings,
                const class Deconvolution& deconvolution,
                const ObservationInfo& observationInfo, size_t majorIterationNr,
                const std::string& commandLine,
                const OutputChannelInfo& channelInfo, bool isModel,
                double startTime);

  WSCFitsWriter(const ImagingTableEntry& entry,
                aocommon::PolarizationEnum polarization, bool isImaginary,
                const Settings& settings,
                const class Deconvolution& deconvolution,
                const ObservationInfo& observationInfo, size_t majorIterationNr,
                const std::string& commandLine,
                const OutputChannelInfo& channelInfo, bool isModel,
                double startTime);

  explicit WSCFitsWriter(FitsReader& templateReader);

  FitsWriter& Writer() { return _writer; }

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

  static ObservationInfo ReadObservationInfo(class FitsReader& reader);

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

  void copyWSCleanKeywords(class FitsReader& reader);

 private:
  FitsWriter _writer;
  std::string _filenamePrefix;
};

#endif
