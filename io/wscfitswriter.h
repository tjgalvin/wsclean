#ifndef WSC_FITS_WRITER_H
#define WSC_FITS_WRITER_H

#include <optional>
#include <string>

#include <aocommon/fits/fitsreader.h>
#include <aocommon/fits/fitswriter.h>
#include <aocommon/image.h>
#include <aocommon/polarization.h>

#include <schaapcommon/facets/facetimage.h>

#include <radler/radler.h>

#include "imagefilename.h"
#include "../main/settings.h"
#include "../structures/observationinfo.h"
#include "../structures/outputchannelinfo.h"
#include "../structures/imagingtable.h"

/**
 * @brief Class to write FITS images with the appropriate generic and WSClean
 * specific keywords
 *
 * The configuration is set in the constructor. This includes a pointing,
 * possibly a shift and a default image size. The underlying FitsWriter object
 * can be obtained by a call to the Writer() function. This can be used to write
 * images of default size from a raw pointer.
 *
 * Images of different sizes can be written using the WriteImage() and
 * WriteFullNameImage() functions. These accept an Image as argument, which can
 * be of a different size than the default size supplied to the constructor.
 *
 * The WriteFullNameImage function also accepts facets, writing the correct size
 * and shift onto the FITS file.
 */
class WSCFitsWriter {
 public:
  WSCFitsWriter(const ImagingTableEntry& entry, bool isImaginary,
                const Settings& settings,
                const std::optional<radler::Radler>& deconvolution,
                const ObservationInfo& observationInfo, double shiftL,
                double shiftM, size_t majorIterationNr,
                const std::string& commandLine,
                const OutputChannelInfo& channelInfo, bool isModel,
                double startTime);

  WSCFitsWriter(const ImagingTableEntry& entry,
                aocommon::PolarizationEnum polarization, bool isImaginary,
                const Settings& settings,
                const std::optional<radler::Radler>& deconvolution,
                const ObservationInfo& observationInfo, double shiftL,
                double shiftM, size_t majorIterationNr,
                const std::string& commandLine,
                const OutputChannelInfo& channelInfo, bool isModel,
                double startTime);

  explicit WSCFitsWriter(aocommon::FitsReader& templateReader);

  explicit WSCFitsWriter(const aocommon::FitsWriter& writer);

  aocommon::FitsWriter& Writer() { return _writer; }
  const aocommon::FitsWriter& Writer() const { return _writer; }

  size_t Width() const { return _writer.Width(); }
  size_t Height() const { return _writer.Height(); }
  void SetImageDimensions(size_t width, size_t height) {
    _writer.SetImageDimensions(width, height);
  }

  void WriteImage(const std::string& suffix, const aocommon::Image& image);

  void WriteFullNameImage(const std::string& fullname,
                          const aocommon::Image& image);

  void WriteFullNameImage(const std::string& fullname,
                          const aocommon::Image& image,
                          const schaapcommon::facets::Facet& facet);

  void WriteFullNameImage(const std::string& fullname,
                          const schaapcommon::facets::FacetImage& facetimage);

  template <typename NumT>
  void WriteUV(const std::string& suffix, const NumT* image);

  /**
   * Restore an elliptical beam using a FFT deconvolution directly from images.
   */
  static void Restore(const class Settings& settings);

  static void RestoreList(const class Settings& settings);

 private:
  void setSettingsKeywords(const Settings& settings,
                           const std::string& commandLine);

  void setGridderConfiguration(const Settings& settings,
                               const ObservationInfo& observationInfo,
                               double shiftL, double shiftM, double startTime);

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
