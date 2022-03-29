#include "wscfitswriter.h"

#include <wscversion.h>

#include "imagefilename.h"

#include "../math/modelrenderer.h"

#include "../model/bbsmodel.h"

#include "../gridding/msgridderbase.h"

WSCFitsWriter::WSCFitsWriter(const ImagingTableEntry& entry, bool isImaginary,
                             const Settings& settings,
                             const std::optional<Deconvolution>& deconvolution,
                             const ObservationInfo& observationInfo,
                             size_t majorIterationNr,
                             const std::string& commandLine,
                             const OutputChannelInfo& channelInfo, bool isModel,
                             double startTime) {
  _filenamePrefix = ImageFilename::GetPrefix(
      settings, entry.polarization, entry.outputChannelIndex,
      entry.outputIntervalIndex, isImaginary);
  setGridderConfiguration(settings, observationInfo, startTime);
  setSettingsKeywords(settings, commandLine);
  setChannelKeywords(entry, entry.polarization, channelInfo);
  setDeconvolutionKeywords(settings);
  if (deconvolution.has_value() && deconvolution->IsInitialized()) {
    setDeconvolutionResultKeywords(deconvolution->IterationNumber(),
                                   majorIterationNr);
  }
  if (isModel) _writer.SetUnit(aocommon::FitsWriter::JanskyPerPixel);
}

WSCFitsWriter::WSCFitsWriter(const ImagingTableEntry& entry,
                             aocommon::PolarizationEnum polarization,
                             bool isImaginary, const Settings& settings,
                             const std::optional<Deconvolution>& deconvolution,
                             const ObservationInfo& observationInfo,
                             size_t majorIterationNr,
                             const std::string& commandLine,
                             const OutputChannelInfo& channelInfo, bool isModel,
                             double startTime) {
  _filenamePrefix =
      ImageFilename::GetPrefix(settings, polarization, entry.outputChannelIndex,
                               entry.outputIntervalIndex, isImaginary);
  setGridderConfiguration(settings, observationInfo, startTime);
  setSettingsKeywords(settings, commandLine);
  setChannelKeywords(entry, polarization, channelInfo);
  setDeconvolutionKeywords(settings);
  if (deconvolution.has_value() && deconvolution->IsInitialized()) {
    setDeconvolutionResultKeywords(deconvolution->IterationNumber(),
                                   majorIterationNr);
  }
  if (isModel) _writer.SetUnit(aocommon::FitsWriter::JanskyPerPixel);
}

WSCFitsWriter::WSCFitsWriter(aocommon::FitsReader& templateReader)
    : _writer(templateReader) {
  copyWSCleanKeywords(templateReader);
}

void WSCFitsWriter::setSettingsKeywords(const Settings& settings,
                                        const std::string& commandLine) {
  _writer.SetOrigin("WSClean", "W-stacking imager written by Andre Offringa");
  _writer.AddHistory(commandLine);
  _writer.SetExtraKeyword("WSCVERSI", WSCLEAN_VERSION_STR);
  _writer.SetExtraKeyword("WSCVDATE", WSCLEAN_VERSION_DATE);
  if (settings.endChannel != 0) {
    _writer.SetExtraKeyword("WSCCHANS", settings.startChannel);
    _writer.SetExtraKeyword("WSCCHANE", settings.endChannel);
  }
  if (settings.endTimestep != 0) {
    _writer.SetExtraKeyword("WSCTIMES", settings.startTimestep);
    _writer.SetExtraKeyword("WSCTIMEE", settings.endTimestep);
  }
  _writer.SetExtraKeyword("WSCFIELD", settings.fieldIds[0]);
}

void WSCFitsWriter::setGridderConfiguration(
    const Settings& settings, const ObservationInfo& observationInfo,
    double startTime) {
  const double ra = observationInfo.phaseCentreRA,
               dec = observationInfo.phaseCentreDec,
               pixelScaleX = settings.pixelScaleX,
               pixelScaleY = settings.pixelScaleY;

  _writer.SetImageDimensions(settings.trimmedImageWidth,
                             settings.trimmedImageHeight, ra, dec, pixelScaleX,
                             pixelScaleY);

  _writer.SetDate(startTime);
  if (observationInfo.hasShiftedPhaseCentre)
    _writer.SetPhaseCentreShift(observationInfo.shiftL, observationInfo.shiftM);
  _writer.SetTelescopeName(observationInfo.telescopeName);
  _writer.SetObserver(observationInfo.observer);
  _writer.SetObjectName(observationInfo.fieldName);

  /* This is the normalization factor that was applied. The factor is useful
   * to undo the normalization for e.g. conversion to Kelvins. */
  _writer.SetExtraKeyword("WSCDATAC", settings.dataColumnName);
  _writer.SetExtraKeyword("WSCWEIGH", settings.weightMode.ToString());
  _writer.SetExtraKeyword("WSCGKRNL", settings.antialiasingKernelSize);
}

void WSCFitsWriter::setDeconvolutionKeywords(const Settings& settings) {
  _writer.SetExtraKeyword("WSCNITER", settings.deconvolutionIterationCount);
  _writer.SetExtraKeyword("WSCTHRES", settings.deconvolutionThreshold);
  _writer.SetExtraKeyword("WSCGAIN", settings.deconvolutionGain);
  _writer.SetExtraKeyword("WSCMGAIN", settings.deconvolutionMGain);
  _writer.SetExtraKeyword("WSCNEGCM", settings.allowNegativeComponents);
  _writer.SetExtraKeyword("WSCNEGST", settings.stopOnNegativeComponents);
}

void WSCFitsWriter::setDeconvolutionResultKeywords(size_t minorIterationNr,
                                                   size_t majorIterationNr) {
  _writer.SetExtraKeyword("WSCMINOR", minorIterationNr);
  _writer.SetExtraKeyword("WSCMAJOR", majorIterationNr);
}

void WSCFitsWriter::setChannelKeywords(const ImagingTableEntry& entry,
                                       aocommon::PolarizationEnum polarization,
                                       const OutputChannelInfo& channelInfo) {
  const double bandStart = entry.bandStartFrequency,
               bandEnd = entry.bandEndFrequency,
               centreFrequency = 0.5 * (bandStart + bandEnd),
               bandwidth = bandEnd - bandStart;
  _writer.SetFrequency(centreFrequency, bandwidth);
  _writer.SetExtraKeyword("WSCIMGWG", channelInfo.weight);
  _writer.SetExtraKeyword("WSCNORMF", channelInfo.normalizationFactor);
  _writer.SetExtraKeyword("WSCNWLAY", channelInfo.wGridSize);
  _writer.SetExtraKeyword("WSCNVIS", channelInfo.visibilityCount);
  _writer.SetExtraKeyword("WSCENVIS", channelInfo.effectiveVisibilityCount);
  _writer.SetExtraKeyword("WSCVWSUM", channelInfo.visibilityWeightSum);
  _writer.SetBeamInfo(channelInfo.beamMaj, channelInfo.beamMin,
                      channelInfo.beamPA);
  _writer.SetPolarization(polarization);
}

void WSCFitsWriter::copyWSCleanKeywords(aocommon::FitsReader& reader) {
  const size_t N_STRKEYWORDS = 4, N_DBLKEYWORDS = 20;
  const char* strKeywords[N_STRKEYWORDS] = {"WSCVERSI", "WSCVDATE", "WSCDATAC",
                                            "WSCWEIGH"};
  const char* dblKeywords[N_DBLKEYWORDS] = {
      "WSCIMGWG", "WSCNWLAY", "WSCGKRNL", "WSCCHANS", "WSCCHANE",
      "WSCTIMES", "WSCTIMEE", "WSCFIELD", "WSCNITER", "WSCNORMF",
      "WSCTHRES", "WSCGAIN",  "WSCMGAIN", "WSCNEGCM", "WSCNEGST",
      "WSCMINOR", "WSCMAJOR", "WSCNVIS",  "WSCENVIS", "WSCVWSUM"};
  for (size_t i = 0; i != N_STRKEYWORDS; ++i)
    _writer.CopyStringKeywordIfExists(reader, strKeywords[i]);
  for (size_t i = 0; i != N_DBLKEYWORDS; ++i)
    _writer.CopyDoubleKeywordIfExists(reader, dblKeywords[i]);
}

template <typename NumT>
void WSCFitsWriter::WriteImage(const std::string& suffix, const NumT* image) {
  std::string name = _filenamePrefix + '-' + suffix;
  _writer.Write(name, image);
}
template void WSCFitsWriter::WriteImage(const std::string& suffix,
                                        const double* image);
template void WSCFitsWriter::WriteImage(const std::string& suffix,
                                        const float* image);

template <typename NumT>
void WSCFitsWriter::WriteUV(const std::string& suffix, const NumT* image) {
  std::string name = _filenamePrefix + '-' + suffix;
  aocommon::FitsWriter::Unit unit = _writer.GetUnit();
  _writer.SetIsUV(true);
  _writer.SetUnit(aocommon::FitsWriter::Jansky);
  _writer.Write(name, image);
  _writer.SetIsUV(false);
  _writer.SetUnit(unit);
}
template void WSCFitsWriter::WriteUV(const std::string& suffix,
                                     const double* image);
template void WSCFitsWriter::WriteUV(const std::string& suffix,
                                     const float* image);

template <typename NumT>
void WSCFitsWriter::WritePSF(const std::string& fullname, const NumT* image) {
  _writer.Write(fullname, image);
}
template void WSCFitsWriter::WritePSF(const std::string& fullname,
                                      const double* image);
template void WSCFitsWriter::WritePSF(const std::string& fullname,
                                      const float* image);

void WSCFitsWriter::Restore(const Settings& settings) {
  aocommon::FitsReader imgReader(settings.restoreInput),
      modReader(settings.restoreModel);
  if (imgReader.ImageWidth() != modReader.ImageWidth() ||
      imgReader.ImageHeight() != modReader.ImageHeight())
    throw std::runtime_error(
        "Image and model images have different dimensions!");
  aocommon::UVector<float> image(imgReader.ImageWidth() *
                                 imgReader.ImageHeight()),
      model(modReader.ImageWidth() * modReader.ImageHeight());
  imgReader.Read(image.data());
  modReader.Read(model.data());

  double beamMaj, beamMin, beamPA;
  if (settings.manualBeamMajorSize != 0.0) {
    beamMaj = settings.manualBeamMajorSize;
    beamMin = settings.manualBeamMinorSize;
    beamPA = settings.manualBeamPA;
  } else {
    beamMaj = imgReader.BeamMajorAxisRad();
    beamMin = imgReader.BeamMinorAxisRad();
    beamPA = imgReader.BeamPositionAngle();
  }

  ModelRenderer::Restore(image.data(), model.data(), imgReader.ImageWidth(),
                         imgReader.ImageHeight(), beamMaj, beamMin, beamPA,
                         imgReader.PixelSizeX(), imgReader.PixelSizeY(),
                         settings.threadCount);

  aocommon::FitsWriter writer(WSCFitsWriter(imgReader).Writer());
  writer.SetBeamInfo(beamMaj, beamMin, beamPA);
  writer.Write(settings.restoreOutput, image.data());
}

void WSCFitsWriter::RestoreList(const Settings& settings) {
  const Model model = BBSModel::Read(settings.restoreModel);
  aocommon::FitsReader imgReader(settings.restoreInput);
  aocommon::UVector<float> image(imgReader.ImageWidth() *
                                 imgReader.ImageHeight());
  imgReader.Read(image.data());

  double beamMaj, beamMin, beamPA;
  if (settings.manualBeamMajorSize != 0.0) {
    beamMaj = settings.manualBeamMajorSize;
    beamMin = settings.manualBeamMinorSize;
    beamPA = settings.manualBeamPA;
  } else {
    beamMaj = imgReader.BeamMajorAxisRad();
    beamMin = imgReader.BeamMinorAxisRad();
    beamPA = imgReader.BeamPositionAngle();
  }

  double frequency = imgReader.Frequency();
  double bandwidth = imgReader.Bandwidth();

  ModelRenderer renderer(imgReader.PhaseCentreRA(), imgReader.PhaseCentreDec(),
                         imgReader.PixelSizeX(), imgReader.PixelSizeY(),
                         imgReader.PhaseCentreDL(), imgReader.PhaseCentreDM());

  renderer.Restore(image.data(), imgReader.ImageWidth(),
                   imgReader.ImageHeight(), model, beamMaj, beamMin, beamPA,
                   frequency - bandwidth * 0.5, frequency + bandwidth * 0.5,
                   aocommon::Polarization::StokesI, settings.threadCount);

  aocommon::FitsWriter writer(WSCFitsWriter(imgReader).Writer());
  writer.SetBeamInfo(beamMaj, beamMin, beamPA);
  writer.Write(settings.restoreOutput, image.data());
}

ObservationInfo WSCFitsWriter::ReadObservationInfo(
    const aocommon::FitsReader& reader) {
  ObservationInfo obsInfo;
  obsInfo.phaseCentreRA = reader.PhaseCentreRA();
  obsInfo.phaseCentreDec = reader.PhaseCentreDec();
  obsInfo.shiftL = reader.PhaseCentreDL();
  obsInfo.shiftM = reader.PhaseCentreDM();
  obsInfo.telescopeName = reader.TelescopeName();
  obsInfo.observer = reader.Observer();
  return obsInfo;
}