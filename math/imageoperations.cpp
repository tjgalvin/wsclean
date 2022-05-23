#include "imageoperations.h"

#include "../io/wscfitswriter.h"

#include "../main/settings.h"

#include "../math/renderer.h"

#include <aocommon/logger.h>
#include <aocommon/fits/fitsreader.h>
#include <aocommon/fits/fitswriter.h>
#include <aocommon/units/angle.h>

#include <schaapcommon/fft/restoreimage.h>
#include <schaapcommon/fitters/gaussianfitter.h>

using aocommon::Logger;
using aocommon::units::Angle;

void ImageOperations::FitBeamSize(const Settings& settings, double& bMaj,
                                  double& bMin, double& bPA, const float* image,
                                  double beamEstimate) {
  schaapcommon::fitters::GaussianFitter beamFitter;
  Logger::Info << "Fitting beam... ";
  Logger::Info.Flush();
  if (settings.circularBeam) {
    bMaj = beamEstimate;
    beamFitter.Fit2DCircularGaussianCentred(image, settings.trimmedImageWidth,
                                            settings.trimmedImageHeight, bMaj,
                                            settings.beamFittingBoxSize);
    bMin = bMaj;
    bPA = 0.0;
  } else {
    beamFitter.Fit2DGaussianCentred(
        image, settings.trimmedImageWidth, settings.trimmedImageHeight,
        beamEstimate, bMaj, bMin, bPA, settings.beamFittingBoxSize);
  }
  bMaj = bMaj * 0.5 * (settings.pixelScaleX + settings.pixelScaleY);
  bMin = bMin * 0.5 * (settings.pixelScaleX + settings.pixelScaleY);
}

void ImageOperations::DetermineBeamSize(const Settings& settings, double& bMaj,
                                        double& bMin, double& bPA,
                                        double& bTheoretical,
                                        const float* image,
                                        double initialEstimate) {
  bTheoretical = initialEstimate;
  if (settings.gaussianTaperBeamSize != 0.0) {
    if (settings.gaussianTaperBeamSize > bTheoretical) {
      bTheoretical = settings.gaussianTaperBeamSize;
      Logger::Debug << "Beam is tapered; using "
                    << Angle::ToNiceString(bTheoretical)
                    << " as initial value in PSF fitting.\n";
    }
  }
  if (settings.manualBeamMajorSize != 0.0) {
    bMaj = settings.manualBeamMajorSize;
    bMin = settings.manualBeamMinorSize;
    bPA = settings.manualBeamPA;
  } else if (settings.fittedBeam) {
    FitBeamSize(
        settings, bMaj, bMin, bPA, image,
        bTheoretical * 2.0 / (settings.pixelScaleX + settings.pixelScaleY));
    Logger::Info << "major=" << Angle::ToNiceString(bMaj)
                 << ", minor=" << Angle::ToNiceString(bMin)
                 << ", PA=" << Angle::ToNiceString(bPA)
                 << ", theoretical=" << Angle::ToNiceString(bTheoretical)
                 << ".\n";
  } else if (settings.theoreticBeam) {
    bMaj = bTheoretical;
    bMin = bTheoretical;
    bPA = 0.0;
    Logger::Info << "Beam size is " << Angle::ToNiceString(bTheoretical)
                 << '\n';
  } else {
    bMaj = std::numeric_limits<double>::quiet_NaN();
    bMin = std::numeric_limits<double>::quiet_NaN();
    bPA = std::numeric_limits<double>::quiet_NaN();
  }
}

void ImageOperations::MakeMFSImage(
    const Settings& settings,
    const std::vector<OutputChannelInfo>& infoPerChannel,
    OutputChannelInfo& mfsInfo, const string& suffix, size_t intervalIndex,
    aocommon::PolarizationEnum pol, bool isImaginary, bool isPSF) {
  double lowestFreq = 0.0, highestFreq = 0.0;
  const size_t size = settings.trimmedImageWidth * settings.trimmedImageHeight;
  aocommon::UVector<float> mfsImage(size, 0.0);
  aocommon::UVector<double> addedImage(size), weightImage(size, 0.0);
  double weightSum = 0.0;
  aocommon::FitsWriter writer;
  for (size_t ch = 0; ch != settings.channelsOut; ++ch) {
    std::string prefixStr =
        isPSF ? ImageFilename::GetPSFPrefix(settings, ch, intervalIndex)
              : ImageFilename::GetPrefix(settings, pol, ch, intervalIndex,
                                         isImaginary);
    const std::string name(prefixStr + '-' + suffix);
    aocommon::FitsReader reader(name);
    if (ch == 0) {
      WSCFitsWriter wscWriter(reader);
      writer = wscWriter.Writer();
      lowestFreq = reader.Frequency() - reader.Bandwidth() * 0.5;
      highestFreq = reader.Frequency() + reader.Bandwidth() * 0.5;
    } else {
      lowestFreq =
          std::min(lowestFreq, reader.Frequency() - reader.Bandwidth() * 0.5);
      highestFreq =
          std::max(highestFreq, reader.Frequency() + reader.Bandwidth() * 0.5);
    }
    const double weight = infoPerChannel[ch].weight;
    weightSum += weight;
    reader.Read(addedImage.data());
    for (size_t i = 0; i != size; ++i) {
      if (std::isfinite(addedImage[i])) {
        mfsImage[i] += addedImage[i] * weight;
        weightImage[i] += weight;
      }
    }
  }
  for (size_t i = 0; i != size; ++i) mfsImage[i] /= weightImage[i];

  if (isPSF) {
    const double pixelScale =
        std::min(settings.pixelScaleX, settings.pixelScaleY);
    const double smallestTheoreticBeamSize =
        std::max(SmallestTheoreticBeamSize(infoPerChannel), pixelScale);

    ImageOperations::DetermineBeamSize(
        settings, mfsInfo.beamMaj, mfsInfo.beamMin, mfsInfo.beamPA,
        mfsInfo.theoreticBeamSize, mfsImage.data(), smallestTheoreticBeamSize);
  }
  if (std::isfinite(mfsInfo.beamMaj))
    writer.SetBeamInfo(mfsInfo.beamMaj, mfsInfo.beamMin, mfsInfo.beamPA);
  else
    writer.SetNoBeamInfo();

  std::string mfsName(ImageFilename::GetMFSPrefix(settings, pol, intervalIndex,
                                                  isImaginary, isPSF) +
                      '-' + suffix);
  Logger::Info << "Writing " << mfsName << "...\n";
  writer.SetFrequency((lowestFreq + highestFreq) * 0.5,
                      highestFreq - lowestFreq);
  writer.SetExtraKeyword("WSCIMGWG", weightSum);
  writer.RemoveExtraKeyword("WSCCHANS");
  writer.RemoveExtraKeyword("WSCCHANE");
  writer.Write(mfsName, mfsImage.data());
}

void ImageOperations::RenderMFSImage(const Settings& settings,
                                     const OutputChannelInfo& mfsInfo,
                                     size_t intervalIndex,
                                     aocommon::PolarizationEnum pol,
                                     bool isImaginary, bool isPBCorrected) {
  const size_t size = settings.trimmedImageWidth * settings.trimmedImageHeight;

  std::string mfsPrefix(ImageFilename::GetMFSPrefix(
      settings, pol, intervalIndex, isImaginary, false));
  std::string postfix = isPBCorrected ? "-pb.fits" : ".fits";
  aocommon::FitsReader residualReader(mfsPrefix + "-residual" + postfix);
  aocommon::FitsReader modelReader(mfsPrefix + "-model" + postfix);
  aocommon::UVector<float> image(size), modelImage(size);
  residualReader.Read(image.data());
  modelReader.Read(modelImage.data());

  double beamMaj = mfsInfo.beamMaj;
  double beamMin, beamPA;
  std::string beamStr;
  if (std::isfinite(beamMaj)) {
    beamMin = mfsInfo.beamMin;
    beamPA = mfsInfo.beamPA;
    beamStr = "(beam=" + Angle::ToNiceString(beamMin) + "-" +
              Angle::ToNiceString(beamMaj) +
              ", PA=" + Angle::ToNiceString(beamPA) + ")";
  } else {
    beamStr = "(beam is neither fitted nor estimated -- using delta scales!)";
    beamMaj = 0.0;
    beamMin = 0.0;
    beamPA = 0.0;
  }
  Logger::Info << "Rendering sources to restored image " + beamStr + "... ";
  Logger::Info.Flush();
  bool hasWarned = false;
  for (float& v : modelImage) {
    if (!std::isfinite(v) || std::fabs(v) > 1.0e9) {
      if (!hasWarned) {
        Logger::Warn
            << "\nWarning: Some beam corrected model components are NaN or of "
               "extremely high flux (> 10^9). These won't\n"
               "be restored on the pb mf image.\n"
               "This can be caused by an undefined or zero beam inside the "
               "FOV. Be sure to check the individual model\n"
               "images to check if the mf image is as expected.\n";
        hasWarned = true;
      }
      v = 0.0;
    }
  }
  schaapcommon::fft::RestoreImage(
      image.data(), modelImage.data(), settings.trimmedImageWidth,
      settings.trimmedImageHeight, beamMaj, beamMin, beamPA,
      settings.pixelScaleX, settings.pixelScaleY, settings.threadCount);
  Logger::Info << "DONE\n";

  Logger::Info << "Writing " << mfsPrefix << "-image" << postfix << "...\n";
  aocommon::FitsWriter imageWriter(residualReader);
  imageWriter.Write(mfsPrefix + "-image" + postfix, image.data());
}
