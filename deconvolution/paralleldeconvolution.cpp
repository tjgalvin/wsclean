#include "paralleldeconvolution.h"

#include "../multiscale/multiscalealgorithm.h"

#include "../io/imagefilename.h"

#include "../main/settings.h"

#include "../math/dijkstrasplitter.h"

#include "../units/fluxdensity.h"

#include "../structures/image.h"
#include "../structures/primarybeam.h"

#include "componentlist.h"

#include <aocommon/parallelfor.h>

ParallelDeconvolution::ParallelDeconvolution(const class Settings& settings)
    : _horImages(0),
      _verImages(0),
      _settings(settings),
      _allocator(nullptr),
      _mask(nullptr),
      _trackPerScaleMasks(false),
      _usePerScaleMasks(false) {}

ParallelDeconvolution::~ParallelDeconvolution() {}

void ParallelDeconvolution::SetAlgorithm(
    std::unique_ptr<class DeconvolutionAlgorithm> algorithm) {
  if (_settings.parallelDeconvolutionMaxSize == 0) {
    _algorithms.resize(1);
    _algorithms.front() = std::move(algorithm);
  } else {
    const size_t width = _settings.trimmedImageWidth,
                 height = _settings.trimmedImageHeight;
    size_t maxSubImageSize = _settings.parallelDeconvolutionMaxSize;
    _horImages = (width + maxSubImageSize - 1) / maxSubImageSize,
    _verImages = (height + maxSubImageSize - 1) / maxSubImageSize;
    _algorithms.resize(_horImages * _verImages);
    _algorithms.front() = std::move(algorithm);
    size_t threadsPerAlg =
        (_settings.parallelDeconvolutionMaxThreads + _algorithms.size() - 1) /
        _algorithms.size();
    _algorithms.front()->SetThreadCount(threadsPerAlg);
    Logger::Debug << "Parallel deconvolution will use " << _algorithms.size()
                  << " subimages.\n";
    for (size_t i = 1; i != _algorithms.size(); ++i)
      _algorithms[i] = _algorithms.front()->Clone();
  }
}

void ParallelDeconvolution::SetRMSFactorImage(Image&& image) {
  if (_settings.parallelDeconvolutionMaxSize == 0)
    _algorithms.front()->SetRMSFactorImage(std::move(image));
  else
    _rmsImage = std::move(image);
}

void ParallelDeconvolution::SetThreshold(double threshold) {
  for (auto& alg : _algorithms) alg->SetThreshold(threshold);
}

void ParallelDeconvolution::SetAutoMaskMode(bool trackPerScaleMasks,
                                            bool usePerScaleMasks) {
  _trackPerScaleMasks = trackPerScaleMasks;
  _usePerScaleMasks = usePerScaleMasks;
  for (auto& alg : _algorithms) {
    class MultiScaleAlgorithm& algorithm =
        static_cast<class MultiScaleAlgorithm&>(*alg);
    algorithm.SetAutoMaskMode(trackPerScaleMasks, usePerScaleMasks);
  }
}

void ParallelDeconvolution::SetCleanMask(const bool* mask) {
  _mask = mask;
  if (_algorithms.size() == 1) _algorithms.front()->SetCleanMask(mask);
}

void ParallelDeconvolution::runSubImage(
    SubImage& subImg, ImageSet& dataImage, ImageSet& modelImage,
    const aocommon::UVector<const float*>& psfImages, double majorIterThreshold,
    bool findPeakOnly, std::mutex* mutex) {
  const size_t width = _settings.trimmedImageWidth,
               height = _settings.trimmedImageHeight;

  std::unique_ptr<ImageSet> subModel, subData;
  {
    std::lock_guard<std::mutex> lock(*mutex);
    subModel = modelImage.Trim(subImg.x, subImg.y, subImg.x + subImg.width,
                               subImg.y + subImg.height, width);
    subData = dataImage.Trim(subImg.x, subImg.y, subImg.x + subImg.width,
                             subImg.y + subImg.height, width);
  }

  // Construct the smaller psfs
  std::vector<Image> subPsfs(psfImages.size());
  aocommon::UVector<const float*> subPsfVector(psfImages.size());
  for (size_t i = 0; i != psfImages.size(); ++i) {
    subPsfs[i] = Image(subImg.width, subImg.height);
    Image::Trim(subPsfs[i].data(), subImg.width, subImg.height, psfImages[i],
                width, height);
    subPsfVector[i] = subPsfs[i].data();
  }
  _algorithms[subImg.index]->SetCleanMask(subImg.mask.data());

  // Construct smaller RMS image if necessary
  if (!_rmsImage.empty()) {
    Image subRmsImage =
        _rmsImage.TrimBox(subImg.x, subImg.y, subImg.width, subImg.height);
    _algorithms[subImg.index]->SetRMSFactorImage(std::move(subRmsImage));
  }

  size_t maxNIter = _algorithms[subImg.index]->MaxNIter();
  if (findPeakOnly)
    _algorithms[subImg.index]->SetMaxNIter(0);
  else
    _algorithms[subImg.index]->SetMajorIterThreshold(majorIterThreshold);

  if (_usePerScaleMasks || _trackPerScaleMasks) {
    std::lock_guard<std::mutex> lock(*mutex);
    MultiScaleAlgorithm& msAlg =
        static_cast<class MultiScaleAlgorithm&>(*_algorithms[subImg.index]);
    // During the first iteration, msAlg will not have scales/masks yet and the
    // nr scales has also not been determined yet.
    if (!_scaleMasks.empty()) {
      // Here we set the scale mask for the multiscale algorithm.
      // The maximum number of scales in the previous iteration can be found by
      // _scaleMasks.size() Not all msAlgs might have used that many scales, so
      // we have to take this into account
      msAlg.SetScaleMaskCount(
          std::max(msAlg.GetScaleMaskCount(), _scaleMasks.size()));
      for (size_t i = 0; i != msAlg.GetScaleMaskCount(); ++i) {
        aocommon::UVector<bool>& output = msAlg.GetScaleMask(i);
        output.assign(subImg.width * subImg.height, false);
        if (i < _scaleMasks.size())
          Image::TrimBox(output.data(), subImg.x, subImg.y, subImg.width,
                         subImg.height, _scaleMasks[i].data(), width, height);
      }
    }
  }

  subImg.peak = _algorithms[subImg.index]->ExecuteMajorIteration(
      *subData, *subModel, subPsfVector, subImg.width, subImg.height,
      subImg.reachedMajorThreshold);

  // Since this was an RMS image specifically for this subimage size, we free it
  // immediately
  _algorithms[subImg.index]->SetRMSFactorImage(Image());

  if (_trackPerScaleMasks) {
    std::lock_guard<std::mutex> lock(*mutex);
    MultiScaleAlgorithm& msAlg =
        static_cast<class MultiScaleAlgorithm&>(*_algorithms[subImg.index]);
    if (_scaleMasks.empty()) {
      _scaleMasks.resize(msAlg.ScaleCount());
      for (aocommon::UVector<bool>& scaleMask : _scaleMasks)
        scaleMask.assign(width * height, false);
    }
    for (size_t i = 0; i != msAlg.ScaleCount(); ++i) {
      const aocommon::UVector<bool>& msMask = msAlg.GetScaleMask(i);
      if (i < _scaleMasks.size())
        Image::CopyMasked(_scaleMasks[i].data(), subImg.x, subImg.y, width,
                          msMask.data(), subImg.width, subImg.height,
                          subImg.mask.data());
    }
  }

  if (_settings.saveSourceList && _settings.useMultiscale) {
    std::lock_guard<std::mutex> lock(*mutex);
    MultiScaleAlgorithm& algorithm =
        static_cast<MultiScaleAlgorithm&>(*_algorithms[subImg.index]);
    if (!_componentList)
      _componentList.reset(new ComponentList(
          width, height, algorithm.ScaleCount(), dataImage.size()));
    _componentList->Add(algorithm.GetComponentList(), subImg.x, subImg.y);
    algorithm.ClearComponentList();
  }

  if (findPeakOnly) _algorithms[subImg.index]->SetMaxNIter(maxNIter);

  if (!findPeakOnly) {
    std::lock_guard<std::mutex> lock(*mutex);
    dataImage.CopyMasked(*subData, subImg.x, subImg.y, width, subImg.width,
                         subImg.height, subImg.mask.data());
    modelImage.CopyMasked(*subModel, subImg.x, subImg.y, width, subImg.width,
                          subImg.height, subImg.mask.data());
  }
}

void ParallelDeconvolution::ExecuteMajorIteration(
    class ImageSet& dataImage, class ImageSet& modelImage,
    const aocommon::UVector<const float*>& psfImages,
    bool& reachedMajorThreshold) {
  const size_t width = _settings.trimmedImageWidth,
               height = _settings.trimmedImageHeight;
  if (_algorithms.size() == 1) {
    ForwardingLogReceiver fwdReceiver;
    _algorithms.front()->SetLogReceiver(fwdReceiver);
    _algorithms.front()->ExecuteMajorIteration(
        dataImage, modelImage, psfImages, width, height, reachedMajorThreshold);
  } else {
    executeParallelRun(dataImage, modelImage, psfImages, reachedMajorThreshold);
  }
}

void ParallelDeconvolution::executeParallelRun(
    class ImageSet& dataImage, class ImageSet& modelImage,
    const aocommon::UVector<const float*>& psfImages,
    bool& reachedMajorThreshold) {
  const size_t width = _settings.trimmedImageWidth,
               height = _settings.trimmedImageHeight,
               avgHSubImageSize = width / _horImages,
               avgVSubImageSize = height / _verImages;

  Image image(width, height), dividingLine(width, height, 0.0);
  aocommon::UVector<bool> largeScratchMask(width * height);
  dataImage.GetLinearIntegrated(image);

  DijkstraSplitter divisor(width, height);

  struct VerticalArea {
    aocommon::UVector<bool> mask;
    size_t x, width;
  };
  std::vector<VerticalArea> verticalAreas(_horImages);

  Logger::Info << "Calculating edge paths...\n";
  aocommon::ParallelFor<size_t> splitLoop(_settings.threadCount);

  // Divide into columns (i.e. construct the vertical lines)
  splitLoop.Run(1, _horImages, [&](size_t divNr, size_t) {
    size_t splitStart = width * divNr / _horImages - avgHSubImageSize / 4,
           splitEnd = width * divNr / _horImages + avgHSubImageSize / 4;
    divisor.DivideVertically(image.data(), dividingLine.data(), splitStart,
                             splitEnd);
  });
  for (size_t divNr = 0; divNr != _horImages; ++divNr) {
    size_t midX = divNr * width / _horImages + avgHSubImageSize / 2;
    VerticalArea& area = verticalAreas[divNr];
    divisor.FloodVerticalArea(dividingLine.data(), midX,
                              largeScratchMask.data(), area.x, area.width);
    area.mask.resize(area.width * height);
    Image::TrimBox(area.mask.data(), area.x, 0, area.width, height,
                   largeScratchMask.data(), width, height);
  }

  // Make the rows (horizontal lines)
  dividingLine = 0.0f;
  splitLoop.Run(1, _verImages, [&](size_t divNr, size_t) {
    size_t splitStart = height * divNr / _verImages - avgVSubImageSize / 4,
           splitEnd = height * divNr / _verImages + avgVSubImageSize / 4;
    divisor.DivideHorizontally(image.data(), dividingLine.data(), splitStart,
                               splitEnd);
  });

  Logger::Info << "Calculating bounding boxes and submasks...\n";

  // Find the bounding boxes and clean masks for each subimage
  aocommon::UVector<bool> mask(width * height);
  std::vector<SubImage> subImages;
  for (size_t y = 0; y != _verImages; ++y) {
    size_t midY = y * height / _verImages + avgVSubImageSize / 2;
    size_t hAreaY, hAreaWidth;
    divisor.FloodHorizontalArea(dividingLine.data(), midY,
                                largeScratchMask.data(), hAreaY, hAreaWidth);

    for (size_t x = 0; x != _horImages; ++x) {
      subImages.emplace_back();
      SubImage& subImage = subImages.back();
      subImage.index = subImages.size() - 1;
      const VerticalArea& vArea = verticalAreas[x];
      divisor.GetBoundingMask(vArea.mask.data(), vArea.x, vArea.width,
                              largeScratchMask.data(), mask.data(), subImage.x,
                              subImage.y, subImage.width, subImage.height);
      Logger::Debug << "Subimage " << subImages.size() << " at (" << subImage.x
                    << "," << subImage.y << ") - ("
                    << subImage.x + subImage.width << ","
                    << subImage.y + subImage.height << ")\n";
      subImage.mask.resize(subImage.width * subImage.height);
      Image::TrimBox(subImage.mask.data(), subImage.x, subImage.y,
                     subImage.width, subImage.height, mask.data(), width,
                     height);

      // If a user mask is active, take the union of that mask with the division
      // mask (note that 'mask' is reused as a scratch space)
      if (_mask != nullptr) {
        Image::TrimBox(mask.data(), subImage.x, subImage.y, subImage.width,
                       subImage.height, _mask, width, height);
        for (size_t i = 0; i != subImage.mask.size(); ++i)
          subImage.mask[i] = subImage.mask[i] && mask[i];
      }
    }
  }
  verticalAreas.clear();

  // Initialize loggers
  std::mutex mutex;
  _logs.Initialize(_horImages, _verImages);
  for (size_t i = 0; i != _algorithms.size(); ++i)
    _algorithms[i]->SetLogReceiver(_logs[i]);

  // Find the starting peak over all subimages
  aocommon::ParallelFor<size_t> loop(_settings.parallelDeconvolutionMaxThreads);
  loop.Run(0, _algorithms.size(), [&](size_t index, size_t) {
    _logs.Activate(index);
    runSubImage(subImages[index], dataImage, modelImage, psfImages, 0.0, true,
                &mutex);
    _logs.Deactivate(index);

    _logs[index].Mute(false);
    _logs[index].Info << "Sub-image " << index << " returned peak position.\n";
    _logs[index].Mute(true);
  });
  double maxValue = 0.0;
  size_t indexOfMax = 0;
  for (SubImage& img : subImages) {
    if (img.peak > maxValue) {
      maxValue = img.peak;
      indexOfMax = img.index;
    }
  }
  Logger::Info << "Subimage " << (indexOfMax + 1) << " has maximum peak of "
               << FluxDensity::ToNiceString(maxValue) << ".\n";
  double mIterThreshold = maxValue * (1.0 - _settings.deconvolutionMGain);

  // Run the deconvolution
  loop.Run(0, _algorithms.size(), [&](size_t index, size_t) {
    _logs.Activate(index);
    runSubImage(subImages[index], dataImage, modelImage, psfImages,
                mIterThreshold, false, &mutex);
    _logs.Deactivate(index);

    _logs[index].Mute(false);
    _logs[index].Info << "Sub-image " << index
                      << " finished its deconvolution iteration.\n";
    _logs[index].Mute(true);
  });

  _rmsImage.reset();

  size_t subImagesFinished = 0;
  reachedMajorThreshold = false;
  bool reachedMaxNIter = false;
  for (SubImage& img : subImages) {
    if (!img.reachedMajorThreshold) ++subImagesFinished;
    if (_algorithms[img.index]->IterationNumber() >=
        _algorithms[img.index]->MaxNIter())
      reachedMaxNIter = true;
  }
  Logger::Info << subImagesFinished << " / " << subImages.size()
               << " sub-images finished";
  reachedMajorThreshold = (subImagesFinished != subImages.size());
  if (reachedMajorThreshold && !reachedMaxNIter)
    Logger::Info << ": Continue next major iteration.\n";
  else if (reachedMajorThreshold && reachedMaxNIter) {
    Logger::Info << ", but nr. of iterations reached at least once: "
                    "Deconvolution finished.\n";
    reachedMajorThreshold = false;
  } else
    Logger::Info << ": Deconvolution finished.\n";
}

void ParallelDeconvolution::SaveSourceList(CachedImageSet& modelImages,
                                           const ImagingTable& table,
                                           long double phaseCentreRA,
                                           long double phaseCentreDec) {
  std::string filename = _settings.prefixName + "-sources.txt";
  if (_settings.useMultiscale) {
    ComponentList* list;
    // If no parallel deconvolution was used, the component list must be
    // retrieved from the deconvolution algorithm.
    if (_algorithms.size() == 1)
      list = &static_cast<MultiScaleAlgorithm*>(_algorithms.front().get())
                  ->GetComponentList();
    else
      list = _componentList.get();
    writeSourceList(*list, filename, phaseCentreRA, phaseCentreDec);
  } else {
    const size_t w = _settings.trimmedImageWidth,
                 h = _settings.trimmedImageHeight;
    ImageSet modelSet(&table, _settings, w, h);
    modelSet.LoadAndAverage(modelImages);
    ComponentList componentList(w, h, modelSet);
    writeSourceList(componentList, filename, phaseCentreRA, phaseCentreDec);
  }
}

void ParallelDeconvolution::correctChannelForPB(
    ComponentList& list, const ImagingTableEntry& entry) const {
  Logger::Debug << "Correcting source list of channel "
                << entry.outputChannelIndex << " for beam\n";
  ImageFilename filename(entry.outputChannelIndex, entry.outputIntervalIndex);
  filename.SetPolarization(entry.polarization);
  PrimaryBeam pb(_settings);
  PrimaryBeamImageSet beam = pb.Load(filename);
  list.CorrectForBeam(beam, entry.outputChannelIndex);
}

void ParallelDeconvolution::SavePBSourceList(CachedImageSet& modelImages,
                                             const ImagingTable& table,
                                             long double phaseCentreRA,
                                             long double phaseCentreDec) const {
  // TODO make this work with subimages
  std::unique_ptr<ComponentList> list;
  const size_t w = _settings.trimmedImageWidth,
               h = _settings.trimmedImageHeight;
  if (_settings.useMultiscale) {
    // If no parallel deconvolution was used, the component list must be
    // retrieved from the deconvolution algorithm.
    if (_algorithms.size() == 1)
      list.reset(new ComponentList(
          static_cast<MultiScaleAlgorithm*>(_algorithms.front().get())
              ->GetComponentList()));
    else
      list.reset(new ComponentList(*_componentList));
  } else {
    ImageSet modelSet(&table, _settings, w, h);
    modelSet.LoadAndAverage(modelImages);
    list.reset(new ComponentList(w, h, modelSet));
  }

  if (_settings.deconvolutionChannelCount == 0 ||
      _settings.deconvolutionChannelCount == table.SquaredGroups().size()) {
    // No beam averaging is required
    for (const ImagingTable::Group& sqGroup : table.SquaredGroups()) {
      correctChannelForPB(*list, *sqGroup.front());
    }
  } else {
    for (size_t ch = 0; ch != _settings.deconvolutionChannelCount; ++ch) {
      Logger::Debug << "Correcting source list of channel " << ch
                    << " for averaged beam\n";
      PrimaryBeamImageSet beamImages = loadAveragePrimaryBeam(ch, table);
      list->CorrectForBeam(beamImages, ch);
    }
  }

  std::string filename = _settings.prefixName + "-sources-pb.txt";
  writeSourceList(*list, filename, phaseCentreRA, phaseCentreDec);
}

void ParallelDeconvolution::writeSourceList(ComponentList& componentList,
                                            const std::string& filename,
                                            long double phaseCentreRA,
                                            long double phaseCentreDec) const {
  if (_settings.useMultiscale) {
    MultiScaleAlgorithm* maxAlgorithm =
        static_cast<MultiScaleAlgorithm*>(_algorithms.front().get());
    for (size_t i = 1; i != _algorithms.size(); ++i) {
      MultiScaleAlgorithm* mAlg =
          static_cast<MultiScaleAlgorithm*>(_algorithms[i].get());
      if (mAlg->ScaleCount() > maxAlgorithm->ScaleCount()) maxAlgorithm = mAlg;
    }
    componentList.Write(filename, *maxAlgorithm, _settings.pixelScaleX,
                        _settings.pixelScaleY, phaseCentreRA, phaseCentreDec);
  } else {
    componentList.WriteSingleScale(filename, *_algorithms.front(),
                                   _settings.pixelScaleX, _settings.pixelScaleY,
                                   phaseCentreRA, phaseCentreDec);
  }
}

PrimaryBeamImageSet ParallelDeconvolution::loadAveragePrimaryBeam(
    size_t imageIndex, const ImagingTable& table) const {
  Logger::Debug << "Averaging beam for deconvolution channel " << imageIndex
                << "\n";

  PrimaryBeamImageSet beamImages;

  Image scratch(_settings.trimmedImageWidth, _settings.trimmedImageHeight);
  size_t deconvolutionChannels = _settings.deconvolutionChannelCount;

  /// TODO : use real weights of images
  size_t count = 0;
  PrimaryBeam pb(_settings);
  const ImagingTable::Groups& squaredGroups = table.SquaredGroups();
  for (size_t sqIndex = 0; sqIndex != squaredGroups.size(); ++sqIndex) {
    size_t curImageIndex =
        (sqIndex * deconvolutionChannels) / squaredGroups.size();
    if (curImageIndex == imageIndex) {
      const ImagingTableEntry& e = *squaredGroups[sqIndex].front();
      Logger::Debug << "Adding beam at " << e.CentralFrequency() * 1e-6
                    << " MHz\n";
      ImageFilename filename(e.outputChannelIndex, e.outputIntervalIndex);

      if (count == 0)
        beamImages = pb.Load(filename);
      else
        beamImages += pb.Load(filename);
      count++;
    }
  }
  beamImages *= (1.0 / double(count));
  return beamImages;
}
