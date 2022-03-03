#include "threadeddeconvolutiontools.h"
#include "multiscaletransforms.h"

#include "../deconvolution/peakfinder.h"
#include "../deconvolution/simpleclean.h"

#include <aocommon/image.h>

#include <memory>

using aocommon::Image;

ThreadedDeconvolutionTools::ThreadedDeconvolutionTools(size_t threadCount)
    : _taskLanes(threadCount),
      _resultLanes(threadCount),
      _threadCount(threadCount) {
  for (size_t i = 0; i != _threadCount; ++i) {
    _taskLanes[i].resize(1);
    _resultLanes[i].resize(1);
    _threadGroup.emplace_back(&ThreadedDeconvolutionTools::threadFunc, this,
                              &_taskLanes[i], &_resultLanes[i]);
  }
}

ThreadedDeconvolutionTools::~ThreadedDeconvolutionTools() {
  for (size_t i = 0; i != _threadCount; ++i) {
    _taskLanes[i].write_end();
  }

  for (std::thread& t : _threadGroup) t.join();
}

void ThreadedDeconvolutionTools::SubtractImage(float* image,
                                               const aocommon::Image& psf,
                                               size_t x, size_t y,
                                               float factor) {
  for (size_t thr = 0; thr != _threadCount; ++thr) {
    auto task = std::make_unique<SubtractionTask>();
    task->image = image;
    task->psf = &psf;
    task->x = x;
    task->y = y;
    task->factor = factor;
    task->startY = psf.Height() * thr / _threadCount;
    task->endY = psf.Height() * (thr + 1) / _threadCount;
    if (thr == _threadCount - 1) {
      (*task)();
    } else {
      _taskLanes[thr].write(std::move(task));
    }
  }
  for (size_t thr = 0; thr != _threadCount - 1; ++thr) {
    std::unique_ptr<ThreadResult> result;
    _resultLanes[thr].read(result);
  }
}

std::unique_ptr<ThreadedDeconvolutionTools::ThreadResult>
ThreadedDeconvolutionTools::SubtractionTask::operator()() {
  SimpleClean::PartialSubtractImage(image, psf->Data(), psf->Width(),
                                    psf->Height(), x, y, factor, startY, endY);
  return std::unique_ptr<ThreadedDeconvolutionTools::ThreadResult>();
}

void ThreadedDeconvolutionTools::FindMultiScalePeak(
    MultiScaleTransforms* msTransforms, const Image& image,
    const aocommon::UVector<float>& scales,
    std::vector<ThreadedDeconvolutionTools::PeakData>& results,
    bool allowNegativeComponents, const bool* mask,
    const std::vector<aocommon::UVector<bool>>& scaleMasks, float borderRatio,
    const Image& rmsFactorImage, bool calculateRMS) {
  size_t imageIndex = 0;
  size_t nextThread = 0;
  size_t resultIndex = 0;

  results.resize(scales.size());

  size_t size = std::min(scales.size(), _threadCount);
  std::vector<Image> imageData;
  std::vector<Image> scratchData;
  imageData.reserve(size);
  scratchData.reserve(size);
  for (size_t i = 0; i != size; ++i) {
    imageData.emplace_back(msTransforms->Width(), msTransforms->Height());
    scratchData.emplace_back(msTransforms->Width(), msTransforms->Height());
  }

  while (imageIndex < scales.size()) {
    std::unique_ptr<FindMultiScalePeakTask> task(new FindMultiScalePeakTask());
    task->msTransforms = msTransforms;
    imageData[nextThread] = image;
    task->image = &imageData[nextThread];
    task->scratch = &scratchData[nextThread];
    task->scale = scales[imageIndex];
    task->allowNegativeComponents = allowNegativeComponents;
    if (scaleMasks.empty())
      task->mask = mask;
    else
      task->mask = scaleMasks[imageIndex].data();
    task->borderRatio = borderRatio;
    task->calculateRMS = calculateRMS;
    task->rmsFactorImage = &rmsFactorImage;
    _taskLanes[nextThread].write(std::move(task));

    ++nextThread;
    if (nextThread == _threadCount) {
      for (size_t thr = 0; thr != nextThread; ++thr) {
        std::unique_ptr<ThreadResult> result;
        _resultLanes[thr].read(result);
        results[resultIndex].normalizedValue =
            static_cast<FindMultiScalePeakResult&>(*result).normalizedValue;
        results[resultIndex].unnormalizedValue =
            static_cast<FindMultiScalePeakResult&>(*result).unnormalizedValue;
        results[resultIndex].x =
            static_cast<FindMultiScalePeakResult&>(*result).x;
        results[resultIndex].y =
            static_cast<FindMultiScalePeakResult&>(*result).y;
        results[resultIndex].rms =
            static_cast<FindMultiScalePeakResult&>(*result).rms;
        ++resultIndex;
      }
      nextThread = 0;
    }
    ++imageIndex;
  }
  for (size_t thr = 0; thr != nextThread; ++thr) {
    std::unique_ptr<ThreadResult> result;
    _resultLanes[thr].read(result);
    results[resultIndex].unnormalizedValue =
        static_cast<FindMultiScalePeakResult&>(*result).unnormalizedValue;
    results[resultIndex].normalizedValue =
        static_cast<FindMultiScalePeakResult&>(*result).normalizedValue;
    results[resultIndex].x = static_cast<FindMultiScalePeakResult&>(*result).x;
    results[resultIndex].y = static_cast<FindMultiScalePeakResult&>(*result).y;
    results[resultIndex].rms =
        static_cast<FindMultiScalePeakResult&>(*result).rms;
    ++resultIndex;
  }
}

std::unique_ptr<ThreadedDeconvolutionTools::ThreadResult>
ThreadedDeconvolutionTools::FindMultiScalePeakTask::operator()() {
  msTransforms->Transform(*image, *scratch, scale);
  const size_t width = msTransforms->Width();
  const size_t height = msTransforms->Height();
  const size_t scaleBorder = size_t(std::ceil(scale * 0.5));
  const size_t horBorderSize =
      std::max<size_t>(std::round(width * borderRatio), scaleBorder);
  const size_t vertBorderSize =
      std::max<size_t>(std::round(height * borderRatio), scaleBorder);
  std::unique_ptr<FindMultiScalePeakResult> result(
      new FindMultiScalePeakResult());
  if (calculateRMS)
    result->rms = RMS(*image, width * height);
  else
    result->rms = -1.0;
  if (rmsFactorImage->Empty()) {
    if (mask == nullptr)
      result->unnormalizedValue = PeakFinder::Find(
          image->Data(), width, height, result->x, result->y,
          allowNegativeComponents, 0, height, horBorderSize, vertBorderSize);
    else
      result->unnormalizedValue =
          PeakFinder::FindWithMask(image->Data(), width, height, result->x,
                                   result->y, allowNegativeComponents, 0,
                                   height, mask, horBorderSize, vertBorderSize);

    result->normalizedValue = result->unnormalizedValue;
  } else {
    for (size_t i = 0; i != rmsFactorImage->Size(); ++i)
      (*scratch)[i] = (*image)[i] * (*rmsFactorImage)[i];

    if (mask == nullptr)
      result->unnormalizedValue = PeakFinder::Find(
          scratch->Data(), width, height, result->x, result->y,
          allowNegativeComponents, 0, height, horBorderSize, vertBorderSize);
    else
      result->unnormalizedValue =
          PeakFinder::FindWithMask(scratch->Data(), width, height, result->x,
                                   result->y, allowNegativeComponents, 0,
                                   height, mask, horBorderSize, vertBorderSize);

    if (result->unnormalizedValue) {
      result->normalizedValue =
          (*result->unnormalizedValue) /
          (*rmsFactorImage)[result->x + result->y * width];
    } else {
      result->normalizedValue.reset();
    }
  }
  return result;
}
