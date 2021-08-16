#ifndef THREADED_DECONVOLUTION_TOOLS_H
#define THREADED_DECONVOLUTION_TOOLS_H

#include "../structures/image.h"

#include <boost/optional/optional.hpp>

#include <aocommon/lane.h>
#include <aocommon/uvector.h>

#include <cmath>
#include <thread>
#include <vector>

class ThreadedDeconvolutionTools {
 public:
  explicit ThreadedDeconvolutionTools(size_t threadCount);
  ~ThreadedDeconvolutionTools();

  struct PeakData {
    boost::optional<float> normalizedValue, unnormalizedValue;
    float rms;
    size_t x, y;
  };

  void SubtractImage(float* image, const float* psf, size_t width,
                     size_t height, size_t x, size_t y, float factor);

  void FindMultiScalePeak(
      class MultiScaleTransforms* msTransforms, const Image& image,
      const aocommon::UVector<float>& scales, std::vector<PeakData>& results,
      bool allowNegativeComponents, const bool* mask,
      const std::vector<aocommon::UVector<bool>>& scaleMasks, float borderRatio,
      const Image& rmsFactorImage, bool calculateRMS);

  static float RMS(const Image& image, size_t n) {
    float result = 0.0;
    for (size_t i = 0; i != n; ++i) result += image[i] * image[i];
    return std::sqrt(result / float(n));
  }

 private:
  struct ThreadResult {};
  struct FindMultiScalePeakResult : public ThreadResult {
    boost::optional<float> unnormalizedValue, normalizedValue;
    float rms;
    size_t x, y;
  };

  struct ThreadTask {
    virtual ThreadResult* operator()() = 0;
    virtual ~ThreadTask() {}
  };
  struct SubtractionTask : public ThreadTask {
    virtual ThreadResult* operator()();

    float* image;
    const float* psf;
    size_t width, height, x, y;
    float factor;
    size_t startY, endY;
  };

  struct FindMultiScalePeakTask : public ThreadTask {
    virtual ThreadResult* operator()();

    class MultiScaleTransforms* msTransforms;
    Image* image;
    Image* scratch;
    float scale;
    bool allowNegativeComponents;
    const bool* mask;
    float borderRatio;
    bool calculateRMS;
    const Image* rmsFactorImage;
  };

  std::vector<aocommon::Lane<ThreadTask*>*> _taskLanes;
  std::vector<aocommon::Lane<ThreadResult*>*> _resultLanes;
  size_t _threadCount;
  std::vector<std::thread> _threadGroup;

  void threadFunc(aocommon::Lane<ThreadTask*>* taskLane,
                  aocommon::Lane<ThreadResult*>* resultLane) {
    ThreadTask* task;
    while (taskLane->read(task)) {
      ThreadResult* result = (*task)();
      resultLane->write(result);
      delete task;
    }
  }
};

#endif
