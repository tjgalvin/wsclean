#ifndef THREADED_DECONVOLUTION_TOOLS_H
#define THREADED_DECONVOLUTION_TOOLS_H

#include "../image.h"

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
    boost::optional<double> normalizedValue, unnormalizedValue;
    double rms;
    size_t x, y;
  };

  void SubtractImage(double* image, const double* psf, size_t width,
                     size_t height, size_t x, size_t y, double factor);

  // This one is for many transforms of the same scale
  void MultiScaleTransform(class MultiScaleTransforms* msTransforms,
                           const aocommon::UVector<double*>& images,
                           double* scratch, double scale);

  // This one is for transform of different scales
  void MultiScaleTransform(class MultiScaleTransforms* msTransforms,
                           const aocommon::UVector<double*>& images,
                           aocommon::UVector<double> scales);

  void FindMultiScalePeak(
      class MultiScaleTransforms* msTransforms, const double* image,
      const aocommon::UVector<double>& scales, std::vector<PeakData>& results,
      bool allowNegativeComponents, const bool* mask,
      const std::vector<aocommon::UVector<bool>>& scaleMasks,
      double borderRatio, const Image& rmsFactorImage, bool calculateRMS);

  static double RMS(const double* image, size_t n) {
    double result = 0.0;
    for (size_t i = 0; i != n; ++i) result += image[i] * image[i];
    return std::sqrt(result / double(n));
  }

 private:
  struct ThreadResult {};
  struct FindMultiScalePeakResult : public ThreadResult {
    boost::optional<double> unnormalizedValue, normalizedValue;
    double rms;
    size_t x, y;
  };

  struct ThreadTask {
    virtual ThreadResult* operator()() = 0;
    virtual ~ThreadTask() {}
  };
  struct SubtractionTask : public ThreadTask {
    virtual ThreadResult* operator()();

    double* image;
    const double* psf;
    size_t width, height, x, y;
    double factor;
    size_t startY, endY;
  };
  struct FinishMultiScaleTransformTask : public ThreadTask {
    virtual ThreadResult* operator()();

    class MultiScaleTransforms* msTransforms;
    double* image;
    double* kernel;
  };
  struct MultiScaleTransformTask : public ThreadTask {
    virtual ThreadResult* operator()();

    class MultiScaleTransforms* msTransforms;
    double* image;
    double* scratch;
    double scale;
  };
  struct FindMultiScalePeakTask : public ThreadTask {
    virtual ThreadResult* operator()();

    class MultiScaleTransforms* msTransforms;
    double* image;
    double* scratch;
    double scale;
    bool allowNegativeComponents;
    const bool* mask;
    double borderRatio;
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
