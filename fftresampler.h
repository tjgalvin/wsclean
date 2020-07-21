#ifndef FFT_RESAMPLE_H
#define FFT_RESAMPLE_H

#include "windowfunction.h"

#include <aocommon/lane.h>
#include <aocommon/uvector.h>

#include <vector>
#include <thread>

#include <fftw3.h>

class FFTResampler {
 private:
  struct Task {
    double *input, *output;
  };

 public:
  FFTResampler(size_t inWidth, size_t inHeight, size_t outWidth,
               size_t outHeight, size_t cpuCount, bool verbose = true);

  ~FFTResampler();

  void AddTask(double* input, double* output) {
    Task task;
    task.input = input;
    task.output = output;
    _tasks.write(task);
  }

  void Start() {
    for (size_t i = 0; i != _tasks.capacity(); ++i) {
      _threads.emplace_back(&FFTResampler::runThread, this);
    }
  }

  void Finish() {
    _tasks.write_end();
    for (std::thread& t : _threads) t.join();
    _threads.clear();
    _tasks.clear();
  }

  void Resample(double* input, double* output) {
    Task task;
    task.input = input;
    task.output = output;
    runSingle(task, false);
  }

  void SingleFT(const double* input, double* realOutput,
                double* imaginaryOutput);

  /**
   * Only to be used with SingleFT (it makes resampling thread unsafe!)
   */
  void SetTukeyWindow(double insetSize, bool correctWindow) {
    _windowFunction = WindowFunction::Tukey;
    _tukeyInsetSize = insetSize;
    _correctWindow = correctWindow;
    _windowRowIn.clear();
    _windowColIn.clear();
    _windowOut.clear();
  }

  void SetWindowFunction(WindowFunction::Type window, bool correctWindow) {
    _windowFunction = window;
    _correctWindow = correctWindow;
    _windowRowIn.clear();
    _windowColIn.clear();
    _windowOut.clear();
  }

 private:
  void runThread();
  void runSingle(const Task& task, bool skipWindow) const;
  void applyWindow(double* data) const;
  void unapplyWindow(double* data) const;
  void makeWindow(aocommon::UVector<double>& data, size_t width) const;
  void makeTukeyWindow(aocommon::UVector<double>& data, size_t width) const;

  size_t _inputWidth, _inputHeight;
  size_t _outputWidth, _outputHeight;
  size_t _fftWidth, _fftHeight;
  WindowFunction::Type _windowFunction;
  double _tukeyInsetSize;
  mutable aocommon::UVector<double> _windowRowIn;
  mutable aocommon::UVector<double> _windowColIn;
  mutable aocommon::UVector<double> _windowOut;
  bool _correctWindow;

  fftw_plan _inToFPlan, _fToOutPlan;

  aocommon::Lane<Task> _tasks;
  std::vector<std::thread> _threads;
  bool _verbose;
};

#endif
