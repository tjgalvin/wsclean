#include "pythondeconvolution.h"

#include <pybind11/attr.h>
#include <pybind11/embed.h>
#include <pybind11/eval.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

struct PyChannel {
  double frequency, weight;
};

class PySpectralFitter {
 public:
  PySpectralFitter(SpectralFitter& fitter) : _fitter(fitter) {}

  pybind11::array_t<double> fit(pybind11::array_t<double> values, size_t x,
                                size_t y) {
    if (values.ndim() != 1)
      throw std::runtime_error(
          "spectral_fitter.fit(): Invalid dimensions of values array");
    if (size_t(values.shape()[0]) != _fitter.NFrequencies())
      throw std::runtime_error(
          "spectral_fitter.fit(): Incorrect size of values array");
    aocommon::UVector<float> vec(_fitter.NFrequencies());
    pybind11::buffer_info info = values.request();
    const unsigned char* buffer = static_cast<const unsigned char*>(info.ptr);
    for (size_t i = 0; i != _fitter.NFrequencies(); ++i) {
      vec[i] = *reinterpret_cast<const double*>(buffer + info.strides[0] * i);
    }
    aocommon::UVector<float> result;
    _fitter.Fit(result, vec.data(), x, y);

    pybind11::buffer_info resultBuf(
        nullptr,  // ask NumPy to allocate
        sizeof(double), pybind11::format_descriptor<double>::value, 1,
        {ptrdiff_t(_fitter.NTerms())}, {sizeof(double)}  // Stride
    );
    pybind11::array_t<double> pyResult(resultBuf);
    std::copy_n(result.data(), _fitter.NTerms(),
                static_cast<double*>(pyResult.request(true).ptr));
    return pyResult;
  }

  pybind11::array_t<double> fit_and_evaluate(pybind11::array_t<double> values,
                                             size_t x, size_t y) {
    if (values.ndim() != 1)
      throw std::runtime_error(
          "spectral_fitter.fit_and_evaluate(): Invalid dimensions of values "
          "array");
    if (size_t(values.shape()[0]) != _fitter.NFrequencies())
      throw std::runtime_error(
          "spectral_fitter.fit_and_evaluate(): Incorrect size of values array");
    aocommon::UVector<float> vec(_fitter.NFrequencies());
    pybind11::buffer_info info = values.request();
    const unsigned char* buffer = static_cast<const unsigned char*>(info.ptr);
    for (size_t i = 0; i != _fitter.NFrequencies(); ++i) {
      vec[i] = *reinterpret_cast<const double*>(buffer + info.strides[0] * i);
    }

    aocommon::UVector<float> fittingScratch;
    _fitter.FitAndEvaluate(vec.data(), x, y, fittingScratch);

    pybind11::buffer_info resultBuf(
        nullptr,  // ask NumPy to allocate
        sizeof(double), pybind11::format_descriptor<double>::value, 1,
        {ptrdiff_t(_fitter.NFrequencies())}, {sizeof(double)}  // Stride
    );
    pybind11::array_t<double> pyResult(resultBuf);
    std::copy_n(vec.data(), _fitter.NFrequencies(),
                static_cast<double*>(pyResult.request(true).ptr));
    return pyResult;
  }

 private:
  SpectralFitter& _fitter;
};

struct PyMetaData {
 public:
  PyMetaData(SpectralFitter& _spectral_fitter)
      : spectral_fitter(_spectral_fitter) {}

  std::vector<PyChannel> channels;
  size_t iteration_number;
  double final_threshold;
  double gain;
  size_t max_iterations;
  double major_iter_threshold;
  double mgain;
  PySpectralFitter spectral_fitter;
  bool square_joined_channels;
};

PythonDeconvolution::PythonDeconvolution(const std::string& filename)
    : _filename(filename), _guard(new pybind11::scoped_interpreter()) {
  pybind11::module main = pybind11::module::import("__main__");
  pybind11::object scope = main.attr("__dict__");
  pybind11::eval_file(_filename, scope);
  _deconvolveFunction = main.attr("deconvolve").cast<pybind11::function>();

  pybind11::class_<PyChannel>(main, "Channel")
      .def_readwrite("frequency", &PyChannel::frequency)
      .def_readwrite("weight", &PyChannel::weight);

  pybind11::class_<PyMetaData>(main, "MetaData")
      .def_readonly("channels", &PyMetaData::channels)
      .def_readonly("final_threshold", &PyMetaData::final_threshold)
      .def_readwrite("iteration_number", &PyMetaData::iteration_number)
      .def_readonly("gain", &PyMetaData::gain)
      .def_readonly("max_iterations", &PyMetaData::max_iterations)
      .def_readonly("major_iter_threshold", &PyMetaData::major_iter_threshold)
      .def_readonly("mgain", &PyMetaData::mgain)
      .def_readonly("spectral_fitter", &PyMetaData::spectral_fitter)
      .def_readonly("square_joined_channels",
                    &PyMetaData::square_joined_channels);

  pybind11::class_<PySpectralFitter>(main, "SpectralFitter")
      .def("fit", &PySpectralFitter::fit)
      .def("fit_and_evaluate", &PySpectralFitter::fit_and_evaluate);
}

void PythonDeconvolution::setBuffer(const ImageSet& imageSet, double* ptr) {
  size_t nFreq = imageSet.NDeconvolutionChannels();
  size_t nPol = imageSet.size() / imageSet.NDeconvolutionChannels();

  for (size_t freq = 0; freq != nFreq; ++freq) {
    for (size_t pol = 0; pol != nPol; ++pol) {
      const aocommon::Image& image = imageSet[freq * nPol + pol];
      std::copy_n(image.Data(), image.Size(), ptr);
      ptr += image.Size();
    }
  }
}

void PythonDeconvolution::getBuffer(ImageSet& imageSet, const double* ptr) {
  size_t nFreq = imageSet.NDeconvolutionChannels();
  size_t nPol = imageSet.size() / imageSet.NDeconvolutionChannels();

  for (size_t freq = 0; freq != nFreq; ++freq) {
    for (size_t pol = 0; pol != nPol; ++pol) {
      const size_t imageIndex = freq * nPol + pol;
      float* img = imageSet.Data(imageIndex);
      const size_t imageSize = imageSet[imageIndex].Size();
      std::copy_n(ptr, imageSize, img);
      ptr += imageSize;
      ;
    }
  }
}

void PythonDeconvolution::setPsf(const std::vector<aocommon::Image>& psfs,
                                 double* pyPtr, size_t width, size_t height) {
  size_t nFreq = psfs.size();

  for (size_t freq = 0; freq != nFreq; ++freq) {
    const float* psf = psfs[freq].Data();

    for (size_t y = 0; y != height; ++y) {
      for (size_t x = 0; x != width; ++x) pyPtr[x] = psf[x];

      pyPtr += width;
      psf += width;
    }
  }
}

float PythonDeconvolution::ExecuteMajorIteration(
    ImageSet& dirtySet, ImageSet& modelSet,
    const std::vector<aocommon::Image>& psfs, bool& reachedMajorThreshold) {
  const size_t width = dirtySet.Width();
  const size_t height = dirtySet.Height();
  size_t nFreq = dirtySet.NDeconvolutionChannels();
  size_t nPol = dirtySet.size() / dirtySet.NDeconvolutionChannels();

  pybind11::object result;

  // A new context block is started to destroy the python data arrays asap
  {
    // Create Residual array
    pybind11::buffer_info residualBuf(
        nullptr,  // ask NumPy to allocate
        sizeof(double), pybind11::format_descriptor<double>::value, 4,
        {ptrdiff_t(nFreq), ptrdiff_t(nPol), ptrdiff_t(height),
         ptrdiff_t(width)},
        {sizeof(double) * width * height * nPol,
         sizeof(double) * width * height, sizeof(double) * width,
         sizeof(double)}  // Strides
    );
    pybind11::array_t<double> pyResiduals(residualBuf);
    setBuffer(dirtySet, static_cast<double*>(pyResiduals.request(true).ptr));

    // Create Model array
    pybind11::buffer_info modelBuf(
        nullptr, sizeof(double), pybind11::format_descriptor<double>::value, 4,
        {ptrdiff_t(nFreq), ptrdiff_t(nPol), ptrdiff_t(height),
         ptrdiff_t(width)},
        {sizeof(double) * width * height * nPol,
         sizeof(double) * width * height, sizeof(double) * width,
         sizeof(double)});
    pybind11::array_t<double> pyModel(modelBuf);
    setBuffer(modelSet, static_cast<double*>(pyModel.request(true).ptr));

    // Create PSF array
    pybind11::buffer_info psfBuf(
        nullptr, sizeof(double), pybind11::format_descriptor<double>::value, 3,
        {ptrdiff_t(nFreq), ptrdiff_t(height), ptrdiff_t(width)},
        {sizeof(double) * width * height, sizeof(double) * width,
         sizeof(double)});
    pybind11::array_t<double> pyPsfs(psfBuf);
    setPsf(psfs, static_cast<double*>(pyPsfs.request(true).ptr), width, height);

    PyMetaData meta(_spectralFitter);
    meta.channels.resize(_spectralFitter.NFrequencies());
    for (size_t i = 0; i != _spectralFitter.NFrequencies(); ++i) {
      meta.channels[i].frequency = _spectralFitter.Frequency(i);
      meta.channels[i].weight = _spectralFitter.Weight(i);
    }
    meta.gain = _gain;
    meta.iteration_number = _iterationNumber;
    meta.major_iter_threshold = _majorIterThreshold;
    meta.max_iterations = _maxIter;
    meta.mgain = _mGain;
    meta.final_threshold = _threshold;

    // Run the python code
    result = _deconvolveFunction(std::move(pyResiduals), std::move(pyModel),
                                 std::move(pyPsfs), &meta);

    _iterationNumber = meta.iteration_number;
  }

  // Extract the results
  pybind11::object resultDict;
  try {
    resultDict = result.cast<pybind11::dict>();
  } catch (std::exception&) {
    throw std::runtime_error(
        "In python deconvolution code: Return value of deconvolve() should be "
        "a dictionary");
  }
  const bool isComplete =
      resultDict.contains("residual") && resultDict.contains("model") &&
      resultDict.contains("level") && resultDict.contains("continue");
  if (!isComplete)
    throw std::runtime_error(
        "In python deconvolution code: Dictionary returned by deconvolve() is "
        "missing items; should have 'residual', 'model', 'level' and "
        "'continue'");
  pybind11::array_t<double> residualRes =
      resultDict["residual"].cast<pybind11::array_t<double>>();
  getBuffer(dirtySet, static_cast<const double*>(residualRes.request().ptr));
  pybind11::array_t<double> modelRes =
      resultDict["model"].cast<pybind11::array_t<double>>();
  getBuffer(modelSet, static_cast<const double*>(modelRes.request().ptr));

  double level = resultDict["level"].cast<double>();
  reachedMajorThreshold = resultDict["continue"].cast<bool>();
  return level;
}
