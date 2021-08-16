#include "moresane.h"

#include <sys/types.h>
#include <sys/wait.h>

#include <unistd.h>

#include <aocommon/fits/fitsreader.h>
#include <aocommon/fits/fitswriter.h>

#include "../math/fftconvolver.h"

#include "../io/logger.h"

void MoreSane::ExecuteMajorIteration(float* dataImage, float* modelImage,
                                     const float* psfImage, size_t width,
                                     size_t height) {
  if (_iterationNumber != 0) {
    Logger::Info << "Convolving model with psf...\n";
    Image preparedPsf(width, height);
    FFTConvolver::PrepareKernel(preparedPsf.data(), psfImage, width, height,
                                _threadCount);
    FFTConvolver::ConvolveSameSize(_fftwManager, modelImage, preparedPsf.data(),
                                   width, height, _threadCount);
    Logger::Info << "Adding model back to residual...\n";
    for (size_t i = 0; i != width * height; ++i) dataImage[i] += modelImage[i];
  }
  std::ostringstream outputStr;
  outputStr << _prefixName << "-tmp-moresaneoutput" << _iterationNumber;
  const std::string dirtyName(_prefixName + "-tmp-moresaneinput-dirty.fits"),
      psfName(_prefixName + "-tmp-moresaneinput-psf.fits"),
      maskName(_prefixName + "-tmp-moresaneinput-mask.fits"),
      outputName(outputStr.str());
  aocommon::FitsWriter writer;
  writer.SetImageDimensions(width, height);
  if (this->_cleanMask != 0) writer.WriteMask(maskName, _cleanMask);
  writer.Write(dirtyName, dataImage);
  writer.Write(psfName, psfImage);

  std::ostringstream commandLine;
  commandLine << "time python \"" << _moresaneLocation << "\" ";
  if (!_allowNegativeComponents) commandLine << "-ep ";
  if (this->_cleanMask != 0) commandLine << "-m \"" << maskName + "\" ";
  if (!_moresaneArguments.empty()) commandLine << _moresaneArguments << ' ';
  commandLine << "\"" << dirtyName << "\" \"" << psfName << "\" \""
              << outputName << '\"';
  if (!_moresaneSigmaLevels.empty()) {
    commandLine << " -sl "
                << _moresaneSigmaLevels[std::min(
                       _iterationNumber, _moresaneSigmaLevels.size() - 1)]
                << " ";
  }

  // TODO should use Application::Run().
  Logger::Info << "Running: " << commandLine.str() << '\n';
  int pid = vfork();
  switch (pid) {
    case -1:  // Error
      throw std::runtime_error(
          "Could not vfork() new process for executing MoreSane");
    case 0:  // Child
      execl("/bin/sh", "sh", "-c", commandLine.str().c_str(), NULL);
      _exit(127);
  }
  // Wait for process to terminate
  int pStatus;
  do {
    int pidReturn;
    do {
      pidReturn = waitpid(pid, &pStatus, 0);
    } while (pidReturn == -1 && errno == EINTR);
  } while (!WIFEXITED(pStatus) && !WIFSIGNALED(pStatus));
  if (WIFEXITED(pStatus)) {
    // all good
    // const int exitStatus = WEXITSTATUS(pStatus);
  } else {
    throw std::runtime_error("MoreSane returned an error");
  }

  aocommon::FitsReader modelReader(outputName + "_model.fits");
  modelReader.Read(modelImage);
  aocommon::FitsReader residualReader(outputName + "_residual.fits");
  residualReader.Read(dataImage);

  unlink(dirtyName.c_str());
  unlink(psfName.c_str());
  unlink(maskName.c_str());
  unlink((outputName + "_model.fits").c_str());
  unlink((outputName + "_residual.fits").c_str());
}

float MoreSane::ExecuteMajorIteration(
    ImageSet& dataImage, ImageSet& modelImage,
    const aocommon::UVector<const float*>& psfImages, size_t width,
    size_t height, bool& reachedMajorThreshold) {
  for (size_t i = 0; i != dataImage.size(); ++i) {
    float* residualData = dataImage[i];
    float* modelData = modelImage[i];
    ExecuteMajorIteration(residualData, modelData,
                          psfImages[dataImage.PSFIndex(i)], width, height);
  }

  ++_iterationNumber;

  reachedMajorThreshold = _iterationNumber < _maxIter;
  return 0.0;
}
