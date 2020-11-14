#include "fftconvolver.h"

#include "../system/fftwmanager.h"

#include <aocommon/uvector.h>

#include <fftw3.h>

#include <complex>
#include <stdexcept>

void FFTConvolver::Convolve(FFTWManager& fftw, float* image, size_t imgWidth,
                            size_t imgHeight, const float* kernel,
                            size_t kernelSize) {
  aocommon::UVector<float> scaledKernel(imgWidth * imgHeight, 0.0);
  PrepareSmallKernel(scaledKernel.data(), imgWidth, imgHeight, kernel,
                     kernelSize);
  ConvolveSameSize(fftw, image, scaledKernel.data(), imgWidth, imgHeight);
}

void FFTConvolver::ReverseAndConvolve(class FFTWManager& fftw, float* image,
                                      size_t imgWidth, size_t imgHeight,
                                      const float* kernel, size_t kernelSize) {
  aocommon::UVector<float> scaledKernel(imgWidth * imgHeight, 0.0);

  PrepareSmallKernel(scaledKernel.data(), imgWidth, imgHeight, kernel,
                     kernelSize);
  ConvolveSameSize(fftw, image, scaledKernel.data(), imgWidth, imgHeight);
}

void FFTConvolver::PrepareSmallKernel(float* dest, size_t imgWidth,
                                      size_t imgHeight, const float* kernel,
                                      size_t kernelSize) {
  if (kernelSize > imgWidth || kernelSize > imgHeight)
    throw std::runtime_error("Kernel size > image dimension");
  const float* kernelIter = kernel;
  for (size_t y = 0; y != kernelSize / 2; ++y) {
    size_t destY = imgHeight - kernelSize / 2 + y;
    size_t firstX = imgWidth - kernelSize / 2;
    float* destIter = &dest[destY * imgWidth + firstX];
    for (size_t x = 0; x != kernelSize / 2; ++x) {
      *destIter = *kernelIter;
      ++kernelIter;
      ++destIter;
    }
    destIter = &dest[destY * imgWidth];
    for (size_t x = kernelSize / 2; x != kernelSize; ++x) {
      *destIter = *kernelIter;
      ++kernelIter;
      ++destIter;
    }
  }
  for (size_t y = kernelSize / 2; y != kernelSize; ++y) {
    size_t firstX = imgWidth - kernelSize / 2;
    float* destIter = &dest[firstX + (y - kernelSize / 2) * imgWidth];
    for (size_t x = 0; x != kernelSize / 2; ++x) {
      *destIter = *kernelIter;
      ++kernelIter;
      ++destIter;
    }
    destIter = &dest[(y - kernelSize / 2) * imgWidth];
    for (size_t x = kernelSize / 2; x != kernelSize; ++x) {
      *destIter = *kernelIter;
      ++kernelIter;
      ++destIter;
    }
  }
}

void FFTConvolver::PrepareKernel(float* dest, const float* source,
                                 size_t imgWidth, size_t imgHeight) {
  const float* sourceIter = source;
  for (size_t y = 0; y != imgHeight / 2; ++y) {
    size_t destY = imgHeight - imgHeight / 2 + y;
    size_t firstX = imgWidth - imgWidth / 2;
    float* destIter = &dest[destY * imgWidth + firstX];
    for (size_t x = 0; x != imgWidth / 2; ++x) {
      *destIter = *sourceIter;
      ++sourceIter;
      ++destIter;
    }
    destIter = &dest[destY * imgWidth];
    for (size_t x = imgWidth / 2; x != imgWidth; ++x) {
      *destIter = *sourceIter;
      ++sourceIter;
      ++destIter;
    }
  }
  for (size_t y = imgHeight / 2; y != imgHeight; ++y) {
    size_t firstX = imgWidth - imgWidth / 2;
    float* destIter = &dest[firstX + (y - imgHeight / 2) * imgWidth];
    for (size_t x = 0; x != imgWidth / 2; ++x) {
      *destIter = *sourceIter;
      ++sourceIter;
      ++destIter;
    }
    destIter = &dest[(y - imgHeight / 2) * imgWidth];
    for (size_t x = imgWidth / 2; x != imgWidth; ++x) {
      *destIter = *sourceIter;
      ++sourceIter;
      ++destIter;
    }
  }
}

void FFTConvolver::ConvolveSameSize(FFTWManager& fftw, float* image,
                                    const float* kernel, size_t imgWidth,
                                    size_t imgHeight) {
  const size_t imgSize = imgWidth * imgHeight;
  const size_t complexSize = (imgWidth / 2 + 1) * imgHeight;
  float* tempData =
      reinterpret_cast<float*>(fftwf_malloc(imgSize * sizeof(float)));
  fftwf_complex* fftImageData = reinterpret_cast<fftwf_complex*>(
      fftwf_malloc(complexSize * sizeof(fftwf_complex)));
  fftwf_complex* fftKernelData = reinterpret_cast<fftwf_complex*>(
      fftwf_malloc(complexSize * sizeof(fftwf_complex)));

  std::unique_lock<std::mutex> lock(fftw.Mutex());
  fftwf_plan inToFPlan = fftwf_plan_dft_r2c_2d(imgHeight, imgWidth, tempData,
                                               fftImageData, FFTW_ESTIMATE);
  fftwf_plan fToOutPlan = fftwf_plan_dft_c2r_2d(
      imgHeight, imgWidth, fftImageData, tempData, FFTW_ESTIMATE);
  lock.unlock();

  fftwf_execute_dft_r2c(inToFPlan, image, fftImageData);

  std::copy_n(kernel, imgSize, tempData);
  fftwf_execute_dft_r2c(inToFPlan, tempData, fftKernelData);

  float fact = 1.0 / imgSize;
  for (size_t i = 0; i != complexSize; ++i)
    reinterpret_cast<std::complex<float>*>(fftImageData)[i] *=
        fact * reinterpret_cast<std::complex<float>*>(fftKernelData)[i];

  fftwf_execute_dft_c2r(fToOutPlan,
                        reinterpret_cast<fftwf_complex*>(fftImageData), image);

  fftwf_free(fftImageData);
  fftwf_free(fftKernelData);
  fftwf_free(tempData);

  lock.lock();
  fftwf_destroy_plan(inToFPlan);
  fftwf_destroy_plan(fToOutPlan);
  lock.unlock();
}

void FFTConvolver::Reverse(float* image, size_t imgWidth, size_t imgHeight) {
  for (size_t y = 0; y != imgHeight / 2; ++y) {
    size_t destY = imgHeight - 1 - y;
    float* sourcePtr = &image[y * imgWidth];
    float* destPtr = &image[destY * imgWidth];
    for (size_t x = 0; x != imgWidth / 2; ++x) {
      size_t destX = imgWidth - 1 - x;
      std::swap(sourcePtr[x], destPtr[destX]);
    }
  }
}
