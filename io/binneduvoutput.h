#ifndef BINNED_UV_OUTPUT_H
#define BINNED_UV_OUTPUT_H

#include <schaapcommon/fft/resampler.h>

#include <aocommon/image.h>
#include <aocommon/fits/fitsreader.h>
#include <aocommon/fits/fitswriter.h>
#include <aocommon/logger.h>

#include <schaapcommon/fft/resampler.h>

#include <fstream>

class BinnedUVOutput {
 public:
  static void Make(const std::string &uvCoveragePrefix,
                   const std::string &dirtyPrefix, double psfLimit = 1e-4) {
    aocommon::FitsReader dirtyReader(dirtyPrefix + "-dirty.fits");
    aocommon::FitsReader psfReader(dirtyPrefix + "-psf.fits");
    aocommon::FitsReader uvCoverageReader(uvCoveragePrefix + "-psf.fits");
    const size_t width = dirtyReader.ImageWidth();
    const size_t height = dirtyReader.ImageHeight();
    aocommon::Image dirty(width, height);
    aocommon::Image psf(width, height);
    aocommon::Image uvCovPsf(width, height);
    aocommon::Image realUV(width, height);
    aocommon::Image imagUV(width, height);
    aocommon::Image realPsfUV(width, height);
    aocommon::Image imagPsfUV(width, height);
    aocommon::Image realPsfUVCoverage(width, height);
    aocommon::Image imagPsfUVCoverage(width, height);
    aocommon::Image binned(width, height);
    dirtyReader.Read(dirty.Data());
    psfReader.Read(psf.Data());
    uvCoverageReader.Read(uvCovPsf.Data());

    double nVis, normF, vwSum;
    if (!dirtyReader.ReadDoubleKeyIfExists("WSCNVIS", nVis))
      throw std::runtime_error("Can't find WSCNVIS keyword in fits file");
    if (!dirtyReader.ReadDoubleKeyIfExists("WSCNORMF", normF))
      throw std::runtime_error("Can't find WSCNORMF keyword in fits file");
    if (!dirtyReader.ReadDoubleKeyIfExists("WSCVWSUM", vwSum))
      throw std::runtime_error("Can't find WSCVWSUM keyword in fits file");

    double nVisUVC, normFUVC, vwSumUVC;
    if (!uvCoverageReader.ReadDoubleKeyIfExists("WSCNVIS", nVisUVC))
      throw std::runtime_error("Can't find WSCNVIS keyword in fits file");
    if (!uvCoverageReader.ReadDoubleKeyIfExists("WSCNORMF", normFUVC))
      throw std::runtime_error("Can't find WSCNORMF keyword in fits file");
    if (!uvCoverageReader.ReadDoubleKeyIfExists("WSCVWSUM", vwSumUVC))
      throw std::runtime_error("Can't find WSCVWSUM keyword in fits file");

    // There are two factors of 2 involved: one coming from
    // SingleFT(), and one from the fact that normF excludes a factor of two.
    dirty *= normF * nVis / (2.0 * sqrt(width * height) * vwSum);
    psf *= normF * nVis / (2.0 * sqrt(width * height) * vwSum);
    uvCovPsf *= normFUVC * nVisUVC / (2.0 * sqrt(width * height) * vwSumUVC);
    schaapcommon::fft::Resampler fft(width, height, width, height, 1);
    fft.SingleFT(dirty.Data(), realUV.Data(), imagUV.Data());
    fft.SingleFT(psf.Data(), realPsfUV.Data(), imagPsfUV.Data());
    fft.SingleFT(uvCovPsf.Data(), realPsfUVCoverage.Data(),
                 imagPsfUVCoverage.Data());
    for (size_t i = 0; i != dirty.Size(); ++i) {
      std::complex<double> psfVal =
          std::complex<double>(realPsfUV[i], imagPsfUV[i]);
      std::complex<double> newVal =
          std::complex<double>(realUV[i], imagUV[i]) / psfVal;
      realUV[i] = newVal.real();
      imagUV[i] = newVal.imag();
    }
    binned = 0.0;
    const double pixelSizeX = dirtyReader.PixelSizeX(),
                 pixelSizeY = dirtyReader.PixelSizeY();
    std::ofstream file(dirtyPrefix + "-binneduvoutput.csv");
    file << "# u (lambda), v (lambda), real (Jy), imaginary (Jy), effective "
            "vis count, weight\n";
    file.precision(15);
    size_t nBins = 0;
    double effVisSum = 0.0;
    for (size_t y = 0; y != height; ++y) {
      const float *realPsfPtr = &realPsfUV[y * width];
      const float *imagPsfPtr = &imagPsfUV[y * width];
      const float *realPsfUVCovPtr = &realPsfUVCoverage[y * width];
      const float *imagPsfUVCovPtr = &imagPsfUVCoverage[y * width];
      float *realPtr = &realUV[y * width];
      float *imagPtr = &imagUV[y * width];
      float *binnedPtr = &binned[y * width];
      double vInLambda =
          double(int(height) / 2 - int(y)) / (pixelSizeY * height);
      for (size_t x = 0; x != width; ++x) {
        double uvCovVal = std::abs(
            std::complex<double>(realPsfUVCovPtr[x], imagPsfUVCovPtr[x]));
        if (uvCovVal >= psfLimit) {
          binnedPtr[x] = std::abs(std::complex<double>(realPtr[x], imagPtr[x]));
          ++nBins;

          double uInLambda =
              double(int(width) / 2 - int(x)) / (pixelSizeX * width);
          double psfVal =
              std::abs(std::complex<double>(realPsfPtr[x], imagPsfPtr[x]));
          file << uInLambda << " " << vInLambda << " " << realPtr[x] << " "
               << imagPtr[x] << " " << psfVal << " " << psfVal * vwSum / nVis
               << '\n';
          effVisSum += psfVal;
        } else {
          binnedPtr[x] = 0.0;
          realPtr[x] = 0.0;
          imagPtr[x] = 0.0;
        }
      }
    }
    aocommon::FitsWriter writer(dirtyReader);
    writer.Write(dirtyPrefix + "-binneduvoutput.fits", binned.Data());
    aocommon::Logger::Info << "UV bins were written. " << nBins
                           << " bins were selected, with an effective nvis of "
                           << effVisSum << " (actual nvis=" << nVis << ").\n";
  }
};

#endif
