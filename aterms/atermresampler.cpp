#include "atermresampler.h"

#include "../fftresampler.h"

#include <aocommon/imagecoordinates.h>

ATermResampler::ATermResampler(
    const ATermBase::CoordinateSystem& coordinateSystem)
    : _coordinateSystem(coordinateSystem),
      _allocatedWidth(coordinateSystem.maxSupport),
      _allocatedHeight(coordinateSystem.maxSupport),
      _downsample(true),
      _window(WindowFunction::RaisedHann),
      _padding(1.0),
      _overrideFitsPhaseCentre(false),
      _overrideRa(0.0),
      _overrideDec(0.0) {}

ATermResampler::~ATermResampler() {}

void ATermResampler::ReadAndResample(FitsReader& reader, size_t fileIndex,
                                     aocommon::UVector<double>& scratch,
                                     aocommon::UVector<double>& output,
                                     double stretchFactor) {
  if (_resampler == nullptr) {
    _resampler.reset(new FFTResampler(_allocatedWidth, _allocatedHeight,
                                      _coordinateSystem.width,
                                      _coordinateSystem.height, 1, false));
    if (_window == WindowFunction::Tukey)
      _resampler->SetTukeyWindow(double(_allocatedWidth) / _padding, false);
    else
      _resampler->SetWindowFunction(_window, true);
  }

  if (_downsample) {
    reader.ReadIndex(output.data(), fileIndex);

    // First, the image is regridded on a smaller image that fits in the kernel
    // support allocated for the aterms
    regrid(reader, scratch.data(), output.data(), stretchFactor);

    // Now, the small image is enlarged so that it matches the kernel size
    _resampler->Resample(scratch.data(), output.data());
  } else {
    scratch.resize(reader.ImageWidth() * reader.ImageHeight());
    reader.ReadIndex(scratch.data(), fileIndex);

    regrid(reader, output.data(), scratch.data(), stretchFactor);
  }
}

void ATermResampler::regrid(const FitsReader& reader, double* dest,
                            const double* source, double stretchFactor) {
  size_t inWidth = reader.ImageWidth(), inHeight = reader.ImageHeight();
  const double inPixelSizeX = reader.PixelSizeX() / stretchFactor,
               inPixelSizeY = reader.PixelSizeY() / stretchFactor,
               inPhaseCentreDL = reader.PhaseCentreDL(),
               inPhaseCentreDM = reader.PhaseCentreDM();
  double inPhaseCentreRA, inPhaseCentreDec;
  if (_overrideFitsPhaseCentre) {
    inPhaseCentreRA = _overrideRa;
    inPhaseCentreDec = _overrideDec;
  } else {
    inPhaseCentreRA = reader.PhaseCentreRA();
    inPhaseCentreDec = reader.PhaseCentreDec();
  }

  size_t outWidth, outHeight;
  if (_downsample) {
    outWidth = _allocatedWidth;
    outHeight = _allocatedHeight;
  } else {
    outWidth = _coordinateSystem.width;
    outHeight = _coordinateSystem.height;
  }
  // The full size is regridded onto the 'Nyquist-sampled' image to remove
  // high-frequency components. atermDL/DM are the pixelsizes of the
  // Nyquist-sampled image.
  double atermDL = _coordinateSystem.dl * _coordinateSystem.width / outWidth;
  double atermDM = _coordinateSystem.dm * _coordinateSystem.height / outHeight;
  /**
   * If phase centra of input and output are the same, i.e. they have the same
   * tangential plane, a few calculations can be saved.
   */
  bool samePlane = inPhaseCentreRA == _coordinateSystem.ra &&
                   inPhaseCentreDec == _coordinateSystem.dec;

  size_t index = 0;
  for (size_t y = 0; y != outWidth; ++y) {
    for (size_t x = 0; x != outWidth; ++x) {
      double l, m;
      aocommon::ImageCoordinates::XYToLM(x, y, atermDL, atermDM, outWidth,
                                         outWidth, l, m);
      l += _coordinateSystem.phaseCentreDL;
      m += _coordinateSystem.phaseCentreDM;
      if (!samePlane) {
        double pixra, pixdec;
        aocommon::ImageCoordinates::LMToRaDec(
            l, m, _coordinateSystem.ra, _coordinateSystem.dec, pixra, pixdec);
        aocommon::ImageCoordinates::RaDecToLM(pixra, pixdec, inPhaseCentreRA,
                                              inPhaseCentreDec, l, m);
      }
      l -= inPhaseCentreDL;
      m -= inPhaseCentreDM;
      int inX, inY;
      aocommon::ImageCoordinates::LMToXY(l, m, inPixelSizeX, inPixelSizeY,
                                         inWidth, inHeight, inX, inY);
      if (inX < 0 || inY < 0 || inX >= int(inWidth) || inY >= int(inHeight))
        dest[index] = 0;
      else {
        dest[index] = source[inX + inY * inWidth];
      }
      ++index;
    }
  }
}
