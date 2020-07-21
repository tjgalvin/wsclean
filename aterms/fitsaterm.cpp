#include "fitsaterm.h"

#include "../wsclean/logger.h"

FitsATerm::FitsATerm(size_t nAntenna, const CoordinateSystem& coordinateSystem)
    : FitsATermBase(nAntenna, coordinateSystem) {}

FitsATerm::~FitsATerm() {}

void FitsATerm::OpenTECFiles(const std::vector<std::string>& filenames) {
  _mode = TECMode;
  _readers.reserve(filenames.size());
  for (const std::string& filename : filenames) {
    _readers.emplace_back(filename, true, true);
    if (_readers.back().NFrequencies() != 1)
      throw std::runtime_error(
          "FITS file for TEC A-terms has multiple frequencies in it");
  }
  initializeFromFiles(_readers);
}

void FitsATerm::OpenDiagGainFiles(const std::vector<std::string>& filenames) {
  _mode = DiagonalMode;
  _readers.reserve(filenames.size());
  for (const std::string& filename : filenames) {
    _readers.emplace_back(filename, true, true);
    if (_readers.back().NMatrixElements() != 4)
      throw std::runtime_error(
          "FITS file for diagonal gains did not have 4 matrix elements in it");
  }
  initializeFromFiles(_readers);
}

bool FitsATerm::Calculate(std::complex<float>* buffer, double time,
                          double frequency, size_t, const double*) {
  size_t timeIndex;
  bool requiresRecalculation;
  bool positionChanged = findFilePosition(buffer, time, frequency, timeIndex,
                                          requiresRecalculation);
  if (!positionChanged)
    return false;
  else {
    if (requiresRecalculation) {
      readImages(buffer, timeIndex, frequency);
      storeInCache(frequency, buffer);
    }
    return true;
  }
}

void FitsATerm::readImages(std::complex<float>* buffer, size_t timeIndex,
                           double frequency) {
  const size_t freqIndex =
      round((frequency - _readers.front().FrequencyDimensionStart()) /
            _readers.front().FrequencyDimensionIncr());
  const size_t imgIndex =
      _timesteps[timeIndex].imgIndex * NFrequencies() + freqIndex;
  FitsReader& reader = _readers[_timesteps[timeIndex].readerIndex];
  _scratchA.resize(Resampler().ScratchASize());
  _scratchB.resize(Resampler().ScratchBSize(reader));
  // TODO do this in parallel. Needs to fix Resampler too, as currently it can't
  // run in parallel when a window is used.
  for (size_t antennaIndex = 0; antennaIndex != NAntenna(); ++antennaIndex) {
    // In case there is only one antenna in the measurement set, copy it
    // to all antennas. This is not very efficient as the single image
    // is still read + resampled NAnt times, but it's a border case.
    size_t antennaFileIndex = antennaIndex;
    if (reader.NAntennas() == 1) antennaFileIndex = 0;
    std::complex<float>* antennaBuffer =
        buffer + antennaIndex * Width() * Height() * 4;

    switch (_mode) {
      case TECMode: {
        // TODO When we are in the same timestep but at a different frequency,
        // it would be possible to skip reading and resampling, and immediately
        // call evaluateTEC() with the "scratch" data still there.
        Resampler().ReadAndResample(reader,
                                    antennaFileIndex + imgIndex * NAntenna(),
                                    _scratchA, _scratchB, 1.0);
        evaluateTEC(antennaBuffer, _scratchB.data(), frequency);
      } break;

      case DiagonalMode: {
        for (size_t p = 0; p != 2; ++p) {
          Resampler().ReadAndResample(
              reader, (antennaFileIndex + imgIndex * NAntenna()) * 4 + p * 2,
              _scratchA, _scratchB, 1.0);
          copyToRealPolarization(antennaBuffer, _scratchB.data(), p * 3);

          Resampler().ReadAndResample(
              reader,
              (antennaFileIndex + imgIndex * NAntenna()) * 4 + p * 2 + 1,
              _scratchA, _scratchB, 1.0);
          copyToImaginaryPolarization(antennaBuffer, _scratchB.data(), p * 3);
        }
        setPolarization(antennaBuffer, 1, std::complex<float>(0.0, 0.0));
        setPolarization(antennaBuffer, 2, std::complex<float>(0.0, 0.0));
      } break;
    }
  }
}

void FitsATerm::evaluateTEC(std::complex<float>* dest, const double* source,
                            double frequency) {
  for (size_t pixel = 0; pixel != Width() * Height(); ++pixel) {
    dest[pixel * 4] =
        std::polar(1.0, source[pixel] * -8.44797245e9 / frequency);
    dest[pixel * 4 + 1] = 0.0;
    dest[pixel * 4 + 2] = 0.0;
    dest[pixel * 4 + 3] = dest[pixel * 4];
  }
}

void FitsATerm::copyToRealPolarization(std::complex<float>* dest,
                                       const double* source, size_t polIndex) {
  dest += polIndex;
  for (size_t i = 0; i != Width() * Height(); ++i) {
    dest[i * 4].real(source[i]);
  }
}

void FitsATerm::copyToImaginaryPolarization(std::complex<float>* dest,
                                            const double* source,
                                            size_t polIndex) {
  dest += polIndex;
  for (size_t i = 0; i != Width() * Height(); ++i) {
    dest[i * 4].imag(source[i]);
  }
}

void FitsATerm::setPolarization(std::complex<float>* dest, size_t polIndex,
                                std::complex<float> value) {
  dest += polIndex;
  for (size_t i = 0; i != Width() * Height(); ++i) {
    dest[i * 4] = value;
  }
}
