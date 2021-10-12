#ifndef CALCULATE_FFT_SIZE_H
#define CALCULATE_FFT_SIZE_H

namespace detail {
inline bool hasLowPrimeFactors(size_t number) {
  while (number > 7) {
    if (number % 2 == 0)
      number /= 2;
    else if (number % 3 == 0)
      number /= 3;
    else if (number % 5 == 0)
      number /= 5;
    else if (number % 7 == 0)
      number /= 7;
    else
      return false;
  }
  return true;
}
}  // namespace detail

inline void CalculateFFTSize(size_t size, double pixelScale, double beamSize,
                             size_t& minimalSize, size_t& optimalSize) {
  double totalSize = double(size) * pixelScale;
  // Calc min res based on Nyquist sampling rate
  minimalSize = size_t(ceil(totalSize * 2 / beamSize));
  optimalSize = minimalSize;
  if (optimalSize % 4 != 0) optimalSize += 4 - (optimalSize % 4);
  while (!detail::hasLowPrimeFactors(optimalSize) && optimalSize < size)
    optimalSize += 4;
  if (optimalSize > size) optimalSize = size;
}

#endif
