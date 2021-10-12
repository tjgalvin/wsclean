#include "msprovider.h"
#include "msreaders/msreader.h"

#include "../io/logger.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/ArrColDesc.h>

#include "../structures/msselection.h"
#include "../structures/multibanddata.h"

namespace {
template <bool add>
void AddOrAssign(std::complex<float>* dest, std::complex<float> source) {
  *dest += source;
}

template <>
void AddOrAssign<false>(std::complex<float>* dest, std::complex<float> source) {
  *dest = source;
}
}  // namespace

MSProvider::~MSProvider() {}

void MSProvider::CopyData(std::complex<float>* dest, size_t startChannel,
                          size_t endChannel,
                          const std::vector<aocommon::PolarizationEnum>& polsIn,
                          const casacore::Array<std::complex<float>>& data,
                          aocommon::PolarizationEnum polOut) {
  const size_t polCount = polsIn.size();
  casacore::Array<std::complex<float>>::const_contiter inPtr =
      data.cbegin() + startChannel * polCount;
  const size_t selectedChannelCount = endChannel - startChannel;

  size_t polIndex;
  if (polOut == aocommon::Polarization::Instrumental) {
    if (polsIn.size() != 4)
      throw std::runtime_error(
          "This mode requires the four polarizations to be present in the "
          "measurement set");
    for (size_t ch = 0; ch != selectedChannelCount * polsIn.size(); ++ch) {
      if (IsCFinite(*inPtr))
        dest[ch] = *inPtr;
      else
        dest[ch] = 0;
      inPtr++;
    }
  } else if (aocommon::Polarization::TypeToIndex(polOut, polsIn, polIndex)) {
    inPtr += polIndex;
    for (size_t ch = 0; ch != selectedChannelCount; ++ch) {
      if (IsCFinite(*inPtr))
        dest[ch] = *inPtr;
      else
        dest[ch] = 0;
      inPtr += polCount;
    }
  } else {
    // Copy the right visibilities with conversion if necessary.
    switch (polOut) {
      case aocommon::Polarization::StokesI: {
        size_t polIndexA = 0, polIndexB = 0;
        bool hasXX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XX, polsIn, polIndexA);
        bool hasYY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YY, polsIn, polIndexB);
        if (!hasXX || !hasYY) {
          bool hasRR = aocommon::Polarization::TypeToIndex(
              aocommon::Polarization::RR, polsIn, polIndexA);
          bool hasLL = aocommon::Polarization::TypeToIndex(
              aocommon::Polarization::LL, polsIn, polIndexB);
          if (!hasRR || !hasLL)
            throw std::runtime_error(
                "Can not form requested polarization (Stokes I) from available "
                "polarizations");
        }

        for (size_t ch = 0; ch != selectedChannelCount; ++ch) {
          inPtr += polIndexA;
          casacore::Complex val = *inPtr;
          inPtr += polIndexB - polIndexA;

          // I = (XX + YY) / 2
          val = (*inPtr + val) * 0.5f;

          if (IsCFinite(val))
            dest[ch] = val;
          else
            dest[ch] = 0.0;

          inPtr += polCount - polIndexB;
        }
      } break;
      case aocommon::Polarization::StokesQ: {
        size_t polIndexA = 0, polIndexB = 0;
        bool hasXX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XX, polsIn, polIndexA);
        bool hasYY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YY, polsIn, polIndexB);
        if (hasXX && hasYY) {
          // Convert to StokesQ from XX and YY
          for (size_t ch = 0; ch != selectedChannelCount; ++ch) {
            inPtr += polIndexA;
            casacore::Complex val = *inPtr;
            inPtr += polIndexB - polIndexA;

            // Q = (XX - YY)/2
            val = (val - *inPtr) * 0.5f;

            if (IsCFinite(val))
              dest[ch] = val;
            else
              dest[ch] = 0.0;

            inPtr += polCount - polIndexB;
          }
        } else {
          // Convert to StokesQ from RR and LL
          bool hasRL = aocommon::Polarization::TypeToIndex(
              aocommon::Polarization::RL, polsIn, polIndexA);
          bool hasLR = aocommon::Polarization::TypeToIndex(
              aocommon::Polarization::LR, polsIn, polIndexB);
          if (!hasRL || !hasLR)
            throw std::runtime_error(
                "Can not form requested polarization (Stokes Q) from available "
                "polarizations");
          for (size_t ch = 0; ch != selectedChannelCount; ++ch) {
            inPtr += polIndexA;
            casacore::Complex val = *inPtr;
            inPtr += polIndexB - polIndexA;

            // Q = (RL + LR)/2
            val = (*inPtr + val) * 0.5f;

            if (IsCFinite(val))
              dest[ch] = val;
            else
              dest[ch] = 0.0;

            inPtr += polCount - polIndexB;
          }
        }
      } break;
      case aocommon::Polarization::StokesU: {
        size_t polIndexA = 0, polIndexB = 0;
        bool hasXY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XY, polsIn, polIndexA);
        bool hasYX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YX, polsIn, polIndexB);
        if (hasXY && hasYX) {
          // Convert to StokesU from XY and YX
          for (size_t ch = 0; ch != selectedChannelCount; ++ch) {
            inPtr += polIndexA;
            casacore::Complex val = *inPtr;
            inPtr += polIndexB - polIndexA;

            // U = (XY + YX)/2
            val = (val + *inPtr) * 0.5f;

            if (IsCFinite(val))
              dest[ch] = val;
            else
              dest[ch] = 0.0;

            inPtr += polCount - polIndexB;
          }
        } else {
          // Convert to StokesU from RR and LL
          bool hasRL = aocommon::Polarization::TypeToIndex(
              aocommon::Polarization::RL, polsIn, polIndexA);
          bool hasLR = aocommon::Polarization::TypeToIndex(
              aocommon::Polarization::LR, polsIn, polIndexB);
          if (!hasRL || !hasLR)
            throw std::runtime_error(
                "Can not form requested polarization (Stokes U) from available "
                "polarizations");
          for (size_t ch = 0; ch != selectedChannelCount; ++ch) {
            inPtr += polIndexA;
            casacore::Complex val = *inPtr;
            inPtr += polIndexB - polIndexA;

            // U = -i (RL - LR)/2
            val = (val - *inPtr) * 0.5f;
            val = casacore::Complex(val.imag(), -val.real());

            if (IsCFinite(val))
              dest[ch] = val;
            else
              dest[ch] = 0.0;

            inPtr += polCount - polIndexB;
          }
        }
      } break;
      case aocommon::Polarization::StokesV: {
        size_t polIndexA = 0, polIndexB = 0;
        bool hasXY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XY, polsIn, polIndexA);
        bool hasYX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YX, polsIn, polIndexB);
        if (hasXY && hasYX) {
          // Convert to StokesV from XX and YY
          for (size_t ch = 0; ch != selectedChannelCount; ++ch) {
            inPtr += polIndexA;
            casacore::Complex val = *inPtr;
            inPtr += polIndexB - polIndexA;

            // V = -i(XY - YX)/2
            val = (val - *inPtr) * 0.5f;
            val = casacore::Complex(val.imag(), -val.real());

            if (IsCFinite(val))
              dest[ch] = val;
            else
              dest[ch] = 0.0;

            inPtr += polCount - polIndexB;
          }
        } else {
          // Convert to StokesV from RR and LL
          bool hasRL = aocommon::Polarization::TypeToIndex(
              aocommon::Polarization::RR, polsIn, polIndexA);
          bool hasLR = aocommon::Polarization::TypeToIndex(
              aocommon::Polarization::LL, polsIn, polIndexB);
          if (!hasRL || !hasLR)
            throw std::runtime_error(
                "Can not form requested polarization (Stokes V) from available "
                "polarizations");
          for (size_t ch = 0; ch != selectedChannelCount; ++ch) {
            inPtr += polIndexA;
            casacore::Complex val = *inPtr;
            inPtr += polIndexB - polIndexA;

            // V = (RR - LL)/2
            val = (val - *inPtr) * 0.5f;

            if (IsCFinite(val))
              dest[ch] = val;
            else
              dest[ch] = 0.0;

            inPtr += polCount - polIndexB;
          }
        }
      } break;
      default:
        throw std::runtime_error(
            "Could not convert ms polarizations to requested polarization");
    }
  }
}

template <typename NumType>
void MSProvider::CopyWeights(
    NumType* dest, size_t startChannel, size_t endChannel,
    const std::vector<aocommon::PolarizationEnum>& polsIn,
    const casacore::Array<std::complex<float>>& data,
    const casacore::Array<float>& weights, const casacore::Array<bool>& flags,
    aocommon::PolarizationEnum polOut) {
  const size_t polCount = polsIn.size();
  casacore::Array<std::complex<float>>::const_contiter inPtr =
      data.cbegin() + startChannel * polCount;
  casacore::Array<float>::const_contiter weightPtr =
      weights.cbegin() + startChannel * polCount;
  casacore::Array<bool>::const_contiter flagPtr =
      flags.cbegin() + startChannel * polCount;
  const size_t selectedChannelCount = endChannel - startChannel;

  size_t polIndex;
  if (polOut == aocommon::Polarization::Instrumental) {
    for (size_t ch = 0; ch != selectedChannelCount * polsIn.size(); ++ch) {
      if (!*flagPtr && IsCFinite(*inPtr))
        // The factor of 4 is to be consistent with StokesI
        // It is for having conjugate visibilities and because IDG doesn't
        // separately count XX and YY visibilities
        dest[ch] = *weightPtr * 4.0f;
      else
        dest[ch] = 0.0f;
      inPtr++;
      weightPtr++;
      flagPtr++;
    }
  } else if (aocommon::Polarization::TypeToIndex(polOut, polsIn, polIndex)) {
    inPtr += polIndex;
    weightPtr += polIndex;
    flagPtr += polIndex;
    for (size_t ch = 0; ch != selectedChannelCount; ++ch) {
      if (!*flagPtr && IsCFinite(*inPtr))
        dest[ch] = *weightPtr;
      else
        dest[ch] = 0.0f;
      inPtr += polCount;
      weightPtr += polCount;
      flagPtr += polCount;
    }
  } else {
    size_t polIndexA = 0, polIndexB = 0;
    switch (polOut) {
      case aocommon::Polarization::StokesI: {
        bool hasXY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XX, polsIn, polIndexA);
        bool hasYX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YY, polsIn, polIndexB);
        if (!hasXY || !hasYX) {
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::RR,
                                              polsIn, polIndexA);
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::LL,
                                              polsIn, polIndexB);
        }
      } break;
      case aocommon::Polarization::StokesQ: {
        bool hasXX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XX, polsIn, polIndexA);
        bool hasYY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YY, polsIn, polIndexB);
        if (!hasXX || !hasYY) {
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::RL,
                                              polsIn, polIndexA);
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::LR,
                                              polsIn, polIndexB);
        }
      } break;
      case aocommon::Polarization::StokesU: {
        bool hasXY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XY, polsIn, polIndexA);
        bool hasYX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YX, polsIn, polIndexB);
        if (!hasXY || !hasYX) {
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::RL,
                                              polsIn, polIndexA);
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::LR,
                                              polsIn, polIndexB);
        }
      } break;
      case aocommon::Polarization::StokesV: {
        bool hasXY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XY, polsIn, polIndexA);
        bool hasYX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YX, polsIn, polIndexB);
        if (!hasXY || !hasYX) {
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::RR,
                                              polsIn, polIndexA);
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::LL,
                                              polsIn, polIndexB);
        }
      } break;
      default:
        throw std::runtime_error(
            "Could not convert ms polarizations to requested polarization");
        break;
    }

    weightPtr += polIndexA;
    inPtr += polIndexA;
    flagPtr += polIndexA;
    for (size_t ch = 0; ch != selectedChannelCount; ++ch) {
      NumType w;
      if (!*flagPtr && IsCFinite(*inPtr))
        w = *weightPtr * 4.0f;
      else
        w = 0.0f;
      inPtr += polIndexB - polIndexA;
      weightPtr += polIndexB - polIndexA;
      flagPtr += polIndexB - polIndexA;
      if (!*flagPtr && IsCFinite(*inPtr))
        w = std::min<NumType>(w, *weightPtr * 4.0f);
      else
        w = 0.0f;
      dest[ch] = w;
      weightPtr += polCount - polIndexB + polIndexA;
      inPtr += polCount - polIndexB + polIndexA;
      flagPtr += polCount - polIndexB + polIndexA;
    }
  }
}

template void MSProvider::CopyWeights<float>(
    float* dest, size_t startChannel, size_t endChannel,
    const std::vector<aocommon::PolarizationEnum>& polsIn,
    const casacore::Array<std::complex<float>>& data,
    const casacore::Array<float>& weights, const casacore::Array<bool>& flags,
    aocommon::PolarizationEnum polOut);

template void MSProvider::CopyWeights<std::complex<float>>(
    std::complex<float>* dest, size_t startChannel, size_t endChannel,
    const std::vector<aocommon::PolarizationEnum>& polsIn,
    const casacore::Array<std::complex<float>>& data,
    const casacore::Array<float>& weights, const casacore::Array<bool>& flags,
    aocommon::PolarizationEnum polOut);

template <bool add>
void MSProvider::ReverseCopyData(
    casacore::Array<std::complex<float>>& dest, size_t startChannel,
    size_t endChannel, const std::vector<aocommon::PolarizationEnum>& polsDest,
    const std::complex<float>* source, aocommon::PolarizationEnum polSource) {
  size_t polCount = polsDest.size();
  const size_t selectedChannelCount = endChannel - startChannel;
  casacore::Array<std::complex<float>>::contiter dataIter =
      dest.cbegin() + startChannel * polCount;

  size_t polIndex;
  if (polSource == aocommon::Polarization::Instrumental) {
    for (size_t chp = 0; chp != selectedChannelCount * polsDest.size(); ++chp) {
      if (std::isfinite(source[chp].real())) {
        AddOrAssign<add>(dataIter, source[chp]);
      }
      dataIter++;
    }
  } else if (aocommon::Polarization::TypeToIndex(polSource, polsDest,
                                                 polIndex)) {
    for (size_t ch = 0; ch != selectedChannelCount; ++ch) {
      if (std::isfinite(source[ch].real())) {
        AddOrAssign<add>(dataIter + polIndex, source[ch]);
      }
      dataIter += polCount;
    }
  } else {
    switch (polSource) {
      case aocommon::Polarization::StokesI: {
        size_t polIndexA = 0, polIndexB = 0;
        bool hasXX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XX, polsDest, polIndexA);
        bool hasYY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YY, polsDest, polIndexB);
        if (!hasXX || !hasYY) {
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::RR,
                                              polsDest, polIndexA);
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::LL,
                                              polsDest, polIndexB);
        }
        for (size_t ch = 0; ch != selectedChannelCount; ++ch) {
          if (std::isfinite(source[ch].real())) {
            AddOrAssign<add>(dataIter + polIndexA,
                             source[ch]);  // XX = I (or rr = I)
            AddOrAssign<add>(dataIter + polIndexB,
                             source[ch]);  // YY = I (or ll = I)
          }
          dataIter += polCount;
        }
      } break;
      case aocommon::Polarization::StokesQ: {
        size_t polIndexA = 0, polIndexB = 0;
        bool hasXX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XX, polsDest, polIndexA);
        bool hasYY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YY, polsDest, polIndexB);
        if (hasXX && hasYY) {
          // StokesQ to linear
          for (size_t ch = 0; ch != selectedChannelCount; ++ch) {
            if (std::isfinite(source[ch].real())) {
              casacore::Complex stokesI =
                  casacore::Complex::value_type(0.5) *
                  (*(dataIter + polIndexB) + *(dataIter + polIndexA));
              AddOrAssign<add>(dataIter + polIndexA,
                               stokesI + source[ch]);  // XX = I + Q
              AddOrAssign<add>(dataIter + polIndexB,
                               stokesI - source[ch]);  // YY = I - Q
            }
            dataIter += polCount;
          }
        } else {
          // StokesQ to circular
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::RL,
                                              polsDest, polIndexA);
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::LR,
                                              polsDest, polIndexB);
          for (size_t ch = 0; ch != selectedChannelCount; ++ch) {
            if (std::isfinite(source[ch].real())) {
              AddOrAssign<add>(dataIter + polIndexA,
                               source[ch]);  // rl = Q + iU (with U still zero)
              AddOrAssign<add>(dataIter + polIndexB,
                               source[ch]);  // lr = Q - iU (with U still zero)
            }
            dataIter += polCount;
          }
        }
      } break;
      case aocommon::Polarization::StokesU: {
        size_t polIndexA = 0, polIndexB = 0;
        bool hasXY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XY, polsDest, polIndexA);
        bool hasYX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YX, polsDest, polIndexB);
        if (hasXY && hasYX) {
          // StokesU to linear
          for (size_t ch = 0; ch != selectedChannelCount; ++ch) {
            if (std::isfinite(source[ch].real())) {
              AddOrAssign<add>(dataIter + polIndexA,
                               source[ch]);  // XY = (U + iV), V still zero
              AddOrAssign<add>(dataIter + polIndexB,
                               source[ch]);  // YX = (U - iV), V still zero
            }
            dataIter += polCount;
          }
        } else {
          // StokesU to circular
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::RL,
                                              polsDest, polIndexA);
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::LR,
                                              polsDest, polIndexB);
          for (size_t ch = 0; ch != selectedChannelCount; ++ch) {
            if (std::isfinite(source[ch].real())) {
              // Q = (RL + LR) / 2
              casacore::Complex stokesQ =
                  casacore::Complex::value_type(0.5) *
                  (*(dataIter + polIndexA) + *(dataIter + polIndexB));
              casacore::Complex iTimesStokesU =
                  casacore::Complex(-source[ch].imag(), source[ch].real());
              AddOrAssign<add>(dataIter + polIndexA,
                               stokesQ + iTimesStokesU);  // rl = Q + iU
              AddOrAssign<add>(dataIter + polIndexB,
                               stokesQ - iTimesStokesU);  // lr = Q - iU
            }
            dataIter += polCount;
          }
        }
      } break;
      case aocommon::Polarization::StokesV: {
        size_t polIndexA = 0, polIndexB = 0;
        bool hasXY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XY, polsDest, polIndexA);
        bool hasYX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YX, polsDest, polIndexB);
        if (hasXY && hasYX) {
          // StokesV to linear
          for (size_t ch = 0; ch != selectedChannelCount; ++ch) {
            if (std::isfinite(source[ch].real())) {
              // U = (YX + XY)/2
              casacore::Complex stokesU =
                  casacore::Complex::value_type(0.5) *
                  (*(dataIter + polIndexB) + *(dataIter + polIndexA));
              casacore::Complex iTimesStokesV =
                  casacore::Complex(-source[ch].imag(), source[ch].real());
              AddOrAssign<add>(dataIter + polIndexA,
                               stokesU + iTimesStokesV);  // XY = (U + iV)
              AddOrAssign<add>(dataIter + polIndexB,
                               stokesU - iTimesStokesV);  // YX = (U - iV)
            }
            dataIter += polCount;
          }
        } else {
          // StokesV to circular
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::RR,
                                              polsDest, polIndexA);
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::LL,
                                              polsDest, polIndexB);
          for (size_t ch = 0; ch != selectedChannelCount; ++ch) {
            if (std::isfinite(source[ch].real())) {
              // I = (RR + LL)/2
              casacore::Complex stokesI =
                  casacore::Complex::value_type(0.5) *
                  (*(dataIter + polIndexA) + *(dataIter + polIndexB));
              AddOrAssign<add>(dataIter + polIndexA,
                               stokesI + source[ch]);  // RR = I + V
              AddOrAssign<add>(dataIter + polIndexB,
                               stokesI - source[ch]);  // LL = I - V
            }
            dataIter += polCount;
          }
        }
      } break;
      default:
        throw std::runtime_error(
            "Can't store polarization in set (not implemented or conversion "
            "not possible)");
    }
  }
}

// Explicit instantiation for true/false
template void MSProvider::ReverseCopyData<true>(
    casacore::Array<std::complex<float>>& dest, size_t startChannel,
    size_t endChannel, const std::vector<aocommon::PolarizationEnum>& polsDest,
    const std::complex<float>* source, aocommon::PolarizationEnum polSource);

template void MSProvider::ReverseCopyData<false>(
    casacore::Array<std::complex<float>>& dest, size_t startChannel,
    size_t endChannel, const std::vector<aocommon::PolarizationEnum>& polsDest,
    const std::complex<float>* source, aocommon::PolarizationEnum polSource);

void MSProvider::ReverseCopyWeights(
    casacore::Array<float>& dest, size_t startChannel, size_t endChannel,
    const std::vector<aocommon::PolarizationEnum>& polsDest,
    const float* source, aocommon::PolarizationEnum polSource) {
  size_t polCount = polsDest.size();
  const size_t selectedChannelCount = endChannel - startChannel;
  casacore::Array<float>::contiter dataIter =
      dest.cbegin() + startChannel * polCount;

  size_t polIndex;
  if (polSource == aocommon::Polarization::Instrumental) {
    std::copy(source, source + selectedChannelCount * polsDest.size(),
              dataIter);
  } else if (aocommon::Polarization::TypeToIndex(polSource, polsDest,
                                                 polIndex)) {
    for (size_t ch = 0; ch != selectedChannelCount; ++ch) {
      *(dataIter + polIndex) = source[ch];
      dataIter += polCount;
    }
  } else {
    switch (polSource) {
      case aocommon::Polarization::StokesI: {
        size_t polIndexA = 0, polIndexB = 0;
        bool hasXX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XX, polsDest, polIndexA);
        bool hasYY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YY, polsDest, polIndexB);
        if (!hasXX || !hasYY) {
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::RR,
                                              polsDest, polIndexA);
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::LL,
                                              polsDest, polIndexB);
        }
        for (size_t ch = 0; ch != selectedChannelCount; ++ch) {
          *(dataIter + polIndexA) = source[ch];  // XX = I (or rr = I)
          *(dataIter + polIndexB) = source[ch];  // YY = I (or ll = I)
          dataIter += polCount;
        }
      } break;
      case aocommon::Polarization::StokesQ:
      case aocommon::Polarization::StokesU:
      case aocommon::Polarization::StokesV:
      default:
        throw std::runtime_error(
            "Can't store weights in measurement set for this combination of "
            "polarizations (not implemented or conversion not possible)");
    }
  }
}

void MSProvider::GetRowRange(casacore::MeasurementSet& ms,
                             const MSSelection& selection, size_t& startRow,
                             size_t& endRow) {
  startRow = 0;
  endRow = ms.nrow();
  if (selection.HasInterval()) {
    Logger::Info << "Determining first and last row index... ";
    Logger::Info.Flush();
    casacore::ROScalarColumn<double> timeColumn(
        ms, casacore::MS::columnName(casacore::MSMainEnums::TIME));
    double time = timeColumn(0);
    size_t timestepIndex = 0;
    for (size_t row = 0; row != ms.nrow(); ++row) {
      if (time != timeColumn(row)) {
        ++timestepIndex;
        if (timestepIndex == selection.IntervalStart()) startRow = row;
        if (timestepIndex == selection.IntervalEnd()) {
          endRow = row;
          break;
        }
        time = timeColumn(row);
      }
    }
    Logger::Info << "DONE (" << startRow << '-' << endRow << ")\n";
  }
}

void MSProvider::GetRowRangeAndIDMap(casacore::MeasurementSet& ms,
                                     const MSSelection& selection,
                                     size_t& startRow, size_t& endRow,
                                     const std::set<size_t>& dataDescIds,
                                     std::vector<size_t>& idToMSRow) {
  startRow = 0;
  endRow = ms.nrow();

  Logger::Info << "Mapping measurement set rows... ";
  Logger::Info.Flush();
  casacore::ArrayColumn<double> uvwColumn(
      ms, casacore::MS::columnName(casacore::MSMainEnums::UVW));
  casacore::ScalarColumn<int> antenna1Column(
      ms, casacore::MS::columnName(casacore::MSMainEnums::ANTENNA1));
  casacore::ScalarColumn<int> antenna2Column(
      ms, casacore::MS::columnName(casacore::MSMainEnums::ANTENNA2));
  casacore::ScalarColumn<int> fieldIdColumn(
      ms, casacore::MS::columnName(casacore::MSMainEnums::FIELD_ID));
  casacore::ScalarColumn<double> timeColumn(
      ms, casacore::MS::columnName(casacore::MSMainEnums::TIME));
  casacore::ScalarColumn<int> dataDescIdColumn(
      ms, ms.columnName(casacore::MSMainEnums::DATA_DESC_ID));
  double time = timeColumn(0);
  size_t timestepIndex = 0;
  bool timeStepSelected =
      !selection.HasInterval() || timestepIndex == selection.IntervalStart();
  for (size_t row = 0; row != ms.nrow(); ++row) {
    if (time != timeColumn(row)) {
      ++timestepIndex;
      if (selection.HasInterval() &&
          timestepIndex == selection.IntervalStart()) {
        startRow = row;
        timeStepSelected = true;
      }
      if (timestepIndex == selection.IntervalEnd()) {
        if (selection.HasInterval()) endRow = row;
        break;
      }
      time = timeColumn(row);
    }
    if (timeStepSelected) {
      const int a1 = antenna1Column(row), a2 = antenna2Column(row),
                fieldId = fieldIdColumn(row),
                dataDescId = dataDescIdColumn(row);
      casacore::Vector<double> uvw = uvwColumn(row);
      std::set<size_t>::const_iterator dataDescIdIter =
          dataDescIds.find(dataDescId);
      if (selection.IsSelected(fieldId, timestepIndex, a1, a2, uvw) &&
          dataDescIdIter != dataDescIds.end())
        idToMSRow.push_back(row);
    }
  }
  Logger::Info << "DONE (" << startRow << '-' << endRow << "; "
               << idToMSRow.size() << " rows)\n";
}

void MSProvider::InitializeModelColumn(casacore::MeasurementSet& ms) {
  casacore::ArrayColumn<casacore::Complex> dataColumn(
      ms, casacore::MS::columnName(casacore::MSMainEnums::DATA));
  if (ms.isColumn(casacore::MSMainEnums::MODEL_DATA)) {
    casacore::ArrayColumn<casacore::Complex> modelColumn(
        ms, casacore::MS::columnName(casacore::MSMainEnums::MODEL_DATA));
    casacore::IPosition dataShape = dataColumn.shape(0);
    bool isDefined = modelColumn.isDefined(0);
    bool isSameShape = false;
    if (isDefined) {
      casacore::IPosition modelShape = modelColumn.shape(0);
      isSameShape = modelShape == dataShape;
    }
    if (!isDefined || !isSameShape) {
      Logger::Warn << "WARNING: Your model column does not have the same shape "
                      "as your data column: resetting MODEL column.\n";
      const casacore::Array<casacore::Complex> zeroArray(
          dataShape, casacore::Complex(0.0, 0.0));
      for (size_t row = 0; row != ms.nrow(); ++row)
        modelColumn.put(row, zeroArray);
    }
  } else {  // if(!_ms.isColumn(casacore::MSMainEnums::MODEL_DATA))
    ms.reopenRW();
    Logger::Info << "Adding model data column... ";
    Logger::Info.Flush();
    casacore::IPosition shape = dataColumn.shape(0);
    casacore::ArrayColumnDesc<casacore::Complex> modelColumnDesc(
        ms.columnName(casacore::MSMainEnums::MODEL_DATA), shape);
    try {
      ms.addColumn(modelColumnDesc, "StandardStMan", true, true);
    } catch (std::exception& e) {
      ms.addColumn(modelColumnDesc, "StandardStMan", false, true);
    }

    const casacore::Array<casacore::Complex> zeroArray(
        shape, casacore::Complex(0.0, 0.0));
    casacore::ArrayColumn<casacore::Complex> modelColumn(
        ms, casacore::MS::columnName(casacore::MSMainEnums::MODEL_DATA));

    for (size_t row = 0; row != ms.nrow(); ++row)
      modelColumn.put(row, zeroArray);

    Logger::Info << "DONE\n";
  }
}

casacore::ArrayColumn<float> MSProvider::InitializeImagingWeightColumn(
    casacore::MeasurementSet& ms) {
  ms.reopenRW();
  casacore::ArrayColumn<casacore::Complex> dataColumn(
      ms, casacore::MS::columnName(casacore::MSMainEnums::DATA));
  if (ms.tableDesc().isColumn("IMAGING_WEIGHT_SPECTRUM")) {
    return casacore::ArrayColumn<float>(ms, "IMAGING_WEIGHT_SPECTRUM");
  } else {
    Logger::Info << "Adding imaging weight spectrum column... ";
    Logger::Info.Flush();
    casacore::IPosition shape = dataColumn.shape(0);
    casacore::ArrayColumnDesc<float> modelColumnDesc("IMAGING_WEIGHT_SPECTRUM",
                                                     shape);
    try {
      ms.addColumn(modelColumnDesc, "StandardStMan", true, true);
    } catch (std::exception& e) {
      ms.addColumn(modelColumnDesc, "StandardStMan", false, true);
    }

    casacore::Array<float> zeroArray(shape);
    for (casacore::Array<float>::contiter i = zeroArray.cbegin();
         i != zeroArray.cend(); ++i)
      *i = 0.0;

    casacore::ArrayColumn<float> imgWColumn(ms, "IMAGING_WEIGHT_SPECTRUM");
    for (size_t row = 0; row != ms.nrow(); ++row)
      imgWColumn.put(row, zeroArray);
    Logger::Info << "DONE\n";
    return imgWColumn;
  }
}

std::vector<aocommon::PolarizationEnum> MSProvider::GetMSPolarizations(
    const casacore::MSPolarization& polTable) {
  std::vector<aocommon::PolarizationEnum> pols;
  casacore::ArrayColumn<int> corrTypeColumn(
      polTable, casacore::MSPolarization::columnName(
                    casacore::MSPolarizationEnums::CORR_TYPE));

  casacore::Array<int> corrTypeVec(corrTypeColumn(0));
  for (casacore::Array<int>::const_contiter p = corrTypeVec.cbegin();
       p != corrTypeVec.cend(); ++p) {
    pols.push_back(aocommon::Polarization::AipsIndexToEnum(*p));
  }

  return pols;
}

void MSProvider::ResetModelColumn() {
  std::unique_ptr<MSReader> msReader = MakeReader();
  SynchronizedMS ms = MS();
  ms->reopenRW();
  const std::vector<std::complex<float>> buffer(NChannels() * NPolarizations(),
                                                {0.0f, 0.0f});
  while (msReader->CurrentRowAvailable()) {
    // Always overwrite
    const bool addToMS = false;
    WriteModel(buffer.data(), addToMS);
    NextOutputRow();
    msReader->NextInputRow();
  }
}

bool MSProvider::OpenWeightSpectrumColumn(
    casacore::MeasurementSet& ms,
    std::unique_ptr<casacore::ROArrayColumn<float>>& weightColumn,
    const casacore::IPosition& dataColumnShape) {
  bool isWeightDefined;
  if (ms.isColumn(casacore::MSMainEnums::WEIGHT_SPECTRUM)) {
    weightColumn.reset(new casacore::ROArrayColumn<float>(
        ms, casacore::MS::columnName(casacore::MSMainEnums::WEIGHT_SPECTRUM)));
    isWeightDefined = weightColumn->isDefined(0);
  } else {
    isWeightDefined = false;
  }
  casacore::Array<float> weightArray(dataColumnShape);
  if (isWeightDefined) {
    casacore::IPosition weightShape = weightColumn->shape(0);
    isWeightDefined = (weightShape == dataColumnShape);
  }
  if (!isWeightDefined) {
    Logger::Warn
        << "WARNING: This measurement set has no or an invalid WEIGHT_SPECTRUM "
           "column; will use less informative WEIGHT column.\n";
    weightColumn.reset();
  }
  return isWeightDefined;
}
