#include "msreader.h"

void MSReader::copyData(std::complex<float>* dest, size_t startChannel,
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
      if (isCFinite(*inPtr))
        dest[ch] = *inPtr;
      else
        dest[ch] = 0;
      inPtr++;
    }
  } else if (aocommon::Polarization::TypeToIndex(polOut, polsIn, polIndex)) {
    inPtr += polIndex;
    for (size_t ch = 0; ch != selectedChannelCount; ++ch) {
      if (isCFinite(*inPtr))
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

          if (isCFinite(val))
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

            if (isCFinite(val))
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

            if (isCFinite(val))
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

            if (isCFinite(val))
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

            if (isCFinite(val))
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

            if (isCFinite(val))
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

            if (isCFinite(val))
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