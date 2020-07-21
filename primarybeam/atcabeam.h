#ifndef ATCA_BEAM_H
#define ATCA_BEAM_H

#include "voltagepattern.h"

#include <aocommon/uvector.h>

#include <string>

class ATCABeam {
 public:
  enum Band {
    ATCA_X,
    ATCA_C,
    ATCA_16,
    ATCA_K,
    ATCA_Q,
    ATCA_W,
    ATCA_L1,
    ATCA_L2,
    ATCA_L3,
    ATCA_S,
    ATCA_C_RI,
    ATCA_Unknown
  };

  static enum Band GetBand(double freqGHz);
  static enum Band GetBand(const std::string& str);

  static VoltagePattern CalculateVoltagePattern(enum Band band);

  static void Calculate(class PrimaryBeamImageSet& beamImages,
                        double pixelScaleX, double pixelScaleY,
                        double phaseCentreRA, double phaseCentreDec,
                        double phaseCentreDL, double phaseCentreDM,
                        double frequencyHz,
                        const VoltagePattern& voltagePattern);

 private:
};

#endif
