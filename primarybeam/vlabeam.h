#ifndef VLA_BEAM_H
#define VLA_BEAM_H

#include <array>
#include <map>
#include <string>

class VLABeam {
 public:
  static std::array<double, 5> GetCoefficients(const std::string& bandName,
                                               double freq);

 private:
  static std::map<int, std::array<double, 5>> getCoefficients();
  static std::map<char, double> getFeedConf();
  static char determineFeed(double freq, double freqCenter = 0.0);
  static void limitFreqForBand(char band, double& freq);
};

#endif
