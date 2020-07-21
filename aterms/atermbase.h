#ifndef ATERM_BASE_H
#define ATERM_BASE_H

#include <complex>

class ATermBase {
 public:
  ATermBase() : _saveATerms(false) {}

  virtual ~ATermBase() {}

  /**
   * Calculate the a-terms for the given time and frequency, for all stations.
   * @param buffer A buffer of size 4 x width x height x nstation, in that
   * order.
   * @param time The time corresponding to the currently gridded visibilities to
   * which the aterm will be applied.
   * @param frequency Frequency of currently gridded visibilities.
   * @param fieldId Field that these visibilities belong to. Different fields
   * might requires different beams.
   * @param uvwInM The UVW of the antennas, referenced to antenna 0.
   * @returns @c True when new aterms are calculated. If these aterms are the
   * same as for the previous call to Calculate(), @c false can be returned and
   * the output buffer does not need to be updated. The gridder will then make
   * sure to use the previous aterms, and not reserve extra memory for it etc.
   */
  virtual bool Calculate(std::complex<float>* buffer, double time,
                         double frequency, size_t fieldId,
                         const double* uvwInM) = 0;

  virtual double AverageUpdateTime() const = 0;

  void StoreATermsEigenvalues(const std::string& filename,
                              const std::complex<float>* buffer,
                              size_t nStations, size_t width, size_t height);

  void StoreATermsReal(const std::string& filename,
                       const std::complex<float>* buffer, size_t nStations,
                       size_t width, size_t height);

  /**
   * Set whether a fits image with the a-terms should be written to disk
   * every time they are calculated.
   * @param saveATerms Fits images are saved when set to true.
   */
  void SetSaveATerms(bool saveATerms, const std::string& prefix) {
    _saveATerms = saveATerms;
    _prefix = prefix;
  }

  struct CoordinateSystem {
    size_t width, height, maxSupport;
    double ra, dec, dl, dm, phaseCentreDL, phaseCentreDM;
  };

 protected:
  void saveATermsIfNecessary(const std::complex<float>* buffer,
                             size_t nStations, size_t width, size_t height) {
    if (_saveATerms) {
      static int index = 0;
      std::ostringstream f;
      StoreATermsEigenvalues(
          _prefix + "-aterm-ev" + std::to_string(index) + ".fits", buffer,
          nStations, width, height);
      StoreATermsReal(
          _prefix + "-aterm-realxx" + std::to_string(index) + ".fits", buffer,
          nStations, width, height);
      ++index;
    }
  }

  bool _saveATerms;
  std::string _prefix;
};

#endif
