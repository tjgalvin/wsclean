#ifndef WSCLEAN_WGRIDDER_BASE_H_
#define WSCLEAN_WGRIDDER_BASE_H_

#include <complex>
#include <cstddef>
#include <span>
#include <vector>

#include <aocommon/banddata.h>
#include "ducc0/wgridder/wgridder.h"

#include "../gridding/gainmode.h"
#include "../gridding/msgridder.h"

namespace wsclean {

class MsGridder;

struct VisibilityCallbackData {
  size_t n_channels;
  const aocommon::BandData &selected_band;
  size_t data_desc_id;
  const std::pair<size_t, size_t> *antennas;
  const std::complex<float> *visibilities;
  const double *uvws;
  const size_t *time_offsets;
  MsGridder *gridder;
  size_t n_antennas;
  const std::complex<float> *parm_response;
  const BeamResponseCacheChunk &beam_response;
};

class WGridderBase {
 public:
  virtual ~WGridderBase() = default;
  virtual size_t ConstantMemoryUsage() const = 0;
  virtual size_t PerVisibilityMemoryUsage() const = 0;
  virtual void InitializeInversion() = 0;
  virtual void AddInversionData(size_t n_rows, size_t n_chan, const double *uvw, const double *freq,
                                const std::complex<float> *vis) = 0;
  virtual void AddIrregularInversionData(std::span<const size_t> band_ids,
                                         std::span<const size_t> n_channels_per_band,
                                         std::span<const double> frequencies,
                                         std::span<const double> uvw,
                                         std::span<const std::complex<float>> vis) = 0;
  virtual void AddInversionDataWithCorrectionCallback(GainMode mode, size_t n_polarizations,
                                                      size_t n_rows, const double *uvws,
                                                      const double *frequencies,
                                                      VisibilityCallbackData &data) = 0;
  virtual void FinalizeImage(double multiplication_factor) = 0;
  virtual std::vector<float> RealImage() = 0;
  virtual void InitializePrediction(const float *image_data) = 0;
  virtual void PredictVisibilities(size_t n_rows, size_t n_channels, const double *uvws,
                                   const double *frequencies,
                                   std::complex<float> *visibilities) const = 0;
};

}  // namespace wsclean

#endif  // WSCLEAN_WGRIDDER_BASE_H_
