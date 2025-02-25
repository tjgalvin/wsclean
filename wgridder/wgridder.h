#ifndef WSCLEAN_WGRIDDER_H_
#define WSCLEAN_WGRIDDER_H_

#include <complex>
#include <cstddef>
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
  const std::pair<size_t, size_t> *antennas;
  const std::complex<float> *visibilities;
  const size_t *time_offsets;
  MsGridder *gridder;
  size_t n_antennas;
  const std::complex<float> *parm_response;
};

class WGridderBase {
 public:
  virtual ~WGridderBase() = default;
  virtual size_t ConstantMemoryUsage() const = 0;
  virtual size_t PerVisibilityMemoryUsage() const = 0;
  virtual void InitializeInversion() = 0;
  virtual void AddInversionData(size_t n_rows, size_t n_chan, const double *uvw,
                                const double *freq,
                                const std::complex<float> *vis) = 0;
  virtual void AddInversionDataWithCorrectionCallback(
      GainMode mode, size_t n_polarizations, size_t n_rows, const double *uvws,
      const double *frequencies, VisibilityCallbackData &data) = 0;
  virtual void FinalizeImage(double multiplication_factor) = 0;
  virtual std::vector<float> RealImage() = 0;
  virtual void InitializePrediction(const float *image_data) = 0;
  virtual void PredictVisibilities(size_t n_rows, size_t n_channels,
                                   const double *uvws,
                                   const double *frequencies,
                                   std::complex<float> *visibilities) const = 0;
};

/* Memory usage of this gridder is:
   between calls:
     width*height*4  between calls (dirty image buffer)
   during gridding/degridding calls:
     width*height*(
       4      +    (dirty image buffer)
       8*2*2  +    (padded complex uv grid)
       8*2*2  )    (second uv grid used during FFT, not really necessary)
     + nvis_unflagged*8 (index arrays, rough guess)
*/
template <typename NumT>
class WGridder final : public WGridderBase {
 private:
  static constexpr double sigma_min = 1.1;
  static constexpr double sigma_max = 2.0;
  size_t width_;
  size_t height_;
  size_t trimmed_width_;
  size_t trimmed_height_;
  size_t n_threads_;
  double pixel_size_x_;
  double pixel_size_y_;
  double l_shift_;
  double m_shift_;
  double epsilon_;
  std::vector<NumT> image_;
  size_t verbosity_;
  bool tuning_;

 public:
  /** Construct a new gridder with given settings.
   * @param epsilon The requested accuracy of the gridding process.
   *   Affects the support of the employed kernel. Useful values
   *   range between 1e-2 and 1e-6 (for single-precision visibilities).
   * @param verbosity The amount of diagnostic output printed
   *   0: no output
   *   1: print short overview for every inversion/prediction
   *   2: print information for every processed w-plane
   */
  WGridder(size_t width, size_t height, size_t trimmed_width,
           size_t trimmed_height, double pixel_size_x, double pixel_size_y,
           double l_shift, double m_shift, size_t n_threads,
           double epsilon = 1e-4, size_t verbosity = 0, bool tuning_ = false);

  WGridder(const WGridder &) = delete;
  WGridder &operator=(const WGridder &) = delete;

  /**
   * @return The constant base memory usage of the object in bytes
   */
  size_t ConstantMemoryUsage() const final;
  /**
   * @return Additional memory required per gridded visibility in bytes.
   */
  size_t PerVisibilityMemoryUsage() const final;

  /**
   * Initialize a new inversion gridding pass. This just
   * intializes the accumulated dirty image with zero.
   */
  void InitializeInversion() final;

  /** Add more data to the current inversion operation.
   * The visibilities will be gridded, and the dirty image
   * will be updated accordingly.
   * visibilities with value 0 will be skipped entirely.
   * @param uvws pointer to n_rows*3 doubles containing UVW in m.
   *        U(row) := uvws[3*row  ]
   *        V(row) := uvws[3*row+1]
   *        W(row) := uvws[3*row+2]
   * @param visibilities pointer to nrow*n_channels complex<float> containing
   * weighted and corrected visibilities: visibility(row, chan) :=
   * vis[row*n_channels + chan]
   */
  void AddInversionData(size_t n_rows, size_t n_channels, const double *uvws,
                        const double *frequencies,
                        const std::complex<float> *visibilities) final;
  /** Equivalent to @ref AddInversionData() but without facet solutions
   * pre-applied and with additional paramaters to allow the creation of a
   * callback that can apply solutions "on the fly" as required
   *
   * It is expected that corrections have already been summed via @ref
   * LoadAndApplyCorrections<ModifierBehaviour::kSum>()
   *
   * @param n_polarizations The number of polarizations per visibility in @ref
   * visibilities
   * @param antennas Pointer to n_rows `std::pair<size_t, size_t>` containing
   * the antenna pair for each row of visibilities.
   * @param uvws Pointer to n_rows*3 doubles containing UVW in meters.
   *        U(row) := uvws[3*row  ]
   *        V(row) := uvws[3*row+1]
   *        W(row) := uvws[3*row+2]
   * @param frequencies Pointer to n_channels doubles containing channel
   * frequencies
   * @param visibilities Pointer to n_rows*n_channels*n_polarizations
   * `complex<float>` containing weighted but uncorrected and not yet collapsed
   * visibilities: visibility(row, chan) := vis[row*n_chan + chan]
   * @param time_offsets Pointer to n_rows `size_t` containing the time offset
   * as calculated by @ref CacheParmResponse() for the corresponding visibility
   * row when applying @ref LoadAndApplyCorrections<ModifierBehaviour::kSum>()
   * on it For further explanation see @ref
   * VisibilityCallbackBuffer::time_offsets_
   * @param gridder Pointer to a gridder that can be called back into in order
   * to apply solutions
   */
  void AddInversionDataWithCorrectionCallback(
      GainMode mode, size_t n_polarizations, size_t n_rows, const double *uvws,
      const double *frequencies, VisibilityCallbackData &data) final;

  /**
   * Finalize inversion once all passes are performed.
   * @param multiplication_factor Apply this factor to all pixels. This can be
   * used to normalize the image for the weighting scheme.
   */
  void FinalizeImage(double multiplication_factor) final;

  /**
   * Get the untrimmed image result of inversion. This is an array of size width
   * x height, and can be indexed with [x + y*width]. It is allowed to change
   * this image, e.g. set the horizon to zero before saving to fits. This call
   * is only valid once @ref FinalizeImage() has been called.
   */
  std::vector<float> RealImage() final;

  /**
   * Initialize gridder for prediction and specify image to predict for.
   * @param image The (untrimmed) model image that is to be predicted for. This
   * is an array of width * height size, index by (x + width*y).
   */
  void InitializePrediction(const float *image_data) final;

  /** Predicts visibilities from the current dirty image.
   * FIXME: how do we indicate flagged visibilities that do not
   *        need to be computed? Some special value on input?
   * @param uvws pointer to n_rows*3 doubles containing UVW in m.
   *        U(row) := uvws[3*row  ]
   *        V(row) := uvws[3*row+1]
   *        W(row) := uvws[3*row+2]
   * @param visibilities pointer to nrow*n_channels complex<float> containing
   * weighted visibilities visibility(row, chan) := vis[row*n_channels + chan]
   */
  void PredictVisibilities(size_t n_rows, size_t n_channels, const double *uvws,
                           const double *frequencies,
                           std::complex<float> *visibilities) const final;

 private:
  /** Internal helper to handle the processing of
   * AddInversionData/AddInversionDataWithCorrectionCallback
   * @tparam TMs Container type for visibilities, must be derived from or
   * compatible with cmav interface
   * @param freq Buffer via which all frequencies for this inversion can be
   * accessed, indexing must be in ascending order
   * @param ms Virtual or actual buffer via which all visibilities for this
   * inversion can be accessed, indexing must be in ascending order
   */
  template <typename Tms>
  void AddInversionMs(size_t n_rows, const double *uvw,
                      const ducc0::cmav<double, 1> &freq, Tms &ms) {
    ducc0::cmav<double, 2> uvw2(uvw, {n_rows, 3});
    ducc0::vmav<NumT, 2> tdirty({trimmed_width_, trimmed_height_});
    ducc0::cmav<float, 2> twgt(nullptr, {0, 0});
    ducc0::cmav<std::uint8_t, 2> tmask(nullptr, {0, 0});
    if (!tuning_)
      ducc0::ms2dirty<NumT, NumT>(uvw2, freq, ms, twgt, tmask, pixel_size_x_,
                                  pixel_size_y_, epsilon_, true, n_threads_,
                                  tdirty, verbosity_, true, false, sigma_min,
                                  sigma_max, -l_shift_, -m_shift_);
    else
      ducc0::ms2dirty_tuning<NumT, NumT>(
          uvw2, freq, ms, twgt, tmask, pixel_size_x_, pixel_size_y_, epsilon_,
          true, n_threads_, tdirty, verbosity_, true, false, sigma_min,
          sigma_max, -l_shift_, -m_shift_);
    for (size_t i = 0; i < trimmed_width_ * trimmed_height_; ++i)
      image_[i] += tdirty.raw(i);
  }

  // Each of the following versions of CreateAndAddInversionMs converts the
  // first parameter into a template argument and then passes the rest of the
  // parameters on to the next call.
  // We add multiple template parameters in a chain to avoid an exponential
  // explosion of boilerplate code.
  // Sorted and numbered in the order that they will be called.
  void CreateAndAddInversionMs(GainMode mode, size_t n_polarizations,
                               size_t n_parms, bool apply_beam,
                               bool apply_forward, bool has_h5_parm,
                               size_t n_rows, const double *uvws,
                               const ducc0::cmav<double, 1> &frequencies,
                               VisibilityCallbackData &data);
  template <GainMode Mode>
  void CreateAndAddInversionMs2(size_t n_polarizations, size_t n_parms,
                                bool apply_beam, bool apply_forward,
                                bool has_h5_parm, size_t n_rows,
                                const double *uvws,
                                const ducc0::cmav<double, 1> &frequencies,
                                VisibilityCallbackData &data);
  template <GainMode Mode, size_t NPolarizations>
  void CreateAndAddInversionMs3(size_t n_parms, bool apply_beam,
                                bool apply_forward, bool has_h5_parm,
                                size_t n_rows, const double *uvws,
                                const ducc0::cmav<double, 1> &frequencies,
                                VisibilityCallbackData &data);
  template <GainMode Mode, size_t NPolarizations, size_t NParms>
  void CreateAndAddInversionMs4(bool apply_beam, bool apply_forward,
                                bool has_h5_parm, size_t n_rows,
                                const double *uvws,
                                const ducc0::cmav<double, 1> &frequencies,
                                VisibilityCallbackData &data);
  template <GainMode Mode, size_t NPolarizations, size_t NParms, bool ApplyBeam>
  void CreateAndAddInversionMs5(bool apply_forward, bool has_h5_parm,
                                size_t n_rows, const double *uvws,
                                const ducc0::cmav<double, 1> &frequencies,
                                VisibilityCallbackData &data);
  template <GainMode Mode, size_t NPolarizations, size_t NParms, bool ApplyBeam,
            bool ApplyForward>
  void CreateAndAddInversionMs6(bool has_h5_parm, size_t n_rows,
                                const double *uvws,
                                const ducc0::cmav<double, 1> &frequencies,
                                VisibilityCallbackData &data);
  // Construct a VisibilityCallbackBuffer object using the remaining paramaters
  // and templatized based on all the paramaters that previous calls in the
  // template chain have parsed.
  // Call AddInversionMs with the constructed callback object.
  template <GainMode Mode, size_t NPolarizations, size_t NParms, bool ApplyBeam,
            bool ApplyForward, bool HasH5Parm>
  void CreateAndAddInversionMs7(size_t n_rows, const double *uvws,
                                const ducc0::cmav<double, 1> &frequencies,
                                VisibilityCallbackData &data);
};

}  // namespace wsclean

#endif  // WSCLEAN_WGRIDDER_H_

/*
Usage scenario:

WGridder gridder(width, height, pixel_size_x, pixel_size_y,
n_threads, 1e-5);
// determine number of visibilities that can be gridded in one go, using
// gridder.memUsage() and information about available memory.
// Making the chunks as large as posssible will improve perfrmance.
// Making chunks compact in w will also help a lot.
gridder.InitializeInversion();
for (auto &chunk: chunks)
  gridder.AddInversionData(...)
auto res = gridder.RealImage();
*/
