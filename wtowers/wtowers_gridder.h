#ifndef WSCLEAN_WTOWERS_GRIDDER_H_
#define WSCLEAN_WTOWERS_GRIDDER_H_

#include <complex>
#include <cstddef>
#include <vector>

#include <ska-sdp-func/grid_data/sdp_gridder_utils.h>

#include "../gridding/msgridder.h"

namespace wsclean {

class MsGridder;

class WTowersGridderBase {
 public:
  virtual ~WTowersGridderBase() = default;
  virtual size_t ConstantMemoryUsage() const = 0;
  virtual size_t PerVisibilityMemoryUsage() const = 0;
  virtual void InitializeInversion() = 0;
  virtual void AddInversionData(size_t n_rows, size_t n_chan, const double *uvw,
                                const double *freq,
                                const std::complex<float> *vis) = 0;
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
class WTowersGridder final : public WTowersGridderBase {
 public:
  /** Construct a new gridder with given settings.
   * @param verbosity The amount of diagnostic output printed
   *   0: no output
   *   1: print information for planes and subgrids
   */
  WTowersGridder(size_t width, size_t height, size_t trimmed_width,
                 size_t trimmed_height, double pixel_size_x,
                 double pixel_size_y, double l_shift, double m_shift,
                 int wtowers_subgrid_size, int wtowers_support,
                 int wtowers_w_support, double wtowers_padding,
                 double wtowers_w_padding, size_t n_threads, double accuracy,
                 double max_abs_w, size_t verbosity = 0);

  WTowersGridder(const WTowersGridder &) = delete;
  WTowersGridder &operator=(const WTowersGridder &) = delete;

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

  /**
   * Print the various computed paramaters with which w-towers will be run to
   * the debug log; size of sub-grid, height of the w-tower etc.
   */
  void LogParameters() const;

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
  // General gridder paramaters shared with other gridders
  size_t width_;
  size_t height_;
  size_t trimmed_width_;
  size_t trimmed_height_;
  size_t n_threads_;
  double pixel_size_x_;
  double pixel_size_y_;
  double l_shift_;
  double m_shift_;
  int subgrid_size_;
  int support_;
  int w_support_;
  double padding_factor_;
  double padding_factor_w_;
  size_t verbosity_;
  std::vector<NumT> image_;

  // W-towers specific paramaters
  struct {
    // Determined algorithmically
    sdp_GridProjection projection;
    double field_of_view_l;
    double field_of_view_m;
    double w_step;
    int w_towers_height;
    int oversampling = 32768;
    int w_oversampling = 32768;
    // User configurable
    double subgrid_frac = 2.0f / 3.0f;
  } wtowers_parameters_;
};

}  // namespace wsclean

#endif  // WSCLEAN_WTOWERS_GRIDDER_SIMPLE_H_
