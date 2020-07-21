#ifndef WGRIDDING_GRIDDER_SIMPLE_H
#define WGRIDDING_GRIDDER_SIMPLE_H

#include <complex>
#include <vector>
#include <cstddef>

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

class WGriddingGridder_Simple {
 private:
  size_t width_, height_, width_t_, height_t_, nthreads_;
  double pixelSizeX_, pixelSizeY_, epsilon_;
  std::vector<float> img;
  size_t verbosity_;

 public:
  /** Construct a new gridder with given settings.
   * @param width The width of the untrimmed image in pixels
   * @param height The height of the untrimmed image in pixels.
   * @param width_t The width of the trimmed image in pixels
   * @param height_t The height of the trimmed image in pixels.
   * @param pixelSizeX The angular width of a pixel in radians.
   * @param pixelSizeY The angular height of a pixel in radians.
   * @param nthreads The number of threads to use
   * @param epsilon The requested accuracy of the gridding process.
   *   Affects the support of the employed kernel. Useful values
   *   range between 1e-2 and 1e-6 (for single-precision visibilities).
   * @param verbosity The amount of diagnostic output printed
   *   0: no output
   *   1: print short overview for every inversion/prediction
   *   2: print information for every processed w-plane
   */
  WGriddingGridder_Simple(size_t width, size_t height, size_t width_t,
                          size_t height_t, double pixelSizeX, double pixelSizeY,
                          size_t nthreads, double epsilon = 1e-4,
                          size_t verbosity = 0);

  WGriddingGridder_Simple(const WGriddingGridder_Simple &) = delete;
  WGriddingGridder_Simple &operator=(const WGriddingGridder_Simple &) = delete;

  /** Returns the base memory requirement of the object, and
   *  the additional overhead per gridded/degridded viibility
   * @param constant The constant base memory usage in bytes
   * @param per_vis additional memory required per gridded visibility in bytes.
   */
  void memUsage(size_t &constant, size_t &per_vis) const;

  /**
   * Initialize a new inversion gridding pass. This just
   * intializes the accumulated dirty image with zero.
   */
  void InitializeInversion();

  /** Add more data to the current inversion operation.
   * The visibilities will be gridded, and the dirty image
   * will be updated accordingly.
   * visibilities with value 0 will be skipped entirely.
   * @param nrows The number of MS rows being passed
   * @param nchan The number of frequency channels
   * @param uvw pointer to nrows*3 doubles containing UVW in m.
   *        U(row) := uvw[3*row  ]
   *        V(row) := uvw[3*row+1]
   *        W(row) := uvw[3*row+2]
   * @param freq pointer to nchan doubles containing channel frequencies
   * @param vis pointer to nrow*nchan complex<float> containing weighted
   *        visibilities
   *        visibility(row, chan) := vis[row*nchan + chan]
   */
  void AddInversionData(size_t nrows, size_t nchan, const double *uvw,
                        const double *freq, const std::complex<float> *vis);

  /**
   * Finalize inversion once all passes are performed.
   * @param multiplicationFactor Apply this factor to all pixels. This can be
   * used to normalize the image for the weighting scheme.
   * @param correctFFTFactor For normal imaging this should be false. It can be
   * set to true when no correcting for weighting is performed, but the FFT
   * factor should be taken out. WSClean uses this for certain values of
   * "-visibility-weighting-mode".
   */
  void FinalizeImage(double multiplicationFactor, bool correctFFTFactor);

  /**
   * Get the untrimmed image result of inversion. This is an array of size width
   * x height, and can be indexed with [x + y*width]. It is allowed to change
   * this image, e.g. set the horizon to zero before saving to fits. This call
   * is only valid once @ref FinalizeImage() has been called.
   */
  std::vector<float> RealImage();

  /**
   * Initialize gridder for prediction and specify image to predict for.
   * @param image The (untrimmed) model image that is to be predicted for. This
   * is an array of width * height size, index by (x + width*y).
   */
  void InitializePrediction(std::vector<float> &&image);

  /** Predicts visibilities from the current dirty image.
   * FIXME: how do we indicate flagged visibilities that do not
   *        need to be computed? Some special value on input?
   * @param nrows The number of MS rows being passed
   * @param nchan The number of frequency channels
   * @param uvw pointer to nrows*3 doubles containing UVW in m.
   *        U(row) := uvw[3*row  ]
   *        V(row) := uvw[3*row+1]
   *        W(row) := uvw[3*row+2]
   * @param freq pointer to nchan doubles containing channel frequencies
   * @param vis pointer to nrow*nchan complex<float> containing weighted
   *        visibilities
   *        visibility(row, chan) := vis[row*nchan + chan]
   */
  void PredictVisibilities(size_t nrows, size_t nchan, const double *uvw,
                           const double *freq, std::complex<float> *vis) const;
};

#endif

/*
Usage scenario:

WGriddingGridder_Simple gridder(width, height, pixelSizeX, pixelSizeY, nthreads,
1e-5);
// determine number of visibilities that can be gridded in one go, using
// gridder.memUsage() and information about available memory.
// Making the chunks as large as posssible will improve perfrmance.
// Making chunks compact in w will also help a lot.
gridder.InitializeInversion();
for (auto &chunk: chunks)
  gridder.AddInversionData(...)
auto res = gridder.RealImage();
*/
