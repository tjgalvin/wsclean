#ifndef MATH_TOPHAT_CONVOLUTION_H_
#define MATH_TOPHAT_CONVOLUTION_H_

#include <aocommon/image.h>

namespace wsclean::tophat_convolution {

/**
 * Perform a convolution with a radial tophat with a given radius (in pixels).
 */
void Convolve(aocommon::Image& input, double radius);

/**
 * Produce an image with a radial tophat function. This means that
 * all pixels outside a specified radius are set to zero, and the other
 * pixels are set to the same non-zero value. The image is normalized such that
 * the full image sums to one.
 *
 * In theory it works for even sizes, but in that case the tophat filter
 * will chose the pixel at (width/2, height/2) as the centre pixel,
 * which means that the filter is not circular/otherwise symmetric.
 * When using an odd size, the tophat filter is circular/symmetric,
 * and that is the recommended usage.
 */
aocommon::Image MakeTopHatImage(size_t width, size_t height, double tophat_radius);

/**
 * The value below which a pixel of a mask that was convolved with Convolve() is
 * considered to be FFT rounding noise instead of a real mask value.
 *
 * Convolve() performs the convolution with FFTs, which leaves rounding noise of
 * either sign in the whole image, including in the regions where the convolved
 * mask should be exactly zero. The smallest value that the convolved mask can
 * legitimately have is one divided by the number of pixels of the tophat
 * kernel, because a single kernel pixel overlapping the edge of the mask is
 * the smallest possible contribution. That number is at most the number of
 * pixels of the square that encloses the kernel, i.e. (2r+1)^2 for a kernel
 * with radius r. The returned threshold is therefore at most half of the
 * smallest legitimate mask value, and no real mask data is removed by it.
 *
 * @param radius The radius that was passed to Convolve(). Using the same
 *        radius for both functions ensures that the threshold is based on the
 *        actual kernel size of the convolution.
 */
float MaskNoiseThreshold(double radius);

/**
 * Remove the FFT rounding noise from a mask that was convolved with Convolve().
 *
 * All values below MaskNoiseThreshold() are set to zero. Since the threshold
 * is positive, this removes every negative value, even if the rounding noise
 * is larger than the threshold for a very large kernel. The mask is therefore
 * non-negative for any radius, so summing the masks of several facets cannot
 * cancel out the real coverage of a pixel. See @ref MaskNoiseThreshold() for
 * the reasoning behind the threshold.
 *
 * @param radius The radius that was passed to Convolve().
 */
void RemoveMaskNoise(aocommon::Image& mask, double radius);

}  // namespace wsclean::tophat_convolution

#endif
