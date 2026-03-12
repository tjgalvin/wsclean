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

}  // namespace wsclean::tophat_convolution

#endif
