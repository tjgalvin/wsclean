#include "simpleclean.h"

#ifdef __SSE__
#define USE_INTRINSICS
#endif

#ifdef USE_INTRINSICS
#include <emmintrin.h>
#include <immintrin.h>
#endif

void SimpleClean::SubtractImage(float *image, const float *psf, size_t width,
                                size_t height, size_t x, size_t y,
                                float factor) {
  size_t startX, startY, endX, endY;
  int offsetX = (int)x - width / 2, offsetY = (int)y - height / 2;

  if (offsetX > 0)
    startX = offsetX;
  else
    startX = 0;

  if (offsetY > 0)
    startY = offsetY;
  else
    startY = 0;

  endX = x + width / 2;
  if (endX > width) endX = width;

  bool isAligned = ((endX - startX) % 2) == 0;
  if (!isAligned) --endX;

  endY = y + height / 2;
  if (endY > height) endY = height;

  for (size_t ypos = startY; ypos != endY; ++ypos) {
    float *imageIter = image + ypos * width + startX;
    const float *psfIter = psf + (ypos - offsetY) * width + startX - offsetX;
    for (size_t xpos = startX; xpos != endX; xpos++) {
      // I've SSE-ified this, but it didn't improve speed at all :-/
      // (Compiler probably already did it)
      *imageIter -= (*psfIter * factor);
      //*(imageIter+1) = *(imageIter+1) - (*(psfIter+1) * factor);
      ++imageIter;
      ++psfIter;
    }
  }
}

void SimpleClean::PartialSubtractImage(float *image, const float *psf,
                                       size_t width, size_t height, size_t x,
                                       size_t y, float factor, size_t startY,
                                       size_t endY) {
  size_t startX, endX;
  int offsetX = (int)x - width / 2, offsetY = (int)y - height / 2;

  if (offsetX > 0)
    startX = offsetX;
  else
    startX = 0;

  if (offsetY > (int)startY) startY = offsetY;

  endX = x + width / 2;
  if (endX > width) endX = width;

  bool isAligned = ((endX - startX) % 2) == 0;
  if (!isAligned) --endX;

  endY = std::min(y + height / 2, endY);

  for (size_t ypos = startY; ypos < endY; ++ypos) {
    float *imageIter = image + ypos * width + startX;
    const float *psfIter = psf + (ypos - offsetY) * width + startX - offsetX;
    for (size_t xpos = startX; xpos != endX; xpos += 2) {
      *imageIter = *imageIter - (*psfIter * factor);
      *(imageIter + 1) = *(imageIter + 1) - (*(psfIter + 1) * factor);
      imageIter += 2;
      psfIter += 2;
    }
    if (!isAligned) *imageIter -= *psfIter * factor;
  }
}

void SimpleClean::PartialSubtractImage(float *image, size_t imgWidth,
                                       size_t /*imgHeight*/, const float *psf,
                                       size_t psfWidth, size_t psfHeight,
                                       size_t x, size_t y, float factor,
                                       size_t startY, size_t endY) {
  size_t startX, endX;
  int offsetX = (int)x - psfWidth / 2, offsetY = (int)y - psfHeight / 2;

  if (offsetX > 0)
    startX = offsetX;
  else
    startX = 0;

  if (offsetY > (int)startY) startY = offsetY;

  endX = std::min(x + psfWidth / 2, imgWidth);

  bool isAligned = ((endX - startX) % 2) == 0;
  if (!isAligned) --endX;

  endY = std::min(y + psfHeight / 2, endY);

  for (size_t ypos = startY; ypos < endY; ++ypos) {
    float *imageIter = image + ypos * imgWidth + startX;
    const float *psfIter = psf + (ypos - offsetY) * psfWidth + startX - offsetX;
    for (size_t xpos = startX; xpos != endX; xpos += 2) {
      *imageIter = *imageIter - (*psfIter * factor);
      *(imageIter + 1) = *(imageIter + 1) - (*(psfIter + 1) * factor);
      imageIter += 2;
      psfIter += 2;
    }
    if (!isAligned) *imageIter -= *psfIter * factor;
  }
}

#if defined __AVX__ && defined USE_INTRINSICS
void SimpleClean::PartialSubtractImageAVX(double *image, size_t imgWidth,
                                          size_t /*imgHeight*/,
                                          const double *psf, size_t psfWidth,
                                          size_t psfHeight, size_t x, size_t y,
                                          double factor, size_t startY,
                                          size_t endY) {
  size_t startX, endX;
  int offsetX = (int)x - psfWidth / 2, offsetY = (int)y - psfHeight / 2;

  if (offsetX > 0)
    startX = offsetX;
  else
    startX = 0;

  if (offsetY > (int)startY) startY = offsetY;

  endX = std::min(x + psfWidth / 2, imgWidth);

  size_t unAlignedCount = (endX - startX) % 4;
  endX -= unAlignedCount;

  endY = std::min(y + psfHeight / 2, endY);

  const __m256d mFactor = _mm256_set1_pd(-factor);
  for (size_t ypos = startY; ypos < endY; ++ypos) {
    double *imageIter = image + ypos * imgWidth + startX;
    const double *psfIter =
        psf + (ypos - offsetY) * psfWidth + startX - offsetX;
    for (size_t xpos = startX; xpos != endX; xpos += 4) {
      __m256d imgVal = _mm256_loadu_pd(imageIter),
              psfVal = _mm256_loadu_pd(psfIter);
#ifdef __FMA__
      _mm256_storeu_pd(imageIter, _mm256_fmadd_pd(psfVal, mFactor, imgVal));
#else
      _mm256_storeu_pd(imageIter,
                       _mm256_add_pd(imgVal, _mm256_mul_pd(psfVal, mFactor)));
#endif
      imageIter += 4;
      psfIter += 4;
    }
    for (size_t xpos = endX; xpos != endX + unAlignedCount; ++xpos) {
      *imageIter -= *psfIter * factor;
      ++imageIter;
      ++psfIter;
    }
  }
}

#endif
