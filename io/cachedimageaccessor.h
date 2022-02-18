#ifndef CACHED_IMAGE_ACCESSOR_H_
#define CACHED_IMAGE_ACCESSOR_H_

#include <aocommon/imageaccessor.h>

#include "cachedimageset.h"

/**
 * @brief ImageAccessor implementation that internally uses CachedImageSet.
 */
class CachedImageAccessor : public aocommon::ImageAccessor {
 public:
  /**
   * @brief Construct a new CachedImageAccessor object
   *
   * @param image_set The CachedImageSet for loading and storing the image.
   * @param polarization Image polarization.
   * @param frequency_index Frequency index of the image.
   * @param is_imaginary False: Image has real values. True: Image has imaginary
   * values.
   */
  CachedImageAccessor(CachedImageSet& image_set,
                      aocommon::PolarizationEnum polarization,
                      size_t frequency_index, bool is_imaginary)
      : image_set_(image_set),
        polarization_(polarization),
        frequency_index_(frequency_index),
        is_imaginary_(is_imaginary) {}

  void Load(aocommon::Image& image) const override {
    image_set_.Load(image.Data(), polarization_, frequency_index_,
                    is_imaginary_);
  }

  void Store(const aocommon::Image& image) override {
    image_set_.Store(image.Data(), polarization_, frequency_index_,
                     is_imaginary_);
  }

  /**
   * Get functions, mainly for testing purposes.
   * @{
   */
  const CachedImageSet& GetImageSet() const { return image_set_; }
  aocommon::PolarizationEnum GetPolarization() const { return polarization_; }
  size_t GetFrequencyIndex() const { return frequency_index_; }
  bool GetIsImaginary() const { return is_imaginary_; }
  /** @} */

 private:
  CachedImageSet& image_set_;
  aocommon::PolarizationEnum polarization_;
  size_t frequency_index_;
  bool is_imaginary_;
};

#endif