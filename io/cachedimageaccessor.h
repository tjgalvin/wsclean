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
        is_imaginary_(is_imaginary),
        facet_id_(0),
        facet_() {}

  /**
   * @brief Construct a new CachedImageAccessor object
   *
   * @param image_set The CachedImageSet for loading and storing the image.
   * @param polarization Image polarization.
   * @param frequency_index Frequency index of the image.
   * @param is_imaginary False: Image has real values. True: Image has imaginary
   * values.
   * @param facet [optional] Facet object, in case the accessor should use a
   * facet of the image. If null, the accessor uses the entire image.
   * @param facet_id The facet index, when the 'facet' argument is valid.
   * Otherwise, this argument is ignored.
   */
  CachedImageAccessor(CachedImageSet& image_set,
                      aocommon::PolarizationEnum polarization,
                      size_t frequency_index, size_t facet_id,
                      std::shared_ptr<const schaapcommon::facets::Facet> facet,
                      bool is_imaginary)
      : image_set_(image_set),
        polarization_(polarization),
        frequency_index_(frequency_index),
        is_imaginary_(is_imaginary),
        facet_id_(facet_id),
        facet_(std::move(facet)) {}

  std::size_t Width() const override {
    return facet_ ? facet_->GetTrimmedBoundingBox().Width()
                  : image_set_.Writer().Width();
  }

  std::size_t Height() const override {
    return facet_ ? facet_->GetTrimmedBoundingBox().Height()
                  : image_set_.Writer().Height();
  }

  void Load(float* data) const override {
    if (facet_) {
      image_set_.LoadFacet(data, polarization_, frequency_index_, facet_id_,
                           facet_, is_imaginary_);
    } else {
      image_set_.Load(data, polarization_, frequency_index_, is_imaginary_);
    }
  }

  void Store(const float* data) override {
    if (facet_) {
      const aocommon::Image image(const_cast<float*>(data), Width(), Height());
      image_set_.StoreFacet(image, polarization_, frequency_index_, facet_id_,
                            facet_, is_imaginary_);
    } else {
      image_set_.Store(data, polarization_, frequency_index_, is_imaginary_);
    }
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
  std::size_t facet_id_;
  std::shared_ptr<const schaapcommon::facets::Facet> facet_;
};

#endif