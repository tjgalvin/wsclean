#ifndef DECONVOLUTION_IMAGE_SET_H
#define DECONVOLUTION_IMAGE_SET_H

#include "deconvolutiontable.h"

#include <aocommon/image.h>
#include <aocommon/uvector.h>

#include <vector>
#include <map>
#include <memory>

class ImageSet {
 public:
  ImageSet(const DeconvolutionTable& table, bool squared_joins,
           const std::set<aocommon::PolarizationEnum>& linked_polarizations,
           size_t width, size_t height);

  /**
   * Constructs an ImageSet with the same settings as an existing ImageSet
   * object, but a different image size.
   *
   * @param image_set An existing ImageSet.
   * @param width The image width for the new ImageSet.
   * @param height The image height for the new ImageSet.
   */
  ImageSet(const ImageSet& image_set, size_t width, size_t height);

  ImageSet(const ImageSet&) = delete;
  ImageSet& operator=(const ImageSet&) = delete;

  aocommon::Image Release(size_t imageIndex) {
    return std::move(_images[imageIndex]);
  }

  void SetImage(size_t imageIndex, aocommon::Image&& data) {
    assert(data.Width() == Width() && data.Height() == Height());
    _images[imageIndex] = std::move(data);
  }

  /**
   * Replace the images of this ImageSet. The images may be of a different size.
   * Both ImageSets are expected to be for the same deconvolution configuration:
   * besides the images and their dimension, no fields are changed.
   */
  void SetImages(ImageSet&& source);

  /**
   * @param use_residual_images: True: Load residual images. False: Load model
   * images.
   */
  void LoadAndAverage(bool use_residual_images);

  std::vector<aocommon::Image> LoadAndAveragePSFs();

  void InterpolateAndStoreModel(const class SpectralFitter& fitter,
                                size_t threadCount);

  void AssignAndStoreResidual();

  /**
   * This function will calculate the integration over all images, squaring
   * images that are in the same squared-image group. For example, with
   * a squared group of [I, Q, ..] and another group [I2, Q2, ...], this
   * will calculate:
   *
   * sqrt(I^2 + Q^2 + ..) + sqrt(I2^2 + Q2^2 ..) + ..
   * ----------------------------------------------
   *           1          +           1          + ..
   *
   * If the 'squared groups' are of size 1, the average of the groups will be
   * returned (i.e., without square-rooting the square).
   *
   * If the squared joining option is set in the provided wsclean settings, the
   * behaviour of this method changes. In that case, it will return the square
   * root of the average squared value:
   *
   *       I^2 + Q^2 + ..  +  I2^2 + Q2^2 ..  + ..
   * sqrt( --------------------------------------- )
   *            1          +        1         + ..
   *
   * These formulae are such that the values will have normal flux values.
   * @param dest Pre-allocated output array that will be filled with the
   * integrated image.
   * @param scratch Pre-allocated scratch space, same size as image.
   */
  void GetSquareIntegrated(aocommon::Image& dest,
                           aocommon::Image& scratch) const {
    if (_squareJoinedChannels)
      getSquareIntegratedWithSquaredChannels(dest);
    else
      getSquareIntegratedWithNormalChannels(dest, scratch);
  }

  /**
   * This function will calculate the 'linear' integration over all images,
   * unless joined channels are requested to be squared. The method will return
   * the weighted average of all images. Normally,
   * @ref GetSquareIntegrated
   * should be used for peak finding, but in case negative values should remain
   * negative, such as with multiscale (otherwise a sidelobe will be fitted with
   * large scales), this function can be used.
   * @param dest Pre-allocated output array that will be filled with the average
   * values.
   */
  void GetLinearIntegrated(aocommon::Image& dest) const {
    if (_squareJoinedChannels)
      getSquareIntegratedWithSquaredChannels(dest);
    else
      getLinearIntegratedWithNormalChannels(dest);
  }

  void GetIntegratedPSF(aocommon::Image& dest,
                        const std::vector<aocommon::Image>& psfs);

  size_t NOriginalChannels() const {
    return _deconvolutionTable.OriginalGroups().size();
  }

  size_t PSFCount() const { return NDeconvolutionChannels(); }

  size_t NDeconvolutionChannels() const {
    return _deconvolutionTable.DeconvolutionGroups().size();
  }

  ImageSet& operator=(float val) {
    for (aocommon::Image& image : _images) image = val;
    return *this;
  }

  /**
   * Exposes image data.
   *
   * ImageSet only exposes a non-const pointer to the image data. When exposing
   * non-const reference to the images themselves, the user could change the
   * image size and violate the invariant that all images have equal sizes.
   * @param index An image index.
   * @return A non-const pointer to the data area for the image.
   */
  float* Data(size_t index) { return _images[index].Data(); }

  /**
   * Exposes the images in the image set.
   *
   * Creating a non-const version of this operator is not desirable. See Data().
   *
   * @param index An image index.
   * @return A const reference to the image with the given index.
   */
  const aocommon::Image& operator[](size_t index) const {
    return _images[index];
  }

  size_t size() const { return _images.size(); }

  size_t PSFIndex(size_t imageIndex) const {
    return _imageIndexToPSFIndex[imageIndex];
  }

  const DeconvolutionTable& Table() const { return _deconvolutionTable; }

  std::unique_ptr<ImageSet> Trim(size_t x1, size_t y1, size_t x2, size_t y2,
                                 size_t oldWidth) const {
    auto p = std::make_unique<ImageSet>(*this, x2 - x1, y2 - y1);
    for (size_t i = 0; i != _images.size(); ++i) {
      copySmallerPart(_images[i], p->_images[i], x1, y1, x2, y2, oldWidth);
    }
    return p;
  }

  /**
   * Like Trim(), but only copies values that are in the mask. All other values
   * are set to zero.
   * @param mask A mask of size (x2-x1) x (y2-y1)
   */
  std::unique_ptr<ImageSet> TrimMasked(size_t x1, size_t y1, size_t x2,
                                       size_t y2, size_t oldWidth,
                                       const bool* mask) const {
    std::unique_ptr<ImageSet> p = Trim(x1, y1, x2, y2, oldWidth);
    for (aocommon::Image& image : p->_images) {
      for (size_t pixel = 0; pixel != image.Size(); ++pixel) {
        if (!mask[pixel]) image[pixel] = 0.0;
      }
    }
    return p;
  }

  void CopyMasked(const ImageSet& fromImageSet, size_t toX, size_t toY,
                  const bool* fromMask) {
    for (size_t i = 0; i != _images.size(); ++i) {
      aocommon::Image::CopyMasked(
          _images[i].Data(), toX, toY, _images[i].Width(),
          fromImageSet._images[i].Data(), fromImageSet._images[i].Width(),
          fromImageSet._images[i].Height(), fromMask);
    }
  }

  /**
   * Place all images in @c from onto the images in this ImageSet at a
   * given position. The dimensions of @c from can be smaller or equal
   * to ones in this.
   */
  void AddSubImage(const ImageSet& from, size_t toX, size_t toY) {
    for (size_t i = 0; i != _images.size(); ++i) {
      aocommon::Image::AddSubImage(_images[i].Data(), toX, toY,
                                   _images[i].Width(), from._images[i].Data(),
                                   from._images[i].Width(),
                                   from._images[i].Height());
    }
  }

  ImageSet& operator*=(float factor) {
    for (aocommon::Image& image : _images) image *= factor;
    return *this;
  }

  ImageSet& operator+=(const ImageSet& other) {
    for (size_t i = 0; i != size(); ++i) _images[i] += other._images[i];
    return *this;
  }

  void FactorAdd(ImageSet& rhs, double factor) {
    for (size_t i = 0; i != size(); ++i)
      _images[i].AddWithFactor(rhs._images[i], factor);
  }

  bool SquareJoinedChannels() const { return _squareJoinedChannels; }

  const std::set<aocommon::PolarizationEnum>& LinkedPolarizations() const {
    return _linkedPolarizations;
  }

  size_t Width() const { return _images.front().Width(); }

  size_t Height() const { return _images.front().Height(); }

  static void CalculateDeconvolutionFrequencies(
      const DeconvolutionTable& groupTable,
      aocommon::UVector<double>& frequencies,
      aocommon::UVector<float>& weights);

 private:
  void initializeIndices();

  void initializePolFactor() {
    const DeconvolutionTable::Group& firstChannelGroup =
        _deconvolutionTable.OriginalGroups().front();
    std::set<aocommon::PolarizationEnum> pols;
    bool all_stokes_without_i = true;
    for (const DeconvolutionTableEntry* entry : firstChannelGroup) {
      if (_linkedPolarizations.empty() ||
          _linkedPolarizations.count(entry->polarization) != 0) {
        if (!aocommon::Polarization::IsStokes(entry->polarization) ||
            entry->polarization == aocommon::Polarization::StokesI)
          all_stokes_without_i = false;
        pols.insert(entry->polarization);
      }
    }
    const bool isDual =
        pols.size() == 2 && aocommon::Polarization::HasDualPolarization(pols);
    const bool isFull =
        pols.size() == 4 &&
        (aocommon::Polarization::HasFullLinearPolarization(pols) ||
         aocommon::Polarization::HasFullCircularPolarization(pols));
    if (all_stokes_without_i)
      _polarizationNormalizationFactor = 1.0 / pols.size();
    else if (isDual || isFull)
      _polarizationNormalizationFactor = 0.5;
    else
      _polarizationNormalizationFactor = 1.0;
  }

  static void copySmallerPart(const aocommon::Image& input,
                              aocommon::Image& output, size_t x1, size_t y1,
                              size_t x2, size_t y2, size_t oldWidth) {
    size_t newWidth = x2 - x1;
    for (size_t y = y1; y != y2; ++y) {
      const float* oldPtr = &input[y * oldWidth];
      float* newPtr = &output[(y - y1) * newWidth];
      for (size_t x = x1; x != x2; ++x) {
        newPtr[x - x1] = oldPtr[x];
      }
    }
  }

  static void copyToLarger(aocommon::Image& to, size_t toX, size_t toY,
                           size_t toWidth, const aocommon::Image& from,
                           size_t fromWidth, size_t fromHeight) {
    for (size_t y = 0; y != fromHeight; ++y) {
      std::copy(from.Data() + y * fromWidth, from.Data() + (y + 1) * fromWidth,
                to.Data() + toX + (toY + y) * toWidth);
    }
  }

  static void copyToLarger(aocommon::Image& to, size_t toX, size_t toY,
                           size_t toWidth, const aocommon::Image& from,
                           size_t fromWidth, size_t fromHeight,
                           const bool* fromMask) {
    for (size_t y = 0; y != fromHeight; ++y) {
      for (size_t x = 0; x != fromWidth; ++x) {
        if (fromMask[y * fromWidth + x])
          to[toX + (toY + y) * toWidth + x] = from[y * fromWidth + x];
      }
    }
  }

  void getSquareIntegratedWithNormalChannels(aocommon::Image& dest,
                                             aocommon::Image& scratch) const;

  void getSquareIntegratedWithSquaredChannels(aocommon::Image& dest) const;

  void getLinearIntegratedWithNormalChannels(aocommon::Image& dest) const;

  const aocommon::Image& entryToImage(
      const DeconvolutionTableEntry& entry) const {
    return _images[_entryIndexToImageIndex[entry.index]];
  }

  std::vector<aocommon::Image> _images;
  // Weight of each deconvolution channels
  aocommon::UVector<float> _weights;
  bool _squareJoinedChannels;
  const DeconvolutionTable& _deconvolutionTable;
  std::vector<size_t> _entryIndexToImageIndex;
  aocommon::UVector<size_t> _imageIndexToPSFIndex;
  float _polarizationNormalizationFactor;
  std::set<aocommon::PolarizationEnum> _linkedPolarizations;
};

#endif  // DECONVOLUTION_IMAGE_SET_H
