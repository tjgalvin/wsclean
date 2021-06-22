#ifndef DECONVOLUTION_IMAGE_SET_H
#define DECONVOLUTION_IMAGE_SET_H

#include <aocommon/uvector.h>

#include "../main/settings.h"

#include "../structures/image.h"
#include "../structures/imagingtable.h"

#include <vector>
#include <map>
#include <memory>

class ImageSet {
 public:
  ImageSet(const ImagingTable* table, const Settings& settings);

  ImageSet(const ImagingTable* table, const Settings& settings, size_t width,
           size_t height);

  void AllocateImages() {
    _images.clear();
    allocateImages();
  }

  void AllocateImages(size_t width, size_t height) {
    _images.clear();
    _width = width;
    _height = height;
    allocateImages();
  }

  Image Release(size_t imageIndex) { return std::move(_images[imageIndex]); }

  void SetImage(size_t imageIndex, Image&& data) {
    _images[imageIndex] = std::move(data);
  }

  bool IsAllocated() const { return _width * _height != 0; }

  void LoadAndAverage(const class CachedImageSet& imageSet);

  void LoadAndAveragePSFs(const class CachedImageSet& psfSet,
                          std::vector<aocommon::UVector<float>>& psfImages,
                          aocommon::PolarizationEnum psfPolarization);

  void InterpolateAndStore(class CachedImageSet& imageSet,
                           const class SpectralFitter& fitter);

  void AssignAndStore(class CachedImageSet& imageSet);

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
  void GetSquareIntegrated(Image& dest, Image& scratch) const {
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
  void GetLinearIntegrated(Image& dest) const {
    if (_squareJoinedChannels)
      getSquareIntegratedWithSquaredChannels(dest);
    else
      getLinearIntegratedWithNormalChannels(dest);
  }

  void GetIntegratedPSF(Image& dest,
                        const aocommon::UVector<const float*>& psfs);

  size_t PSFCount() const { return _channelsInDeconvolution; }

  size_t ChannelsInDeconvolution() const { return _channelsInDeconvolution; }

  ImageSet& operator=(float val) {
    for (Image& image : _images) image = val;
    return *this;
  }

  float* operator[](size_t index) { return _images[index].data(); }

  const float* operator[](size_t index) const { return _images[index].data(); }

  size_t size() const { return _images.size(); }

  size_t PSFIndex(size_t imageIndex) const {
    return _imageIndexToPSFIndex[imageIndex];
  }

  const ImagingTable& Table() const { return _imagingTable; }

  std::unique_ptr<ImageSet> Trim(size_t x1, size_t y1, size_t x2, size_t y2,
                                 size_t oldWidth) const {
    std::unique_ptr<ImageSet> p(
        new ImageSet(&_imagingTable, _settings, x2 - x1, y2 - y1));
    for (size_t i = 0; i != _images.size(); ++i) {
      copySmallerPart(_images[i], p->_images[i], x1, y1, x2, y2, oldWidth);
    }
    return p;
  }

  void Copy(const ImageSet& from, size_t toX, size_t toY, size_t toWidth,
            size_t fromWidth, size_t fromHeight) {
    for (size_t i = 0; i != _images.size(); ++i) {
      copyToLarger(_images[i], toX, toY, toWidth, from._images[i], fromWidth,
                   fromHeight);
    }
  }

  void CopyMasked(const ImageSet& from, size_t toX, size_t toY, size_t toWidth,
                  size_t fromWidth, size_t fromHeight, const bool* fromMask) {
    for (size_t i = 0; i != _images.size(); ++i) {
      copyToLarger(_images[i], toX, toY, toWidth, from._images[i], fromWidth,
                   fromHeight, fromMask);
    }
  }

  ImageSet& operator*=(float factor) {
    for (Image& image : _images) image *= factor;
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

  const class Settings& Settings() const { return _settings; }

  static void CalculateDeconvolutionFrequencies(
      const ImagingTable& groupTable, aocommon::UVector<double>& frequencies,
      aocommon::UVector<float>& weights, size_t nDeconvolutionChannels);

 private:
  ImageSet(const ImageSet&) = delete;
  ImageSet& operator=(const ImageSet&) = delete;

  void allocateImages() {
    for (Image& img : _images) {
      img = Image(_width, _height);
    }
  }

  void assignMultiply(Image& lhs, const Image& rhs, float factor) const {
    for (size_t i = 0; i != _width * _height; ++i) lhs[i] = rhs[i] * factor;
  }

  void squareRootMultiply(Image& image, float factor) const {
    for (size_t i = 0; i != _width * _height; ++i)
      image[i] = std::sqrt(image[i]) * factor;
  }

  void initializeIndices();

  void initializePolFactor() {
    const ImagingTable::Group& firstChannelGroup =
        _imagingTable.SquaredGroups().front();
    std::set<aocommon::PolarizationEnum> pols;
    for (const ImagingTable::EntryPtr& entry : firstChannelGroup) {
      if (_linkedPolarizations.empty() ||
          _linkedPolarizations.count(entry->polarization) != 0) {
        pols.insert(entry->polarization);
      }
    }
    bool isDual =
        pols.size() == 2 && aocommon::Polarization::HasDualPolarization(pols);
    bool isFull = pols.size() == 4 &&
                  (aocommon::Polarization::HasFullLinearPolarization(pols) ||
                   aocommon::Polarization::HasFullCircularPolarization(pols));
    if (isDual || isFull)
      _polarizationNormalizationFactor = 0.5;
    else
      _polarizationNormalizationFactor = 1.0;
  }

  static void copySmallerPart(const Image& input, Image& output, size_t x1,
                              size_t y1, size_t x2, size_t y2,
                              size_t oldWidth) {
    size_t newWidth = x2 - x1;
    for (size_t y = y1; y != y2; ++y) {
      const float* oldPtr = &input[y * oldWidth];
      float* newPtr = &output[(y - y1) * newWidth];
      for (size_t x = x1; x != x2; ++x) {
        newPtr[x - x1] = oldPtr[x];
      }
    }
  }

  static void copyToLarger(Image& to, size_t toX, size_t toY, size_t toWidth,
                           const Image& from, size_t fromWidth,
                           size_t fromHeight) {
    for (size_t y = 0; y != fromHeight; ++y) {
      std::copy(from.data() + y * fromWidth, from.data() + (y + 1) * fromWidth,
                to.data() + toX + (toY + y) * toWidth);
    }
  }

  static void copyToLarger(Image& to, size_t toX, size_t toY, size_t toWidth,
                           const Image& from, size_t fromWidth,
                           size_t fromHeight, const bool* fromMask) {
    for (size_t y = 0; y != fromHeight; ++y) {
      for (size_t x = 0; x != fromWidth; ++x) {
        if (fromMask[y * fromWidth + x])
          to[toX + (toY + y) * toWidth + x] = from[y * fromWidth + x];
      }
    }
  }

  void directStore(class CachedImageSet& imageSet);

  void getSquareIntegratedWithNormalChannels(Image& dest, Image& scratch) const;

  void getSquareIntegratedWithSquaredChannels(Image& dest) const;

  void getLinearIntegratedWithNormalChannels(Image& dest) const;

  size_t channelToSqIndex(size_t channel) const {
    // Calculate reverse of
    // (outChannel*_channelsInDeconvolution)/_imagingTable.SquaredGroups().size();
    size_t fromFloor = channel * _imagingTable.SquaredGroups().size() /
                       _channelsInDeconvolution;
    while (fromFloor * _channelsInDeconvolution /
               _imagingTable.SquaredGroups().size() !=
           channel)
      ++fromFloor;
    return fromFloor;
  }

  const Image& entryToImage(const ImagingTable::EntryPtr& entry) const {
    size_t imageIndex = _entryIndexToImageIndex.find(entry->index)->second;
    return _images[imageIndex];
  }

  std::vector<Image> _images;
  size_t _width, _height, _channelsInDeconvolution;
  // Weight of each deconvolution channels
  aocommon::UVector<float> _weights;
  bool _squareJoinedChannels;
  const ImagingTable& _imagingTable;
  std::map<size_t, size_t> _entryIndexToImageIndex;
  aocommon::UVector<size_t> _imageIndexToPSFIndex;
  float _polarizationNormalizationFactor;
  std::set<aocommon::PolarizationEnum> _linkedPolarizations;
  const class Settings& _settings;
};

#endif  // DECONVOLUTION_IMAGE_SET_H
