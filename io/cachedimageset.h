#ifndef CACHED_IMAGE_SET_H
#define CACHED_IMAGE_SET_H

#include <aocommon/fits/fitsreader.h>
#include <aocommon/fits/fitswriter.h>
#include <aocommon/image.h>
#include <aocommon/logger.h>

#include <schaapcommon/facets/facet.h>
#include <schaapcommon/facets/facetimage.h>

#include "../io/fitsiomutex.h"
#include "../io/wscfitswriter.h"

#include <iomanip>
#include <mutex>
#include <set>
#include <string.h>

namespace wsclean {

/**
 * Temporary image cache used while stitching and during intermediate imaging
 * steps.
 *
 * Concurrent Load()/Store() operations are supported after the cache has been
 * initialized. Initialize() reconfigures the cache state and may therefore not
 * be called concurrently with other methods.
 */
class CachedImageSet {
 public:
  CachedImageSet() : _polCount(0), _freqCount(0), _facetCount(0), _image() {}

  ~CachedImageSet() {
    for (const std::string& filename : _storedNames) std::remove(filename.c_str());
  }

  CachedImageSet(const CachedImageSet& source) = delete;
  CachedImageSet& operator=(const CachedImageSet& source) = delete;

  void Initialize(const WSCFitsWriter& writer, size_t polCount, size_t freqCount, size_t facetCount,
                  const std::string& prefix) {
    std::lock_guard<std::mutex> lock(mutex_);
    _writer = writer;
    _polCount = polCount;
    _freqCount = freqCount;
    _facetCount = facetCount;
    _prefix = prefix;
    _image.Reset();
  }

  void Initialize(const aocommon::FitsWriter& writer, size_t polCount, size_t freqCount,
                  size_t facetCount, const std::string& prefix) {
    Initialize(WSCFitsWriter(writer), polCount, freqCount, facetCount, prefix);
  }

  void SetWSCFitsWriter(const WSCFitsWriter& writer) {
    std::lock_guard<std::mutex> lock(mutex_);
    _writer = writer;
  }

  std::pair<size_t, size_t> GetImageDimensions() const {
    std::lock_guard<std::mutex> lock(mutex_);
    CheckWriterIsSetLocked();
    if (_writer->Width() == 0 || _writer->Height() == 0) {
      throw std::runtime_error("Writer is not set.");
    }
    return {_writer->Width(), _writer->Height()};
  }

  void SetImageDimensions(size_t width, size_t height) {
    std::lock_guard<std::mutex> lock(mutex_);
    CheckWriterIsSetLocked();
    _writer->Writer().SetImageDimensions(width, height);
  }

  template <typename NumT>
  void Load(NumT* image, aocommon::PolarizationEnum polarization, size_t freqIndex,
            bool isImaginary) const {
    if (IsSingleImageMode()) {
      const std::pair<size_t, size_t> image_dimensions = GetImageDimensions();
      std::lock_guard<std::mutex> lock(mutex_);
      assert(!isImaginary);
      if (_image.Empty()) {
        throw std::runtime_error("Loading image before store");
      }
      std::copy(_image.Data(), _image.Data() + image_dimensions.first * image_dimensions.second,
                image);
      return;
    }
    const std::string filename = name(polarization, freqIndex, isImaginary);
    aocommon::Logger::Debug << "Loading " << filename << '\n';
    {
      std::lock_guard<std::mutex> fits_lock(io::FitsIOMutex());
      aocommon::FitsReader reader(filename);
      reader.Read(image);
    }
  }

  template <typename NumT>
  void LoadFacet(NumT* image, aocommon::PolarizationEnum polarization, size_t freqIndex,
                 size_t facetIndex, const std::shared_ptr<const schaapcommon::facets::Facet>& facet,
                 bool isImaginary) const {
    if (!facet) {
      Load<NumT>(image, polarization, freqIndex, isImaginary);
    } else {
      // No need to check width and height here, because we do not cache
      const std::string filename = FacetFilename(polarization, freqIndex, facetIndex, isImaginary);
      aocommon::Logger::Debug << "Loading " << filename << '\n';
      std::lock_guard<std::mutex> fits_lock(io::FitsIOMutex());
      aocommon::FitsReader reader(filename);
      reader.Read(image);
    }
  }

  /**
   * @brief Store image from raw pointer
   *
   * @tparam NumT Numeric type, template parameter, however if there is only
   * one polarization and and frequency in cache then the data is stored in
   * memory in an Image object which is of type float. So in that case NumT
   * should be float.
   * @param image pointer to data of type NumT
   */
  template <typename NumT>
  void Store(const NumT* image, aocommon::PolarizationEnum polarization, size_t freqIndex,
             bool isImaginary) {
    WSCFitsWriter writer = GetWriterCopy();
    Store(image, writer, polarization, freqIndex, isImaginary);
  }

  // These overloads accept a caller-owned writer copy for hot paths where each
  // worker thread already keeps its own WSCFitsWriter instance.
  template <typename NumT>
  void Store(const NumT* image, WSCFitsWriter& writer, aocommon::PolarizationEnum polarization,
             size_t freqIndex, bool isImaginary) {
    if (IsSingleImageMode()) {
      assert(!isImaginary);
      const size_t width = writer.Width();
      const size_t height = writer.Height();
      std::lock_guard<std::mutex> lock(mutex_);
      if (_image.Empty()) {
        _image = aocommon::Image(width, height);
      }
      std::copy(image, image + width * height, _image.Data());
    } else {
      const std::string filename = name(polarization, freqIndex, isImaginary);
      aocommon::Logger::Debug << "Storing " << filename << '\n';
      {
        std::lock_guard<std::mutex> fits_lock(io::FitsIOMutex());
        writer.Writer().Write(filename, image);
      }
      AddStoredName(filename);
    }
  }

  /**
   * @brief Store an Image object.
   */
  void Store(const aocommon::Image& image, aocommon::PolarizationEnum polarization,
             size_t freqIndex, bool isImaginary) {
    WSCFitsWriter writer = GetWriterCopy();
    Store(image, writer, polarization, freqIndex, isImaginary);
  }

  void Store(const aocommon::Image& image, WSCFitsWriter& writer,
             aocommon::PolarizationEnum polarization, size_t freqIndex, bool isImaginary) {
    if (IsSingleImageMode()) {
      assert(!isImaginary);
      const size_t width = writer.Width();
      const size_t height = writer.Height();
      std::lock_guard<std::mutex> lock(mutex_);
      if (_image.Empty()) {
        _image = aocommon::Image(width, height);
      }
      std::copy(image.Data(), image.Data() + width * height, _image.Data());
    } else {
      const std::string filename = name(polarization, freqIndex, isImaginary);
      aocommon::Logger::Debug << "Storing " << filename << '\n';
      writer.WriteFullNameImage(filename, image);
      AddStoredName(filename);
    }
  }

  /**
   * @brief Store a facet image
   *
   * main Store is used if no facet is given,
   * otherwise facet parameter is used to shift the coordinate system
   */
  void StoreFacet(const aocommon::Image& image, aocommon::PolarizationEnum polarization,
                  size_t freqIndex, size_t facetIndex,
                  const std::shared_ptr<const schaapcommon::facets::Facet>& facet,
                  bool isImaginary) {
    WSCFitsWriter writer = GetWriterCopy();
    StoreFacet(image, writer, polarization, freqIndex, facetIndex, facet, isImaginary);
  }

  void StoreFacet(const aocommon::Image& image, WSCFitsWriter& writer,
                  aocommon::PolarizationEnum polarization, size_t freqIndex, size_t facetIndex,
                  const std::shared_ptr<const schaapcommon::facets::Facet>& facet,
                  bool isImaginary) {
    if (!facet) {
      // If _facetCount 0, use the main "Store" as is
      Store(image, writer, polarization, freqIndex, isImaginary);
    } else {
      const std::string filename = FacetFilename(polarization, freqIndex, facetIndex, isImaginary);
      aocommon::Logger::Debug << "Storing " << filename << '\n';
      writer.WriteFullNameImage(filename, image, *facet);
      AddStoredName(filename);
    }
  }

  /**
   * @brief Store a facet image
   *
   * @param facetimage contains both image data and facet metadata
   */
  void StoreFacet(const schaapcommon::facets::FacetImage& facetimage,
                  aocommon::PolarizationEnum polarization, size_t freqIndex, size_t facetIndex,
                  bool isImaginary) {
    WSCFitsWriter writer = GetWriterCopy();
    StoreFacet(facetimage, writer, polarization, freqIndex, facetIndex, isImaginary);
  }

  void StoreFacet(const schaapcommon::facets::FacetImage& facetimage, WSCFitsWriter& writer,
                  aocommon::PolarizationEnum polarization, size_t freqIndex, size_t facetIndex,
                  bool isImaginary) {
    const std::string filename = FacetFilename(polarization, freqIndex, facetIndex, isImaginary);
    aocommon::Logger::Debug << "Storing " << filename << '\n';
    writer.WriteFullNameImage(filename, facetimage);
    AddStoredName(filename);
  }

  /**
   * A CachedImageSet is empty as long as Store() has not been called.
   */
  bool Empty() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return _storedNames.empty() && _image.Empty();
  }

  /**
   * @return The filenames of the temporarily stored files, for testing only.
   */
  std::set<std::string> GetStoredNames() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return _storedNames;
  }

 private:
  void CheckWriterIsSetLocked() const {
    if (!_writer) throw std::runtime_error("Writer is not set.");
  }

  WSCFitsWriter GetWriterCopy() const {
    std::lock_guard<std::mutex> lock(mutex_);
    CheckWriterIsSetLocked();
    return *_writer;
  }

  bool IsSingleImageMode() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return _polCount == 1 && _freqCount == 1 && _facetCount == 0;
  }

  void AddStoredName(const std::string& filename) {
    std::lock_guard<std::mutex> lock(mutex_);
    _storedNames.insert(filename);
  }

  std::string name(aocommon::PolarizationEnum polarization, size_t freqIndex,
                   bool isImaginary) const {
    std::lock_guard<std::mutex> lock(mutex_);
    return nameTrunkLocked(polarization, freqIndex, isImaginary) + "-tmp.fits";
  }

  std::string FacetFilename(aocommon::PolarizationEnum polarization, size_t freqIndex,
                            size_t facetIndex, bool isImaginary) const {
    std::lock_guard<std::mutex> lock(mutex_);
    std::ostringstream str;
    str << nameTrunkLocked(polarization, freqIndex, isImaginary);
    if (_facetCount > 0) {
      str << "-f";
      str << std::setw(4) << std::setfill('0') << facetIndex;
    }
    str << "-tmp.fits";
    return str.str();
  }

  std::string nameTrunkLocked(aocommon::PolarizationEnum polarization, size_t freq_index,
                              bool is_imaginary) const {
    const std::string pol_string =
        (_polCount == 1) ? "" : '-' + aocommon::Polarization::TypeToShortString(polarization);
    if (_freqCount == 1) {
      if (is_imaginary)
        return _prefix + pol_string + "i";
      else
        return _prefix + pol_string;
    } else {
      std::ostringstream str;
      str << _prefix + pol_string;
      if (is_imaginary) str << 'i';
      str << '-';
      str << std::setw(4) << std::setfill('0') << freq_index;
      return str.str();
    }
  }

  std::optional<WSCFitsWriter> _writer;
  size_t _polCount, _freqCount, _facetCount;
  std::string _prefix;

  aocommon::Image _image;
  std::set<std::string> _storedNames;
  mutable std::mutex mutex_;
};

}  // namespace wsclean

#endif
