#ifndef CACHED_IMAGE_SET_H
#define CACHED_IMAGE_SET_H

#include "logger.h"
#include "../structures/image.h"

#include <aocommon/fits/fitsreader.h>
#include <aocommon/fits/fitswriter.h>

#include <schaapcommon/facets/facet.h>

#include <string.h>
#include <set>
#include <iomanip>

using aocommon::FitsReader;
using aocommon::FitsWriter;
using schaapcommon::facets::Facet;

class CachedImageSet {
 public:
  CachedImageSet() : _polCount(0), _freqCount(0), _facetCount(0), _image() {}

  ~CachedImageSet() {
    for (const std::string& filename : _storedNames)
      std::remove(filename.c_str());
  }

  CachedImageSet(const CachedImageSet& source) = delete;
  CachedImageSet& operator=(const CachedImageSet& source) = delete;

  void Initialize(const FitsWriter& writer, size_t polCount, size_t freqCount,
                  size_t facetCount, const std::string& prefix) {
    _writer = writer;
    _polCount = polCount;
    _freqCount = freqCount;
    _facetCount = facetCount;
    _prefix = prefix;
    _image.reset();
  }

  void SetFitsWriter(const FitsWriter& writer) { _writer = writer; }

  template <typename NumT>
  void Load(NumT* image, aocommon::PolarizationEnum polarization,
            size_t freqIndex, bool isImaginary) const {
    if (_writer.Width() == 0 || _writer.Height() == 0)
      throw std::runtime_error("Writer is not set.");
    Logger::Debug << "Loading " << name(polarization, freqIndex, isImaginary)
                  << '\n';
    if (_polCount == 1 && _freqCount == 1 && _facetCount == 0)
      if (_image.empty())
        throw std::runtime_error("Loading image before store");
      else
        std::copy(_image.data(),
                  _image.data() + _writer.Width() * _writer.Height(), image);
    else {
      FitsReader reader(name(polarization, freqIndex, isImaginary));
      reader.Read(image);
    }
  }

  template <typename NumT>
  void LoadFacet(NumT* image, aocommon::PolarizationEnum polarization,
                 size_t freqIndex, size_t facetIndex,
                 const std::shared_ptr<const Facet>& facet,
                 bool isImaginary) const {
    if (!facet) {
      Load<NumT>(image, polarization, freqIndex, isImaginary);
    } else {
      // No need to check width and height here, because we do not cache
      std::string filename =
          nameFacet(polarization, freqIndex, facetIndex, isImaginary);
      Logger::Debug << "Loading " << filename << '\n';
      FitsReader reader(filename);
      reader.Read(image);
    }
  }

  template <typename NumT>
  void Store(const NumT* image, aocommon::PolarizationEnum polarization,
             size_t freqIndex, bool isImaginary) {
    if (_writer.Width() == 0 || _writer.Height() == 0)
      throw std::runtime_error("Writer is not set.");
    Logger::Debug << "Storing " << name(polarization, freqIndex, isImaginary)
                  << '\n';
    if (_polCount == 1 && _freqCount == 1 && _facetCount == 0) {
      if (_image.empty()) {
        _image = Image(_writer.Width(), _writer.Height());
      }
      std::copy(image, image + _writer.Width() * _writer.Height(),
                _image.data());
    } else {
      std::string filename = name(polarization, freqIndex, isImaginary);
      _writer.Write(filename, image);
      _storedNames.insert(filename);
    }
  }

  template <typename NumT>
  void StoreFacet(const NumT* image, aocommon::PolarizationEnum polarization,
                  size_t freqIndex, size_t facetIndex,
                  const std::shared_ptr<const Facet>& facet, bool isImaginary) {
    if (!facet) {
      // If _facetCount 0, use the main "Store" as is
      Store<NumT>(image, polarization, freqIndex, isImaginary);
    } else {
      std::string filename =
          nameFacet(polarization, freqIndex, facetIndex, isImaginary);
      Logger::Debug << "Storing " << filename << '\n';

      // Initialize FacetWriter, use the trimmed facet width and
      // height as dimensions. The image argument that is fed into
      // the Write() function should have the same size.
      FitsWriter facetWriter;
      facetWriter.SetImageDimensions(facet->GetTrimmedBoundingBox().Width(),
                                     facet->GetTrimmedBoundingBox().Height());
      facetWriter.Write(filename, image);
      _storedNames.insert(filename);
    }
  }

  /**
   * @return The filenames of the temporarily stored files, for testing only.
   */
  const std::set<std::string>& GetStoredNames() const { return _storedNames; };

 private:
  std::string name(aocommon::PolarizationEnum polarization, size_t freqIndex,
                   bool isImaginary) const {
    return nameTrunk(polarization, freqIndex, isImaginary) + "-tmp.fits";
  }

  std::string nameFacet(aocommon::PolarizationEnum polarization,
                        size_t freqIndex, size_t facetIndex,
                        bool isImaginary) const {
    std::ostringstream str;
    str << nameTrunk(polarization, freqIndex, isImaginary);
    if (_facetCount > 0) {
      str << "-f";
      str << std::setw(4) << std::setfill('0') << facetIndex;
    }
    str << "-tmp.fits";
    return str.str();
  }

  std::string nameTrunk(aocommon::PolarizationEnum polarization,
                        size_t freqIndex, bool isImaginary) const {
    if (_freqCount == 1) {
      if (isImaginary)
        return _prefix + '-' +
               aocommon::Polarization::TypeToShortString(polarization) + "i";
      else
        return _prefix + '-' +
               aocommon::Polarization::TypeToShortString(polarization);
    } else {
      std::ostringstream str;
      str << _prefix + '-' +
                 aocommon::Polarization::TypeToShortString(polarization);
      if (isImaginary) str << 'i';
      str << '-';
      str << std::setw(4) << std::setfill('0') << freqIndex;
      return str.str();
    }
  }

  FitsWriter _writer;
  size_t _polCount, _freqCount, _facetCount;
  std::string _prefix;

  Image _image;
  std::set<std::string> _storedNames;
};

#endif
