#ifndef COMPONENT_LIST_H
#define COMPONENT_LIST_H

#include "imageset.h"

#include "../structures/image.h"

#include <aocommon/uvector.h>

#include <vector>

class ComponentList {
 public:
  /**
   * Constructor for single-scale clean
   */
  ComponentList(size_t width, size_t height, ImageSet& imageSet)
      : _width(width),
        _height(height),
        _nFrequencies(imageSet.size()),
        _componentsAddedSinceLastMerge(0),
        _maxComponentsBeforeMerge(100000),
        _listPerScale(1) {
    loadFromImageSet(imageSet, 0);
  }

  /**
   * Constructor for multi-scale clean
   */
  ComponentList(size_t width, size_t height, size_t nScales,
                size_t nFrequencies)
      : _width(width),
        _height(height),
        _nFrequencies(nFrequencies),
        _componentsAddedSinceLastMerge(0),
        _maxComponentsBeforeMerge(100000),
        _listPerScale(nScales) {}

  void Add(size_t x, size_t y, size_t scaleIndex, const float* values) {
    _listPerScale[scaleIndex].values.push_back(values, values + _nFrequencies);
    _listPerScale[scaleIndex].positions.emplace_back(x, y);
    ++_componentsAddedSinceLastMerge;
    if (_componentsAddedSinceLastMerge >= _maxComponentsBeforeMerge)
      MergeDuplicates();
  }

  void Add(const ComponentList& other, int offsetX, int offsetY) {
    if (other._nFrequencies != _nFrequencies)
      throw std::runtime_error(
          "Add(ComponentList...) called with incorrect frequency count");
    if (other.NScales() > NScales()) SetNScales(other.NScales());
    for (size_t scale = 0; scale != other.NScales(); ++scale) {
      const ScaleList& list = other._listPerScale[scale];
      for (size_t i = 0; i != list.positions.size(); ++i) {
        Add(list.positions[i].x + offsetX, list.positions[i].y + offsetY, scale,
            &list.values[i * _nFrequencies]);
      }
    }
  }

  void Write(const std::string& filename,
             const class MultiScaleAlgorithm& multiscale,
             long double pixelScaleX, long double pixelScaleY,
             long double phaseCentreRA, long double phaseCentreDec);

  void WriteSingleScale(const std::string& filename,
                        const class DeconvolutionAlgorithm& algorithm,
                        long double pixelScaleX, long double pixelScaleY,
                        long double phaseCentreRA, long double phaseCentreDec);

  void MergeDuplicates() {
    for (size_t scaleIndex = 0; scaleIndex != _listPerScale.size();
         ++scaleIndex) {
      mergeDuplicates(scaleIndex);
    }
    _componentsAddedSinceLastMerge = 0;
  }

  void Clear() {
    for (ScaleList& list : _listPerScale) {
      list.positions.clear();
      list.values.clear();
    }
  }

  size_t Width() const { return _width; }
  size_t Height() const { return _height; }

  size_t ComponentCount(size_t scaleIndex) const {
    return _listPerScale[scaleIndex].positions.size();
  }

  void GetComponent(size_t scaleIndex, size_t index, size_t& x, size_t& y,
                    float* values) const {
    x = _listPerScale[scaleIndex].positions[index].x;
    y = _listPerScale[scaleIndex].positions[index].y;
    for (size_t f = 0; f != _nFrequencies; ++f)
      values[f] = _listPerScale[scaleIndex].values[index * _nFrequencies + f];
  }

  void CorrectForBeam(class PrimaryBeamImageSet& beam, size_t channel);

  size_t NScales() const { return _listPerScale.size(); }

  void SetNScales(size_t nScales) { _listPerScale.resize(nScales); }

 private:
  struct Position {
    Position(size_t _x, size_t _y) : x(_x), y(_y) {}
    size_t x, y;
  };
  struct ScaleList {
    /**
     * This list contains nFrequencies values for each
     * component, such that _positions[i] corresponds with the values
     * starting at _values[i * _nFrequencies].
     */
    aocommon::UVector<float> values;
    aocommon::UVector<Position> positions;
  };

  void write(const std::string& filename,
             const class DeconvolutionAlgorithm& algorithm,
             const aocommon::UVector<double>& scaleSizes,
             long double pixelScaleX, long double pixelScaleY,
             long double phaseCentreRA, long double phaseCentreDec);

  void loadFromImageSet(ImageSet& imageSet, size_t scaleIndex);

  void mergeDuplicates(size_t scaleIndex) {
    ScaleList& list = _listPerScale[scaleIndex];
    aocommon::UVector<float> newValues;
    aocommon::UVector<Position> newPositions;

    std::vector<Image> images(_nFrequencies);
    for (Image& image : images) image = Image(_width, _height, 0.0);
    size_t valueIndex = 0;
    for (size_t index = 0; index != list.positions.size(); ++index) {
      size_t position =
          list.positions[index].x + list.positions[index].y * _width;
      for (size_t frequency = 0; frequency != _nFrequencies; ++frequency) {
        images[frequency][position] += list.values[valueIndex];
        valueIndex++;
      }
    }

    list.values.clear();
    list.positions.clear();

    for (size_t imageIndex = 0; imageIndex != images.size(); ++imageIndex) {
      Image& image = images[imageIndex];
      size_t posIndex = 0;
      for (size_t y = 0; y != _height; ++y) {
        for (size_t x = 0; x != _width; ++x) {
          if (image[posIndex] != 0.0) {
            for (size_t i = 0; i != images.size(); ++i) {
              newValues.push_back(images[i][posIndex]);
              images[i][posIndex] = 0.0;
            }
            newPositions.emplace_back(x, y);
          }
          ++posIndex;
        }
      }
    }
    std::swap(_listPerScale[scaleIndex].values, newValues);
    std::swap(_listPerScale[scaleIndex].positions, newPositions);
  }
  const size_t _width, _height;
  size_t _nFrequencies;
  size_t _componentsAddedSinceLastMerge;
  size_t _maxComponentsBeforeMerge;
  std::vector<ScaleList> _listPerScale;
};

#endif
