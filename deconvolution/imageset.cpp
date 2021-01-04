#include "imageset.h"
#include "spectralfitter.h"

#include "../math/nlplfitter.h"

#include "../io/cachedimageset.h"
#include "../io/logger.h"

#include "../structures/primarybeam.h"
#include "../structures/primarybeamimageset.h"

#include <aocommon/staticfor.h>

ImageSet::ImageSet(const ImagingTable* table, const class Settings& settings)
    : _images(),
      _width(0),
      _height(0),
      _channelsInDeconvolution((settings.deconvolutionChannelCount == 0)
                                   ? table->SquaredGroupCount()
                                   : settings.deconvolutionChannelCount),
      _squareJoinedChannels(settings.squaredJoins),
      _imagingTable(*table),
      _imageIndexToPSFIndex(),
      _linkedPolarizations(settings.linkedPolarizations),
      _settings(settings) {
  size_t nPol = table->GetSquaredGroup(0).EntryCount();
  size_t nImages = nPol * _channelsInDeconvolution;
  _images.resize(nImages);
  _imageIndexToPSFIndex.resize(nImages);

  initializePolFactor();
  initializeIndices();
  CalculateDeconvolutionFrequencies(*table, _frequencies, _weights,
                                    _channelsInDeconvolution);
}

ImageSet::ImageSet(const ImagingTable* table, const class Settings& settings,
                   size_t width, size_t height)
    : _images(),
      _width(width),
      _height(height),
      _channelsInDeconvolution((settings.deconvolutionChannelCount == 0)
                                   ? table->SquaredGroupCount()
                                   : settings.deconvolutionChannelCount),
      _squareJoinedChannels(settings.squaredJoins),
      _imagingTable(*table),
      _imageIndexToPSFIndex(),
      _linkedPolarizations(settings.linkedPolarizations),
      _settings(settings) {
  size_t nPol = table->GetSquaredGroup(0).EntryCount();
  size_t nImages = nPol * _channelsInDeconvolution;
  _images.resize(nImages);
  _imageIndexToPSFIndex.resize(nImages);

  initializePolFactor();
  initializeIndices();
  CalculateDeconvolutionFrequencies(*table, _frequencies, _weights,
                                    _channelsInDeconvolution);
  allocateImages();
}

void ImageSet::initializeIndices() {
  size_t lastDeconvolutionChannel = 0;
  size_t deconvolutionChannelStartIndex = 0, lastOutChannel = 0;
  size_t imgIndex = 0;
  for (const ImagingTableEntry& entry : _imagingTable) {
    size_t outChannel = entry.outputChannelIndex;
    size_t chIndex = (outChannel * _channelsInDeconvolution) /
                     _imagingTable.SquaredGroupCount();
    if (outChannel != lastOutChannel && chIndex == lastDeconvolutionChannel) {
      // New output channel maps to an earlier deconvolution channel:
      // start at index of previous deconvolution channel
      imgIndex = deconvolutionChannelStartIndex;
    }
    if (chIndex != lastDeconvolutionChannel) {
      deconvolutionChannelStartIndex = imgIndex;
    }
    _tableIndexToImageIndex.insert(std::make_pair(entry.index, imgIndex));
    lastOutChannel = outChannel;
    lastDeconvolutionChannel = chIndex;
    ++imgIndex;
  }
  for (size_t channel = 0; channel != _channelsInDeconvolution; ++channel) {
    size_t sqIndex = channelToSqIndex(channel);
    ImagingTable subTable = _imagingTable.GetSquaredGroup(sqIndex);
    for (size_t eIndex = 0; eIndex != subTable.EntryCount(); ++eIndex) {
      const ImagingTableEntry& entry = subTable[eIndex];
      size_t imageIndex = _tableIndexToImageIndex.find(entry.index)->second;
      _imageIndexToPSFIndex[imageIndex] = channel;
    }
  }
}

void ImageSet::LoadAndAverage(CachedImageSet& imageSet) {
  for (size_t i = 0; i != _images.size(); ++i) _images[i] = 0.0;

  ImageF scratch(_width, _height);

  /// TODO : use real weights of images
  aocommon::UVector<size_t> weights(_images.size(), 0.0);
  size_t imgIndex = 0;
  for (size_t sqIndex = 0; sqIndex != _imagingTable.SquaredGroupCount();
       ++sqIndex) {
    // The next loop iterates over the polarizations. The logic in the next loop
    // makes sure that images of the same polarizations and that belong to the
    // same deconvolution channel are averaged together.
    size_t imgIndexForChannel = imgIndex;
    ImagingTable subTable = _imagingTable.GetSquaredGroup(sqIndex);
    for (size_t eIndex = 0; eIndex != subTable.EntryCount(); ++eIndex) {
      const ImagingTableEntry& e = subTable[eIndex];
      for (size_t i = 0; i != e.imageCount; ++i) {
        imageSet.Load(scratch.data(), e.polarization, e.outputChannelIndex,
                      i == 1);
        _images[imgIndex] += scratch;
        weights[imgIndex]++;
        ++imgIndex;
      }
    }
    size_t thisChannelIndex = (sqIndex * _channelsInDeconvolution) /
                              _imagingTable.SquaredGroupCount();
    size_t nextChannelIndex = ((sqIndex + 1) * _channelsInDeconvolution) /
                              _imagingTable.SquaredGroupCount();
    // If the next loaded image belongs to the same deconvolution channel as the
    // previously loaded, they need to be averaged together.
    if (thisChannelIndex == nextChannelIndex) imgIndex = imgIndexForChannel;
  }

  for (size_t i = 0; i != _images.size(); ++i)
    _images[i] *= 1.0 / double(weights[i]);
}

void ImageSet::LoadAndAveragePSFs(
    CachedImageSet& psfSet, std::vector<aocommon::UVector<float>>& psfImages,
    aocommon::PolarizationEnum psfPolarization) {
  for (size_t chIndex = 0; chIndex != _channelsInDeconvolution; ++chIndex)
    psfImages[chIndex].assign(_width * _height, 0.0);

  ImageF scratch(_width, _height);

  /// TODO : use real weights of images
  aocommon::UVector<size_t> weights(_channelsInDeconvolution, 0.0);
  for (size_t sqIndex = 0; sqIndex != _imagingTable.SquaredGroupCount();
       ++sqIndex) {
    size_t chIndex = (sqIndex * _channelsInDeconvolution) /
                     _imagingTable.SquaredGroupCount();
    ImagingTable subTable = _imagingTable.GetSquaredGroup(sqIndex);
    const ImagingTableEntry& e = subTable.Front();
    psfSet.Load(scratch.data(), psfPolarization, e.outputChannelIndex, 0);
    for (size_t i = 0; i != _width * _height; ++i)
      psfImages[chIndex][i] += scratch[i];
    weights[chIndex]++;
  }

  for (size_t chIndex = 0; chIndex != ChannelsInDeconvolution(); ++chIndex) {
    for (size_t i = 0; i != _width * _height; ++i) {
      psfImages[chIndex][i] *= 1.0 / double(weights[chIndex]);
    }
  }
}

void ImageSet::InterpolateAndStore(CachedImageSet& imageSet,
                                   const SpectralFitter& fitter) {
  if (_channelsInDeconvolution == _imagingTable.SquaredGroupCount()) {
    directStore(imageSet);
  } else {
    Logger::Info << "Interpolating from " << _channelsInDeconvolution << " to "
                 << _imagingTable.SquaredGroupCount() << " channels...\n";

    // TODO should use spectralimagefitter to do the interpolation of images;
    // here we should just unpack the data structure

    // The following loop will make an 'image' with at each pixel
    // the terms of the fit. By doing this first, it is not necessary
    // to have all channel images in memory at the same time.
    // TODO: this assumes that polarizations are not joined!
    size_t nTerms = fitter.NTerms();
    aocommon::UVector<float> termsImage(_width * _height * nTerms);
    aocommon::StaticFor<size_t> loop(_settings.threadCount);
    loop.Run(0, _width * _height, [&](size_t pxStart, size_t pxEnd) {
      aocommon::UVector<float> spectralPixel(_channelsInDeconvolution);
      aocommon::UVector<float> termsPixel(nTerms);
      for (size_t px = pxStart; px != pxEnd; ++px) {
        bool isZero = true;
        for (size_t s = 0; s != _images.size(); ++s) {
          float value = _images[s][px];
          spectralPixel[s] = value;
          isZero = isZero && (value == 0.0);
        }
        float* termsPtr = &termsImage[px * nTerms];
        // Skip fitting if it is zero; most of model images will be zero, so
        // this can save a lot of time.
        if (isZero) {
          for (float* p = termsPtr; p != termsPtr + nTerms; ++p) *p = 0.0;
        } else {
          fitter.Fit(termsPixel, spectralPixel.data());
          for (size_t i = 0; i != nTerms; ++i) termsPtr[i] = termsPixel[i];
        }
      }
    });

    // Now that we know the fit for each pixel, evaluate the function for each
    // pixel of each output channel.
    Image scratch(_width, _height);
    size_t imgIndex = 0;
    for (size_t eIndex = 0; eIndex != _imagingTable.EntryCount(); ++eIndex) {
      const ImagingTableEntry& e = _imagingTable[eIndex];
      double freq = e.CentralFrequency();
      loop.Run(0, _width * _height, [&](size_t pxStart, size_t pxEnd) {
        aocommon::UVector<float> termsPixel(nTerms);
        for (size_t px = pxStart; px != pxEnd; ++px) {
          const float* termsPtr = &termsImage[px * nTerms];
          for (size_t i = 0; i != nTerms; ++i) termsPixel[i] = termsPtr[i];
          scratch[px] = fitter.Evaluate(termsPixel, freq);
        }
      });

      imageSet.Store(scratch.data(), e.polarization, e.outputChannelIndex,
                     false);
      ++imgIndex;
    }
  }
}

void ImageSet::AssignAndStore(CachedImageSet& imageSet) {
  if (_channelsInDeconvolution == _imagingTable.SquaredGroupCount()) {
    directStore(imageSet);
  } else {
    Logger::Info << "Assigning from " << _channelsInDeconvolution << " to "
                 << _imagingTable.SquaredGroupCount() << " channels...\n";
    size_t imgIndex = 0;
    for (size_t sqIndex = 0; sqIndex != _imagingTable.SquaredGroupCount();
         ++sqIndex) {
      size_t imgIndexForChannel = imgIndex;
      ImagingTable subTable = _imagingTable.GetSquaredGroup(sqIndex);
      for (size_t eIndex = 0; eIndex != subTable.EntryCount(); ++eIndex) {
        const ImagingTableEntry& e = subTable[eIndex];
        for (size_t i = 0; i != e.imageCount; ++i) {
          imageSet.Store(_images[imgIndex].data(), e.polarization,
                         e.outputChannelIndex, i == 1);
          ++imgIndex;
        }
      }
      size_t thisChannelIndex = (sqIndex * _channelsInDeconvolution) /
                                _imagingTable.SquaredGroupCount();
      size_t nextChannelIndex = ((sqIndex + 1) * _channelsInDeconvolution) /
                                _imagingTable.SquaredGroupCount();
      if (thisChannelIndex == nextChannelIndex) imgIndex = imgIndexForChannel;
    }
  }
}

void ImageSet::directStore(CachedImageSet& imageSet) {
  size_t imgIndex = 0;
  for (size_t i = 0; i != _imagingTable.EntryCount(); ++i) {
    const ImagingTableEntry& e = _imagingTable[i];
    for (size_t i = 0; i != e.imageCount; ++i) {
      imageSet.Store(_images[imgIndex].data(), e.polarization,
                     e.outputChannelIndex, i == 1);
      ++imgIndex;
    }
  }
}

void ImageSet::getSquareIntegratedWithNormalChannels(ImageF& dest,
                                                     ImageF& scratch) const {
  // In case only one frequency channel is used, we do not have to use
  // 'scratch', which saves copying and normalizing the data.
  if (_channelsInDeconvolution == 1) {
    ImagingTable subTable = _imagingTable.GetSquaredGroup(0);
    if (subTable.EntryCount() == 1) {
      const ImagingTableEntry& entry = subTable[0];
      size_t imageIndex = _tableIndexToImageIndex.find(entry.index)->second;
      dest = _images[imageIndex];
    } else {
      const bool useAllPolarizations = _linkedPolarizations.empty();
      for (size_t eIndex = 0; eIndex != subTable.EntryCount(); ++eIndex) {
        const ImagingTableEntry& entry = subTable[eIndex];
        if (useAllPolarizations ||
            _linkedPolarizations.count(entry.polarization) != 0) {
          size_t imageIndex = _tableIndexToImageIndex.find(entry.index)->second;
          if (eIndex == 0) {
            dest = _images[0];
            dest.Square();
          } else {
            dest.AddSquared(_images[imageIndex]);
          }
        }
      }
      squareRootMultiply(dest, std::sqrt(_polarizationNormalizationFactor));
    }
  } else {
    double weightSum = 0.0;
    for (size_t chIndex = 0; chIndex != _channelsInDeconvolution; ++chIndex) {
      size_t sqIndex = channelToSqIndex(chIndex);
      ImagingTable subTable = _imagingTable.GetSquaredGroup(sqIndex);
      const double groupWeight = _weights[chIndex];
      weightSum += groupWeight;
      if (subTable.EntryCount() == 1) {
        const ImagingTableEntry& entry = subTable[0];
        size_t imageIndex = _tableIndexToImageIndex.find(entry.index)->second;
        scratch = _images[imageIndex];
      } else {
        const bool useAllPolarizations = _linkedPolarizations.empty();
        for (size_t eIndex = 0; eIndex != subTable.EntryCount(); ++eIndex) {
          const ImagingTableEntry& entry = subTable[eIndex];
          if (useAllPolarizations ||
              _linkedPolarizations.count(entry.polarization) != 0) {
            size_t imageIndex =
                _tableIndexToImageIndex.find(entry.index)->second;
            if (eIndex == 0) {
              scratch = _images[imageIndex];
              scratch.Square();
            } else {
              scratch.AddSquared(_images[imageIndex]);
            }
          }
        }
        scratch.Sqrt();
      }

      if (chIndex == 0)
        assignMultiply(dest, scratch, groupWeight);
      else
        dest.AddWithFactor(scratch, groupWeight);
    }
    if (_channelsInDeconvolution > 0)
      dest *= std::sqrt(_polarizationNormalizationFactor) / weightSum;
    else
      dest = 0.0;
  }
}

void ImageSet::getSquareIntegratedWithSquaredChannels(ImageF& dest) const {
  bool isFirst = true;
  const bool useAllPolarizations = _linkedPolarizations.empty();
  for (size_t channel = 0; channel != _channelsInDeconvolution; ++channel) {
    // TODO this should be weighted
    size_t sqIndex = channelToSqIndex(channel);
    ImagingTable subTable = _imagingTable.GetSquaredGroup(sqIndex);
    for (size_t eIndex = 0; eIndex != subTable.EntryCount(); ++eIndex) {
      const ImagingTableEntry& entry = subTable[eIndex];
      if (useAllPolarizations ||
          _linkedPolarizations.count(entry.polarization) != 0) {
        size_t imageIndex = _tableIndexToImageIndex.find(entry.index)->second;
        if (isFirst) {
          dest = _images[imageIndex];
          dest.Square();
          isFirst = false;
        } else {
          dest.AddSquared(_images[imageIndex]);
        }
      }
    }
  }
  double factor = _channelsInDeconvolution > 0
                      ? sqrt(_polarizationNormalizationFactor) /
                            double(_channelsInDeconvolution)
                      : 0.0;
  squareRootMultiply(dest, factor);
}

void ImageSet::getLinearIntegratedWithNormalChannels(ImageF& dest) const {
  const bool useAllPolarizations = _linkedPolarizations.empty();
  if (_channelsInDeconvolution == 1 &&
      _imagingTable.GetSquaredGroup(0).EntryCount() == 1) {
    ImagingTable subTable = _imagingTable.GetSquaredGroup(0);
    const ImagingTableEntry& entry = subTable[0];
    size_t imageIndex = _tableIndexToImageIndex.find(entry.index)->second;
    dest = _images[imageIndex];
  } else {
    bool isFirst = true;
    double weightSum = 0.0;
    for (size_t channel = 0; channel != _channelsInDeconvolution; ++channel) {
      size_t sqIndex = channelToSqIndex(channel);
      ImagingTable subTable = _imagingTable.GetSquaredGroup(sqIndex);
      const double groupWeight = _weights[channel];
      weightSum += groupWeight;
      for (size_t eIndex = 0; eIndex != subTable.EntryCount(); ++eIndex) {
        const ImagingTableEntry& entry = subTable[eIndex];
        if (useAllPolarizations ||
            _linkedPolarizations.count(entry.polarization) != 0) {
          size_t imageIndex = _tableIndexToImageIndex.find(entry.index)->second;
          if (isFirst) {
            assignMultiply(dest, _images[imageIndex], groupWeight);
            isFirst = false;
          } else {
            dest.AddWithFactor(_images[imageIndex], groupWeight);
          }
        }
      }
    }
    if (weightSum > 0.0)
      dest *= _polarizationNormalizationFactor / weightSum;
    else
      dest = 0.0;
  }
}

void ImageSet::CalculateDeconvolutionFrequencies(
    const ImagingTable& groupTable, aocommon::UVector<double>& frequencies,
    aocommon::UVector<float>& weights, size_t nDeconvolutionChannels) {
  size_t nInputChannels = groupTable.SquaredGroupCount();
  if (nDeconvolutionChannels == 0) nDeconvolutionChannels = nInputChannels;
  frequencies.assign(nDeconvolutionChannels, 0.0);
  weights.assign(nDeconvolutionChannels, 0.0);
  aocommon::UVector<double> weightSums(nDeconvolutionChannels, 0);
  for (size_t i = 0; i != nInputChannels; ++i) {
    const ImagingTableEntry& entry = groupTable.GetSquaredGroup(i)[0];
    double freq = entry.CentralFrequency(), weight = entry.imageWeight;
    size_t deconvolutionChannel = i * nDeconvolutionChannels / nInputChannels;
    frequencies[deconvolutionChannel] += freq * weight;
    weights[deconvolutionChannel] += weight;
  }
  for (size_t i = 0; i != nDeconvolutionChannels; ++i)
    frequencies[i] /= weights[i];
}

void ImageSet::GetIntegratedPSF(ImageF& dest,
                                const aocommon::UVector<const float*>& psfs) {
  // TODO should use weighting!
  std::copy_n(psfs[0], _width * _height, dest.data());
  for (size_t img = 1; img != PSFCount(); ++img) {
    for (size_t i = 0; i != _width * _height; ++i) dest[img] += psfs[img][i];
  }
  if (PSFCount() != 1) {
    for (size_t i = 0; i != _width * _height; ++i)
      dest[i] *= 1.0 / float(PSFCount());
  }
}
