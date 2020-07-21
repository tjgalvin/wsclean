#include "subminorloop.h"

#include "../deconvolution/spectralfitter.h"
#include "../deconvolution/componentlist.h"

#include "../fftconvolver.h"
#include "../image.h"

#include "../wsclean/logger.h"

template <bool AllowNegatives>
size_t SubMinorModel::GetMaxComponent(double* scratch, double& maxValue) const {
  _residual->GetLinearIntegrated(scratch);
  if (!_rmsFactorImage.empty()) {
    for (size_t i = 0; i != size(); ++i) scratch[i] *= _rmsFactorImage[i];
  }
  size_t maxComponent = 0;
  maxValue = scratch[0];
  for (size_t i = 0; i != size(); ++i) {
    double value;
    if (AllowNegatives)
      value = std::fabs(scratch[i]);
    else
      value = scratch[i];
    if (value > maxValue) {
      maxComponent = i;
      maxValue = value;
    }
  }
  maxValue = scratch[maxComponent];  // If it was negative, make sure a negative
                                     // value is returned
  return maxComponent;
}

boost::optional<double> SubMinorLoop::Run(
    ImageSet& convolvedResidual,
    const aocommon::UVector<const double*>& doubleConvolvedPsfs) {
  _subMinorModel = SubMinorModel(_width, _height);

  findPeakPositions(convolvedResidual);

  _subMinorModel.MakeSets(convolvedResidual);
  if (!_rmsFactorImage.empty())
    _subMinorModel.MakeRMSFactorImage(_rmsFactorImage);
  _logReceiver.Debug << "Number of components selected > " << _threshold << ": "
                     << _subMinorModel.size() << '\n';

  if (_subMinorModel.size() == 0) return boost::optional<double>();

  aocommon::UVector<double> scratch(_subMinorModel.size());
  double maxValue;
  size_t maxComponent = _subMinorModel.GetMaxComponent(
      scratch.data(), maxValue, _allowNegativeComponents);

  while (std::fabs(maxValue) > _threshold &&
         _currentIteration < _maxIterations &&
         (!_stopOnNegativeComponent || maxValue >= 0.0)) {
    aocommon::UVector<double> componentValues(_subMinorModel.Residual().size());
    for (size_t imgIndex = 0; imgIndex != _subMinorModel.Residual().size();
         ++imgIndex)
      componentValues[imgIndex] =
          _subMinorModel.Residual()[imgIndex][maxComponent] * _gain;
    _fluxCleaned += maxValue * _gain;

    if (_fitter) _fitter->FitAndEvaluate(componentValues.data());

    for (size_t imgIndex = 0; imgIndex != _subMinorModel.Model().size();
         ++imgIndex)
      _subMinorModel.Model()[imgIndex][maxComponent] +=
          componentValues[imgIndex];

    size_t x = _subMinorModel.X(maxComponent),
           y = _subMinorModel.Y(maxComponent);
    /*
      Commented out because even in verbose mode this is a bit too verbose, but
    useful in case divergence occurs: _logReceiver.Debug << x << ", " << y << "
    " << maxValue << " -> "; for(size_t imgIndex=0;
    imgIndex!=_clarkModel.Model().size(); ++imgIndex) _logReceiver.Debug <<
    componentValues[imgIndex] << ' '; _logReceiver.Debug << '\n';
    */
    for (size_t imgIndex = 0; imgIndex != _subMinorModel.Residual().size();
         ++imgIndex) {
      double* image = _subMinorModel.Residual()[imgIndex];
      const double* psf =
          doubleConvolvedPsfs[_subMinorModel.Residual().PSFIndex(imgIndex)];
      double psfFactor = componentValues[imgIndex];
      for (size_t px = 0; px != _subMinorModel.size(); ++px) {
        int psfX = _subMinorModel.X(px) - x + _width / 2;
        int psfY = _subMinorModel.Y(px) - y + _height / 2;
        if (psfX >= 0 && psfX < int(_width) && psfY >= 0 && psfY < int(_height))
          image[px] -= psf[psfX + psfY * _width] * psfFactor;
      }
    }

    maxComponent = _subMinorModel.GetMaxComponent(scratch.data(), maxValue,
                                                  _allowNegativeComponents);
    ++_currentIteration;
  }
  return maxValue;
}

void SubMinorModel::MakeSets(const ImageSet& residualSet) {
  _residual.reset(
      new ImageSet(&residualSet.Table(), residualSet.Settings(), size(), 1));
  _model.reset(
      new ImageSet(&residualSet.Table(), residualSet.Settings(), size(), 1));
  for (size_t imgIndex = 0; imgIndex != _model->size(); ++imgIndex) {
    std::fill((*_model)[imgIndex], (*_model)[imgIndex] + size(), 0.0);

    const double* sourceResidual = residualSet[imgIndex];
    double* destResidual = (*_residual)[imgIndex];
    for (size_t pxIndex = 0; pxIndex != size(); ++pxIndex) {
      size_t srcIndex =
          _positions[pxIndex].second * _width + _positions[pxIndex].first;
      destResidual[pxIndex] = sourceResidual[srcIndex];
    }
  }
}

void SubMinorModel::MakeRMSFactorImage(Image& rmsFactorImage) {
  _rmsFactorImage = Image(size(), 1);
  for (size_t pxIndex = 0; pxIndex != size(); ++pxIndex) {
    size_t srcIndex =
        _positions[pxIndex].second * _width + _positions[pxIndex].first;
    _rmsFactorImage[pxIndex] = rmsFactorImage[srcIndex];
  }
}

void SubMinorLoop::findPeakPositions(ImageSet& convolvedResidual) {
  Image integratedScratch(_width, _height);
  convolvedResidual.GetLinearIntegrated(integratedScratch.data());

  if (!_rmsFactorImage.empty()) {
    integratedScratch *= _rmsFactorImage;
  }

  const size_t xiStart = _horizontalBorder,
               xiEnd = std::max<long>(xiStart, _width - _horizontalBorder),
               yiStart = _verticalBorder,
               yiEnd = std::max<long>(yiStart, _height - _verticalBorder);

  if (_mask) {
    for (size_t y = yiStart; y != yiEnd; ++y) {
      const bool* maskPtr = _mask + y * _width;
      double* imagePtr = integratedScratch.data() + y * _width;
      for (size_t x = xiStart; x != xiEnd; ++x) {
        double value;
        if (_allowNegativeComponents)
          value = fabs(imagePtr[x]);
        else
          value = imagePtr[x];
        if (value >= _threshold && maskPtr[x]) _subMinorModel.AddPosition(x, y);
      }
    }
  } else {
    for (size_t y = yiStart; y != yiEnd; ++y) {
      double* imagePtr = integratedScratch.data() + y * _width;
      for (size_t x = xiStart; x != xiEnd; ++x) {
        double value;
        if (_allowNegativeComponents)
          value = fabs(imagePtr[x]);
        else
          value = imagePtr[x];
        if (value >= _threshold) _subMinorModel.AddPosition(x, y);
      }
    }
  }
}

void SubMinorLoop::GetFullIndividualModel(size_t imageIndex,
                                          double* individualModelImg) const {
  std::fill(individualModelImg, individualModelImg + _width * _height, 0.0);
  const double* data = _subMinorModel.Model()[imageIndex];
  for (size_t px = 0; px != _subMinorModel.size(); ++px) {
    individualModelImg[_subMinorModel.FullIndex(px)] = data[px];
  }
}

void SubMinorLoop::CorrectResidualDirty(
    class FFTWManager& fftw, double* scratchA, double* scratchB,
    double* scratchC, size_t imageIndex, double* residual,
    const double* singleConvolvedPsf) const {
  // Get padded kernel in scratchB
  Image::Untrim(scratchA, _paddedWidth, _paddedHeight, singleConvolvedPsf,
                _width, _height);
  FFTConvolver::PrepareKernel(scratchB, scratchA, _paddedWidth, _paddedHeight);

  // Get padded model image in scratchA
  GetFullIndividualModel(imageIndex, scratchC);
  Image::Untrim(scratchA, _paddedWidth, _paddedHeight, scratchC, _width,
                _height);

  // Convolve and store in scratchA
  FFTConvolver::ConvolveSameSize(fftw, scratchA, scratchB, _paddedWidth,
                                 _paddedHeight);

  // Trim the result into scratchC
  Image::Trim(scratchC, _width, _height, scratchA, _paddedWidth, _paddedHeight);

  for (size_t i = 0; i != _width * _height; ++i) residual[i] -= scratchC[i];
}

void SubMinorLoop::UpdateAutoMask(bool* mask) const {
  for (size_t imageIndex = 0; imageIndex != _subMinorModel.Model().size();
       ++imageIndex) {
    const double* image = _subMinorModel.Model()[imageIndex];
    for (size_t px = 0; px != _subMinorModel.size(); ++px) {
      if (image[px] != 0.0) mask[_subMinorModel.FullIndex(px)] = true;
    }
  }
}

void SubMinorLoop::UpdateComponentList(class ComponentList& list,
                                       size_t scaleIndex) const {
  aocommon::UVector<double> values(_subMinorModel.Model().size());
  for (size_t px = 0; px != _subMinorModel.size(); ++px) {
    bool isNonZero = false;
    for (size_t imageIndex = 0; imageIndex != _subMinorModel.Model().size();
         ++imageIndex) {
      values[imageIndex] = _subMinorModel.Model()[imageIndex][px];
      if (values[imageIndex] != 0.0) isNonZero = true;
    }
    if (isNonZero) {
      size_t posIndex = _subMinorModel.FullIndex(px);
      size_t x = posIndex % _width, y = posIndex / _width;
      list.Add(x, y, scaleIndex, values.data());
    }
  }
}
