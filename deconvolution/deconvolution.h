#ifndef DECONVOLUTION_H
#define DECONVOLUTION_H

#include "paralleldeconvolution.h"

#include <aocommon/polarization.h>
#include <aocommon/uvector.h>

#include <cstring>

class Deconvolution {
 public:
  explicit Deconvolution(const class Settings& settings);
  ~Deconvolution();

  void Perform(const class ImagingTable& groupTable,
               bool& reachedMajorThreshold, size_t majorIterationNr);

  void InitializeDeconvolutionAlgorithm(
      const ImagingTable& groupTable,
      aocommon::PolarizationEnum psfPolarization, double beamSize,
      size_t threadCount);

  void InitializeImages(class CachedImageSet& residuals, CachedImageSet& models,
                        CachedImageSet& psfs) {
    _residualImages = &residuals;
    _modelImages = &models;
    _psfImages = &psfs;
  }

  void FreeDeconvolutionAlgorithms() {
    _parallelDeconvolution.FreeDeconvolutionAlgorithms();
  }

  class DeconvolutionAlgorithm& GetAlgorithm() {
    return _parallelDeconvolution.FirstAlgorithm();
  }
  const DeconvolutionAlgorithm& GetAlgorithm() const {
    return _parallelDeconvolution.FirstAlgorithm();
  }

  bool IsInitialized() const { return _parallelDeconvolution.IsInitialized(); }

  void SaveSourceList(const class ImagingTable& table,
                      long double phaseCentreRA, long double phaseCentreDec) {
    _parallelDeconvolution.SaveSourceList(*_modelImages, table, phaseCentreRA,
                                          phaseCentreDec);
  }

  void SavePBSourceList(const class ImagingTable& table,
                        long double phaseCentreRA, long double phaseCentreDec) {
    _parallelDeconvolution.SavePBSourceList(*_modelImages, table, phaseCentreRA,
                                            phaseCentreDec);
  }

 private:
  void correctChannelForPB(class ComponentList& list,
                           const struct ImagingTableEntry& entry) const;

  void readMask(const ImagingTable& groupTable);

  const class Settings& _settings;

  ParallelDeconvolution _parallelDeconvolution;

  aocommon::UVector<bool> _cleanMask;

  bool _autoMaskIsFinished;
  aocommon::UVector<double> _channelFrequencies;
  aocommon::UVector<float> _channelWeights;
  std::set<aocommon::PolarizationEnum> _polarizations;
  aocommon::PolarizationEnum _psfPolarization;
  size_t _imgWidth, _imgHeight;
  CachedImageSet *_psfImages, *_modelImages, *_residualImages;
  aocommon::UVector<bool> _autoMask;
  double _beamSize, _pixelScaleX, _pixelScaleY;
};

#endif
