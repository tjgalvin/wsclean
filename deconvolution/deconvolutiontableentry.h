#ifndef WSCLEAN_DECONVOLUTION_TABLE_ENTRY_H_
#define WSCLEAN_DECONVOLUTION_TABLE_ENTRY_H_

#include <memory>
#include <vector>

#include <aocommon/imageaccessor.h>
#include <aocommon/polarization.h>

struct DeconvolutionTableEntry {
  double CentralFrequency() const {
    return 0.5 * (bandStartFrequency + bandEndFrequency);
  }

  /**
   * Index of the entry in its DeconvolutionTable.
   */
  size_t index = 0;

  /**
   * Note that mses might have overlapping frequencies.
   */
  double bandStartFrequency = 0.0;
  double bandEndFrequency = 0.0;

  aocommon::PolarizationEnum polarization = aocommon::Polarization::StokesI;

  size_t outputChannelIndex = 0;
  size_t outputIntervalIndex = 0;

  /**
   * The group of entries with equal squaredDeconvolutionIndex should be
   * 'joinedly' deconvolved by adding their squared flux density values
   * together. Normally, all the polarizations from a single (output)channel /
   * timestep form such a group.
   */
  size_t squaredDeconvolutionIndex = 0;

  /**
   * A number that scales with the estimated inverse-variance of the image. It
   * can be used when averaging images or fitting functions through the images
   * to get the optimal sensitivity. It is set after the first inversion.
   */
  double imageWeight = 0.0;

  /**
   * Image accessor for the PSF image for this entry. This accessor is only used
   * for the first entry of each channel group.
   */
  std::unique_ptr<aocommon::ImageAccessor> psf_accessor;

  /**
   * Image accessor for the model image for this entry.
   */
  std::unique_ptr<aocommon::ImageAccessor> model_accessor;

  /**
   * Image accessor for the residual image for this entry.
   */
  std::unique_ptr<aocommon::ImageAccessor> residual_accessor;
};

#endif
