#include "componentlistwriter.h"

#include "../main/settings.h"

#include "../multiscale/multiscalealgorithm.h"

#include "../structures/primarybeam.h"

#include <aocommon/image.h>
#include <aocommon/logger.h>

using aocommon::Logger;

void ComponentListWriter::SaveSourceList(const Deconvolution& deconvolution,
                                         long double phase_centre_ra,
                                         long double phase_centre_dec) const {
  const std::string filename = settings_.prefixName + "-sources.txt";
  ComponentList list = deconvolution.GetComponentList();
  const DeconvolutionAlgorithm& deconvolution_algorithm =
      deconvolution.MaxScaleCountAlgorithm();
  WriteSourceList(list, deconvolution_algorithm, filename, phase_centre_ra,
                  phase_centre_dec);
}

void ComponentListWriter::SavePbCorrectedSourceList(
    const Deconvolution& deconvolution, long double phase_centre_ra,
    long double phase_centre_dec) const {
  const std::string filename = settings_.prefixName + "-sources-pb.txt";
  ComponentList list = deconvolution.GetComponentList();

  if (settings_.deconvolutionChannelCount == 0 ||
      settings_.deconvolutionChannelCount ==
          deconvolution_table_->OriginalGroups().size()) {
    // No beam averaging is required
    for (const DeconvolutionTable::Group& channel_group :
         deconvolution_table_->OriginalGroups()) {
      CorrectChannelForPrimaryBeam(list, *channel_group.front());
    }
  } else {
    for (size_t ch = 0; ch != settings_.deconvolutionChannelCount; ++ch) {
      Logger::Debug << "Correcting source list of channel " << ch
                    << " for averaged beam\n";
      PrimaryBeamImageSet beam_images = LoadAveragePrimaryBeam(ch);
      beam_images.CorrectComponentList(list, ch);
    }
  }

  const DeconvolutionAlgorithm& deconvolution_algorithm =
      deconvolution.MaxScaleCountAlgorithm();
  WriteSourceList(list, deconvolution_algorithm, filename, phase_centre_ra,
                  phase_centre_dec);
}

void ComponentListWriter::CorrectChannelForPrimaryBeam(
    ComponentList& list, const DeconvolutionTableEntry& entry) const {
  Logger::Debug << "Correcting source list of channel "
                << entry.original_channel_index << " for beam\n";
  ImageFilename filename(entry.original_channel_index,
                         entry.original_interval_index);
  filename.SetPolarization(entry.polarization);
  PrimaryBeam beam(settings_);
  PrimaryBeamImageSet beam_images = beam.Load(filename);
  beam_images.CorrectComponentList(list, entry.original_channel_index);
}

PrimaryBeamImageSet ComponentListWriter::LoadAveragePrimaryBeam(
    size_t image_index) const {
  Logger::Debug << "Averaging beam for deconvolution channel " << image_index
                << "\n";

  PrimaryBeamImageSet beam_images;

  aocommon::Image scratch(settings_.trimmedImageWidth,
                          settings_.trimmedImageHeight);
  size_t deconvolution_channels = settings_.deconvolutionChannelCount;

  /// TODO : use real weights of images
  size_t count = 0;
  PrimaryBeam beam(settings_);
  const std::vector<DeconvolutionTable::Group>& channel_groups =
      deconvolution_table_->OriginalGroups();
  for (size_t channel_index = 0; channel_index != channel_groups.size();
       ++channel_index) {
    size_t current_image_index =
        (channel_index * deconvolution_channels) / channel_groups.size();
    if (current_image_index == image_index) {
      const DeconvolutionTableEntry& e = *channel_groups[channel_index].front();
      Logger::Debug << "Adding beam at " << e.CentralFrequency() * 1e-6
                    << " MHz\n";
      ImageFilename filename(e.original_channel_index,
                             e.original_interval_index);

      if (count == 0) {
        beam_images = beam.Load(filename);
      } else {
        beam_images += beam.Load(filename);
      }
      count++;
    }
  }
  beam_images *= 1.0 / count;
  return beam_images;
}

void ComponentListWriter::WriteSourceList(
    const ComponentList& list,
    const DeconvolutionAlgorithm& deconvolution_algorithm,
    const std::string& filename, long double phase_centre_ra,
    long double phase_centre_dec) const {
  if (const auto* multiscale_algorithm =
          dynamic_cast<const MultiScaleAlgorithm*>(&deconvolution_algorithm)) {
    list.Write(filename, *multiscale_algorithm, settings_.pixelScaleX,
               settings_.pixelScaleY, phase_centre_ra, phase_centre_dec);
  } else {
    list.WriteSingleScale(filename, deconvolution_algorithm,
                          settings_.pixelScaleX, settings_.pixelScaleY,
                          phase_centre_ra, phase_centre_dec);
  }
}
