#include "componentlistwriter.h"

#include "../main/settings.h"

#include "../structures/primarybeam.h"

#include <aocommon/image.h>
#include <aocommon/logger.h>

#include <radler/component_list.h>

using aocommon::Logger;

void ComponentListWriter::SaveSourceList(const radler::Radler& deconvolution,
                                         long double phase_centre_ra,
                                         long double phase_centre_dec) const {
  const std::string filename = settings_.prefixName + "-sources.txt";
  radler::ComponentList list = deconvolution.GetComponentList();
  list.WriteSources(deconvolution, filename, settings_.pixelScaleX,
                    settings_.pixelScaleY, phase_centre_ra, phase_centre_dec);
}

void ComponentListWriter::SavePbCorrectedSourceList(
    const radler::Radler& deconvolution, long double phase_centre_ra,
    long double phase_centre_dec) const {
  const std::string filename = settings_.prefixName + "-sources-pb.txt";
  radler::ComponentList list = deconvolution.GetComponentList();

  if (settings_.deconvolutionChannelCount == 0 ||
      settings_.deconvolutionChannelCount ==
          deconvolution_table_->OriginalGroups().size()) {
    // No beam averaging is required
    for (const radler::DeconvolutionTable::Group& channel_group :
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

  list.WriteSources(deconvolution, filename, settings_.pixelScaleX,
                    settings_.pixelScaleY, phase_centre_ra, phase_centre_dec);
}

void ComponentListWriter::CorrectChannelForPrimaryBeam(
    radler::ComponentList& list,
    const radler::DeconvolutionTableEntry& entry) const {
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
  const std::vector<radler::DeconvolutionTable::Group>& channel_groups =
      deconvolution_table_->OriginalGroups();
  for (size_t channel_index = 0; channel_index != channel_groups.size();
       ++channel_index) {
    size_t current_image_index =
        (channel_index * deconvolution_channels) / channel_groups.size();
    if (current_image_index == image_index) {
      const radler::DeconvolutionTableEntry& e =
          *channel_groups[channel_index].front();
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
