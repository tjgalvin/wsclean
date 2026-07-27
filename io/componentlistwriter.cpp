#include "componentlistwriter.h"

#include "../main/settings.h"

#include "../structures/primarybeam.h"
#include "../structures/imagingtable.h"

#include <aocommon/image.h>
#include <aocommon/logger.h>

#include <schaapcommon/facets/facet.h>

#include <radler/component_list.h>

using aocommon::Logger;

namespace wsclean {
namespace {

std::vector<std::vector<double>> AverageCorrectionFactors(
    const std::vector<std::vector<double>>& correction_factors, size_t output_channels) {
  const size_t input_channels = correction_factors.size();
  assert(input_channels >= output_channels);
  std::vector<std::vector<double>> result(output_channels);
  std::vector<size_t> counts(output_channels, 0);

  // TODO use weights
  for (size_t input_index = 0; input_index != input_channels; ++input_index) {
    const size_t output_index = (input_index * output_channels) / input_channels;
    const size_t n_facets = correction_factors[input_index].size();
    result[output_index].resize(n_facets);
    counts[output_index]++;
    for (size_t facet_index = 0; facet_index != n_facets; ++facet_index) {
      result[output_index][facet_index] += correction_factors[input_index][facet_index];
    }
  }
  for (size_t output_index = 0; output_index != output_channels; ++output_index) {
    for (double& facet_factor : result[output_index]) {
      facet_factor /= counts[output_index];
    }
  }
  return result;
}

}  // namespace

void ComponentListWriter::SaveSourceList(const radler::Radler& deconvolution,
                                         long double phase_centre_ra, long double phase_centre_dec,
                                         long double l_shift, long double m_shift) const {
  const std::string filename = settings_.prefixName + "-sources.txt";
  radler::ComponentList list = deconvolution.GetComponentList();
  list.WriteSources(deconvolution, filename, settings_.pixelScaleX, settings_.pixelScaleY,
                    phase_centre_ra, phase_centre_dec, l_shift, m_shift);
}

void ComponentListWriter::SavePbCorrectedSourceList(const radler::Radler& deconvolution,
                                                    long double phase_centre_ra,
                                                    long double phase_centre_dec,
                                                    long double l_shift,
                                                    long double m_shift) const {
  const std::string filename = settings_.prefixName + "-sources-pb.txt";
  radler::ComponentList list = deconvolution.GetComponentList();

  if (settings_.deconvolutionChannelCount == 0 ||
      settings_.deconvolutionChannelCount == deconvolution_table_->OriginalGroups().size()) {
    // No beam averaging is required
    for (const radler::WorkTable::Group& channel_group : deconvolution_table_->OriginalGroups()) {
      CorrectChannelForPrimaryBeam(list, *channel_group.front());
    }
  } else {
    for (size_t ch = 0; ch != settings_.deconvolutionChannelCount; ++ch) {
      Logger::Debug << "Correcting source list of channel " << ch << " for averaged beam\n";
      PrimaryBeamImageSet beam_images = LoadAveragePrimaryBeam(ch);
      beam_images.CorrectComponentList(list, ch);
    }
  }

  list.WriteSources(deconvolution, filename, settings_.pixelScaleX, settings_.pixelScaleY,
                    phase_centre_ra, phase_centre_dec, l_shift, m_shift);
}

std::vector<std::vector<double>> GetFacetPbCorrectionMap(
    std::span<OutputChannelInfo> channel_infos) {
  std::vector<std::vector<double>> result;
  for (const OutputChannelInfo& channel_info : channel_infos) {
    std::vector<double>& channel_values = result.emplace_back();
    channel_values.reserve(channel_info.averageFacetCorrection.size());
    for (const AverageCorrection& correction : channel_info.averageFacetCorrection) {
      const double stokes_i = correction.GetStokesIValue();
      channel_values.emplace_back(stokes_i != 0.0 ? 1.0 / stokes_i : 0.0);
    }
  }
  return result;
}

void ComponentListWriter::SaveFacetPbCorrectedSourceList(
    const radler::Radler& deconvolution, const std::vector<std::vector<double>>& correction_factors,
    const ImagingTable& imaging_table, long double phase_centre_ra, long double phase_centre_dec,
    long double l_shift, long double m_shift) const {
  const std::string filename = settings_.prefixName + "-sources-fpb.txt";
  radler::ComponentList component_list = deconvolution.GetComponentList();

  assert(component_list.NFrequencies() > 0);
  assert(correction_factors.size() >= component_list.NFrequencies());

  std::vector<std::vector<double>> averaged_factors;
  if (correction_factors.size() > component_list.NFrequencies()) {
    averaged_factors = AverageCorrectionFactors(correction_factors, component_list.NFrequencies());
  } else {
    averaged_factors = correction_factors;
  }

  // Make a list of the facets
  const size_t n_facets = imaging_table.FacetCount();
  std::vector<const schaapcommon::facets::Facet*> facets;
  facets.reserve(n_facets);
  for (size_t facet_index = 0; facet_index != n_facets; ++facet_index) {
    const ImagingTable table = imaging_table.GetFacet(facet_index);
    facets.emplace_back(table.Front().facet.get());
  }

  // Correct all sources for their facet correction
  const size_t n_channels = averaged_factors.size();
  for (size_t scale_index = 0; scale_index != component_list.NScales(); ++scale_index) {
    const aocommon::UVector<radler::ComponentList::Position>& positions =
        component_list.GetPositions(scale_index);
    size_t previous_facet = 0;
    for (size_t position_index = 0; position_index != positions.size(); ++position_index) {
      const schaapcommon::facets::PixelPosition pos(positions[position_index].x,
                                                    positions[position_index].y);
      const size_t facet_index = schaapcommon::facets::FindFacet(facets, pos, previous_facet);
      if (facet_index < facets.size()) {
        previous_facet = facet_index;

        for (size_t channel_index = 0; channel_index != n_channels; ++channel_index) {
          const double correction_factor = std::sqrt(averaged_factors[channel_index][facet_index]);
          component_list.MultiplyScaleComponent(scale_index, position_index, channel_index,
                                                correction_factor);
        }
      } else {
        for (size_t channel_index = 0; channel_index != n_channels; ++channel_index) {
          // A source could be outside a facet if it is e.g. in the feathering boundary. In this
          // case, it is best to remove the source as we have no way of correcting it to a
          // physical flux density.
          const double correction_factor = 0.0;
          component_list.MultiplyScaleComponent(scale_index, position_index, channel_index,
                                                correction_factor);
        }
      }
    }
  }

  component_list.WriteSources(deconvolution, filename, settings_.pixelScaleX, settings_.pixelScaleY,
                              phase_centre_ra, phase_centre_dec, l_shift, m_shift);
}

void ComponentListWriter::CorrectChannelForPrimaryBeam(radler::ComponentList& list,
                                                       const radler::WorkTableEntry& entry) const {
  Logger::Debug << "Correcting source list of channel " << entry.original_channel_index
                << " for beam\n";
  ImageFilename filename(entry.original_channel_index, entry.original_interval_index);
  filename.SetPolarization(entry.polarization);
  PrimaryBeam beam(settings_);
  PrimaryBeamImageSet beam_images = beam.LoadStokesI(filename);
  beam_images.CorrectComponentList(list, entry.original_channel_index);
}

PrimaryBeamImageSet ComponentListWriter::LoadAveragePrimaryBeam(size_t image_index) const {
  Logger::Debug << "Averaging beam for deconvolution channel " << image_index << "\n";

  PrimaryBeamImageSet beam_images;

  aocommon::Image scratch(settings_.trimmedImageWidth, settings_.trimmedImageHeight);
  size_t deconvolution_channels = settings_.deconvolutionChannelCount;

  /// TODO : use real weights of images
  size_t count = 0;
  PrimaryBeam beam(settings_);
  const std::vector<radler::WorkTable::Group>& channel_groups =
      deconvolution_table_->OriginalGroups();
  for (size_t channel_index = 0; channel_index != channel_groups.size(); ++channel_index) {
    size_t current_image_index = (channel_index * deconvolution_channels) / channel_groups.size();
    if (current_image_index == image_index) {
      const radler::WorkTableEntry& e = *channel_groups[channel_index].front();
      Logger::Debug << "Adding beam at " << e.CentralFrequency() * 1e-6 << " MHz\n";
      ImageFilename filename(e.original_channel_index, e.original_interval_index);

      if (count == 0) {
        beam_images = beam.LoadStokesI(filename);
      } else {
        beam_images += beam.LoadStokesI(filename);
      }
      count++;
    }
  }
  beam_images *= 1.0 / count;
  return beam_images;
}

}  // namespace wsclean
