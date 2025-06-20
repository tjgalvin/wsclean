#include "imageweightinitializer.h"

#include <aocommon/logger.h>

#include "../io/imagefilename.h"
#include "../msproviders/contiguousms.h"

using aocommon::Logger;
using schaapcommon::reordering::MSSelection;

namespace wsclean {

std::shared_ptr<ImageWeights> ImageWeightInitializer::Initialize(
    const ImagingTableEntry& entry, const std::vector<MsListItem>& ms_list,
    ImageWeightCache& cache) const {
  if (settings_.mfWeighting) {
    return cache.GetMFWeights();
  } else {
    std::shared_ptr<ImageWeights> weights =
        cache.Get(ms_list, entry.outputChannelIndex, entry.outputIntervalIndex);
    if (settings_.isWeightImageSaved) {
      const std::string prefix = ImageFilename::GetPSFPrefix(
          settings_, entry.outputChannelIndex, entry.outputIntervalIndex);
      weights->Save(prefix + "-weights.fits");
    }
    return weights;
  }
}

void ImageWeightInitializer::GridMfReorderedBand(
    size_t data_desc_id, const aocommon::MultiBandData& bands,
    ImageWeights& weights, size_t ms_index,
    const ImagingTableEntry& entry) const {
  MSSelection part_selection(global_selection_);
  const bool has_selection =
      SelectMsChannels(part_selection, bands, data_desc_id, entry);
  if (has_selection) {
    const ImagingTableEntry::MSInfo& entry_ms_info = entry.msData[ms_index];
    const aocommon::PolarizationEnum pol =
        settings_.GetProviderPolarization(entry.polarization);
    ReorderedMsProvider ms_provider(reordered_ms_handles_[ms_index],
                                    entry_ms_info.bands[data_desc_id].partIndex,
                                    pol, data_desc_id);
    aocommon::BandData selected_band(bands[data_desc_id]);
    if (part_selection.HasChannelRange()) {
      selected_band =
          aocommon::BandData(selected_band, part_selection.ChannelRangeStart(),
                             part_selection.ChannelRangeEnd());
    }
    weights.Grid(ms_provider, selected_band);
  }
}

void ImageWeightInitializer::GridMfContiguousBand(size_t filename_index,
                                                  size_t data_desc_id,
                                                  ImageWeights& weights) const {
  const aocommon::PolarizationEnum pol =
      settings_.GetProviderPolarization(*settings_.polarizations.begin());
  ContiguousMS msProvider(settings_.filenames[filename_index],
                          settings_.dataColumnName, settings_.modelColumnName,
                          settings_.modelStorageManager, global_selection_, pol,
                          data_desc_id, settings_.UseMpi());
  aocommon::BandData selected_band = ms_bands_[filename_index][data_desc_id];
  if (global_selection_.HasChannelRange())
    selected_band =
        aocommon::BandData(selected_band, global_selection_.ChannelRangeStart(),
                           global_selection_.ChannelRangeEnd());
  weights.Grid(msProvider, selected_band);
}

void ImageWeightInitializer::InitializeMf(const ImagingTable& imaging_table,
                                          ImageWeightCache& cache) {
  Logger::Info << "Precalculating MF weights for "
               << settings_.weightMode.ToString() << " weighting...\n";
  std::unique_ptr<ImageWeights> weights = cache.MakeEmptyWeights();
  if (settings_.doReorder) {
    for (const ImagingTable::Group& group : imaging_table.SquaredGroups()) {
      const ImagingTableEntry& entry = *group.front();

      for (size_t ms_index = 0; ms_index != settings_.filenames.size();
           ++ms_index) {
        const aocommon::MultiBandData& band_data = ms_bands_[ms_index];

        for (size_t data_desc_id = 0;
             data_desc_id != band_data.HighestDataDescId() + 1;
             ++data_desc_id) {
          if (band_data.HasDataDescId(data_desc_id)) {
            const size_t band_index = band_data.GetBandIndex(data_desc_id);

            if (settings_.IsBandSelected(band_index)) {
              GridMfReorderedBand(data_desc_id, band_data, *weights, ms_index,
                                  entry);
            }
          }
        }
      }
    }
  } else {
    for (size_t i = 0; i != settings_.filenames.size(); ++i) {
      for (size_t d = 0; d != ms_bands_[i].HighestDataDescId() + 1; ++d) {
        if (ms_bands_[i].HasDataDescId(d)) {
          GridMfContiguousBand(i, d, *weights);
        }
      }
    }
  }
  weights->FinishGridding();
  cache.SetMFWeights(std::move(weights));
  if (settings_.isWeightImageSaved)
    cache.GetMFWeights()->Save(settings_.prefixName + "-weights.fits");
}

}  // namespace wsclean
