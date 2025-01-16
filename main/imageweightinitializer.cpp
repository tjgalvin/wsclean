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
        const ImagingTableEntry::MSInfo& entry_ms_info = entry.msData[ms_index];
        const aocommon::MultiBandData& band_data = ms_bands_[ms_index];

        for (size_t data_description_id = 0;
             data_description_id != band_data.DataDescCount();
             ++data_description_id) {
          const size_t band_index = band_data.GetBandIndex(data_description_id);

          if (settings_.IsBandSelected(band_index)) {
            MSSelection part_selection(global_selection_);
            const bool has_selection = SelectMsChannels(
                part_selection, band_data, data_description_id, entry);
            if (has_selection) {
              const aocommon::PolarizationEnum pol =
                  settings_.GetProviderPolarization(entry.polarization);
              ReorderedMsProvider msProvider(
                  reordered_ms_handles_[ms_index],
                  entry_ms_info.bands[data_description_id].partIndex, pol,
                  data_description_id);
              aocommon::BandData selected_band(band_data[data_description_id]);
              if (part_selection.HasChannelRange()) {
                selected_band = aocommon::BandData(
                    selected_band, part_selection.ChannelRangeStart(),
                    part_selection.ChannelRangeEnd());
              }
              weights->Grid(msProvider, selected_band);
            }
          }
        }
      }
    }
  } else {
    for (size_t i = 0; i != settings_.filenames.size(); ++i) {
      for (size_t d = 0; d != ms_bands_[i].DataDescCount(); ++d) {
        const aocommon::PolarizationEnum pol =
            settings_.GetProviderPolarization(*settings_.polarizations.begin());
        ContiguousMS msProvider(
            settings_.filenames[i], settings_.dataColumnName,
            settings_.modelColumnName, settings_.modelStorageManager,
            global_selection_, pol, d, settings_.UseMpi());
        aocommon::BandData selected_band = ms_bands_[i][d];
        if (global_selection_.HasChannelRange())
          selected_band = aocommon::BandData(
              selected_band, global_selection_.ChannelRangeStart(),
              global_selection_.ChannelRangeEnd());
        weights->Grid(msProvider, selected_band);
        Logger::Info << '.';
        Logger::Info.Flush();
      }
    }
  }
  weights->FinishGridding();
  cache.SetMFWeights(std::move(weights));
  if (settings_.isWeightImageSaved)
    cache.GetMFWeights()->Save(settings_.prefixName + "-weights.fits");
}

}  // namespace wsclean
