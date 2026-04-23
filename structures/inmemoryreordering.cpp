#include "inmemoryreordering.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <aocommon/banddata.h>

#include <schaapcommon/reordering/reordering.h>

#include "../interface/inmemoryms.h"
#include "../main/settings.h"
#include "../msproviders/operations.h"

using aocommon::Polarization;
using aocommon::PolarizationEnum;
using aocommon::VectorMap;

using schaapcommon::reordering::ChannelRange;
using schaapcommon::reordering::ExtractData;
using schaapcommon::reordering::ExtractWeights;

namespace wsclean {
namespace {

std::shared_ptr<InMemoryInfo> GetInMemoryInfo(InMemoryMs& data) {
  std::shared_ptr<InMemoryInfo> meta_data = std::make_shared<InMemoryInfo>();
  meta_data->antenna_names = std::move(data.antenna_names);
  meta_data->has_frequency_bda = data.has_frequency_bda;
  meta_data->interval = data.interval;
  meta_data->observation_info.fieldName = data.field_name;
  meta_data->observation_info.observer = data.observer;
  meta_data->observation_info.phaseCentreDec = data.phase_centre_dec;
  meta_data->observation_info.phaseCentreRA = data.phase_centre_ra;
  meta_data->observation_info.telescopeName = data.telescope_name;
  meta_data->start_time = data.start_time;
  return meta_data;
}

}  // namespace

InMemoryHandle ReorderInMemory(InMemoryMs&& data,
                               const std::vector<VectorMap<ChannelRange>>& channels,
                               const std::string& data_column_name,
                               const std::string& model_column_name, bool set_model_to_zero,
                               const Settings& settings) {
  InMemoryHandle result;

  // Make a list of the independent output polarizations. Note that this list
  // may have Instrumental or DiagonalInstrumental values in them that expand
  // to multiple correlations per output part.
  std::set<PolarizationEnum> output_polarizations;
  for (PolarizationEnum p : settings.polarizations) {
    output_polarizations.emplace(settings.GetProviderPolarization(p));
  }

  std::shared_ptr<InMemoryInfo> info = GetInMemoryInfo(data);

  const size_t n_provider_correlations =
      Polarization::GetVisibilityCount(*output_polarizations.begin());
  const std::set<PolarizationEnum> input_polarizations(data.polarizations.begin(),
                                                       data.polarizations.end());
  const size_t max_n_flags = n_provider_correlations * data.bands.MaxBandChannels();
  const std::unique_ptr<bool[]> flag_buffer(std::make_unique<bool[]>(max_n_flags));

  std::shared_ptr<std::vector<MSProvider::MetaData>> row_meta_data =
      std::make_shared<std::vector<MSProvider::MetaData>>();

  const std::vector<aocommon::MultiBandData> bands_per_part =
      MakeSelectedBands(data.bands, channels);

  for (size_t output_channel = 0; output_channel != channels.size(); ++output_channel) {
    for (aocommon::PolarizationEnum out_polarization : output_polarizations) {
      InMemoryPart& part_data = result.GetPart(out_polarization, output_channel);

      part_data.info = info;
      part_data.meta_data = row_meta_data;
      part_data.is_regular = bands_per_part[output_channel].BandCount() == 1;
      part_data.polarization = out_polarization;
      part_data.selected_bands = bands_per_part[output_channel];
      part_data.description = "Output channel " + std::to_string(output_channel) +
                              ", polarization " + Polarization::TypeToFullString(out_polarization);
    }
  }

  for (const InMemoryRow& row : data.rows) {
    MSProvider::MetaData& meta_data = row_meta_data->emplace_back();
    meta_data.antenna1 = row.antenna1;
    meta_data.antenna2 = row.antenna2;
    meta_data.data_desc_id = row.data_desc_id;
    meta_data.field_id = row.field_id;
    meta_data.time = row.time;
    meta_data.u_in_m = row.uvw[0];
    meta_data.v_in_m = row.uvw[1];
    meta_data.w_in_m = row.uvw[2];

    for (size_t output_channel = 0; output_channel != channels.size(); ++output_channel) {
      const VectorMap<ChannelRange>& ranges = channels[output_channel];
      if (ContainsDataDescId(ranges, row.data_desc_id)) {
        for (aocommon::PolarizationEnum out_polarization : output_polarizations) {
          InMemoryPart& part_data = result.GetPart(out_polarization, output_channel);
          InMemoryPartRow& part_row = part_data.rows.emplace_back();

          const size_t part_start_ch = ranges[row.data_desc_id].start;
          const size_t part_end_ch = ranges[row.data_desc_id].end;
          const size_t n_channels = part_end_ch - part_start_ch;
          part_row.data_model_weights = InMemoryPartRow::DataModelWeightsArray(
              n_provider_correlations * n_channels, aocommon::kSkipInitialization);
          if (set_model_to_zero) {
            std::fill(part_row.Model().begin(), part_row.Model().end(), 0.0f);
          }
          ExtractData(part_row.Data().data(), part_start_ch, part_end_ch, input_polarizations,
                      row.data.data(), out_polarization);
          ExtractWeights(part_row.Weights().data(), part_start_ch, part_end_ch, input_polarizations,
                         row.data.data(), row.weights.data(), flag_buffer.get(), out_polarization);
        }
      }
    }
  }

  return result;
}

}  // namespace wsclean
