#ifndef WSCLEAN_MSPROVIDERS_IN_MEMORY_PROVIDER_H_
#define WSCLEAN_MSPROVIDERS_IN_MEMORY_PROVIDER_H_

#include "msprovider.h"

#include "../structures/inmemorypart.h"

namespace wsclean {

class InMemoryProvider final : public MSProvider {
 public:
  InMemoryProvider(InMemoryPart& data)
      : data_(data), output_iterator_(data_.rows.begin()) {}

  std::string PartDescription() const final { return data_.description; }

  void NextOutputRow() final { ++output_iterator_; }

  void ResetWritePosition() final { output_iterator_ = data_.rows.begin(); }

  void WriteModel(const std::complex<float>* buffer, bool add) final {
    std::vector<std::complex<float>>& data = output_iterator_->model_data;
    if (add) {
      for (size_t i = 0; i != data.size(); ++i) {
        data[i] += buffer[i];
      }
    } else {
      std::copy_n(buffer, data.size(), data.begin());
    }
  }

  void ReopenRW() final {}

  double StartTime() { return data_.info->start_time; }

  aocommon::PolarizationEnum Polarization() final { return data_.polarization; }

  size_t NMaxChannels() final { return data_.selected_bands.MaxBandChannels(); }

  bool IsRegular() const final { return data_.is_regular; }

  bool HasFrequencyBda() const final { return data_.info->has_frequency_bda; }

  size_t NAntennas() final { return data_.info->antenna_names.size(); }

  size_t NPolarizations() final {
    return aocommon::Polarization::GetVisibilityCount(data_.polarization);
  }

  size_t NRows() final { return data_.rows.size(); }

  double Interval() final { return data_.info->interval; }

  ObservationInfo GetObservationInfo() final {
    return data_.info->observation_info;
  }

  std::vector<std::string> GetAntennaNames() final {
    return data_.info->antenna_names;
  }

  const aocommon::MultiBandData& SelectedBands() final {
    return data_.selected_bands;
  }

  std::unique_ptr<MSReader> MakeReader() final;

 public:
  InMemoryPart& data_;
  std::vector<InMemoryPartRow>::iterator output_iterator_;
};

}  // namespace wsclean

#endif
