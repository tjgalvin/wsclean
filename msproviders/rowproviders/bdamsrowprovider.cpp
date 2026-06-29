#include "bdamsrowprovider.h"

#include <casacore/tables/Tables/TableRecord.h>

#include <cassert>

#include "../operations.h"

using schaapcommon::reordering::MSSelection;

namespace wsclean {

static std::optional<casacore::ArrayColumn<casacore::Complex>> GetModel(
    const casacore::MeasurementSet& ms, const std::string& model_column_name, bool require_model) {
  if (require_model) return casacore::ArrayColumn<casacore::Complex>(ms, model_column_name);

  return {};
}

BdaMsRowProvider::BdaMsRowProvider(const casacore::MeasurementSet& ms, const MSSelection& selection,
                                   const aocommon::VectorSet<size_t>& selected_data_description_ids,
                                   const std::string& data_column_name,
                                   const std::string& model_column_name, bool require_model)
    : MsRowProviderBase(ms, selection, data_column_name, model_column_name),
      max_bda_interval_(GetBdaMaxTimeInterval(Ms())),
      selected_data_description_ids_(selected_data_description_ids),
      weight_(Ms()),
      model_(GetModel(Ms(), model_column_name, require_model)) {
  if (Selection().HasInterval() || Selection().EvenOrOddTimesteps() != MSSelection::kAllTimesteps)
    throw std::runtime_error("An interval selection isn't supported for a BDA measurement set.");

  // Make a sorted index of the table
  sort_data_.reserve(ms.nrow());
  for (size_t row_index = 0; row_index != ms.nrow(); ++row_index) {
    sort_data_.emplace_back(Columns().time(row_index), row_index);
  }
  std::sort(sort_data_.begin(), sort_data_.end());

  data_buffer_offset_ = load_position_;
  output_position_ = 0;

  FillBuffer();

  if (!data_buffer_.empty()) {
    // If the current output row is not selected, move to the first
    // selected one
    const SortDataEntry& entry = sort_data_[output_position_];
    const size_t buffer_position = entry.row_index - data_buffer_offset_;
    if (!IsRowSelected(data_buffer_[buffer_position])) NextRow();
  }
}

void BdaMsRowProvider::FillBuffer() {
  // Fill the data_buffer_ such that it spans twice the max bda interval
  double earliest_time = std::numeric_limits<double>::max();
  double latest_time = std::numeric_limits<double>::lowest();
  if (!data_buffer_.empty()) {
    auto [min, max] = std::minmax_element(data_buffer_.begin(), data_buffer_.end());
    earliest_time = min->time;
    latest_time = max->time;
  }
  while (load_position_ != sort_data_.size() &&
         latest_time - earliest_time < 2.0 * max_bda_interval_) {
    RowData& new_row = data_buffer_.emplace_back();
    LoadRow(new_row, load_position_);
    earliest_time = std::min(earliest_time, new_row.time);
    latest_time = std::max(latest_time, new_row.time);
    ++load_position_;
  }
}

void BdaMsRowProvider::MakeAvailable(size_t row_index) {
  const SortDataEntry& entry = sort_data_[row_index];
  assert(entry.row_index >= data_buffer_offset_);
  const size_t index_position = entry.row_index - data_buffer_offset_;
  if (index_position >= data_buffer_.size()) {
    // Shift the buffer; remove rows that are old enough to not be necessary anymore
    size_t new_first = 0;
    while (new_first != data_buffer_.size() &&
           entry.time - data_buffer_[new_first].time > max_bda_interval_)
      ++new_first;
    data_buffer_.erase(data_buffer_.begin(), data_buffer_.begin() + new_first);
    data_buffer_offset_ += new_first;

    FillBuffer();
  }
}

void BdaMsRowProvider::NextRow() {
  assert(output_position_ != sort_data_.size());
  bool selected;
  do {
    ++output_position_;
    if (output_position_ == sort_data_.size()) return;
    MakeAvailable(output_position_);
    const SortDataEntry& entry = sort_data_[output_position_];
    const size_t index_position = entry.row_index - data_buffer_offset_;
    selected = IsRowSelected(data_buffer_[index_position]);
  } while (!selected);
}

BdaMsRowProvider::RowData& BdaMsRowProvider::CurrentRow() {
  assert(output_position_ != sort_data_.size());
  const SortDataEntry& entry = sort_data_[output_position_];
  assert(entry.row_index >= data_buffer_offset_);
  const size_t index_position = entry.row_index - data_buffer_offset_;

  return data_buffer_[index_position];
}

void BdaMsRowProvider::ReadData(DataArray& data, FlagArray& flags, WeightArray& weights, double& u,
                                double& v, double& w, uint32_t& data_description_id,
                                uint32_t& antenna_1, uint32_t& antenna_2, uint32_t& field_id,
                                double& time) {
  RowData& row_data = CurrentRow();
  u = row_data.uvw[0];
  v = row_data.uvw[1];
  w = row_data.uvw[2];
  data_description_id = row_data.data_description_id;
  antenna_1 = row_data.antenna_1;
  antenna_2 = row_data.antenna_2;
  field_id = row_data.field_id;
  time = row_data.time;
  casacore::swap(data, row_data.casa.data);
  casacore::swap(flags, row_data.casa.flags);
  casacore::swap(weights, row_data.casa.weights);
}

void BdaMsRowProvider::ReadModel(DataArray& model) {
  RowData& row_data = CurrentRow();
  casacore::swap(model, row_data.casa.model);
}

void BdaMsRowProvider::LoadRow(RowData& data, size_t row_index) {
  MsColumns& columns = Columns();
  data.antenna_1 = columns.antenna_1(row_index);
  data.antenna_2 = columns.antenna_2(row_index);
  data.field_id = columns.field_id(row_index);
  data.data_description_id = columns.data_description_id(row_index);
  columns.uvw.get(row_index, uvw_);
  data.uvw[0] = uvw_(0);
  data.uvw[1] = uvw_(1);
  data.uvw[2] = uvw_(2);
  data.row_index = row_index;
  data.time = columns.time(row_index);
  columns.data.get(row_index, data.casa.data, true);
  columns.flag.get(row_index, data.casa.flags, true);
  weight_.ReadData(data.casa.weights, row_index, data.casa.data.shape());
  if (model_) model_->get(row_index, data.casa.model, true);
}

bool BdaMsRowProvider::IsRowSelected(const RowData& data) const {
  if (!selected_data_description_ids_.Contains(data.data_description_id)) {
    return false;
  }

  return Selection().IsSelected(data.field_id, /*timestep=*/-1, data.antenna_1, data.antenna_2,
                                data.uvw.data());
}

BdaMsRowProvider::WriteInfo BdaMsRowProvider::GetWriteInfo() {
  const RowData& row_data = CurrentRow();
  return WriteInfo{.row_index = row_data.row_index,
                   .current_progress = CurrentProgress(),
                   .data_desc_id = row_data.data_description_id};
}

}  // namespace wsclean
