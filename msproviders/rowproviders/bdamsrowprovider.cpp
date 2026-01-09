#include "bdamsrowprovider.h"

#include <casacore/tables/Tables/TableRecord.h>

#include <cassert>

using schaapcommon::reordering::MSSelection;

namespace wsclean {

static std::optional<casacore::ArrayColumn<casacore::Complex>> GetModel(
    const casacore::MeasurementSet& ms, const std::string& model_column_name,
    bool require_model) {
  if (require_model)
    return casacore::ArrayColumn<casacore::Complex>(ms, model_column_name);

  return {};
}

BdaMsRowProvider::BdaMsRowProvider(
    const casacore::MeasurementSet& ms, const MSSelection& selection,
    const std::map<size_t, size_t>& selected_data_description_ids,
    const std::string& data_column_name, const std::string& model_column_name,
    bool require_model)
    : MsRowProviderBase(ms, selection, data_column_name, model_column_name),
      selected_data_description_ids_(selected_data_description_ids),
      weight_(Ms()),
      model_(GetModel(Ms(), model_column_name, require_model)),
      current_row_(BeginRow()),
      last_read_data_(Columns().time(BeginRow())),
      max_bda_interval_(GetBdaMaxTimeInterval(Ms())) {
  if (Selection().HasInterval() ||
      Selection().EvenOrOddTimesteps() != MSSelection::kAllTimesteps)
    throw std::runtime_error(
        "An interval selection isn't supported for a BDA measurement set.");

  // Initializes last_read_data_.
  if (current_row_ != EndRow() && !LoadCurrentRow()) {
    MoveToNextSelectedRow();
  }
  // Fill the queue until we have a full time span in the queue.
  while (current_row_ != EndRow() && TimeSpanInQueue() <= max_bda_interval_) {
    queue_.emplace(last_read_data_);
    MoveToNextSelectedRow();
  }

  if (AtEnd())
    throw std::runtime_error(
        "The measurement set contains no data for the current selection.");
}

inline double BdaMsRowProvider::TimeSpanInQueue() const {
  if (queue_.empty()) return 0.0;
  return last_read_data_.time - queue_.top().time;
}

void BdaMsRowProvider::NextRow() {
  queue_.pop();
  while (current_row_ != EndRow() && TimeSpanInQueue() <= max_bda_interval_) {
    queue_.push(last_read_data_);
    MoveToNextSelectedRow();
  }
}

void BdaMsRowProvider::MoveToNextSelectedRow() {
  do {
    ++current_row_;
  } while (current_row_ != EndRow() && !LoadCurrentRow());
}

void BdaMsRowProvider::ReadData(DataArray& data, FlagArray& flags,
                                WeightArray& weights, double& u, double& v,
                                double& w, uint32_t& data_description_id,
                                uint32_t& antenna_1, uint32_t& antenna_2,
                                uint32_t& field_id, double& time) {
  assert(!queue_.empty());
  const Data& row_data = queue_.top();
  u = row_data.uvw[0];
  v = row_data.uvw[1];
  w = row_data.uvw[2];
  data_description_id = row_data.data_description_id;
  MsColumns& columns = Columns();
  const size_t row_index = row_data.row_index;
  columns.data.get(row_index, data, true);
  columns.flag.get(row_index, flags, true);
  antenna_1 = columns.antenna_1(row_index);
  antenna_2 = columns.antenna_2(row_index);
  field_id = columns.field_id(row_index);
  time = row_data.time;
  weight_.ReadData(weights, row_index, data.shape());
}

void BdaMsRowProvider::ReadModel(DataArray& model) {
  assert(!queue_.empty());
  model_->get(queue_.top().row_index, model, true);
}

bool BdaMsRowProvider::LoadCurrentRow() {
  MsColumns& columns = Columns();
  last_read_data_.antenna_1 = columns.antenna_1(current_row_);
  last_read_data_.antenna_2 = columns.antenna_2(current_row_);
  last_read_data_.field_id = columns.field_id(current_row_);
  last_read_data_.data_description_id =
      columns.data_description_id(current_row_);
  last_read_data_.time = columns.time(current_row_);
  columns.uvw.get(current_row_, uvw_);
  last_read_data_.uvw[0] = uvw_(0);
  last_read_data_.uvw[1] = uvw_(1);
  last_read_data_.uvw[2] = uvw_(2);
  last_read_data_.row_index = current_row_;

  return IsRowSelected(last_read_data_);
}

bool BdaMsRowProvider::IsRowSelected(const Data& data) const {
  if (!selected_data_description_ids_.contains(data.data_description_id)) {
    return false;
  }

  return Selection().IsSelected(data.field_id, /*timestep=*/-1, data.antenna_1,
                                data.antenna_2, data.uvw.data());
}

}  // namespace wsclean
