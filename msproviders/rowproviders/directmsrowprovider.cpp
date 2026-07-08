#include "directmsrowprovider.h"

namespace wsclean {

void DirectMSRowProvider::ReadData(MSRowProvider::DataArray& data, MSRowProvider::FlagArray& flags,
                                   WeightArray& weights, double& u, double& v, double& w,
                                   uint32_t& dataDescId, uint32_t& antenna1, uint32_t& antenna2,
                                   uint32_t& fieldId, double& time, float& interval) {
  u = _currentUVWArray(0);
  v = _currentUVWArray(1);
  w = _currentUVWArray(2);
  dataDescId = _currentDataDescId;
  MsColumns& columns = Columns();
  columns.data.get(_currentRow, data, true);
  columns.flag.get(_currentRow, flags, true);
  antenna1 = columns.antenna_1(_currentRow);
  antenna2 = columns.antenna_2(_currentRow);
  fieldId = columns.field_id(_currentRow);
  time = columns.time(_currentRow);
  interval = columns.interval(_currentRow);

  getCurrentWeights(weights, data.shape());
}

void DirectMSRowProvider::ReadModel(MSRowProvider::DataArray& model) {
  _modelColumn->get(_currentRow, model, true);
}

DirectMSRowProvider::WriteInfo DirectMSRowProvider::GetWriteInfo() {
  return WriteInfo{.row_index = _currentRow,
                   .current_progress = CurrentProgress(),
                   .data_desc_id = static_cast<uint32_t>(_currentDataDescId)};
}

}  // namespace wsclean
