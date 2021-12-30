#include "bdamsrowprovider.h"

#include "msprovider.h"

#include <aocommon/throwruntimeerror.h>

#include <casacore/tables/Tables/TableRecord.h>

#include <cassert>

const char* BdaMsRowProvider::kBDAFactorsTable = "BDA_FACTORS";

static boost::optional<casacore::ArrayColumn<casacore::Complex>> GetModel(
    const casacore::MeasurementSet& ms, bool require_model) {
  if (require_model)
    return casacore::ArrayColumn<casacore::Complex>(
        ms, casacore::MS::columnName(casacore::MSMainEnums::MODEL_DATA));

  return {};
}

BdaMsRowProvider::BdaMsRowProvider(
    const casacore::MeasurementSet& ms, const MSSelection& selection,
    const std::map<size_t, size_t>& selected_data_description_ids,
    const std::string& data_column_name, bool require_model)
    : MsRowProviderBase(ms, selection, data_column_name),
      selected_data_description_ids_(selected_data_description_ids),
      weight_(Ms()),
      model_(GetModel(Ms(), require_model)),
      current_row_(BeginRow()),
      data_(Columns().time(BeginRow())) {
  if (!MsHasBdaData(Ms()))
    aocommon::ThrowRuntimeError("A BDA measurement set requires a ",
                                kBDAFactorsTable, " table.");

  if (Selection().HasInterval() ||
      Selection().EvenOrOddTimesteps() != MSSelection::AllTimesteps)
    aocommon::ThrowRuntimeError(
        "An interval selection isn't supported for a BDA measurement set.");

  // Initializes the uninitalized field of data_.
  if (current_row_ != EndRow() && !LoadCurrentRow()) {
    NextRow();
  }

  if (current_row_ == EndRow())
    aocommon::ThrowRuntimeError(
        "The measurement set contains no data for the current selection.");
}

void BdaMsRowProvider::NextRow() {
  do {
    ++current_row_;
  } while (current_row_ != EndRow() && !LoadCurrentRow());
}

void BdaMsRowProvider::ReadData(DataArray& data, FlagArray& flags,
                                WeightArray& weights, double& u, double& v,
                                double& w, uint32_t& data_description_id,
                                uint32_t& antenna_1, uint32_t& antenna_2,
                                uint32_t& field_id, double& time) {
  assert(data_.uvw.size() == 3 &&
         "Invalid UVW dimensions used or an empty data set.");

  u = data_.uvw(0);
  v = data_.uvw(1);
  w = data_.uvw(2);
  data_description_id = data_.data_description_id;
  MsColumns& columns = Columns();
  columns.data.get(current_row_, data, true);
  columns.flag.get(current_row_, flags, true);
  antenna_1 = columns.antenna_1(current_row_);
  antenna_2 = columns.antenna_2(current_row_);
  field_id = columns.field_id(current_row_);
  time = data_.time;
  weight_.ReadData(weights, current_row_, data.shape());
}

void BdaMsRowProvider::ReadModel(DataArray& model) {
  model_->get(current_row_, model, true);
}

template <class Container, class T>
static bool Contains(const Container& container, const T& value) {
  return container.find(value) != container.end();
}

bool BdaMsRowProvider::LoadCurrentRow() {
  MsColumns& columns = Columns();
  data_.antenna_1 = columns.antenna_1(current_row_);
  data_.antenna_2 = columns.antenna_2(current_row_);
  data_.field_id = columns.field_id(current_row_);
  data_.data_description_id = columns.data_description_id(current_row_);
  data_.time = columns.time(current_row_);
  data_.uvw = columns.uvw(current_row_);

  return IsCurrentRowSelected();
}

bool BdaMsRowProvider::IsCurrentRowSelected() const {
  if (!Contains(selected_data_description_ids_, data_.data_description_id)) {
    return false;
  }

  return Selection().IsSelected(data_.field_id, /*timestep=*/-1,
                                data_.antenna_1, data_.antenna_2, data_.uvw);
}
