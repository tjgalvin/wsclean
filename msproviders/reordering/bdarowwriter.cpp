#include "bdarowwriter.h"

#include <cassert>
#include <exception>

#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include "../operations.h"

namespace wsclean {

BdaRowWriter::BdaRowWriter(casacore::MeasurementSet& ms,
                           casacore::ArrayColumn<std::complex<float>>& data_column,
                           casacore::ArrayColumn<std::complex<float>>& model_column)
    : max_bda_interval_(GetBdaMaxTimeInterval(ms)),
      data_column_(data_column),
      model_column_(model_column),
      time_column_(ms, casacore::MS::columnName(casacore::MS::TIME)) {}

BdaRowWriter::~BdaRowWriter() {
  if (std::uncaught_exceptions() == 0) {
    while (!rows_.empty()) {
      FillGapAndWrite(rows_.top());
    }

    FillGap(data_column_.nrow());
  }
}

void BdaRowWriter::FillGap(size_t up_to_row_nr) {
  while (current_row_ < up_to_row_nr) {
    scratch_.resize(data_column_.shape(current_row_));
    scratch_ = 0.0f;
    model_column_.put(current_row_, scratch_);
    ++current_row_;
  }
  assert(up_to_row_nr == current_row_);
}

void BdaRowWriter::FillGapAndWrite(const Row& row) {
  FillGap(row.row_nr);
  model_column_.put(row.row_nr, row.data);
  rows_.pop();
  ++current_row_;
}

void BdaRowWriter::Write(size_t row_nr, const casacore::Array<std::complex<float>>& data) {
  if (row_nr == current_row_) {
    model_column_.put(row_nr, data);
    ++current_row_;
    return;
  }

  const double time = time_column_(row_nr);
  rows_.emplace(row_nr, time, data.copy());

  const Row* top = &rows_.top();
  while (top->time - time > max_bda_interval_) {
    // The first item in the list is so much older than the new value that we can be sure no even
    // older values will arrive
    FillGapAndWrite(*top);

    if (rows_.empty()) break;
    top = &rows_.top();
  }
}

}  // namespace wsclean
