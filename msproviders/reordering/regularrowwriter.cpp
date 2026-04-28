#include "regularrowwriter.h"

#include <exception>

#include <casacore/tables/Tables/ArrayColumn.h>

namespace wsclean {

RegularRowWriter::~RegularRowWriter() {
  if (std::uncaught_exceptions() == 0) {
    // Fill any final unselected rows in the measurement set
    FillGap(data_column_.nrow());
  }
}

void RegularRowWriter::FillGap(size_t up_to_row_nr) {
  while (current_row_ < up_to_row_nr) {
    scratch_.resize(data_column_.shape(current_row_));
    scratch_ = 0.0f;
    model_column_.put(current_row_, scratch_);
    ++current_row_;
  }
  assert(up_to_row_nr == current_row_);
}

void RegularRowWriter::Write(size_t row_nr, const casacore::Array<std::complex<float>>& data) {
  assert(current_row_ <= row_nr);
  FillGap(row_nr);
  model_column_.put(row_nr, data);
  ++current_row_;
}

}  // namespace wsclean
