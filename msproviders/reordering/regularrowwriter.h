#ifndef WSCLEAN_MSPROVIDERS_REORDERING_STANDARD_ROW_WRITER_H_
#define WSCLEAN_MSPROVIDERS_REORDERING_STANDARD_ROW_WRITER_H_

#include <queue>

#include <casacore/casa/Arrays/Array.h>

#include "rowwriter.h"

namespace casacore {
template <typename T>
class ArrayColumn;
class MeasurementSet;
}  // namespace casacore

namespace wsclean {

class RegularRowWriter final : public RowWriter {
 public:
  RegularRowWriter(casacore::ArrayColumn<std::complex<float>>& data_column,
                   casacore::ArrayColumn<std::complex<float>>& model_column)
      : data_column_(data_column), model_column_(model_column) {}

  ~RegularRowWriter() override;

  void Write(size_t row_nr, const casacore::Array<std::complex<float>>& data) override;

 private:
  /**
   * Fill any unselected rows in the measurement set
   */
  void FillGap(size_t up_to_row_nr);

  casacore::ArrayColumn<std::complex<float>>& data_column_;
  casacore::ArrayColumn<std::complex<float>>& model_column_;
  size_t current_row_ = 0;
  // Used as scratch space and declared in the class to allow reuse without reallocation.
  casacore::Array<std::complex<float>> scratch_;
};

}  // namespace wsclean

#endif
