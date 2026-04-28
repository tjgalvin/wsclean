#ifndef WSCLEAN_MSPROVIDERS_REORDERING_BDA_ROW_WRITER_H_
#define WSCLEAN_MSPROVIDERS_REORDERING_BDA_ROW_WRITER_H_

#include <queue>

#include <casacore/casa/Arrays/Array.h>
#include <casacore/tables/Tables/ScalarColumn.h>

#include "rowwriter.h"

namespace casacore {
template <typename T>
class ArrayColumn;
class MeasurementSet;
}  // namespace casacore

namespace wsclean {

/**
 * This class makes sure that the measurement set is written sequentially. With time BDA, the
 * reordering performs on the fly sorting in time, which makes the row nrs no longer
 * sequentially. We reverse this operation here again to write sequentially by row nr instead
 * of time. Sorting rows sequentially is particularly important for storage managers like Sisco.
 */
class BdaRowWriter final : public RowWriter {
 public:
  BdaRowWriter(casacore::MeasurementSet& ms,
               casacore::ArrayColumn<std::complex<float>>& data_column,
               casacore::ArrayColumn<std::complex<float>>& model_column);
  ~BdaRowWriter() override;

  void Write(size_t row_nr, const casacore::Array<std::complex<float>>& data) override;

 private:
  struct Row {
    Row(size_t _row_nr, double _time, casacore::Array<std::complex<float>> _data)
        : row_nr(_row_nr), time(_time), data(_data) {}
    Row(Row&&) = default;
    Row& operator=(Row&& rhs) {
      // This asignment operator is necessary because casacore's Array class does
      // not properly support r-value semantics.
      row_nr = rhs.row_nr;
      time = rhs.time;
      casacore::swap(data, rhs.data);
      return *this;
    }

    size_t row_nr;
    double time;
    casacore::Array<std::complex<float>> data;

    /**
     * Returns row_nr > rhs.row_nr ; order is reversed such that earliest
     * rows are placed on top of the priority queue.
     */
    bool operator<(const Row& rhs) const { return row_nr > rhs.row_nr; }
  };

  void FillGapAndWrite(const Row& row);
  /**
   * Fill any unselected rows in the measurement set
   */
  void FillGap(size_t up_to_row_nr);

  std::priority_queue<Row> rows_;

  double max_bda_interval_ = 0.0;
  casacore::ArrayColumn<std::complex<float>>& data_column_;
  casacore::ArrayColumn<std::complex<float>>& model_column_;
  casacore::ScalarColumn<double> time_column_;
  size_t current_row_ = 0;
  double current_time_ = 0.0;
  casacore::Array<std::complex<float>> scratch_;
};

}  // namespace wsclean

#endif
