#ifndef MSPROVIDERS_BDAMSROWPROVIDER_H
#define MSPROVIDERS_BDAMSROWPROVIDER_H

#include "msrowproviderbase.h"

#include <array>
#include <cstdint>
#include <optional>
#include <vector>

#include <aocommon/vectorset.h>

#include <casacore/casa/Arrays/Vector.h>

#include "msweightcolumn.h"

namespace wsclean {

/**
 * This row provider sorts the rows in time and maintains the reading order (rows are only accessed
 * in increasing row nr), which may be important for compressing storage managers.
 *
 * It works as follows. The entire TIME column is read at the start, and is sorted in combination
 * with the row number. Then, a buffer is filled that spans at least twice the BDA order uncertainty
 * (rows aren't completely sorted on time when using BDA).  This buffer is used until a row is
 * requested that is not inside the buffer. At that point, half the buffer is removed from memory
 * (i.e. the maximum BDA time distance), and the buffer is filled with another span of that same
 * size.
 *
 * If we read twice the max bda interval, we know that the first part is guaranteed to be complete,
 * but the second part may not be (as a sample maybe off by the max bda interval). After processing
 * the first part, the second part becomes the first part and another max bda interval is read, such
 * that the first part is complete.
 *
 * If it is ever a problem that all time values are stored in memory, it would be possible to change
 * the algorithm to store only chunks of time values that correspond with the buffer. This makes the
 * algorithm more complex, so has not yet been done.
 *
 * An initial implementation used a priority queue, but it was found that the priority queue would
 * hold thousands of row indicides at one time, which was not very fast (probably because of the
 * random memory access it causes for each inserted/removed time-row index pair).
 */
class BdaMsRowProvider final : public MsRowProviderBase {
 public:
  /**
   * @pre MsHasBdaData(ms)
   * @pre !selection.HasInterval()
   */
  explicit BdaMsRowProvider(const casacore::MeasurementSet& ms,
                            const schaapcommon::reordering::MSSelection& selection,
                            const aocommon::VectorSet<size_t>& selected_data_description_ids,
                            const std::string& data_column_name,
                            const std::string& model_column_name, bool require_model);

  bool AtEnd() const override { return output_position_ == sort_data_.size(); }

  void NextRow() override;

  void OutputStatistics() const override {}

  void ReadData(DataArray& data, FlagArray& flags, WeightArray& weights, double& u, double& v,
                double& w, uint32_t& data_description_id, uint32_t& antenna_1, uint32_t& antenna_2,
                uint32_t& field_id, double& time) override;

  void ReadModel(DataArray& model) override;

  size_t CurrentProgress() const override { return output_position_ - BeginRow(); }

  WriteInfo GetWriteInfo() final override;

 private:
  struct CasaCoreData {
    DataArray data;
    FlagArray flags;
    WeightArray weights;
    DataArray model;

    CasaCoreData() noexcept = default;
    CasaCoreData(CasaCoreData& source) = delete;
    CasaCoreData& operator=(CasaCoreData& rhs) = delete;

    CasaCoreData(CasaCoreData&& source) noexcept {
      casacore::swap(data, source.data);
      casacore::swap(flags, source.flags);
      casacore::swap(weights, source.weights);
      casacore::swap(model, source.model);
    }

    CasaCoreData& operator=(CasaCoreData&& rhs) noexcept {
      casacore::swap(data, rhs.data);
      casacore::swap(flags, rhs.flags);
      casacore::swap(weights, rhs.weights);
      casacore::swap(model, rhs.model);
      return *this;
    }
  };

  struct RowData {
    size_t row_index;
    double time;
    int antenna_1;
    int antenna_2;
    int field_id;
    uint32_t data_description_id;
    std::array<double, 3> uvw;
    CasaCoreData casa;

    friend bool operator<(const RowData& lhs, const RowData& rhs) noexcept {
      return lhs.time < rhs.time;
    }
  };

  struct SortDataEntry {
    double time;
    size_t row_index;

    friend bool operator<(const SortDataEntry& lhs, const SortDataEntry& rhs) noexcept {
      return lhs.time < rhs.time;
    }
  };

  /**
   * Load the data for the given row index from disk.
   */
  void LoadRow(RowData& data, size_t row_index);

  bool IsRowSelected(const RowData& data) const;

  /**
   * Makes sure the data_buffer_ is filled such that it spans two
   * times the max BDA interval. It will use and update load_position_.
   */
  void FillBuffer();

  /**
   * Will make sure that the buffer contains the specified row. If the row
   * is already part of the output buffer, it will do nothing. Otherwise,
   * it will first remove all old items from the buffer, and insert new
   * rows to the buffer such that the buffer spans twice the max BDA interval.
   */
  void MakeAvailable(size_t row_index);

  /**
   * Points to the next row to be read from the measurement set. It is an
   * absolute index.
   */
  size_t load_position_ = 0;

  /**
   * Index into sort_data_, pointing to the next entry that is
   * to be returned in a call to ReadData()/ReadModel(). It always
   * points to a selected row, unless it is at the end.
   *
   * It does not correspond to a row number; it needs to be added
   * to data_buffer_offset_ to make it a row number.
   */
  size_t output_position_ = 0;

  /**
   * Offset that needs to be added to the data_buffer to map it to row numbers.
   * data_buffer_[0] contains information about the row with index data_buffer_offset_.
   */
  size_t data_buffer_offset_ = 0;

  /**
   * Stores a range of the full data. It is ordered by row number.
   */
  std::vector<RowData> data_buffer_;

  /**
   * Stores the sorted time values and maps these to row indices.
   */
  std::vector<SortDataEntry> sort_data_;

  /**
   * Returns the row data for the row pointed to by load_position_.
   */
  RowData& CurrentRow();

  double max_bda_interval_ = 0.0;
  aocommon::VectorSet<size_t> selected_data_description_ids_;
  MsWeightColumn weight_;
  std::optional<casacore::ArrayColumn<casacore::Complex>> model_;
  casacore::Vector<double> uvw_;
};

}  // namespace wsclean

#endif  // MSPROVIDERS_BDAMSROWPROVIDER_H
