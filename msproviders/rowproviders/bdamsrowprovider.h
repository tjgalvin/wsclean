#ifndef MSPROVIDERS_BDAMSROWPROVIDER_H
#define MSPROVIDERS_BDAMSROWPROVIDER_H

#include "msrowproviderbase.h"
#include "msweightcolumn.h"

#include <optional>
#include <queue>

namespace wsclean {

/** A MSRowProvider providing the selected rows in a BDA measurement set. */
class BdaMsRowProvider final : public MsRowProviderBase {
 public:
  /** The table name for the BDA factors. */
  static inline const std::string kBdaFactorsTable = "BDA_FACTORS";
  static inline const std::string kBdaTimeAxisTable = "BDA_TIME_AXIS";
  static inline const std::string kBdaMaxTimeIntervalColumn =
      "MAX_TIME_INTERVAL";
  /**
   * @pre MsHasBdaData(ms)
   * @pre !selection.HasInterval()
   */
  explicit BdaMsRowProvider(
      const casacore::MeasurementSet& ms,
      const schaapcommon::reordering::MSSelection& selection,
      const std::map<size_t, size_t>& selected_data_description_ids,
      const std::string& data_column_name, const std::string& model_column_name,
      bool require_model);

  bool AtEnd() const override {
    return current_row_ == EndRow() && queue_.empty();
  }

  void NextRow() override;

  void OutputStatistics() const override {}

  void ReadData(DataArray& data, FlagArray& flags, WeightArray& weights,
                double& u, double& v, double& w, uint32_t& data_description_id,
                uint32_t& antenna_1, uint32_t& antenna_2, uint32_t& field_id,
                double& time) override;

  void ReadModel(DataArray& model) override;

  size_t CurrentProgress() const override {
    return current_row_ - BeginRow() - queue_.size();
  }

 private:
  struct Data {
    explicit Data(double t) : time(t) {}
    int antenna_1;
    int antenna_2;
    int field_id;
    double time;
    uint32_t data_description_id;
    size_t row_index;
    std::array<double, 3> uvw;

    inline friend bool operator>(const Data& lhs, const Data& rhs) {
      return lhs.time > rhs.time;
    }
  };

  std::map<size_t, size_t> selected_data_description_ids_;
  MsWeightColumn weight_;
  std::optional<casacore::ArrayColumn<casacore::Complex>> model_;
  size_t current_row_;
  /**
   * BDA might not be fully time-ordered, because small baselines might be
   * averaged more, causing a later time (centroid) to appear before a longer
   * baselines's time. We could sort the entire measurement set, but this is
   * costly. Instead, we keep a sorted queue with row data, and keep this filled
   * with as many rows such that it can never happen that the top item is not
   * the next to be processed row.
   *
   * This condition is certainly true when we read from the measurement set a
   * row with a time that is 'BDA_TIME_AXIS.MAX_TIME_INTERVAL' later than the
   * top of the queue.
   */
  std::priority_queue<Data, std::vector<Data>, std::greater<Data>> queue_;
  Data last_read_data_;
  double max_bda_interval_ = 0.0;
  casacore::Vector<double> uvw_;

  void MoveToNextSelectedRow();
  double TimeSpanInQueue() const;

  /**
   * Loads the data for the current row.
   *
   * @returns Whether the current row is selected.
   */
  bool LoadCurrentRow();
  bool IsRowSelected(const Data& data) const;
};

}  // namespace wsclean

#endif  // MSPROVIDERS_BDAMSROWPROVIDER_H
