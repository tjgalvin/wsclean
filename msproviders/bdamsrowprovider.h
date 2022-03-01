#ifndef MSPROVIDERS_BDAMSROWPROVIDER_H
#define MSPROVIDERS_BDAMSROWPROVIDER_H

#include "msrowproviderbase.h"
#include "msweightcolumn.h"

#include <optional>

/** A MSRowProvider providing the selected rows in a BDA measurement set. */
class BdaMsRowProvider final : public MsRowProviderBase {
 public:
  /** The table name for the BDA factors. */
  static const char* kBDAFactorsTable;

  /**
   * @pre MsHasBdaData(ms)
   * @pre !selection.HasInterval()
   */
  explicit BdaMsRowProvider(
      const casacore::MeasurementSet& ms, const MSSelection& selection,
      const std::map<size_t, size_t>& selected_data_description_ids,
      const std::string& data_column_name, bool require_model);

  bool AtEnd() const override { return current_row_ == EndRow(); }

  void NextRow() override;

  void OutputStatistics() const override {}

  void ReadData(DataArray& data, FlagArray& flags, WeightArray& weights,
                double& u, double& v, double& w, uint32_t& data_description_id,
                uint32_t& antenna_1, uint32_t& antenna_2, uint32_t& field_id,
                double& time) override;

  void ReadModel(DataArray& model) override;

  size_t CurrentProgress() const override { return current_row_ - BeginRow(); }

 private:
  std::map<size_t, size_t> selected_data_description_ids_;
  MsWeightColumn weight_;
  std::optional<casacore::ArrayColumn<casacore::Complex>> model_;
  size_t current_row_;

  struct Data {
    explicit Data(double t) : time(t) {}
    int antenna_1;
    int antenna_2;
    int field_id;
    double time;
    uint32_t data_description_id;
    casacore::Vector<double> uvw;
  };

  Data data_;

  /**
   * Loads the data for the current row.
   *
   * @returns Whether the current row is selected.
   */
  bool LoadCurrentRow();
  bool IsCurrentRowSelected() const;
};

#endif  // MSPROVIDERS_BDAMSROWPROVIDER_H
