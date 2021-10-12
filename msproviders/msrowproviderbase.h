#ifndef MSPROVIDERS_MSROWPROVIDER_H
#define MSPROVIDERS_MSROWPROVIDER_H

#include "../structures/msselection.h"
#include "mscolumns.h"

#include <casacore/casa/Arrays/Array.h>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <complex>
#include <memory>
#include <map>
#include <string>

/** The abstract base of the row providers. */
class MsRowProviderBase {
 public:
  using DataArray = casacore::Array<std::complex<float>>;
  using WeightArray = casacore::Array<float>;
  using FlagArray = casacore::Array<bool>;

  explicit MsRowProviderBase(const casacore::MeasurementSet& ms,
                             const MSSelection& selection,
                             const std::string& data_column_name);

  virtual ~MsRowProviderBase() = default;

  virtual bool AtEnd() const = 0;

  virtual void NextRow() = 0;

  /**
   * Read one row of data. The data, flags and weights arrays are resized.
   */
  virtual void ReadData(DataArray& data, FlagArray& flags, WeightArray& weights,
                        double& u, double& v, double& w, uint32_t& dataDescId,
                        uint32_t& antenna1, uint32_t& antenna2,
                        uint32_t& fieldId, double& time) = 0;

  /**
   * Read one model row of data. The model array is resized.
   */
  virtual void ReadModel(DataArray& model) = 0;

  virtual void OutputStatistics() const = 0;

  double StartTime() const {
    return columns_.epoch_as_time(begin_row_).getValue().get();
  }
  double EndTime() const {
    return columns_.epoch_as_time(end_row_ - 1).getValue().get();
  }

  size_t BeginRow() const { return begin_row_; }
  size_t EndRow() const { return end_row_; }

  virtual size_t CurrentProgress() const = 0;
  size_t TotalProgress() const { return EndRow() - BeginRow(); }

  casacore::MeasurementSet& Ms() { return ms_; }

  const MSSelection& Selection() const { return selection_; }

  MsColumns& Columns() { return columns_; }

 private:
  casacore::MeasurementSet ms_;
  MSSelection selection_;
  MsColumns columns_;
  /**
   * Index of the beginning of the selected rows.
   *
   * This is an index to the rows in the original ms.
   */
  size_t begin_row_{0};
  /**
   * Index of the end of the selected rows.
   *
   * This is an index to the rows in the original ms.
   */
  size_t end_row_{0};
};

/**
 * Create a row provider for a measurement set.
 *
 * Create a row provider for a measurement set based on real data and not one
 * of the "dummy" providers. If the measurement set contains the BDA tables
 * created by DP3 it creates a BDA row provider, else a direct row provider is
 * created.
 *
 * @todo AST-630 Return BdaMsRowProvider when BDA tables are present.
 */
[[nodiscard]] std::unique_ptr<MsRowProviderBase> MakeMsRowProvider(
    const std::string& ms_name, const MSSelection& selection,
    const std::map<size_t, size_t>& selected_data_description_ids,
    const std::string& data_column_name, bool require_model);

#endif
