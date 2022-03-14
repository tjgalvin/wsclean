#ifndef MS_ROW_PROVIDER_H
#define MS_ROW_PROVIDER_H

#include "msrowproviderbase.h"

#include <casacore/casa/Arrays/Array.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/ArrayColumn.h>

#include <cstring>
#include <map>
#include <memory>
#include <vector>

/**
 * An MSRowProvider provides the selected rows of a data set.
 */
class MSRowProvider : public MsRowProviderBase {
 public:
  MSRowProvider(const string& msPath, const MSSelection& selection,
                const std::map<size_t, size_t>& selectedDataDescIds,
                const std::string& dataColumnName, bool requireModel);

  explicit MSRowProvider(
      const casacore::MeasurementSet& ms, const MSSelection& selection,
      const std::map<size_t, size_t>& selected_data_description_ids,
      const std::string& data_column_name, bool require_model);

  bool AtEnd() const override { return _currentRow == EndRow(); }

  void NextRow() override;

  void OutputStatistics() const override {}

  size_t StartTimestep() const { return _startTimestep; }
  size_t EndTimestep() const { return _endTimestep; }

  size_t CurrentProgress() const final { return _currentRow - BeginRow(); }

 protected:
  std::unique_ptr<casacore::ArrayColumn<float>> _weightSpectrumColumn;
  std::unique_ptr<casacore::ArrayColumn<float>> _weightScalarColumn;
  std::unique_ptr<casacore::ArrayColumn<casacore::Complex>> _modelColumn;

  void getCurrentWeights(WeightArray& weights,
                         const casacore::IPosition& shape);
  const std::map<size_t, size_t>& selectedDataDescIds() const {
    return _selectedDataDescIds;
  }
  bool requireModel() const { return _requireModel; }

  bool isCurrentRowSelected(int fieldId, int a1, int a2) const;

  size_t _currentRow;
  size_t _currentTimestep;
  double _currentTime;
  casacore::Vector<double> _currentUVWArray;
  size_t _currentDataDescId;

 private:
  std::map<size_t, size_t> _selectedDataDescIds;
  bool _msHasWeights;
  bool _requireModel;
  casacore::Array<float> _scratchWeightScalarArray;

  size_t _startTimestep, _endTimestep;

  void Initialize();
};

#endif
