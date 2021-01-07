#ifndef MS_ROW_PROVIDER_H
#define MS_ROW_PROVIDER_H

#include "../structures/msselection.h"

#include <casacore/casa/Arrays/Array.h>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>

#include <cstring>
#include <map>
#include <memory>
#include <vector>

/**
 * An MSRowProvider provides the selected rows of a data set.
 */
class MSRowProvider {
 public:
  MSRowProvider(const string& msPath, const MSSelection& selection,
                const std::map<size_t, size_t>& selectedDataDescIds,
                const std::string& dataColumnName, bool requireModel);
  virtual ~MSRowProvider() {}

  typedef casacore::Array<std::complex<float>> DataArray;
  typedef casacore::Array<float> WeightArray;
  typedef casacore::Array<bool> FlagArray;

  virtual bool AtEnd() const { return _currentRow == _endRow; }

  virtual void NextRow();

  virtual void ReadData(DataArray& data, FlagArray& flags, WeightArray& weights,
                        double& u, double& v, double& w, uint32_t& dataDescId,
                        uint32_t& antenna1, uint32_t& antenna2,
                        uint32_t& fieldId, double& time) = 0;

  virtual void ReadModel(DataArray& model) = 0;

  virtual void OutputStatistics() const {}

  casacore::MeasurementSet& MS() { return _ms; }
  casacore::IPosition DataShape() const { return _dataColumn.shape(0); }
  double StartTime() const {
    return _timeEpochColumn(_startRow).getValue().get();
  }
  double EndTime() const {
    return _timeEpochColumn(_endRow - 1).getValue().get();
  }
  size_t StartTimestep() const { return _startTimestep; }
  size_t EndTimestep() const { return _endTimestep; }

  size_t CurrentProgress() const { return _currentRow - _startRow; }
  size_t TotalProgress() const { return _endRow - _startRow; }

 protected:
  casacore::MeasurementSet _ms;
  casacore::ScalarColumn<int> _antenna1Column;
  casacore::ScalarColumn<int> _antenna2Column;
  casacore::ScalarColumn<int> _fieldIdColumn;
  casacore::ScalarColumn<double> _timeColumn;
  casacore::MEpoch::ScalarColumn _timeEpochColumn;
  casacore::ArrayColumn<double> _uvwColumn;
  std::unique_ptr<casacore::ArrayColumn<float>> _weightSpectrumColumn;
  std::unique_ptr<casacore::ArrayColumn<float>> _weightScalarColumn;
  casacore::ArrayColumn<casacore::Complex> _dataColumn;
  casacore::ArrayColumn<bool> _flagColumn;
  casacore::ScalarColumn<int> _dataDescIdColumn;
  std::unique_ptr<casacore::ArrayColumn<casacore::Complex>> _modelColumn;

  void getCurrentWeights(WeightArray& weights);
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
  MSSelection _selection;
  bool _msHasWeights;
  bool _requireModel;
  casacore::Array<float> _scratchWeightScalarArray;

  size_t _startRow, _endRow;
  size_t _startTimestep, _endTimestep;
};

#endif
