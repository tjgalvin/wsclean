#ifndef CONTIGUOUSMS_H
#define CONTIGUOUSMS_H

#include "msprovider.h"

#include "../structures/msselection.h"
#include "../structures/multibanddata.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>

#include <memory>

class ContiguousMS : public MSProvider {
 public:
  ContiguousMS(const string& msPath, const std::string& dataColumnName,
               const MSSelection& selection, aocommon::PolarizationEnum polOut,
               size_t dataDescIndex);

  ContiguousMS(const ContiguousMS&) = delete;

  ContiguousMS& operator=(const ContiguousMS&) = delete;

  SynchronizedMS MS() final override { return _ms; }

  const std::string& DataColumnName() final override { return _dataColumnName; }

  size_t RowId() const final override { return _currentRowId; }

  bool CurrentRowAvailable() final override;

  void NextInputRow() final override;

  void NextOutputRow() final override;

  void Reset() final override;

  void ReadMeta(double& u, double& v, double& w,
                size_t& dataDescId) final override;

  void ReadMeta(MetaData& metaData) final override;

  void ReadData(std::complex<float>* buffer) final override;

  void ReadModel(std::complex<float>* buffer) final override;

  void WriteModel(const std::complex<float>* buffer,
                  bool addToMS) final override;

  void ReadWeights(std::complex<float>* buffer) final override;

  void ReadWeights(float* buffer) final override;

  void WriteImagingWeights(const float* buffer) final override;

  void ReopenRW() final override { _ms->reopenRW(); }

  double StartTime() final override;

  void MakeIdToMSRowMapping(std::vector<size_t>& idToMSRow) final override;

  aocommon::PolarizationEnum Polarization() final override { return _polOut; }

  size_t NChannels() override;

  size_t NPolarizations() override;

  size_t NAntennas() override { return _nAntenna; }

 private:
  void open();

  size_t _currentInputRow;
  size_t _currentInputTimestep;
  double _currentInputTime;
  size_t _currentOutputRow;
  size_t _currentOutputTimestep;
  double _currentOutputTime;
  size_t _currentRowId;
  const int _dataDescId;
  size_t _nAntenna;
  bool _isDataRead, _isModelRead, _isWeightRead;
  bool _isModelColumnPrepared;
  size_t _startRow, _endRow;
  std::vector<size_t> _idToMSRow;
  std::vector<aocommon::PolarizationEnum> _inputPolarizations;
  MSSelection _selection;
  aocommon::PolarizationEnum _polOut;
  std::string _msPath;
  SynchronizedMS _ms;
  MultiBandData _bandData;
  bool _msHasWeightSpectrum;

  casacore::ScalarColumn<int> _antenna1Column, _antenna2Column, _fieldIdColumn,
      _dataDescIdColumn;
  casacore::ScalarColumn<double> _timeColumn;
  casacore::ArrayColumn<double> _uvwColumn;
  std::unique_ptr<casacore::ArrayColumn<float>> _weightSpectrumColumn;
  std::unique_ptr<casacore::ArrayColumn<float>> _weightScalarColumn;
  std::string _dataColumnName;
  casacore::ArrayColumn<casacore::Complex> _dataColumn;
  casacore::ArrayColumn<bool> _flagColumn;
  std::unique_ptr<casacore::ArrayColumn<casacore::Complex>> _modelColumn;
  std::unique_ptr<casacore::ArrayColumn<float>> _imagingWeightsColumn;

  casacore::Array<std::complex<float>> _dataArray, _modelArray;
  casacore::Array<float> _weightSpectrumArray, _weightScalarArray,
      _imagingWeightSpectrumArray;
  casacore::Array<bool> _flagArray;

  void prepareModelColumn();
  void readData() {
    if (!_isDataRead) {
      _dataColumn.get(_currentInputRow, _dataArray);
      _isDataRead = true;
    }
  }
  void readWeights() {
    if (!_isWeightRead) {
      _flagColumn.get(_currentInputRow, _flagArray);
      if (_msHasWeightSpectrum)
        _weightSpectrumColumn->get(_currentInputRow, _weightSpectrumArray);
      else {
        _weightScalarColumn->get(_currentInputRow, _weightScalarArray);
        expandScalarWeights(_weightScalarArray, _weightSpectrumArray);
      }
      _isWeightRead = true;
    }
  }
  void readModel() {
    if (!_isModelRead) {
      _modelColumn->get(_currentInputRow, _modelArray);
      _isModelRead = true;
    }
  }
};

#endif
