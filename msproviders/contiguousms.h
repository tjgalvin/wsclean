#ifndef CONTIGUOUSMS_H
#define CONTIGUOUSMS_H

#include "msprovider.h"

#include "../structures/msselection.h"

#include <aocommon/multibanddata.h>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>

#include <memory>

class ContiguousMSReader;

class ContiguousMS final : public MSProvider {
  friend class ContiguousMSReader;

 public:
  ContiguousMS(const string& msPath, const std::string& dataColumnName,
               const MSSelection& selection, aocommon::PolarizationEnum polOut,
               size_t dataDescIndex, bool useMPI);
  virtual ~ContiguousMS(){};

  ContiguousMS(const ContiguousMS&) = delete;

  ContiguousMS& operator=(const ContiguousMS&) = delete;

  std::unique_ptr<MSReader> MakeReader() override;

  SynchronizedMS MS() override { return _ms; }

  const std::string& DataColumnName() override { return _dataColumnName; }

  void NextOutputRow() override;

  void ResetWritePosition() override;

  void WriteModel(const std::complex<float>* buffer, bool addToMS) override;

  void ReopenRW() override { _ms->reopenRW(); }

  double StartTime() override;

  void MakeIdToMSRowMapping(std::vector<size_t>& idToMSRow) override;

  aocommon::PolarizationEnum Polarization() override { return _polOut; }

  size_t DataDescId() override { return _dataDescId; }

  size_t NChannels() override;

  size_t NPolarizations() override;

  size_t NAntennas() override { return _nAntenna; }

 private:
  void open();

  size_t _currentOutputRow;
  size_t _currentOutputTimestep;
  double _currentOutputTime;
  const int _dataDescId;
  const bool _useMPI;
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
  aocommon::MultiBandData _bandData;
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
  casacore::ArrayColumn<casacore::Complex> _modelColumn;

  casacore::Array<std::complex<float>> _dataArray, _modelArray;
  casacore::Array<float> _weightSpectrumArray, _weightScalarArray,
      _imagingWeightSpectrumArray;
  casacore::Array<bool> _flagArray;

  void prepareModelColumn();
};

#endif
