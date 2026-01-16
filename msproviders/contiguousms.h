#ifndef CONTIGUOUSMS_H
#define CONTIGUOUSMS_H

#include "msprovider.h"

#include "../structures/msselection.h"

#include <aocommon/multibanddata.h>
#include <aocommon/vectormap.h>

#include <schaapcommon/reordering/channelrange.h>
#include <schaapcommon/reordering/msselection.h>
#include <schaapcommon/reordering/storagemanagertype.h>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>

#include <memory>

namespace wsclean {

class ContiguousMSReader;

class ContiguousMS final : public MSProvider {
  friend class ContiguousMSReader;

 public:
  ContiguousMS(
      const string& msPath, const std::string& dataColumnName,
      const std::string& modelColumnName,
      schaapcommon::reordering::StorageManagerType modelStorageManager,
      const schaapcommon::reordering::MSSelection& selection,
      const aocommon::VectorMap<schaapcommon::reordering::ChannelRange>&
          channel_selection,
      aocommon::PolarizationEnum polOut, bool useMPI);

  ContiguousMS(const ContiguousMS&) = delete;

  ContiguousMS& operator=(const ContiguousMS&) = delete;

  std::unique_ptr<MSReader> MakeReader() override;

  void NextOutputRow() override;

  void ResetWritePosition() override;

  void WriteModel(const std::complex<float>* buffer, bool addToMS) override;

  void ReopenRW() override { _ms->reopenRW(); }

  double StartTime() override;

  aocommon::PolarizationEnum Polarization() override {
    return _outputPolarization;
  }

  size_t NMaxChannels() override;

  bool IsRegular() const override;

  bool HasFrequencyBda() const override { return has_frequency_bda_; }

  std::string PartDescription() const final { return _msPath; }

  size_t NPolarizations() override;

  size_t NAntennas() override { return _nAntenna; }

  size_t NRows() final { return _endRow - _startRow; }

  const aocommon::MultiBandData& SelectedBands() final {
    return selected_bands_;
  }

  double Interval() final;

  ObservationInfo GetObservationInfo() final;

  std::vector<std::string> GetAntennaNames() final;

  std::optional<SynchronizedMS> MsIfAvailable() final { return _ms; }

 private:
  void open();

  size_t _currentOutputRow = 0;
  size_t _currentOutputTimestep = 0;
  double _currentOutputTime = 0.0;
  bool _useMPI;
  size_t _nAntenna = 0;
  bool _isDataRead, _isModelRead, _isWeightRead;
  bool _isModelColumnPrepared = false;
  size_t _startRow, _endRow;
  aocommon::VectorMap<std::set<aocommon::PolarizationEnum>>
      input_polarizations_;
  schaapcommon::reordering::MSSelection selection_;
  aocommon::VectorMap<schaapcommon::reordering::ChannelRange> channel_ranges_;
  bool has_frequency_bda_ = false;
  aocommon::PolarizationEnum _outputPolarization;
  std::string _msPath;
  SynchronizedMS _ms;
  aocommon::MultiBandData original_bands_;
  aocommon::MultiBandData selected_bands_;
  bool _msHasWeightSpectrum;

  casacore::ScalarColumn<int> _antenna1Column, _antenna2Column, _fieldIdColumn,
      _dataDescIdColumn;
  casacore::ScalarColumn<double> _timeColumn;
  casacore::ArrayColumn<double> _uvwColumn;
  std::unique_ptr<casacore::ArrayColumn<float>> _weightSpectrumColumn;
  std::unique_ptr<casacore::ArrayColumn<float>> _weightScalarColumn;
  std::string _dataColumnName;
  std::string _modelColumnName;
  casacore::ArrayColumn<casacore::Complex> _dataColumn;
  casacore::ArrayColumn<bool> _flagColumn;
  casacore::ArrayColumn<casacore::Complex> _modelColumn;
  schaapcommon::reordering::StorageManagerType _modelStorageManager;

  casacore::Array<std::complex<float>> _dataArray, _modelArray;
  casacore::Array<float> _weightSpectrumArray, _weightScalarArray,
      _imagingWeightSpectrumArray;
  casacore::Array<bool> _flagArray;

  void prepareModelColumn();
};

}  // namespace wsclean

#endif
