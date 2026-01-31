#include "contiguousms.h"

#include <memory>

#include <aocommon/logger.h>

#include <schaapcommon/reordering/reordering.h>

#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>
#include <casacore/tables/Tables/TableLocker.h>

#include "operations.h"
#include "msreaders/contiguousmsreader.h"

using schaapcommon::reordering::ChannelRange;
using schaapcommon::reordering::ContainsDataDescId;
using schaapcommon::reordering::MSSelection;
using schaapcommon::reordering::StorageManagerType;

namespace wsclean {
namespace {
void OutputOpeningMessage(
    const std::string& ms_path,
    const aocommon::VectorMap<ChannelRange>& channel_ranges) {
  std::vector<size_t> ids;
  for (const ChannelRange& range : channel_ranges) {
    if (!range.Empty()) ids.emplace_back(range.data_desc_id);
  }
  aocommon::Logger::Info << "Opening " << ms_path << ", ";
  if (ids.size() == 1) {
    aocommon::Logger::Info << "spw " << ids.front();
  } else {
    aocommon::Logger::Info << "multiple spws";
  }
  aocommon::Logger::Info << " with contiguous MS reader.\n";
}
}  // namespace

ContiguousMS::ContiguousMS(
    const string& msPath, const std::string& dataColumnName,
    const std::string& modelColumnName, StorageManagerType modelStorageManager,
    const MSSelection& selection,
    const aocommon::VectorMap<ChannelRange>& channel_ranges,
    aocommon::PolarizationEnum outputPolarization, bool useMPI)
    : _useMPI(useMPI),
      selection_(selection),
      channel_ranges_(channel_ranges),
      _outputPolarization(outputPolarization),
      _msPath(msPath),
      _dataColumnName(dataColumnName),
      _modelColumnName(modelColumnName),
      _modelStorageManager(modelStorageManager) {
  OutputOpeningMessage(_msPath, channel_ranges_);
  open();
}

void ContiguousMS::open() {
  _ms = SynchronizedMS(_msPath, _useMPI ? casacore::TableLock::UserNoReadLocking
                                        : casacore::TableLock::DefaultLocking);
  has_frequency_bda_ = wsclean::HasFrequencyBda(*_ms);

  _antenna1Column = casacore::ScalarColumn<int>(
      *_ms, casacore::MS::columnName(casacore::MSMainEnums::ANTENNA1));
  _antenna2Column = casacore::ScalarColumn<int>(
      *_ms, casacore::MS::columnName(casacore::MSMainEnums::ANTENNA2));
  _fieldIdColumn = casacore::ScalarColumn<int>(
      *_ms, casacore::MS::columnName(casacore::MSMainEnums::FIELD_ID));
  _dataDescIdColumn = casacore::ScalarColumn<int>(
      *_ms, casacore::MS::columnName(casacore::MSMainEnums::DATA_DESC_ID));
  _timeColumn = casacore::ScalarColumn<double>(
      *_ms, casacore::MS::columnName(casacore::MSMainEnums::TIME));
  _uvwColumn = casacore::ArrayColumn<double>(
      *_ms, casacore::MS::columnName(casacore::MSMainEnums::UVW));
  _dataColumn = casacore::ArrayColumn<casacore::Complex>(*_ms, _dataColumnName);
  _flagColumn = casacore::ArrayColumn<bool>(
      *_ms, casacore::MS::columnName(casacore::MSMainEnums::FLAG));

  const casacore::IPosition shape(_dataColumn.shape(0));
  _dataArray = casacore::Array<std::complex<float>>(shape);
  _weightSpectrumArray = casacore::Array<float>(shape);
  _imagingWeightSpectrumArray = casacore::Array<float>(shape);
  _flagArray = casacore::Array<bool>(shape);
  original_bands_ = aocommon::MultiBandData(*_ms);
  selected_bands_ = MakeSelectedPartBands(original_bands_, channel_ranges_);
  if (original_bands_.BandCount() > 1) {
    throw std::runtime_error(
        "This set contains multiple spws, and can therefore not be opened "
        "directly due to possible synchronization issues between spws. You can "
        "force reordering of the measurement by adding -reorder to the command "
        "line.");
  }

  for (const ChannelRange& range : channel_ranges_) {
    const size_t data_desc_id = range.data_desc_id;
    input_polarizations_.AlwaysEmplace(data_desc_id,
                                       GetMSPolarizations(data_desc_id, *_ms));
  }

  _nAntenna = _ms->antenna().nrow();

  _msHasWeightSpectrum = OpenWeightSpectrumColumn(*_ms, _weightSpectrumColumn);
  if (!_msHasWeightSpectrum) {
    casacore::IPosition scalarShape(1, shape[0]);
    _weightScalarArray = casacore::Array<float>(scalarShape);
    _weightScalarColumn.reset(new casacore::ArrayColumn<float>(
        *_ms, casacore::MS::columnName(casacore::MSMainEnums::WEIGHT)));
  }

  GetRowRange(*_ms, selection_, _startRow, _endRow);
  ResetWritePosition();
}

std::unique_ptr<MSReader> ContiguousMS::MakeReader() {
  std::unique_ptr<MSReader> reader(new ContiguousMSReader(this));
  return reader;
}

void ContiguousMS::ResetWritePosition() {
  _currentOutputRow = _startRow - 1;
  _currentOutputTime = 0.0;
  // TODO: something similar needed in the ContiguousMSReader class?
  if (selection_.HasInterval()) {
    _currentOutputTimestep = selection_.IntervalStart() - 1;
  } else {
    _currentOutputTimestep = -1;
  }
  NextOutputRow();
}

void ContiguousMS::NextOutputRow() {
  int fieldId, a1, a2, dataDescId;
  casacore::Vector<double> uvw;
  do {
    ++_currentOutputRow;
    if (_currentOutputRow >= _endRow) return;

    fieldId = _fieldIdColumn(_currentOutputRow);
    a1 = _antenna1Column(_currentOutputRow);
    a2 = _antenna2Column(_currentOutputRow);
    uvw = _uvwColumn(_currentOutputRow);
    dataDescId = _dataDescIdColumn(_currentOutputRow);
    if (_currentOutputTime != _timeColumn(_currentOutputRow)) {
      ++_currentOutputTimestep;
      _currentOutputTime = _timeColumn(_currentOutputRow);
    }
  } while (!selection_.IsSelected(fieldId, _currentOutputTimestep, a1, a2,
                                  uvw.data()) ||
           !ContainsDataDescId(channel_ranges_, dataDescId));
}

double ContiguousMS::StartTime() {
  return casacore::MEpoch::ScalarColumn(
             *_ms, casacore::MS::columnName(casacore::MS::TIME))(_startRow)
      .getValue()
      .get();
}

size_t ContiguousMS::NMaxChannels() {
  return selected_bands_.MaxBandChannels();
}

size_t ContiguousMS::NPolarizations() {
  return aocommon::Polarization::GetVisibilityCount(_outputPolarization);
}

void ContiguousMS::prepareModelColumn() {
  InitializeModelColumn(*_ms, _dataColumnName, _modelColumnName,
                        _modelStorageManager);
  _modelColumn =
      casacore::ArrayColumn<casacore::Complex>(*_ms, _modelColumnName);
  const casacore::IPosition shape(_modelColumn.shape(0));
  _modelArray = casacore::Array<std::complex<float>>(shape);
  _isModelColumnPrepared = true;
}

void ContiguousMS::WriteModel(const std::complex<float>* buffer, bool addToMS) {
  std::unique_ptr<casacore::TableLocker> lock;
  if (_useMPI) {
    // When using different MPI processes, automatic casacore locks do not work.
    // -> Use UserNoReadLocking with explicit write locks.
    lock = std::make_unique<casacore::TableLocker>(*_ms,
                                                   casacore::FileLocker::Write);

    // This resync() is required for synchronizing different MPI processes.
    _ms->resync();
  }

  if (!_isModelColumnPrepared) prepareModelColumn();

  const size_t data_desc_id = _dataDescIdColumn(_currentOutputRow);
  const size_t start_channel = channel_ranges_[data_desc_id].start;
  const size_t end_channel = channel_ranges_[data_desc_id].end;

  _modelColumn.get(_currentOutputRow, _modelArray);
  if (addToMS) {
    schaapcommon::reordering::StoreData<true>(
        _modelArray.data(), start_channel, end_channel,
        input_polarizations_[data_desc_id], buffer, _outputPolarization);
  } else {
    schaapcommon::reordering::StoreData<false>(
        _modelArray.data(), start_channel, end_channel,
        input_polarizations_[data_desc_id], buffer, _outputPolarization);
  }
  _modelColumn.put(_currentOutputRow, _modelArray);
}

bool ContiguousMS::IsRegular() const {
  auto IsNotEmpty = [](const ChannelRange& range) { return !range.Empty(); };

  const size_t n_data_desc_ids = channel_ranges_.CountKeysIf(IsNotEmpty);
  return n_data_desc_ids <= 1;
}

double ContiguousMS::Interval() { return GetMsInterval(*_ms); }

ObservationInfo ContiguousMS::GetObservationInfo() {
  return ReadObservationInfo(*_ms, selection_.FieldIds().front());
}

std::vector<std::string> ContiguousMS::GetAntennaNames() {
  return wsclean::GetAntennaNames(_ms->antenna());
}

}  // namespace wsclean
