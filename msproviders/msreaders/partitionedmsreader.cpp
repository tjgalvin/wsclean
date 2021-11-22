#include "partitionedmsreader.h"
#include "../partitionedms.h"

PartitionedMSReader::PartitionedMSReader(PartitionedMS* partitionedMS)
    : MSReader(partitionedMS),
      _currentInputRow(0),
      _readPtrIsOk(true),
      _metaPtrIsOk(true),
      _weightPtrIsOk(true) {
  _metaFile.open(PartitionedMS::getMetaFilename(
                     partitionedMS->_handle._data->_msPath,
                     partitionedMS->_handle._data->_temporaryDirectory,
                     partitionedMS->_partHeader.dataDescId),
                 std::ios::in);
  std::vector<char> msPath(partitionedMS->_metaHeader.filenameLength + 1,
                           char(0));
  // meta and data header were read in PartitionedMS constructor
  _metaFile.seekg(sizeof(PartitionedMS::MetaHeader), std::ios::beg);
  _metaFile.read(msPath.data(), partitionedMS->_metaHeader.filenameLength);
  std::string partPrefix = PartitionedMS::getPartPrefix(
      msPath.data(), partitionedMS->_partIndex, partitionedMS->_polarization,
      partitionedMS->_partHeader.dataDescId,
      partitionedMS->_handle._data->_temporaryDirectory);
  _dataFile.open(partPrefix + ".tmp", std::ios::in);
  if (!_dataFile.good())
    throw std::runtime_error("Error opening temporary data file in '" +
                             partPrefix + ".tmp'");
  _dataFile.seekg(sizeof(PartitionedMS::PartHeader), std::ios::beg);

  _weightFile.open(partPrefix + "-w.tmp", std::ios::in);
  if (!_weightFile.good())
    throw std::runtime_error("Error opening temporary data weight file '" +
                             partPrefix + "-w.tmp'");
}

bool PartitionedMSReader::CurrentRowAvailable() {
  const PartitionedMS& partitionedms =
      static_cast<const PartitionedMS&>(*_msProvider);
  return _currentInputRow < partitionedms._metaHeader.selectedRowCount;
}

void PartitionedMSReader::NextInputRow() {
  const PartitionedMS& partitionedms =
      static_cast<const PartitionedMS&>(*_msProvider);

  ++_currentInputRow;
  if (_currentInputRow < partitionedms._metaHeader.selectedRowCount) {
    if (_readPtrIsOk)
      _dataFile.seekg(partitionedms._partHeader.channelCount *
                          partitionedms._polarizationCountInFile *
                          sizeof(std::complex<float>),
                      std::ios::cur);
    else
      _readPtrIsOk = true;

    if (_metaPtrIsOk)
      _metaFile.seekg(PartitionedMS::MetaRecord::BINARY_SIZE, std::ios::cur);
    else
      _metaPtrIsOk = true;

    if (_weightPtrIsOk)
      _weightFile.seekg(partitionedms._partHeader.channelCount *
                            partitionedms._polarizationCountInFile *
                            sizeof(float),
                        std::ios::cur);
    _weightPtrIsOk = true;
  }
}

void PartitionedMSReader::ReadMeta(double& u, double& v, double& w) {
  if (!_metaPtrIsOk)
    _metaFile.seekg(-PartitionedMS::MetaRecord::BINARY_SIZE, std::ios::cur);
  _metaPtrIsOk = false;

  PartitionedMS::MetaRecord record;
  record.read(_metaFile);
  u = record.u;
  v = record.v;
  w = record.w;
}

void PartitionedMSReader::ReadMeta(MSProvider::MetaData& metaData) {
  if (!_metaPtrIsOk)
    _metaFile.seekg(-PartitionedMS::MetaRecord::BINARY_SIZE, std::ios::cur);
  _metaPtrIsOk = false;

  PartitionedMS::MetaRecord record;
  record.read(_metaFile);
  metaData.uInM = record.u;
  metaData.vInM = record.v;
  metaData.wInM = record.w;
  metaData.fieldId = record.fieldId;
  metaData.antenna1 = record.antenna1;
  metaData.antenna2 = record.antenna2;
  metaData.time = record.time;
}

void PartitionedMSReader::ReadData(std::complex<float>* buffer) {
  const PartitionedMS& partitionedms =
      static_cast<const PartitionedMS&>(*_msProvider);

  if (!_readPtrIsOk) {
    _dataFile.seekg(-partitionedms._partHeader.channelCount *
                        partitionedms._polarizationCountInFile *
                        sizeof(std::complex<float>),
                    std::ios::cur);
  }
#ifndef NDEBUG
  size_t pos = size_t(_dataFile.tellg()) - sizeof(PartitionedMS::PartHeader);
  if (pos != _currentInputRow * partitionedms._partHeader.channelCount *
                 partitionedms._polarizationCountInFile *
                 sizeof(std::complex<float>)) {
    std::ostringstream s;
    s << "Not on right pos: " << pos << " instead of "
      << _currentInputRow * partitionedms._partHeader.channelCount *
             partitionedms._polarizationCountInFile *
             sizeof(std::complex<float>)
      << " (row "
      << (pos / (partitionedms._partHeader.channelCount *
                 partitionedms._polarizationCountInFile *
                 sizeof(std::complex<float>)))
      << " instead of " << _currentInputRow << ")";
    throw std::runtime_error(s.str());
  }
#endif
  _dataFile.read(reinterpret_cast<char*>(buffer),
                 partitionedms._partHeader.channelCount *
                     partitionedms._polarizationCountInFile *
                     sizeof(std::complex<float>));
  _readPtrIsOk = false;
}

void PartitionedMSReader::ReadModel(std::complex<float>* buffer) {
  const PartitionedMS& partitionedms =
      static_cast<const PartitionedMS&>(*_msProvider);

#ifndef NDEBUG
  if (!partitionedms._partHeader.hasModel)
    throw std::runtime_error("Partitioned MS initialized without model");
#endif
  size_t rowLength = partitionedms._partHeader.channelCount *
                     partitionedms._polarizationCountInFile *
                     sizeof(std::complex<float>);
  memcpy(reinterpret_cast<char*>(buffer),
         partitionedms._modelFileMap + rowLength * _currentInputRow, rowLength);
}

void PartitionedMSReader::ReadWeights(std::complex<float>* buffer) {
  const PartitionedMS& partitionedms =
      static_cast<const PartitionedMS&>(*_msProvider);

  if (!_weightPtrIsOk)
    _weightFile.seekg(-partitionedms._partHeader.channelCount * sizeof(float),
                      std::ios::cur);
  float* displacedBuffer = reinterpret_cast<float*>(buffer) +
                           partitionedms._partHeader.channelCount *
                               partitionedms._polarizationCountInFile;
  _weightFile.read(reinterpret_cast<char*>(displacedBuffer),
                   partitionedms._partHeader.channelCount *
                       partitionedms._polarizationCountInFile * sizeof(float));
  _weightPtrIsOk = false;
  MSProvider::CopyRealToComplex(buffer, displacedBuffer,
                                partitionedms._partHeader.channelCount *
                                    partitionedms._polarizationCountInFile);
}

void PartitionedMSReader::ReadWeights(float* buffer) {
  const PartitionedMS& partitionedms =
      static_cast<const PartitionedMS&>(*_msProvider);

  if (!_weightPtrIsOk)
    _weightFile.seekg(-partitionedms._partHeader.channelCount *
                          partitionedms._polarizationCountInFile *
                          sizeof(float),
                      std::ios::cur);
  _weightFile.read(reinterpret_cast<char*>(buffer),
                   partitionedms._partHeader.channelCount *
                       partitionedms._polarizationCountInFile * sizeof(float));
  _weightPtrIsOk = false;
}

void PartitionedMSReader::WriteImagingWeights(const float* buffer) {
  const PartitionedMS& partitionedms =
      static_cast<const PartitionedMS&>(*_msProvider);

  if (_imagingWeightsFile == nullptr) {
    std::string partPrefix = PartitionedMS::getPartPrefix(
        partitionedms._handle._data->_msPath, partitionedms._partIndex,
        partitionedms._polarization, partitionedms._partHeader.dataDescId,
        partitionedms._handle._data->_temporaryDirectory);
    _imagingWeightsFile.reset(
        new std::fstream(partPrefix + "-imgw.tmp",
                         std::ios::in | std::ios::out | std::ios::binary));
  }
  const size_t nVis = partitionedms._partHeader.channelCount *
                      partitionedms._polarizationCountInFile;
  _imagingWeightBuffer.resize(nVis);
  const size_t chunkSize = nVis * sizeof(float);
  _imagingWeightsFile->seekg(chunkSize * _currentInputRow, std::ios::beg);
  _imagingWeightsFile->read(
      reinterpret_cast<char*>(_imagingWeightBuffer.data()),
      partitionedms._partHeader.channelCount *
          partitionedms._polarizationCountInFile * sizeof(float));
  for (size_t i = 0; i != partitionedms._partHeader.channelCount *
                              partitionedms._polarizationCountInFile;
       ++i) {
    if (std::isfinite(buffer[i])) _imagingWeightBuffer[i] = buffer[i];
  }
  _imagingWeightsFile->seekp(chunkSize * _currentInputRow, std::ios::beg);
  _imagingWeightsFile->write(
      reinterpret_cast<const char*>(_imagingWeightBuffer.data()),
      partitionedms._partHeader.channelCount *
          partitionedms._polarizationCountInFile * sizeof(float));
}
