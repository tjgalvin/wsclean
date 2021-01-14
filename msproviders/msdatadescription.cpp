#include "msdatadescription.h"

#include "contiguousms.h"

#include "../main/settings.h"

#include <aocommon/io/serialostream.h>
#include <aocommon/io/serialistream.h>

std::unique_ptr<MSProvider> MSDataDescription::GetProvider() const {
  if (_isPartitioned)
    return std::unique_ptr<MSProvider>(new PartitionedMS(
        _partitionHandle, _partIndex, _polarization, _dataDescId));
  else
    return std::unique_ptr<MSProvider>(new ContiguousMS(
        _filename, _dataColumnName, _selection, _polarization, _dataDescId));
}

void MSDataDescription::Serialize(aocommon::SerialOStream& stream) const {
  stream.Bool(_isPartitioned).UInt16(_polarization).UInt32(_dataDescId);
  _selection.Serialize(stream);
  stream.String(_filename).String(_dataColumnName);
  _partitionHandle.Serialize(stream);
  stream.UInt64(_partIndex);
}

std::unique_ptr<MSDataDescription> MSDataDescription::Unserialize(
    aocommon::SerialIStream& stream) {
  std::unique_ptr<MSDataDescription> mdd(new MSDataDescription());
  stream.Bool(mdd->_isPartitioned);
  mdd->_polarization = (aocommon::PolarizationEnum)stream.UInt16();
  stream.UInt32(mdd->_dataDescId);
  mdd->_selection.Unserialize(stream);
  stream.String(mdd->_filename).String(mdd->_dataColumnName);
  mdd->_partitionHandle.Unserialize(stream);
  stream.UInt64(mdd->_partIndex);
  return mdd;
}
