#include "msdatadescription.h"

#include "contiguousms.h"

#include "../main/settings.h"

#include <aocommon/io/serialostream.h>
#include <aocommon/io/serialistream.h>

#include <boost/make_unique.hpp>

#include <cassert>

std::unique_ptr<MSProvider> MSDataDescription::GetProvider() const {
  if (_isPartitioned)
    return boost::make_unique<PartitionedMS>(_partitionHandle, _partIndex,
                                             _polarization, _dataDescId);
  else
    return boost::make_unique<ContiguousMS>(_filename, _dataColumnName,
                                            _selection, _polarization,
                                            _dataDescId, _useMPI);
}

void MSDataDescription::Serialize(aocommon::SerialOStream& stream) const {
  // Serialization is only used with MPI.
  assert(_useMPI);
  stream.Bool(_isPartitioned)
      .UInt16(_polarization)
      .UInt32(_dataDescId)
      .Object(_selection)
      .String(_filename)
      .String(_dataColumnName)
      .Object(_partitionHandle)
      .UInt64(_partIndex);
}

std::unique_ptr<MSDataDescription> MSDataDescription::Unserialize(
    aocommon::SerialIStream& stream) {
  std::unique_ptr<MSDataDescription> mdd(new MSDataDescription());
  stream.Bool(mdd->_isPartitioned)
      .UInt16(mdd->_polarization)
      .UInt32(mdd->_dataDescId)
      .Object(mdd->_selection)
      .String(mdd->_filename)
      .String(mdd->_dataColumnName)
      .Object(mdd->_partitionHandle)
      .UInt64(mdd->_partIndex);
  mdd->_useMPI = true;  // Serialization only happens with MPI.
  return mdd;
}
