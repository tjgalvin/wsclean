#include "msdatadescription.h"

#include "contiguousms.h"

#include "../main/settings.h"
#include "reorderedmsprovider.h"

#include <aocommon/io/serialostream.h>
#include <aocommon/io/serialistream.h>

#include <cassert>
#include <memory>

using schaapcommon::reordering::StorageManagerType;

namespace wsclean {

std::unique_ptr<MSProvider> MSDataDescription::GetProvider() const {
  if (_isReordered)
    return std::make_unique<ReorderedMsProvider>(_reorderedHandle, _partIndex,
                                                 _polarization);
  else
    return std::make_unique<ContiguousMS>(
        _filename, _dataColumnName, _modelColumnName, _modelStorageManager,
        _selection, _channel_selection, _polarization, _useMPI);
}

void MSDataDescription::Serialize(aocommon::SerialOStream& stream) const {
  // Serialization is only used with MPI.
  assert(_useMPI);
  stream.Bool(_isReordered)
      .UInt16(_polarization)
      .String(_filename)
      .String(_dataColumnName)
      .String(_modelColumnName)
      .UInt32(static_cast<unsigned>(_modelStorageManager))
      .Object(_selection)
      .UInt64(_channel_selection.Size());
  for (const schaapcommon::reordering::ChannelRange& range :
       _channel_selection) {
    stream.UInt64(range.data_desc_id).UInt64(range.start).UInt64(range.end);
  }
  stream.Object(_reorderedHandle).UInt64(_partIndex);
}

std::unique_ptr<MSDataDescription> MSDataDescription::Unserialize(
    aocommon::SerialIStream& stream) {
  std::unique_ptr<MSDataDescription> mdd(new MSDataDescription());
  stream.Bool(mdd->_isReordered)
      .UInt16(mdd->_polarization)
      .String(mdd->_filename)
      .String(mdd->_dataColumnName)
      .String(mdd->_modelColumnName);
  mdd->_modelStorageManager = static_cast<StorageManagerType>(stream.UInt32());
  stream.Object(mdd->_selection);
  const size_t n = stream.UInt64();
  for (size_t i = 0; i != n; ++i) {
    schaapcommon::reordering::ChannelRange& range =
        mdd->_channel_selection.EmplaceBack();
    stream.UInt64(range.data_desc_id).UInt64(range.start).UInt64(range.end);
  }
  stream.Object(mdd->_reorderedHandle).UInt64(mdd->_partIndex);
  mdd->_useMPI = true;  // Serialization only happens with MPI.
  return mdd;
}

}  // namespace wsclean
