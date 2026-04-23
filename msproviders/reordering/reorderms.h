#ifndef WSCLEAN_MSPROVIDERS_REORDERMS_H_
#define WSCLEAN_MSPROVIDERS_REORDERMS_H_

#include <string>
#include <vector>

#include <aocommon/vectormap.h>

#include <schaapcommon/reordering/storagemanagertype.h>

#include "reorderedhandle.h"

namespace schaapcommon::reordering {
struct ChannelRange;
class MSSelection;
}  // namespace schaapcommon::reordering

namespace wsclean {

class Settings;

ReorderedHandle ReorderMS(
    const std::string& ms_path,
    const std::vector<aocommon::VectorMap<schaapcommon::reordering::ChannelRange>>& channels,
    const schaapcommon::reordering::MSSelection& selection, const std::string& data_column_name,
    const std::string& model_column_name,
    schaapcommon::reordering::StorageManagerType model_storage_manager, bool include_model,
    bool initial_model_required, const Settings& settings);

}  // namespace wsclean

#endif
