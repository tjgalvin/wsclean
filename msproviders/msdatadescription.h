#ifndef MS_DATA_DESCRIPTION_H
#define MS_DATA_DESCRIPTION_H

#include "msprovider.h"
#include "reorderedmsprovider.h"

#include <aocommon/io/serialstreamfwd.h>
#include <schaapcommon/reordering/storagemanagertype.h>

#include <memory>

namespace wsclean {

/**
 * This class contains all the information necessary to open
 * a dataset. In particular, it provides all the information
 * to create an MSProvider object.
 *
 * For distributed computations, an object of this class can
 * be transferred to another node, and thereby provide all the
 * information to that node for reading the data.
 */
class MSDataDescription {
 public:
  static std::unique_ptr<MSDataDescription> ForContiguous(
      const std::string& filename, const std::string& dataColumnName,
      const std::string& modelColumnName,
      schaapcommon::reordering::StorageManagerType modelStorageManager,
      const schaapcommon::reordering::MSSelection& selection,
      aocommon::PolarizationEnum polarization, size_t dataDescId, bool useMPI) {
    std::unique_ptr<MSDataDescription> mdd(new MSDataDescription());
    mdd->_isReordered = false;
    mdd->_useMPI = useMPI;
    mdd->_polarization = polarization;
    mdd->_dataDescId = dataDescId;
    mdd->_selection = selection;
    mdd->_filename = filename;
    mdd->_dataColumnName = dataColumnName;
    mdd->_modelColumnName = modelColumnName;
    mdd->_modelStorageManager = modelStorageManager;
    return mdd;
  }

  static std::unique_ptr<MSDataDescription> ForReordered(
      ReorderedMsProvider::ReorderedHandle reorderedHandle,
      const schaapcommon::reordering::MSSelection& selection, size_t partIndex,
      aocommon::PolarizationEnum polarization, size_t dataDescId, bool useMPI) {
    std::unique_ptr<MSDataDescription> mdd(new MSDataDescription());
    mdd->_isReordered = true;
    mdd->_useMPI = useMPI;
    mdd->_polarization = polarization;
    mdd->_dataDescId = dataDescId;
    mdd->_selection = selection;
    mdd->_reorderedHandle = std::move(reorderedHandle);
    mdd->_partIndex = partIndex;
    return mdd;
  }

  std::unique_ptr<MSProvider> GetProvider() const;

  /**
   * A MSSelection object that identifies the data range of the
   * measurement set that is selected. This includes separating
   * channels caused by e.g. -channels-out and -channel-range.
   */
  const schaapcommon::reordering::MSSelection& Selection() const {
    return _selection;
  }

  size_t DataDescId() const { return _dataDescId; }

  void Serialize(aocommon::SerialOStream& stream) const;
  static std::unique_ptr<MSDataDescription> Unserialize(
      aocommon::SerialIStream& stream);

 private:
  MSDataDescription(){};

  // Common
  bool _isReordered;
  bool _useMPI;
  aocommon::PolarizationEnum _polarization;
  size_t _dataDescId;
  schaapcommon::reordering::MSSelection _selection;

  // Contiguous
  std::string _filename;
  std::string _dataColumnName;
  std::string _modelColumnName;
  schaapcommon::reordering::StorageManagerType _modelStorageManager;

  // Reordered
  ReorderedMsProvider::ReorderedHandle _reorderedHandle;
  size_t _partIndex;
};

}  // namespace wsclean

#endif
