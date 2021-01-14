#ifndef MS_DATA_DESCRIPTION_H
#define MS_DATA_DESCRIPTION_H

#include "msprovider.h"
#include "partitionedms.h"

#include <aocommon/io/serialstreamfwd.h>

#include <memory>

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
      const MSSelection& selection, aocommon::PolarizationEnum polarization,
      size_t dataDescId) {
    std::unique_ptr<MSDataDescription> mdd(new MSDataDescription());
    mdd->_isPartitioned = false;
    mdd->_polarization = polarization;
    mdd->_dataDescId = dataDescId;
    mdd->_selection = selection;
    mdd->_filename = filename;
    mdd->_dataColumnName = dataColumnName;
    return mdd;
  }

  static std::unique_ptr<MSDataDescription> ForPartitioned(
      PartitionedMS::Handle partitionHandle, const MSSelection& selection,
      size_t partIndex, aocommon::PolarizationEnum polarization,
      size_t dataDescId) {
    std::unique_ptr<MSDataDescription> mdd(new MSDataDescription());
    mdd->_isPartitioned = true;
    mdd->_polarization = polarization;
    mdd->_dataDescId = dataDescId;
    mdd->_selection = selection;
    mdd->_partitionHandle = std::move(partitionHandle);
    mdd->_partIndex = partIndex;
    return mdd;
  }

  std::unique_ptr<MSProvider> GetProvider() const;

  const MSSelection& Selection() const { return _selection; }

  void Serialize(aocommon::SerialOStream& stream) const;
  static std::unique_ptr<MSDataDescription> Unserialize(
      aocommon::SerialIStream& stream);

 private:
  MSDataDescription(){};

  // Common
  bool _isPartitioned;
  aocommon::PolarizationEnum _polarization;
  size_t _dataDescId;
  MSSelection _selection;

  // Contiguous
  std::string _filename;
  std::string _dataColumnName;

  // Partitioned
  PartitionedMS::Handle _partitionHandle;
  size_t _partIndex;
};

#endif
