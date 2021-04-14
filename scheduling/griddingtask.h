#ifndef GRIDDING_TASK_H
#define GRIDDING_TASK_H

#include <aocommon/io/serialstreamfwd.h>
#include <aocommon/polarization.h>

#include "../structures/image.h"
#include "../structures/imageweights.h"
#include "../structures/observationinfo.h"

#include "../msproviders/msdatadescription.h"

#include "metadatacache.h"

namespace schaapcommon {
namespace facets {
class Facet;
}
}  // namespace schaapcommon

class GriddingTask {
 public:
  enum Operation { Invert, Predict } operation;
  bool imagePSF;
  bool subtractModel;
  aocommon::PolarizationEnum polarization;
  bool verbose;
  std::unique_ptr<MetaDataCache> cache;
  bool storeImagingWeights;

  std::shared_ptr<ImageWeights> imageWeights;
  std::vector<std::unique_ptr<MSDataDescription>> msList;

  // For prediction
  bool addToModel;
  Image modelImageReal;
  Image modelImageImaginary;
  ObservationInfo observationInfo;

  const schaapcommon::facets::Facet* facet = nullptr;
  size_t facetIndex;

  void Serialize(aocommon::SerialOStream& stream) const;
  void Unserialize(aocommon::SerialIStream& stream);
};

#endif
