#ifndef GRIDDING_TASK_H
#define GRIDDING_TASK_H

#include <aocommon/io/serialstreamfwd.h>
#include <aocommon/image.h>
#include <aocommon/polarization.h>

#include "../structures/imageweights.h"
#include "../structures/observationinfo.h"

#include "../msproviders/msdatadescription.h"

#include "metadatacache.h"

class AverageBeam;

namespace schaapcommon {
namespace facets {
class Facet;
}
}  // namespace schaapcommon

class GriddingTask {
 public:
  GriddingTask();
  GriddingTask(const GriddingTask&) = delete;
  GriddingTask(GriddingTask&& source) noexcept;
  ~GriddingTask() noexcept;
  GriddingTask& operator=(const GriddingTask& source) = delete;
  GriddingTask& operator=(GriddingTask&& source) noexcept;

  enum Operation { Invert, Predict } operation;
  bool imagePSF;
  bool subtractModel;
  aocommon::PolarizationEnum polarization;
  bool verbose;
  std::unique_ptr<MetaDataCache> cache;
  std::unique_ptr<AverageBeam> averageBeam;
  bool storeImagingWeights;

  std::shared_ptr<ImageWeights> imageWeights;
  std::vector<std::unique_ptr<MSDataDescription>> msList;

  /**
   * Images for prediction. See the documentation of
   * @ref GriddingResult::images for an explanation of why this is a vector.
   */
  std::vector<aocommon::Image> modelImages;
  ObservationInfo observationInfo;
  double shiftL;
  double shiftM;

  std::shared_ptr<schaapcommon::facets::Facet> facet;
  size_t facetIndex;
  size_t facetGroupIndex;

  void Serialize(aocommon::SerialOStream& stream) const;
  void Unserialize(aocommon::SerialIStream& stream);
};

#endif
