#ifndef GRIDDING_TASK_H
#define GRIDDING_TASK_H

#include <aocommon/polarization.h>
#include "../image.h"
#include "../imageweights.h"

#include "../msproviders/msdatadescription.h"

#include "metadatacache.h"

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

  void Serialize(class SerialOStream& stream) const;
  void Unserialize(class SerialIStream& stream);
};

#endif
