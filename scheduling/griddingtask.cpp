#include "griddingtask.h"

#include <aocommon/io/serialostream.h>
#include <aocommon/io/serialistream.h>

void GriddingTask::Serialize(aocommon::SerialOStream& stream) const {
  stream.UInt32(operation)
      .Bool(imagePSF)
      .Bool(subtractModel)
      .UInt32(polarization)
      .Bool(verbose)
      .Ptr(cache)
      .Bool(storeImagingWeights)
      .Ptr(imageWeights);

  // msList
  stream.UInt64(msList.size());
  for (const std::unique_ptr<MSDataDescription>& dataDesc : msList)
    dataDesc->Serialize(stream);

  stream.Bool(addToModel);
  modelImageReal.Serialize(stream);
  modelImageImaginary.Serialize(stream);
  observationInfo.Serialize(stream);
}

void GriddingTask::Unserialize(aocommon::SerialIStream& stream) {
  operation = (Operation)stream.UInt32();
  stream.Bool(imagePSF).Bool(subtractModel);
  polarization = (aocommon::PolarizationEnum)stream.UInt32();
  stream.Bool(verbose).Ptr(cache).Bool(storeImagingWeights).Ptr(imageWeights);

  // msList
  msList.resize(stream.UInt64());
  for (std::unique_ptr<MSDataDescription>& dataDesc : msList)
    dataDesc = MSDataDescription::Unserialize(stream);

  addToModel = stream.Bool();
  modelImageReal.Unserialize(stream);
  modelImageImaginary.Unserialize(stream);
  observationInfo.Unserialize(stream);
}
