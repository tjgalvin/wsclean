#include "griddingtask.h"

#include <limits>

#include "../idg/averagebeam.h"

#include <aocommon/io/serialostream.h>
#include <aocommon/io/serialistream.h>
#include <schaapcommon/facets/facet.h>

namespace wsclean {

void GriddingTask::Serialize(aocommon::SerialOStream& stream) const {
  stream.UInt32(unique_id)
      .UInt32(operation)
      .Bool(imagePSF)
      .Bool(subtractModel)
      .UInt32(polarization)
      .Bool(isFirstTask)
      .Bool(storeImagingWeights)
      .Ptr(imageWeights)
      .Object(observationInfo)
      .ObjectVector(facets)
      .UInt64(facetGroupIndex);

  // msList
  stream.UInt64(msList.size());
  for (const MsListItem& item : msList) item.Serialize(stream);
}

void GriddingTask::Unserialize(aocommon::SerialIStream& stream) {
  stream.UInt32(unique_id);
  operation = static_cast<Operation>(stream.UInt32());
  stream.Bool(imagePSF).Bool(subtractModel);
  polarization = static_cast<aocommon::PolarizationEnum>(stream.UInt32());
  stream.Bool(isFirstTask)
      .Bool(storeImagingWeights)
      .Ptr(imageWeights)
      .Object(observationInfo)
      .ObjectVector(facets)
      .UInt64(facetGroupIndex);

  assert(!facets.empty());

  // msList
  msList.resize(stream.UInt64());
  for (MsListItem& item : msList) item.Unserialize(stream);
}

GriddingTask::FacetData::FacetData(
    size_t _index, double _l_shift, double _m_shift,
    std::unique_ptr<MetaDataCache> _cache,
    std::unique_ptr<AverageBeam> _average_beam,
    const std::shared_ptr<schaapcommon::facets::Facet>& _facet,
    std::vector<aocommon::Image>&& _model_images)
    : index{_index},
      l_shift{_l_shift},
      m_shift{_m_shift},
      cache{std::move(_cache)},
      averageBeam{std::move(_average_beam)},
      facet{_facet},
      modelImages{std::move(_model_images)} {}

void GriddingTask::FacetData::Serialize(aocommon::SerialOStream& stream) const {
  stream.ObjectVector(modelImages)
      .UInt64(index)
      .Double(l_shift)
      .Double(m_shift)
      .Ptr(cache)
      .Ptr(averageBeam)
      .Ptr(facet);
}

void GriddingTask::FacetData::Unserialize(aocommon::SerialIStream& stream) {
  stream.ObjectVector(modelImages)
      .UInt64(index)
      .Double(l_shift)
      .Double(m_shift)
      .Ptr(cache)
      .Ptr(averageBeam)
      .Ptr(facet);
}

}  // namespace wsclean
