#ifndef GRIDDING_TASK_H
#define GRIDDING_TASK_H

#include <aocommon/io/serialstreamfwd.h>
#include <aocommon/image.h>
#include <aocommon/polarization.h>

#include "../structures/imageweights.h"
#include "../structures/mslistitem.h"
#include "../structures/observationinfo.h"
#include "../system/completionsignal.h"

#include "metadatacache.h"

namespace schaapcommon {
namespace facets {
class Facet;
}
}  // namespace schaapcommon

namespace wsclean {

class AverageBeam;

class GriddingTask {
 public:
  GriddingTask() = default;
  GriddingTask(const GriddingTask&) = delete;
  GriddingTask(GriddingTask&& source) noexcept = default;
  ~GriddingTask() noexcept = default;
  GriddingTask& operator=(const GriddingTask& source) = delete;
  GriddingTask& operator=(GriddingTask&& source) noexcept = default;

  uint32_t unique_id;
  enum Operation {
    Invert,
    Predict,
    Wait  // Consumes a thread from the task pool until the task it is
          // associated with is complete.
  } operation;
  bool imagePSF;       // Only for invert tasks.
  bool subtractModel;  // Only for invert tasks.
  aocommon::PolarizationEnum polarization;
  bool isFirstTask;
  bool storeImagingWeights;  // Only for invert tasks.

  std::shared_ptr<ImageWeights> imageWeights;
  std::vector<MsListItem> msList;

  ObservationInfo observationInfo;

  struct FacetData {
    /// The default constructor is required since deserializing an ObjectVector
    /// first creates empty objects and then calls Unserialize on those objects.
    FacetData() = default;

    explicit FacetData(
        size_t _index, double _l_shift, double _m_shift,
        std::unique_ptr<MetaDataCache> _cache,
        std::unique_ptr<AverageBeam> _average_beam,
        const std::shared_ptr<schaapcommon::facets::Facet>& _facet,
        std::vector<aocommon::Image>&& _model_images);

    void Serialize(aocommon::SerialOStream& stream) const;
    void Unserialize(aocommon::SerialIStream& stream);

    size_t index;    ///< Index of the facet, between zero and n_facets.
    double l_shift;  ///< l_shift, adjusted to the center of the facet.
    double m_shift;  ///< m_shift, adjusted to the center of the facet.
    std::unique_ptr<MetaDataCache> cache;
    std::unique_ptr<AverageBeam> averageBeam;
    /// The facet itself. If null, faceting is disabled.
    std::shared_ptr<schaapcommon::facets::Facet> facet;

    /// Images for prediction. The documentation of
    /// @ref GriddingResult::FacetData::images explains why it is a vector.
    std::vector<aocommon::Image> modelImages;
  };

  /// 'facets' always contains at least one element.
  /// When faceting is disabled, there is always a single element with
  /// a null 'facet' pointer. That 'facet' covers the entire image.
  std::vector<FacetData> facets;
  size_t facetGroupIndex;

  /// This variable should not be serialized. Only the main MPI
  /// node uses it to select a target node according to the
  /// channel-to-node option in the command line settings.
  size_t outputChannelIndex;

  void Serialize(aocommon::SerialOStream& stream) const;
  void Unserialize(aocommon::SerialIStream& stream);

  // Below members are used only when doing shared reads
  // Stop scheduler resource blockers from closing until task is complete
  std::shared_ptr<CompletionSignal> lock_excess_scheduler_tasks_;
  // Number of parallel gridders that we are allowed to launch. NB! This is not
  // necessarilly the same as the global number.
  size_t num_parallel_gridders_ = 1;
};

}  // namespace wsclean

#endif
