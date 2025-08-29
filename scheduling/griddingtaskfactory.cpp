#include "griddingtaskfactory.h"

#include <limits>

#include "../idg/averagebeam.h"
#include "../io/imagefilename.h"

namespace wsclean {

void GriddingTaskFactory::AddFacet(
    std::vector<GriddingTask>& tasks, const ImagingTableEntry& entry,
    std::unique_ptr<AverageBeam> average_beam,
    std::vector<aocommon::Image>&& model_images) {
  assert(!tasks.empty());  // AddFacet adds a facet to an existing task.

  if (!tasks.front().facets.empty()) {
    // Both the first facet data and the entry should contain a facet, in the
    // same facet group.
    assert(tasks.front().facets.front().facet);
    assert(entry.facet);
    assert(tasks.front().facetGroupIndex == entry.facetGroupIndex);
  }

  double l_shift = l_shift_;
  double m_shift = m_shift_;

  if (entry.facet) {
    l_shift -= entry.centreShiftX * settings_.pixelScaleX;
    m_shift += entry.centreShiftY * settings_.pixelScaleY;
  } else {
    assert(entry.facetIndex == 0);
  }
  tasks.back().facets.emplace_back(entry.facetIndex, l_shift, m_shift,
                                   std::move(metadata_cache_[entry.index]),
                                   std::move(average_beam), entry.facet,
                                   std::move(model_images));
}

aocommon::PolarizationEnum GriddingTaskFactory::DeterminePolarization(
    const ImagingTableEntry& entry) const {
  if (settings_.gridderType == GridderType::IDG &&
      settings_.polarizations.size() != 1)
    return aocommon::Polarization::FullStokes;
  else
    return entry.polarization;
}

GriddingTask GriddingTaskFactory::CreateBase(
    const ImagingTableEntry& entry, ImageWeightCache& image_weight_cache,
    bool is_first_task) {
  GriddingTask task;

  assert(entry.index <= std::numeric_limits<decltype(task.unique_id)>::max());
  task.unique_id = entry.index;
  task.observationInfo = observation_info_;
  task.isFirstTask = is_first_task;
  task.facetGroupIndex = entry.facetGroupIndex;
  task.outputChannelIndex = entry.outputChannelIndex;

  task.msList = ms_helper_.InitializeMsList(entry);
  task.imageWeights =
      InitializeWeights(settings_, entry, task.msList, image_weight_cache);

  return task;
}

std::vector<GriddingTask> GriddingTaskFactory::CreatePsfTasks(
    const ImagingTable::Group& facet_group,
    ImageWeightCache& image_weight_cache, bool combine_facets,
    bool is_first_task) {
  const bool store_imaging_weights = settings_.writeImagingWeightSpectrumColumn;

  std::vector<GriddingTask> tasks;
  tasks.reserve(combine_facets ? 1 : facet_group.size());

  for (const std::shared_ptr<ImagingTableEntry>& entry : facet_group) {
    // During PSF imaging, the average beam will never exist, so it is not
    // necessary to set the average beam in the task.
    if (tasks.empty() || !combine_facets) {
      GriddingTask task = CreateBase(*entry, image_weight_cache, is_first_task);
      task.operation = GriddingTask::Invert;
      task.imagePSF = true;
      task.polarization = entry->polarization;
      task.subtractModel = false;
      task.storeImagingWeights = store_imaging_weights;
      tasks.push_back(std::move(task));
    }
    AddFacet(tasks, *entry);

    is_first_task = false;
  }

  return tasks;
}

std::vector<GriddingTask> GriddingTaskFactory::CreateInvertTasks(
    const ImagingTable::Group& facet_group,
    ImageWeightCache& image_weight_cache, bool combine_facets,
    bool is_first_task, bool is_first_inversion,
    std::vector<std::unique_ptr<AverageBeam>>&& average_beams) {
  assert(average_beams.empty() || average_beams.size() == facet_group.size());

  std::vector<GriddingTask> tasks;
  tasks.reserve(combine_facets ? 1 : facet_group.size());

  for (std::size_t i = 0; i < facet_group.size(); ++i) {
    const ImagingTableEntry& entry = *facet_group[i];

    if (tasks.empty() || !combine_facets) {
      GriddingTask task = CreateBase(entry, image_weight_cache, is_first_task);

      task.operation = GriddingTask::Invert;
      task.imagePSF = false;
      task.polarization = DeterminePolarization(entry);
      task.subtractModel = !is_first_inversion || settings_.subtractModel ||
                           settings_.continuedRun;
      task.storeImagingWeights =
          is_first_inversion && settings_.writeImagingWeightSpectrumColumn;
      tasks.push_back(std::move(task));
    }

    std::unique_ptr<AverageBeam> average_beam;
    if (!average_beams.empty()) average_beam = std::move(average_beams[i]);
    AddFacet(tasks, entry, std::move(average_beam));

    is_first_task = false;
  }

  return tasks;
}

std::vector<GriddingTask> GriddingTaskFactory::CreatePredictTasks(
    const ImagingTable::Group& facet_group,
    ImageWeightCache& image_weight_cache, bool combine_facets,
    std::vector<std::vector<aocommon::Image>>&& model_images,
    std::vector<std::unique_ptr<AverageBeam>>&& average_beams) {
  assert(model_images.size() == facet_group.size());
  assert(average_beams.empty() || average_beams.size() == facet_group.size());

  std::vector<GriddingTask> tasks;
  tasks.reserve(combine_facets ? 1 : facet_group.size());

  for (std::size_t i = 0; i < facet_group.size(); ++i) {
    const ImagingTableEntry& entry = *facet_group[i];

    if (tasks.empty() || !combine_facets) {
      GriddingTask task = CreateBase(entry, image_weight_cache, false);
      task.operation = GriddingTask::Predict;
      task.imagePSF = false;
      task.polarization = DeterminePolarization(entry);
      task.subtractModel = false;
      task.storeImagingWeights = false;
      tasks.push_back(std::move(task));
    }

    std::unique_ptr<AverageBeam> average_beam;
    if (!average_beams.empty()) average_beam = std::move(average_beams[i]);
    AddFacet(tasks, entry, std::move(average_beam), std::move(model_images[i]));
  }

  return tasks;
}

}  // namespace wsclean
