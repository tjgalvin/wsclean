#include "../../scheduling/griddingtaskfactory.h"

#include <boost/test/unit_test.hpp>

#include <schaapcommon/facets/facet.h>

#include "../../idg/averagebeam.h"

using schaapcommon::facets::Facet;
using schaapcommon::reordering::MSSelection;

namespace wsclean {

namespace {
constexpr double kPixelScale = 0.0042;
constexpr size_t kImageSize = 142;
constexpr float kImageValue = 1.42;
constexpr double kLShift = 0.42;
constexpr double kMShift = 0.1;
constexpr bool kCombineFacets = true;
constexpr bool kIsFirstTask = true;
constexpr bool kIsFirstInversion = true;

std::shared_ptr<Facet> CreateFacet() {
  Facet::InitializationData initialization_data{kPixelScale, kPixelScale,
                                                kImageSize, kImageSize};
  schaapcommon::facets::BoundingBox box;
  return std::make_shared<Facet>(initialization_data, box);
}

struct FactoryFixture {
  FactoryFixture()
      : settings(),
        global_selection(),
        ms_bands(),
        ms_helper(settings, global_selection, ms_bands),
        image_weight_initializer(settings, global_selection, ms_bands,
                                 ms_helper.GetReorderedMsHandles()),
        observation_info{42.0, 6.0, "test_name", "test_observer", "test_field"},
        image_weight_cache(WeightMode(WeightClass::Natural), kImageSize,
                           kImageSize, kPixelScale, kPixelScale, 0.0, 1.0, 0.0,
                           0, false),
        group(2),
        average_beams(group.size()),
        average_beam_pointers(group.size()),
        model_images(group.size()),
        metadata_cache_pointers(group.size()),
        factory(ms_helper, image_weight_initializer, observation_info, kLShift,
                kMShift, group.size()) {
    settings.pixelScaleX = kPixelScale;
    settings.pixelScaleY = kPixelScale;
    settings.writeImagingWeightSpectrumColumn = true;

    // Create two entries with identical facet group index:
    group[0] = std::make_shared<ImagingTableEntry>();
    group[0]->index = 0;
    group[0]->polarization = aocommon::PolarizationEnum::LR;
    group[0]->facetIndex = 4;
    group[0]->facetGroupIndex = 42;
    group[0]->centreShiftX = 10;
    group[0]->centreShiftY = 20;
    group[0]->facet = CreateFacet();

    group[1] = std::make_shared<ImagingTableEntry>();
    group[1]->index = 1;
    group[1]->polarization = aocommon::PolarizationEnum::XY;
    group[1]->facetIndex = 2;
    group[1]->facetGroupIndex = group[0]->facetGroupIndex;
    group[1]->centreShiftX = 30;
    group[1]->centreShiftY = 40;
    group[1]->facet = CreateFacet();

    for (size_t index = 0; index < group.size(); ++index) {
      // Create meta data cache for each entry.
      ImagingTableEntry entry;
      entry.index = index;
      auto cache = std::make_unique<MetaDataCache>();
      metadata_cache_pointers[index] = cache.get();
      factory.SetMetaDataCacheEntry(*group[index], std::move(cache));

      // Create an AverageBeam object for each entry.
      average_beams[index] = std::make_unique<AverageBeam>();
      average_beam_pointers[index] = average_beams[index].get();

      model_images[index].emplace_back(kImageSize, kImageSize,
                                       kImageValue + index);
    }
  }

  /// Performs common checks for the main task fields.
  void CheckCommonTaskData(const GriddingTask& task,
                           const ImagingTableEntry& entry) const {
    BOOST_TEST(task.polarization == entry.polarization);
    BOOST_TEST(task.imageWeights);
    BOOST_TEST(task.msList.empty());
    BOOST_CHECK(task.observationInfo == observation_info);
    BOOST_TEST(task.facetGroupIndex == entry.facetGroupIndex);
  }

  /// Performs common checks for facet data in a task.
  void CheckCommonFacetData(const GriddingTask::FacetData& facet_data,
                            const ImagingTableEntry& entry) const {
    BOOST_TEST(facet_data.index == entry.facetIndex);
    BOOST_CHECK_CLOSE(facet_data.l_shift,
                      kLShift - entry.centreShiftX * kPixelScale, 1.0e-7);
    BOOST_CHECK_CLOSE(facet_data.m_shift,
                      kMShift + entry.centreShiftY * kPixelScale, 1.0e-7);
    BOOST_TEST(facet_data.cache.get() == metadata_cache_pointers[entry.index]);
    BOOST_TEST(facet_data.facet == entry.facet);
  }

  void CheckPsfTask(const GriddingTask& task, ImagingTable::Group group,
                    bool is_first_task) const {
    CheckCommonTaskData(task, *group.front());
    BOOST_TEST(task.operation == GriddingTask::Invert);
    BOOST_TEST(task.imagePSF == true);
    BOOST_TEST(task.subtractModel == false);
    BOOST_TEST(task.isFirstTask == is_first_task);
    BOOST_TEST(task.storeImagingWeights ==
               settings.writeImagingWeightSpectrumColumn);

    BOOST_REQUIRE(task.facets.size() == group.size());
    for (size_t i = 0; i < group.size(); ++i) {
      const GriddingTask::FacetData& facet_data = task.facets[i];
      CheckCommonFacetData(facet_data, *group[i]);
      BOOST_TEST(facet_data.modelImages.empty());
      BOOST_TEST(!facet_data.averageBeam);
    }
  };

  void CheckInvertTask(const GriddingTask& task, ImagingTable::Group group,
                       bool is_first_task) const {
    CheckCommonTaskData(task, *group.front());
    BOOST_TEST(task.operation == GriddingTask::Invert);
    BOOST_TEST(task.imagePSF == false);
    BOOST_TEST(task.subtractModel == false);
    BOOST_TEST(task.isFirstTask == is_first_task);
    BOOST_TEST(task.storeImagingWeights ==
               settings.writeImagingWeightSpectrumColumn);

    BOOST_TEST_REQUIRE(task.facets.size() == group.size());
    for (size_t i = 0; i < group.size(); ++i) {
      const GriddingTask::FacetData& facet_data = task.facets[i];
      CheckCommonFacetData(facet_data, *group[i]);
      BOOST_TEST(facet_data.modelImages.empty());
      BOOST_TEST(facet_data.averageBeam.get() ==
                 average_beam_pointers[group[i]->index]);
    }
  }

  void CheckPredictTask(const GriddingTask& task,
                        ImagingTable::Group group) const {
    CheckCommonTaskData(task, *group.front());
    BOOST_TEST(task.operation == GriddingTask::Predict);
    BOOST_TEST(!task.isFirstTask);

    BOOST_TEST_REQUIRE(task.facets.size() == group.size());
    for (size_t i = 0; i < group.size(); ++i) {
      const GriddingTask::FacetData& facet_data = task.facets[i];
      CheckCommonFacetData(facet_data, *group[i]);
      const size_t entry_index = group[i]->index;
      BOOST_TEST(facet_data.modelImages == model_images[entry_index]);
      BOOST_TEST(facet_data.averageBeam.get() ==
                 average_beam_pointers[entry_index]);
    }
  }

  // For initializing the Initializer.
  Settings settings;
  MSSelection global_selection;
  std::vector<aocommon::MultiBandData> ms_bands;

  // For initializing a GriddingTaskFactory:
  MsHelper ms_helper;
  ImageWeightInitializer image_weight_initializer;
  ObservationInfo observation_info;

  // Data for the Create* calls in the tests.
  ImageWeightCache image_weight_cache;

  ImagingTable::Group group;

  std::vector<std::unique_ptr<AverageBeam>> average_beams;
  std::vector<AverageBeam*> average_beam_pointers;

  std::vector<std::vector<aocommon::Image>> model_images;

  // Gets pointers to the cache items in the factory.
  std::vector<MetaDataCache*> metadata_cache_pointers;

  // The factory object, which is the object under test.
  GriddingTaskFactory factory;
};
}  // namespace

BOOST_AUTO_TEST_SUITE(gridding_task_factory)

BOOST_FIXTURE_TEST_CASE(get_meta_data_cache, FactoryFixture) {
  BOOST_REQUIRE(factory.GetMetaDataCache().size() == group.size());
  for (size_t i = 0; i < group.size(); ++i) {
    BOOST_TEST(factory.GetMetaDataCache()[i].get() ==
               metadata_cache_pointers[i]);
  }
}

BOOST_FIXTURE_TEST_CASE(psf_separate_tasks, FactoryFixture) {
  const std::vector<GriddingTask> separate_tasks = factory.CreatePsfTasks(
      group, image_weight_cache, !kCombineFacets, kIsFirstTask);

  BOOST_REQUIRE(separate_tasks.size() == group.size());
  CheckPsfTask(separate_tasks[0], {group[0]}, kIsFirstTask);
  CheckPsfTask(separate_tasks[1], {group[1]}, !kIsFirstTask);
}

BOOST_FIXTURE_TEST_CASE(psf_combined_facets, FactoryFixture) {
  const std::vector<GriddingTask> combined_tasks = factory.CreatePsfTasks(
      group, image_weight_cache, kCombineFacets, kIsFirstTask);

  BOOST_REQUIRE(combined_tasks.size() == 1);
  CheckPsfTask(combined_tasks[0], group, kIsFirstTask);
}

BOOST_FIXTURE_TEST_CASE(invert_separate_tasks, FactoryFixture) {
  const std::vector<GriddingTask> separate_tasks = factory.CreateInvertTasks(
      group, image_weight_cache, !kCombineFacets, kIsFirstTask,
      kIsFirstInversion, std::move(average_beams));

  BOOST_REQUIRE(separate_tasks.size() == group.size());
  CheckInvertTask(separate_tasks[0], {group[0]}, kIsFirstTask);
  CheckInvertTask(separate_tasks[1], {group[1]}, !kIsFirstTask);
}

BOOST_FIXTURE_TEST_CASE(invert_combined_facets, FactoryFixture) {
  const std::vector<GriddingTask> combined_tasks = factory.CreateInvertTasks(
      group, image_weight_cache, kCombineFacets, kIsFirstTask,
      kIsFirstInversion, std::move(average_beams));

  BOOST_REQUIRE(combined_tasks.size() == 1);
  CheckInvertTask(combined_tasks[0], group, kIsFirstTask);
}

BOOST_FIXTURE_TEST_CASE(predict_separate_tasks, FactoryFixture) {
  const std::vector<GriddingTask> separate_tasks = factory.CreatePredictTasks(
      group, image_weight_cache, !kCombineFacets,
      std::vector<std::vector<aocommon::Image>>{model_images},
      std::move(average_beams));

  BOOST_REQUIRE(separate_tasks.size() == group.size());
  CheckPredictTask(separate_tasks[0], {group[0]});
  CheckPredictTask(separate_tasks[1], {group[1]});
}

BOOST_FIXTURE_TEST_CASE(predict_combined_facets, FactoryFixture) {
  const std::vector<GriddingTask> combined_tasks = factory.CreatePredictTasks(
      group, image_weight_cache, kCombineFacets,
      std::vector<std::vector<aocommon::Image>>{model_images},
      std::move(average_beams));

  BOOST_REQUIRE(combined_tasks.size() == 1);
  CheckPredictTask(combined_tasks[0], group);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace wsclean
