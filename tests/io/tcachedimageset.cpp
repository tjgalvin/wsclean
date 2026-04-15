#include "../../io/cachedimageaccessor.h"
#include "../../io/cachedimageset.h"
#include "../../main/stitchingmemory.h"
#include <radler/image_set.h>

#include <aocommon/fits/fitswriter.h>
#include <aocommon/image.h>
#include <aocommon/polarization.h>
#include <aocommon/uvector.h>

#include <schaapcommon/facets/facet.h>
#include <schaapcommon/facets/facetimage.h>
#include <schaapcommon/fitters/spectralfitter.h>

#include <boost/test/unit_test.hpp>

#include <atomic>
#include <chrono>
#include <cmath>
#include <limits>
#include <memory>
#include <thread>

using aocommon::FitsWriter;
using aocommon::Image;
using aocommon::PolarizationEnum;
using radler::ImageSet;
using schaapcommon::facets::Facet;

namespace wsclean {

namespace {

class DummyImageAccessor : public aocommon::ImageAccessor {
 public:
  DummyImageAccessor() {}
  ~DummyImageAccessor() override {}

  void Load(float*) const override { BOOST_FAIL("Unexpected ImageAccessor::Load() call"); }

  void Store(const float*) override { BOOST_FAIL("Unexpected ImageAccessor::Store() call"); }

  std::size_t Width() const override {
    BOOST_FAIL("Unexpected ImageAccessor::Width() call");
    return 0;
  }

  std::size_t Height() const override {
    BOOST_FAIL("Unexpected ImageAccessor::Height() call");
    return 0;
  }
};

}  // namespace

class ImageSetFixtureBase {
 public:
  ImageSetFixtureBase() {}

  void initTable(size_t n_original_channels, size_t n_deconvolution_channels) {
    std::vector<radler::PsfOffset> psf_offsets;
    table_ = std::make_unique<radler::WorkTable>(psf_offsets, n_original_channels,
                                                 n_deconvolution_channels);
  }

  void addToImageSet(size_t outChannel, PolarizationEnum pol, size_t frequencyMHz,
                     double imageWeight = 1.0) {
    auto e = std::make_unique<radler::WorkTableEntry>();
    e->original_channel_index = outChannel;
    e->polarization = pol;
    e->band_start_frequency = frequencyMHz;
    e->band_end_frequency = frequencyMHz;
    e->image_weight = imageWeight;
    e->model_accessor = std::make_unique<CachedImageAccessor>(cSet_, pol, outChannel, false);
    e->residual_accessor = std::make_unique<DummyImageAccessor>();
    table_->AddEntry(std::move(e));
  }

  std::unique_ptr<radler::WorkTable> table_;
  FitsWriter writer_;
  CachedImageSet cSet_;
};

template <size_t NDeconvolutionChannels>
class ImageSetFixture : public ImageSetFixtureBase {
 public:
  ImageSetFixture() : image_(4, 0.0) {
    initTable(2, NDeconvolutionChannels);
    writer_.SetImageDimensions(2, 2);
    cSet_.Initialize(writer_, 2, 2, 0, "wsctest");
    addToImageSet(0, aocommon::Polarization::XX, 100);
    addToImageSet(0, aocommon::Polarization::YY, 100);
    addToImageSet(1, aocommon::Polarization::XX, 200);
    addToImageSet(1, aocommon::Polarization::YY, 200);
    image_[0] = 2.0;
    cSet_.Store(image_.data(), aocommon::Polarization::XX, 0, false);
    image_[0] = -1.0;
    cSet_.Store(image_.data(), aocommon::Polarization::YY, 0, false);
    image_[0] = 20.0;
    cSet_.Store(image_.data(), aocommon::Polarization::XX, 1, false);
    image_[0] = -10.0;
    cSet_.Store(image_.data(), aocommon::Polarization::YY, 1, false);
  }

  aocommon::UVector<double> image_;
};

BOOST_AUTO_TEST_SUITE(cachedimageset)

BOOST_FIXTURE_TEST_CASE(load, ImageSetFixture<1>) {
  cSet_.Load(image_.data(), aocommon::Polarization::XX, 1, false);
  BOOST_CHECK_EQUAL(image_[0], 20.0);
  cSet_.Load(image_.data(), aocommon::Polarization::YY, 1, false);
  BOOST_CHECK_EQUAL(image_[0], -10.0);
  cSet_.Load(image_.data(), aocommon::Polarization::XX, 0, false);
  BOOST_CHECK_EQUAL(image_[0], 2.0);
  cSet_.Load(image_.data(), aocommon::Polarization::YY, 0, false);
  BOOST_CHECK_EQUAL(image_[0], -1.0);
}

BOOST_FIXTURE_TEST_CASE(load_and_average, ImageSetFixtureBase) {
  // Almost equivalent to radler::timageset.cc::load_and_average
  // test, except now the CachedImageSet class is used.
  initTable(6, 2);
  const size_t nPol = 2;
  const PolarizationEnum pols[nPol] = {PolarizationEnum::XX, PolarizationEnum::YY};
  const size_t width = 7;
  const size_t height = 9;
  writer_.SetImageDimensions(width, height);
  const std::vector<double> weights{4.0, 4.0, 0.0, 0.0, 1.0, 1.0};
  cSet_.Initialize(writer_, 4, 6, 0, "imagesettest");
  Image storedImage(width, height);
  for (size_t ch = 0; ch != table_->OriginalGroups().size(); ++ch) {
    for (size_t p = 0; p != nPol; ++p) {
      size_t index = ch * nPol + p;
      addToImageSet(ch, pols[p], 100 + ch, weights[ch]);

      storedImage = (1 << index);  // assign the entire image to 2^index
      cSet_.Store(storedImage.Data(), pols[p], ch, false);
    }
  }
  const std::set<PolarizationEnum> kLinkedPolarizations{aocommon::Polarization::XX,
                                                        aocommon::Polarization::YY};

  ImageSet imageSet(*table_, false, kLinkedPolarizations, width, height);
  imageSet.LoadAndAverage(false);
  // The first image has all values set to 2^0, the second image 2^1, etc...
  // The XX polarizations of deconvolution channel 1 consists of
  // images 0, 2 and 4. These have been weighted with 4, 4, 0:
  BOOST_CHECK_CLOSE_FRACTION(imageSet[0 * nPol + 0][0], double(1 * 4 + 4 * 4 + 16 * 0) / 8.0, 1e-6);
  // The YY polarizations consists of images 1, 3 and 5, weights 4, 4, 0:
  BOOST_CHECK_CLOSE_FRACTION(imageSet[0 * nPol + 1][0], double(2 * 4 + 8 * 4 + 32 * 0) / 8.0, 1e-6);
  // The XX polarizations of deconvolution channel 2 consists of images 6, 8 and
  // 10 Weights 0, 1, 1
  BOOST_CHECK_CLOSE_FRACTION(imageSet[1 * nPol + 0][0], double(64 * 0 + 256 * 1 + 1024 * 1) / 2.0,
                             1e-6);
  // YY: images 7, 9, 10, weights 0, 1, 1
  BOOST_CHECK_CLOSE_FRACTION(imageSet[1 * nPol + 1][0], double(128 * 0 + 512 * 1 + 2048 * 1) / 2.0,
                             1e-6);

  // The total linear integrated sum should be a complete
  // weighting of all input channels
  Image linearIntegrated(width, height);
  imageSet.GetLinearIntegrated(linearIntegrated);
  BOOST_CHECK_CLOSE_FRACTION(linearIntegrated[0],
                             double(1 * 4 + 4 * 4 + 16 * 0 + 2 * 4 + 8 * 4 + 32 * 0 + 64 * 0 +
                                    256 * 1 + 1024 * 1 + 128 * 0 + 512 * 1 + 2048 * 1) /
                                 20.0,
                             1e-6);
}

// NOTE: only check whether point x,y collides with the bounding box of a facet
size_t FacetCollision(const std::vector<std::shared_ptr<Facet>>& facets, int x, int y) {
  size_t facetCollision = std::numeric_limits<size_t>::max();
  for (size_t f = 0; f < facets.size(); ++f) {
    schaapcommon::facets::BoundingBox bbox(facets[f]->GetPixels(), 1, false);
    if (y >= bbox.Min().y && y < bbox.Max().y && x >= bbox.Min().x && x < bbox.Max().x) {
      facetCollision = f;
      break;
    }
  }
  return facetCollision;
}

BOOST_AUTO_TEST_CASE(store_and_load_facet) {
  std::string prefix = "facettest";
  size_t image_width = 8;
  size_t image_height = 8;
  double dl_dm = 0.0125;

  aocommon::FitsWriter writer;
  writer.SetImageDimensions(image_width, image_height, dl_dm, dl_dm);

  // Make two 4x4 facets
  std::vector<schaapcommon::facets::Coord> coords0{
      {0.05, -0.05}, {0.0, -0.05}, {0.0, 0.0}, {0.05, 0.0}};

  // Second facet (i=1) is mirrored in origin
  std::vector<schaapcommon::facets::Coord> coords1;
  for (schaapcommon::facets::Coord coordinate : coords0) {
    coords1.emplace_back(-1.0 * coordinate.ra, -1.0 * coordinate.dec);
  }

  Facet::InitializationData facet_data(writer.PixelSizeX(), writer.PixelSizeY(), writer.Width(),
                                       writer.Height());
  facet_data.phase_centre.ra = writer.RA();
  facet_data.phase_centre.dec = writer.Dec();
  facet_data.l_shift = writer.LShift();
  facet_data.m_shift = writer.MShift();
  // The bounding box is padded such that it is partially outside the main image
  facet_data.padding = 1.5;

  std::vector<std::shared_ptr<Facet>> facets{
      std::make_shared<Facet>(facet_data, std::move(coords0)),
      std::make_shared<Facet>(facet_data, std::move(coords1))};
  std::vector<aocommon::Image> facet_images;
  facet_images.reserve(facets.size());

  for (size_t i = 0; i < facets.size(); ++i) {
    facet_images.emplace_back(facets[i]->GetTrimmedBoundingBox().Width(),
                              facets[i]->GetTrimmedBoundingBox().Height(),
                              static_cast<float>(i + 1));
  }

  CachedImageSet cSet;
  cSet.Initialize(writer, 2, 1, facets.size(), prefix);

  std::vector<aocommon::PolarizationEnum> polarizations{aocommon::Polarization::XX,
                                                        aocommon::Polarization::YY};

  for (const auto& polarization : polarizations) {
    for (size_t facet_idx = 0; facet_idx < facets.size(); ++facet_idx) {
      cSet.StoreFacet(facet_images[facet_idx], polarization, 1, facet_idx, facets[facet_idx],
                      false);
    }
  }

  // Retrieve the cached tmp filenames
  const std::set<std::string> stored_name_set = cSet.GetStoredNames();
  std::vector<std::string> storedNames(stored_name_set.begin(), stored_name_set.end());

  schaapcommon::facets::FacetImage imageStorage(image_width, image_height, 1);
  for (size_t pol_idx = 0; pol_idx < polarizations.size(); ++pol_idx) {
    aocommon::Image imageMain(image_width, image_height, 0.0f);
    for (size_t facet_idx = 0; facet_idx < facets.size(); ++facet_idx) {
      // Offset in file list
      size_t offset = pol_idx * facets.size() + facet_idx;
      imageStorage.SetFacet(*facets[facet_idx], true);
      BOOST_CHECK_EQUAL(storedNames[offset],
                        prefix + "-" +
                            aocommon::Polarization::TypeToShortString(polarizations[pol_idx]) +
                            "-f000" + std::to_string(facet_idx) + "-tmp.fits");

      size_t num_facet_pixels = facets[facet_idx]->GetTrimmedBoundingBox().Width() *
                                facets[facet_idx]->GetTrimmedBoundingBox().Height();
      cSet.LoadFacet(imageStorage.Data(0), polarizations[pol_idx], 1, facet_idx, facets[facet_idx],
                     false);
      imageStorage.AddToImage({imageMain.Data()});
      BOOST_CHECK_EQUAL_COLLECTIONS(imageStorage.Data(0), imageStorage.Data(0) + num_facet_pixels,
                                    facet_images[facet_idx].begin(), facet_images[facet_idx].end());
    }

    // Check whether data in imageMain is correct
    for (size_t y = 0; y != image_height; ++y) {
      for (size_t x = 0; x != image_width; ++x) {
        size_t offset = y * image_width + x;
        size_t facetCollision = FacetCollision(facets, static_cast<int>(x), static_cast<int>(y));
        BOOST_CHECK_EQUAL(imageMain[offset], (facetCollision == std::numeric_limits<size_t>::max())
                                                 ? 0.0f
                                                 : facetCollision + 1.0f);
      }
    }
    cSet.Store(imageMain.Data(), polarizations[pol_idx], 1, false);
    const std::set<std::string> stored_names = cSet.GetStoredNames();
    BOOST_CHECK(
        stored_names.find(prefix + "-" +
                          aocommon::Polarization::TypeToShortString(polarizations[pol_idx]) +
                          "-tmp.fits") != stored_names.end());
  }
}

BOOST_AUTO_TEST_CASE(parallel_store_with_local_writer) {
  aocommon::FitsWriter writer;
  writer.SetImageDimensions(2, 2);

  CachedImageSet cSet;
  cSet.Initialize(writer, 1, 4, 0, "parallel-store-test");

  std::vector<std::thread> threads;
  for (size_t channel = 0; channel != 4; ++channel) {
    threads.emplace_back([&, channel]() {
      WSCFitsWriter local_writer(writer);
      Image image(2, 2, static_cast<float>(channel + 1));
      cSet.Store(image, local_writer, aocommon::Polarization::StokesI, channel, false);
    });
  }

  for (std::thread& thread : threads) thread.join();

  BOOST_TEST(cSet.GetStoredNames().size() == 4u);
  for (size_t channel = 0; channel != 4; ++channel) {
    Image image(2, 2);
    cSet.Load(image.Data(), aocommon::Polarization::StokesI, channel, false);
    BOOST_TEST(image[0] == static_cast<float>(channel + 1));
  }
}

BOOST_AUTO_TEST_CASE(parallel_load_facet_reads_consistent_data) {
  const std::string kPrefix = "parallel-load-facet-test";
  constexpr size_t kImageWidth = 8;
  constexpr size_t kImageHeight = 8;
  constexpr double kDlDm = 0.0125;

  aocommon::FitsWriter writer;
  writer.SetImageDimensions(kImageWidth, kImageHeight, kDlDm, kDlDm);

  std::vector<schaapcommon::facets::Coord> coords{
      {0.05, -0.05}, {0.0, -0.05}, {0.0, 0.0}, {0.05, 0.0}};
  Facet::InitializationData facet_data(writer.PixelSizeX(), writer.PixelSizeY(), writer.Width(),
                                       writer.Height());
  facet_data.phase_centre.ra = writer.RA();
  facet_data.phase_centre.dec = writer.Dec();
  facet_data.l_shift = writer.LShift();
  facet_data.m_shift = writer.MShift();
  facet_data.padding = 1.5;

  auto facet = std::make_shared<Facet>(facet_data, std::move(coords));
  Image stored(facet->GetTrimmedBoundingBox().Width(), facet->GetTrimmedBoundingBox().Height(),
               7.0f);

  CachedImageSet cSet;
  cSet.Initialize(writer, 1, 1, 1, kPrefix);
  cSet.StoreFacet(stored, aocommon::Polarization::StokesI, 0, 0, facet, false);

  std::atomic<bool> mismatch{false};
  std::vector<std::thread> threads;
  for (size_t i = 0; i != 4; ++i) {
    threads.emplace_back([&]() {
      Image loaded(stored.Width(), stored.Height(), 0.0f);
      for (size_t iter = 0; iter != 16; ++iter) {
        cSet.LoadFacet(loaded.Data(), aocommon::Polarization::StokesI, 0, 0, facet, false);
        if (!std::equal(loaded.begin(), loaded.end(), stored.begin(), stored.end())) {
          mismatch = true;
        }
      }
    });
  }

  for (std::thread& thread : threads) thread.join();
  BOOST_TEST(mismatch.load() == false);
}

BOOST_AUTO_TEST_CASE(local_writer_single_image_mode) {
  aocommon::FitsWriter writer;
  writer.SetImageDimensions(2, 2);

  CachedImageSet cSet;
  cSet.Initialize(writer, 1, 1, 0, "single-image-mode-test");

  WSCFitsWriter local_writer(writer);
  Image stored(2, 2, 3.5f);
  cSet.Store(stored, local_writer, aocommon::Polarization::StokesI, 0, false);

  Image loaded(2, 2);
  cSet.Load(loaded.Data(), aocommon::Polarization::StokesI, 0, false);
  BOOST_TEST(loaded[0] == 3.5f);
  BOOST_TEST(cSet.Empty() == false);
}

BOOST_AUTO_TEST_CASE(stitch_memory_estimation) {
  const detail::StitchMemoryModel baseline{256,   128,   2,    StitchImageType::kOther,
                                           false, false, false};
  const size_t image_slot_bytes = 256u * 128u * sizeof(float);
  const size_t baseline_bytes = detail::EstimateStitchJobBytes(baseline);

  BOOST_TEST(baseline_bytes == (2u * image_slot_bytes * 5u + 3u) / 4u);

  detail::StitchMemoryModel feathering = baseline;
  feathering.has_feathering = true;
  BOOST_TEST(detail::EstimateStitchJobBytes(feathering) > baseline_bytes);

  detail::StitchMemoryModel matrix = baseline;
  matrix.apply_matrix = true;
  BOOST_TEST(detail::EstimateStitchJobBytes(matrix) > detail::EstimateStitchJobBytes(feathering));

  detail::StitchMemoryModel psf = baseline;
  psf.type = StitchImageType::kPSF;
  psf.save_uv = true;
  BOOST_TEST(detail::EstimateStitchJobBytes(psf) > baseline_bytes);

  detail::StitchMemoryModel facet_pb_model = baseline;
  facet_pb_model.type = StitchImageType::kFacetPbModel;
  facet_pb_model.has_feathering = false;
  BOOST_TEST(detail::EstimateStitchJobBytes(facet_pb_model) <=
             detail::EstimateStitchJobBytes(feathering));

  const size_t job_slots =
      detail::CalculateRequiredStitchMemorySlots(baseline_bytes, image_slot_bytes);
  BOOST_TEST(job_slots >= 2u);
}

BOOST_AUTO_TEST_CASE(stitch_memory_gate_blocks_until_release) {
  detail::MemoryGate gate(2);
  detail::MemoryReservation reservation = gate.Acquire(2);

  std::atomic<bool> acquired{false};
  std::thread waiter([&]() {
    detail::MemoryReservation delayed_reservation = gate.Acquire(1);
    acquired.store(true);
  });

  std::this_thread::sleep_for(std::chrono::milliseconds(100));
  BOOST_TEST(acquired.load() == false);

  reservation.Release();
  waiter.join();
  BOOST_TEST(acquired.load());
  BOOST_TEST(gate.AvailableSlots() == gate.TotalSlots());
}

BOOST_AUTO_TEST_CASE(stitch_memory_gate_allows_small_jobs_to_pass_large_waiter) {
  detail::MemoryGate gate(3);
  detail::MemoryReservation reservation = gate.Acquire(1);

  std::atomic<bool> large_acquired{false};
  std::atomic<bool> small_acquired{false};

  std::thread large_waiter([&]() {
    detail::MemoryReservation delayed_reservation = gate.Acquire(3);
    large_acquired.store(true);
  });
  std::thread small_waiter([&]() {
    detail::MemoryReservation delayed_reservation = gate.Acquire(1);
    small_acquired.store(true);
  });

  for (size_t i = 0; i != 20 && !small_acquired.load(); ++i) {
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
  }
  BOOST_TEST(small_acquired.load());
  BOOST_TEST(large_acquired.load() == false);

  reservation.Release();
  small_waiter.join();
  large_waiter.join();
  BOOST_TEST(large_acquired.load());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace wsclean
