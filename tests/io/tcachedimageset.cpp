#include "../../io/cachedimageaccessor.h"
#include "../../io/cachedimageset.h"
#include <radler/image_set.h>

#include <aocommon/fits/fitswriter.h>
#include <aocommon/image.h>
#include <aocommon/polarization.h>

#include <schaapcommon/fitters/spectralfitter.h>

#include <boost/test/unit_test.hpp>

#include <memory>

using aocommon::FitsWriter;
using aocommon::Image;
using aocommon::PolarizationEnum;
using radler::ImageSet;

namespace {

class DummyImageAccessor : public aocommon::ImageAccessor {
 public:
  DummyImageAccessor() {}
  ~DummyImageAccessor() override {}

  void Load(float*) const override {
    BOOST_FAIL("Unexpected ImageAccessor::Load() call");
  }

  void Store(const float*) override {
    BOOST_FAIL("Unexpected ImageAccessor::Store() call");
  }

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
    table_ = std::make_unique<radler::WorkTable>(
        psf_offsets, n_original_channels, n_deconvolution_channels);
  }

  void addToImageSet(size_t outChannel, PolarizationEnum pol,
                     size_t frequencyMHz, double imageWeight = 1.0) {
    auto e = std::make_unique<radler::WorkTableEntry>();
    e->original_channel_index = outChannel;
    e->polarization = pol;
    e->band_start_frequency = frequencyMHz;
    e->band_end_frequency = frequencyMHz;
    e->image_weight = imageWeight;
    e->model_accessor =
        std::make_unique<CachedImageAccessor>(cSet_, pol, outChannel, false);
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

BOOST_AUTO_TEST_SUITE(imageset)

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

BOOST_FIXTURE_TEST_CASE(loadAndAverage, ImageSetFixture<1>) {
  ImageSet dset(*table_, false, {}, 2, 2);
  dset.LoadAndAverage(false);
  BOOST_CHECK_CLOSE_FRACTION(dset[0][0], 0.5 * (2.0 + 20.0), 1e-8);
  BOOST_CHECK_CLOSE_FRACTION(dset[1][0], 0.5 * (-1.0 - 10.0), 1e-8);
}

BOOST_FIXTURE_TEST_CASE(interpolateAndStore, ImageSetFixture<2>) {
  ImageSet dset(*table_, false, {}, 2, 2);
  schaapcommon::fitters::SpectralFitter fitter(
      schaapcommon::fitters::SpectralFittingMode::kNoFitting, 2);
  dset.LoadAndAverage(false);
  dset.InterpolateAndStoreModel(fitter, 1);
  BOOST_CHECK_CLOSE_FRACTION(dset[0][0], 2.0, 1e-8);
  BOOST_CHECK_CLOSE_FRACTION(dset[1][0], -1.0, 1e-8);
}

BOOST_FIXTURE_TEST_CASE(load_and_average, ImageSetFixtureBase) {
  // Almost equivalent to radler::timageset.cc::load_and_average
  // test, except now the CachedImageSet class is used.
  initTable(6, 2);
  const size_t nPol = 2;
  const PolarizationEnum pols[nPol] = {PolarizationEnum::XX,
                                       PolarizationEnum::YY};
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
  const std::set<PolarizationEnum> kLinkedPolarizations{
      aocommon::Polarization::XX, aocommon::Polarization::YY};

  ImageSet imageSet(*table_, false, kLinkedPolarizations, width, height);
  imageSet.LoadAndAverage(false);
  // The first image has all values set to 2^0, the second image 2^1, etc...
  // The XX polarizations of deconvolution channel 1 consists of
  // images 0, 2 and 4. These have been weighted with 4, 4, 0:
  BOOST_CHECK_CLOSE_FRACTION(imageSet[0 * nPol + 0][0],
                             double(1 * 4 + 4 * 4 + 16 * 0) / 8.0, 1e-6);
  // The YY polarizations consists of images 1, 3 and 5, weights 4, 4, 0:
  BOOST_CHECK_CLOSE_FRACTION(imageSet[0 * nPol + 1][0],
                             double(2 * 4 + 8 * 4 + 32 * 0) / 8.0, 1e-6);
  // The XX polarizations of deconvolution channel 2 consists of images 6, 8 and
  // 10 Weights 0, 1, 1
  BOOST_CHECK_CLOSE_FRACTION(imageSet[1 * nPol + 0][0],
                             double(64 * 0 + 256 * 1 + 1024 * 1) / 2.0, 1e-6);
  // YY: images 7, 9, 10, weights 0, 1, 1
  BOOST_CHECK_CLOSE_FRACTION(imageSet[1 * nPol + 1][0],
                             double(128 * 0 + 512 * 1 + 2048 * 1) / 2.0, 1e-6);

  // The total linear integrated sum should be a complete
  // weighting of all input channels
  Image linearIntegrated(width, height);
  imageSet.GetLinearIntegrated(linearIntegrated);
  BOOST_CHECK_CLOSE_FRACTION(
      linearIntegrated[0],
      double(1 * 4 + 4 * 4 + 16 * 0 + 2 * 4 + 8 * 4 + 32 * 0 + 64 * 0 +
             256 * 1 + 1024 * 1 + 128 * 0 + 512 * 1 + 2048 * 1) /
          20.0,
      1e-6);
}

BOOST_AUTO_TEST_SUITE_END()
