#include "../deconvolution/imageset.h"
#include "../deconvolution/spectralfitter.h"
#include "../io/cachedimageset.h"

#include <aocommon/image.h>
#include <aocommon/polarization.h>

#include <boost/test/unit_test.hpp>

using aocommon::Image;
using aocommon::PolarizationEnum;

struct ImageSetFixtureBase {
  void addToImageSet(ImagingTable& table, size_t index, size_t joinedGroup,
                     size_t outChannel, size_t squaredIndex,
                     PolarizationEnum pol, size_t frequencyMHz,
                     double imageWeight = 1.0) {
    std::unique_ptr<ImagingTableEntry> e(new ImagingTableEntry());
    e->index = index;
    e->joinedGroupIndex = joinedGroup;
    e->outputChannelIndex = outChannel;
    e->squaredDeconvolutionIndex = squaredIndex;
    e->polarization = pol;
    e->lowestFrequency = frequencyMHz;
    e->highestFrequency = frequencyMHz;
    e->bandStartFrequency = frequencyMHz;
    e->bandEndFrequency = frequencyMHz;
    e->imageCount = 1;
    e->imageWeight = imageWeight;
    table.AddEntry(std::move(e));
  }

  void checkLinearValue(size_t index, float value, const ImageSet& dset) {
    Image dest(2, 2, 1.0);
    dset.GetLinearIntegrated(dest);
    BOOST_CHECK_CLOSE_FRACTION(dest[index], value, 1e-6);
  }

  void checkSquaredValue(size_t index, float value, const ImageSet& dset) {
    Image dest(2, 2, 1.0), scratch(2, 2);
    dset.GetSquareIntegrated(dest, scratch);
    BOOST_CHECK_CLOSE_FRACTION(dest[index], value, 1e-6);
  }

  ImagingTable table;
  Settings settings;
};

struct ImageSetFixture : public ImageSetFixtureBase {
  ImageSetFixture() {
    settings.deconvolutionChannelCount = 1;
    settings.squaredJoins = false;
    settings.linkedPolarizations = std::set<PolarizationEnum>();

    addToImageSet(table, 0, 0, 0, 0, aocommon::Polarization::XX, 100);
    addToImageSet(table, 1, 0, 0, 0, aocommon::Polarization::YY, 100);
    addToImageSet(table, 2, 0, 1, 1, aocommon::Polarization::XX, 200);
    addToImageSet(table, 3, 0, 1, 1, aocommon::Polarization::YY, 200);
    table.Update();
  }
};

BOOST_FIXTURE_TEST_SUITE(imageset, ImageSetFixture)

BOOST_AUTO_TEST_CASE(squaredGroupCount) {
  BOOST_CHECK_EQUAL(table.SquaredGroups().size(), 2u);
}

BOOST_AUTO_TEST_CASE(entryCount) { BOOST_CHECK_EQUAL(table.EntryCount(), 4u); }

BOOST_AUTO_TEST_CASE(entriesInGroup) {
  BOOST_CHECK_EQUAL(table.SquaredGroups().front().size(), 2u);
}

BOOST_AUTO_TEST_CASE(psfCount1) {
  Settings settings;
  settings.deconvolutionChannelCount = 1;
  settings.squaredJoins = false;
  settings.linkedPolarizations = std::set<PolarizationEnum>();
  ImageSet dset(&table, settings, 2, 2);
  BOOST_CHECK_EQUAL(dset.PSFCount(), 1u);
}

BOOST_AUTO_TEST_CASE(psfCount2) {
  Settings settings;
  settings.deconvolutionChannelCount = 2;
  settings.squaredJoins = false;
  settings.linkedPolarizations = std::set<PolarizationEnum>();
  ImageSet dset(&table, settings, 2, 2);
  BOOST_CHECK_EQUAL(dset.PSFCount(), 2u);
}

struct AdvImageSetFixture : public ImageSetFixture {
  FitsWriter writer;
  CachedImageSet cSet;
  aocommon::UVector<double> image;

  AdvImageSetFixture() : image(4, 0.0) {
    writer.SetImageDimensions(2, 2);
    cSet.Initialize(writer, 2, 2, 0, "wsctest");
    image[0] = 2.0;
    cSet.Store(image.data(), aocommon::Polarization::XX, 0, false);
    image[0] = -1.0;
    cSet.Store(image.data(), aocommon::Polarization::YY, 0, false);
    image[0] = 20.0;
    cSet.Store(image.data(), aocommon::Polarization::XX, 1, false);
    image[0] = -10.0;
    cSet.Store(image.data(), aocommon::Polarization::YY, 1, false);
  }
};

BOOST_FIXTURE_TEST_CASE(load, AdvImageSetFixture) {
  cSet.Load(image.data(), aocommon::Polarization::XX, 1, false);
  BOOST_CHECK_EQUAL(image[0], 20.0);
}

BOOST_FIXTURE_TEST_CASE(loadAndAverage, AdvImageSetFixture) {
  ImageSet dset(&table, settings, 2, 2);
  dset.LoadAndAverage(cSet);
  BOOST_CHECK_CLOSE_FRACTION(dset[0][0], 0.5 * (2.0 + 20.0), 1e-8);
  BOOST_CHECK_CLOSE_FRACTION(dset[1][0], 0.5 * (-1.0 - 10.0), 1e-8);
}

BOOST_FIXTURE_TEST_CASE(interpolateAndStore, AdvImageSetFixture) {
  settings.deconvolutionChannelCount = 2;
  ImageSet dset(&table, settings, 2, 2);
  SpectralFitter fitter(SpectralFittingMode::NoFitting, 2);
  dset.LoadAndAverage(cSet);
  dset.InterpolateAndStore(cSet, fitter);
  BOOST_CHECK_CLOSE_FRACTION(dset[0][0], 2.0, 1e-8);
  BOOST_CHECK_CLOSE_FRACTION(dset[1][0], -1.0, 1e-8);
}

BOOST_FIXTURE_TEST_CASE(xxNormalization, ImageSetFixtureBase) {
  settings.linkedPolarizations =
      std::set<PolarizationEnum>{aocommon::Polarization::XX};
  addToImageSet(table, 0, 0, 0, 0, aocommon::Polarization::XX, 100);
  table.Update();
  ImageSet dset(&table, settings, 2, 2);
  dset = 0.0;
  dset[0][1] = 5.0;
  checkLinearValue(1, 5.0, dset);
  checkSquaredValue(1, 5.0, dset);
}

BOOST_FIXTURE_TEST_CASE(iNormalization, ImageSetFixtureBase) {
  addToImageSet(table, 0, 0, 0, 0, aocommon::Polarization::StokesI, 100);
  table.Update();
  settings.linkedPolarizations =
      std::set<PolarizationEnum>{aocommon::Polarization::StokesI};
  ImageSet dset(&table, settings, 2, 2);
  dset = 0.0;
  dset[0][2] = 6.0;
  checkLinearValue(2, 6.0, dset);
  checkSquaredValue(2, 6.0, dset);
}

BOOST_FIXTURE_TEST_CASE(i_2channel_Normalization, ImageSetFixtureBase) {
  addToImageSet(table, 0, 0, 0, 0, aocommon::Polarization::StokesI, 100);
  addToImageSet(table, 1, 0, 1, 1, aocommon::Polarization::StokesI, 200);
  table.Update();
  settings.deconvolutionChannelCount = 2;
  settings.linkedPolarizations =
      std::set<PolarizationEnum>{aocommon::Polarization::StokesI};
  std::set<PolarizationEnum> pols{aocommon::Polarization::StokesI};
  ImageSet dset(&table, settings, 2, 2);
  dset = 0.0;
  dset[0][0] = 12.0;
  dset[1][0] = 13.0;
  checkLinearValue(0, 12.5, dset);
  checkSquaredValue(0, 12.5, dset);
}

BOOST_FIXTURE_TEST_CASE(i_2channel_NaNs, ImageSetFixtureBase) {
  addToImageSet(table, 0, 0, 0, 0, aocommon::Polarization::StokesI, 100, 0.0);
  addToImageSet(table, 1, 0, 1, 1, aocommon::Polarization::StokesI, 200, 1.0);
  table.Update();
  settings.deconvolutionChannelCount = 2;
  settings.linkedPolarizations =
      std::set<PolarizationEnum>{aocommon::Polarization::StokesI};
  std::set<PolarizationEnum> pols{aocommon::Polarization::StokesI};
  ImageSet dset(&table, settings, 2, 2);
  dset = 0.0;
  dset[0][0] = std::numeric_limits<float>::quiet_NaN();
  dset[1][0] = 42.0;
  checkLinearValue(0, 42.0f, dset);
  checkSquaredValue(0, 42.0f, dset);
}

BOOST_FIXTURE_TEST_CASE(xxyyNormalization, ImageSetFixtureBase) {
  addToImageSet(table, 0, 0, 0, 0, aocommon::Polarization::XX, 100);
  addToImageSet(table, 1, 0, 0, 0, aocommon::Polarization::YY, 100);
  table.Update();
  settings.linkedPolarizations = std::set<PolarizationEnum>{
      aocommon::Polarization::XX, aocommon::Polarization::YY};
  ImageSet dset(&table, settings, 2, 2);
  dset = 0.0;
  dset[0][3] = 7.0;
  dset[1][3] = 8.0;
  checkLinearValue(3, 7.5, dset);
  dset[0][3] = -7.0;
  checkSquaredValue(3, std::sqrt((7.0 * 7.0 + 8.0 * 8.0) * 0.5), dset);
}

BOOST_FIXTURE_TEST_CASE(iqNormalization, ImageSetFixtureBase) {
  addToImageSet(table, 0, 0, 0, 0, aocommon::Polarization::StokesI, 100);
  addToImageSet(table, 1, 0, 0, 0, aocommon::Polarization::StokesQ, 100);
  table.Update();
  settings.linkedPolarizations = std::set<PolarizationEnum>{
      aocommon::Polarization::StokesI, aocommon::Polarization::StokesQ};
  ImageSet dset(&table, settings, 2, 2);
  dset = 0.0;
  dset[0][0] = 6.0;
  dset[1][0] = -1.0;
  checkLinearValue(0, 5.0, dset);
  checkSquaredValue(0, std::sqrt(6.0 * 6.0 + -1.0 * -1.0), dset);
}

BOOST_FIXTURE_TEST_CASE(linkedINormalization, ImageSetFixtureBase) {
  addToImageSet(table, 0, 0, 0, 0, aocommon::Polarization::StokesI, 100);
  addToImageSet(table, 1, 0, 0, 0, aocommon::Polarization::StokesQ, 100);
  table.Update();
  settings.linkedPolarizations =
      std::set<PolarizationEnum>{aocommon::Polarization::StokesI};
  ImageSet dset(&table, settings, 2, 2);
  dset = 0.0;
  dset[0][0] = 3.0;
  dset[1][0] = -1.0;
  checkLinearValue(0, 3.0, dset);
  checkSquaredValue(0, 3.0, dset);
}

BOOST_FIXTURE_TEST_CASE(iquvNormalization, ImageSetFixtureBase) {
  addToImageSet(table, 0, 0, 0, 0, aocommon::Polarization::StokesI, 100);
  addToImageSet(table, 1, 0, 0, 0, aocommon::Polarization::StokesQ, 100);
  addToImageSet(table, 2, 0, 0, 0, aocommon::Polarization::StokesU, 100);
  addToImageSet(table, 3, 0, 0, 0, aocommon::Polarization::StokesV, 100);
  table.Update();
  settings.linkedPolarizations = std::set<PolarizationEnum>{
      aocommon::Polarization::StokesI, aocommon::Polarization::StokesQ,
      aocommon::Polarization::StokesU, aocommon::Polarization::StokesV};
  ImageSet dset(&table, settings, 2, 2);
  dset = 0.0;
  dset[0][0] = 9.0;
  dset[1][0] = 0.2;
  dset[2][0] = 0.2;
  dset[3][0] = 0.2;
  checkLinearValue(0, 9.6, dset);
  checkSquaredValue(0, std::sqrt(9.0 * 9.0 + 3.0 * 0.2 * 0.2), dset);
}

BOOST_FIXTURE_TEST_CASE(xx_xy_yx_yyNormalization, ImageSetFixtureBase) {
  addToImageSet(table, 0, 0, 0, 0, aocommon::Polarization::XX, 100);
  addToImageSet(table, 1, 0, 0, 0, aocommon::Polarization::XY, 100);
  addToImageSet(table, 2, 0, 0, 0, aocommon::Polarization::YX, 100);
  addToImageSet(table, 3, 0, 0, 0, aocommon::Polarization::YY, 100);
  table.Update();
  settings.linkedPolarizations = std::set<PolarizationEnum>{
      aocommon::Polarization::XX, aocommon::Polarization::XY,
      aocommon::Polarization::YX, aocommon::Polarization::YY};
  ImageSet dset(&table, settings, 2, 2);
  dset = 0.0;
  dset[0][1] = 10.0;
  dset[1][1] = 0.25;
  dset[2][1] = 0.25;
  dset[3][1] = 10.0;
  checkLinearValue(1, 10.25f, dset);
  checkSquaredValue(
      1, std::sqrt((10.0f * 10.0f * 2.0f + 0.25f * 0.25f * 2.0f) * 0.5f), dset);
}

BOOST_FIXTURE_TEST_CASE(xx_xy_yx_yy_2channel_Normalization,
                        ImageSetFixtureBase) {
  addToImageSet(table, 0, 0, 0, 0, aocommon::Polarization::XX, 100);
  addToImageSet(table, 1, 0, 0, 0, aocommon::Polarization::XY, 100);
  addToImageSet(table, 2, 0, 0, 0, aocommon::Polarization::YX, 100);
  addToImageSet(table, 3, 0, 0, 0, aocommon::Polarization::YY, 100);
  addToImageSet(table, 4, 0, 1, 1, aocommon::Polarization::XX, 200);
  addToImageSet(table, 5, 0, 1, 1, aocommon::Polarization::XY, 200);
  addToImageSet(table, 6, 0, 1, 1, aocommon::Polarization::YX, 200);
  addToImageSet(table, 7, 0, 1, 1, aocommon::Polarization::YY, 200);
  table.Update();
  settings.deconvolutionChannelCount = 2;
  settings.linkedPolarizations = std::set<PolarizationEnum>{
      aocommon::Polarization::XX, aocommon::Polarization::XY,
      aocommon::Polarization::YX, aocommon::Polarization::YY};
  ImageSet dset(&table, settings, 2, 2);
  dset = 0.0;
  dset[0][2] = 5.0;
  dset[1][2] = 0.1;
  dset[2][2] = 0.2;
  dset[3][2] = 6.0;
  dset[4][2] = 7.0;
  dset[5][2] = 0.3;
  dset[6][2] = 0.4;
  dset[7][2] = 8.0;
  double sqVal1 = 0.0, sqVal2 = 0.0;
  for (size_t i = 0; i != 4; ++i) {
    sqVal1 += dset[i][2] * dset[i][2];
    sqVal2 += dset[i + 4][2] * dset[i + 4][2];
  }
  checkLinearValue(2, 27.0 * 0.25, dset);
  checkSquaredValue(
      2, (std::sqrt(sqVal1 * 0.5) + std::sqrt(sqVal2 * 0.5)) * 0.5, dset);
}

BOOST_FIXTURE_TEST_CASE(linked_xx_yy_2channel_Normalization,
                        ImageSetFixtureBase) {
  addToImageSet(table, 0, 0, 0, 0, aocommon::Polarization::XX, 100);
  addToImageSet(table, 1, 0, 0, 0, aocommon::Polarization::XY, 100);
  addToImageSet(table, 2, 0, 0, 0, aocommon::Polarization::YX, 100);
  addToImageSet(table, 3, 0, 0, 0, aocommon::Polarization::YY, 100);
  addToImageSet(table, 4, 0, 1, 1, aocommon::Polarization::XX, 200);
  addToImageSet(table, 5, 0, 1, 1, aocommon::Polarization::XY, 200);
  addToImageSet(table, 6, 0, 1, 1, aocommon::Polarization::YX, 200);
  addToImageSet(table, 7, 0, 1, 1, aocommon::Polarization::YY, 200);
  table.Update();
  settings.deconvolutionChannelCount = 2;
  settings.linkedPolarizations = std::set<PolarizationEnum>{
      aocommon::Polarization::XX, aocommon::Polarization::YY};
  ImageSet dset(&table, settings, 2, 2);
  dset = 0.0;
  dset[0][2] = 7.5;
  dset[1][2] = 0.1;
  dset[2][2] = -0.2;
  dset[3][2] = 6.5;
  dset[4][2] = 8.5;
  dset[5][2] = 0.3;
  dset[6][2] = -0.4;
  dset[7][2] = 9.5;
  double sqVal1 = dset[0][2] * dset[0][2] + dset[3][2] * dset[3][2],
         sqVal2 = dset[4][2] * dset[4][2] + dset[7][2] * dset[7][2];
  checkLinearValue(2, 32.0 * 0.25, dset);
  checkSquaredValue(
      2, (std::sqrt(sqVal1 * 0.5) + std::sqrt(sqVal2 * 0.5)) * 0.5, dset);
}

BOOST_FIXTURE_TEST_CASE(linked_xx_2channel_Normalization, ImageSetFixtureBase) {
  addToImageSet(table, 0, 0, 0, 0, aocommon::Polarization::XX, 100);
  addToImageSet(table, 1, 0, 0, 0, aocommon::Polarization::XY, 100);
  addToImageSet(table, 2, 0, 0, 0, aocommon::Polarization::YX, 100);
  addToImageSet(table, 3, 0, 0, 0, aocommon::Polarization::YY, 100);
  addToImageSet(table, 4, 0, 1, 1, aocommon::Polarization::XX, 200);
  addToImageSet(table, 5, 0, 1, 1, aocommon::Polarization::XY, 200);
  addToImageSet(table, 6, 0, 1, 1, aocommon::Polarization::YX, 200);
  addToImageSet(table, 7, 0, 1, 1, aocommon::Polarization::YY, 200);
  table.Update();
  settings.deconvolutionChannelCount = 2;
  settings.linkedPolarizations =
      std::set<PolarizationEnum>{aocommon::Polarization::XX};
  ImageSet dset(&table, settings, 2, 2);
  dset = 0.0;
  dset[0][2] = 7.5;
  dset[1][2] = 0.1;
  dset[2][2] = -0.2;
  dset[3][2] = 6.5;
  dset[4][2] = 8.5;
  dset[5][2] = 0.3;
  dset[6][2] = -0.4;
  dset[7][2] = 9.5;
  double sqVal1 = dset[0][2] * dset[0][2], sqVal2 = dset[4][2] * dset[4][2];
  checkLinearValue(2, 32.0 * 0.25, dset);
  checkSquaredValue(2, (std::sqrt(sqVal1) + std::sqrt(sqVal2)) * 0.5, dset);
}

BOOST_FIXTURE_TEST_CASE(deconvchannels_normalization, ImageSetFixtureBase) {
  addToImageSet(table, 0, 0, 0, 0, aocommon::Polarization::StokesI, 100, 1);
  addToImageSet(table, 1, 0, 1, 1, aocommon::Polarization::StokesI, 200, 1);
  addToImageSet(table, 2, 0, 2, 2, aocommon::Polarization::StokesI, 300, 2);
  addToImageSet(table, 3, 0, 3, 3, aocommon::Polarization::StokesI, 400, 2);
  table.Update();
  settings.deconvolutionChannelCount = 2;
  settings.linkedPolarizations =
      std::set<PolarizationEnum>{aocommon::Polarization::StokesI};
  ImageSet dset(&table, settings, 2, 2);
  dset = 0.0;
  dset[0][0] = 10.0;
  dset[1][0] = 13.0;
  checkLinearValue(0, 12.0, dset);
  checkSquaredValue(0, 12.0, dset);
}

BOOST_FIXTURE_TEST_CASE(deconvchannels_zeroweight, ImageSetFixtureBase) {
  addToImageSet(table, 0, 0, 0, 0, aocommon::Polarization::StokesI, 100, 1);
  addToImageSet(table, 1, 0, 1, 1, aocommon::Polarization::StokesI, 200, 0);
  addToImageSet(table, 2, 0, 2, 2, aocommon::Polarization::StokesI, 300, 2);
  addToImageSet(table, 3, 0, 3, 3, aocommon::Polarization::StokesI, 400, 2);
  table.Update();
  settings.deconvolutionChannelCount = 2;
  settings.linkedPolarizations =
      std::set<PolarizationEnum>{aocommon::Polarization::StokesI};
  ImageSet dset(&table, settings, 2, 2);
  dset = 0.0;
  dset[0][0] = 10.0;
  dset[1][0] = 5.0;
  checkLinearValue(0, 6.0, dset);
  checkSquaredValue(0, 6.0, dset);
}

BOOST_FIXTURE_TEST_CASE(deconvchannels_divisor, ImageSetFixtureBase) {
  for (size_t ch = 0; ch != 16; ++ch)
    addToImageSet(table, ch, 0, ch, ch, aocommon::Polarization::StokesI,
                  100 + ch, 1);
  table.Update();
  settings.deconvolutionChannelCount = 3;
  settings.linkedPolarizations =
      std::set<PolarizationEnum>{aocommon::Polarization::StokesI};
  ImageSet dset(&table, settings, 2, 2);
  dset = 0.0;
  for (size_t ch = 0; ch != 3; ++ch) dset[ch][0] = 7.0;
  checkLinearValue(0, 7.0, dset);
  checkSquaredValue(0, 7.0, dset);

  BOOST_CHECK_EQUAL(dset.PSFIndex(0), 0u);
  BOOST_CHECK_EQUAL(dset.PSFIndex(1), 1u);
  BOOST_CHECK_EQUAL(dset.PSFIndex(2), 2u);
}

BOOST_FIXTURE_TEST_CASE(psfindex, ImageSetFixtureBase) {
  for (size_t ch = 0; ch != 4; ++ch) {
    addToImageSet(table, ch * 4 + 0, 0, ch, ch, aocommon::Polarization::XX, 100,
                  1);
    addToImageSet(table, ch * 4 + 1, 0, ch, ch, aocommon::Polarization::XY, 200,
                  0);
    addToImageSet(table, ch * 4 + 2, 0, ch, ch, aocommon::Polarization::YX, 300,
                  2);
    addToImageSet(table, ch * 4 + 3, 0, ch, ch, aocommon::Polarization::YY, 400,
                  2);
  }
  table.Update();
  settings.deconvolutionChannelCount = 2;
  settings.linkedPolarizations = std::set<PolarizationEnum>{
      aocommon::Polarization::XX, aocommon::Polarization::XY,
      aocommon::Polarization::YX, aocommon::Polarization::YY};
  ImageSet dset(&table, settings, 2, 2);

  BOOST_CHECK_EQUAL(dset.PSFIndex(0), 0u);
  BOOST_CHECK_EQUAL(dset.PSFIndex(1), 0u);
  BOOST_CHECK_EQUAL(dset.PSFIndex(2), 0u);
  BOOST_CHECK_EQUAL(dset.PSFIndex(3), 0u);
  BOOST_CHECK_EQUAL(dset.PSFIndex(4), 1u);
  BOOST_CHECK_EQUAL(dset.PSFIndex(5), 1u);
  BOOST_CHECK_EQUAL(dset.PSFIndex(6), 1u);
  BOOST_CHECK_EQUAL(dset.PSFIndex(7), 1u);
}

BOOST_FIXTURE_TEST_CASE(load_and_average, ImageSetFixtureBase) {
  const size_t nPol = 2;
  const size_t nFreq = 6;
  const PolarizationEnum pols[nPol] = {PolarizationEnum::XX,
                                       PolarizationEnum::YY};
  const size_t width = 7;
  const size_t height = 9;
  FitsWriter writer;
  writer.SetImageDimensions(width, height);
  CachedImageSet cachedImgs;
  const std::vector<double> weights{4.0, 4.0, 0.0, 0.0, 1.0, 1.0};
  cachedImgs.Initialize(writer, 4, 6, 0, "imagesettest");
  Image storedImage(width, height);
  for (size_t ch = 0; ch != nFreq; ++ch) {
    for (size_t p = 0; p != nPol; ++p) {
      size_t index = ch * nPol + p;
      addToImageSet(table, index, 0, ch, ch, pols[p], 100 + ch, weights[ch]);

      storedImage = (1 << index);  // assign the entire image to 2^index
      cachedImgs.Store(storedImage.Data(), pols[p], ch, false);
    }
  }
  settings.linkedPolarizations = std::set<PolarizationEnum>{
      aocommon::Polarization::XX, aocommon::Polarization::YY};
  settings.deconvolutionChannelCount = 2;
  table.Update();

  ImageSet imageSet(&table, settings, width, height);
  imageSet.LoadAndAverage(cachedImgs);
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
