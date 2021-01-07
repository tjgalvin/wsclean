#include "../../math/dijkstrasplitter.h"

#include <boost/test/unit_test.hpp>

#include "../../structures/image.h"

#include <random>

BOOST_AUTO_TEST_SUITE(dijkstra_splitter)

ImageF MakeImage(size_t width, const std::string& str) {
  const size_t height = str.size() / width;
  BOOST_CHECK_EQUAL(width * height, str.size());
  ImageF image(width, height);

  for (size_t y = 0; y != height; ++y) {
    for (size_t x = 0; x != width; ++x) {
      if (str[x + y * width] == 'X')
        image[x + y * width] = 0.1f;
      else
        image[x + y * width] = 10.0f;
    }
  }
  return image;
}

std::string PathStr(const ImageF& pathImage) {
  std::ostringstream str;
  for (size_t y = 0; y != pathImage.Height(); ++y) {
    for (size_t x = 0; x != pathImage.Width(); ++x) {
      if (pathImage[y * pathImage.Width() + x] == 0.0f)
        str << ' ';
      else
        str << 'X';
    }
    str << '\n';
  }
  return str.str();
}

std::string PathStr(const aocommon::UVector<bool>& mask, size_t width) {
  size_t height = mask.size() / width;
  BOOST_CHECK_EQUAL(mask.size(), width * height);
  std::ostringstream str;
  for (size_t y = 0; y != height; ++y) {
    for (size_t x = 0; x != width; ++x) {
      if (mask[y * width + x])
        str << 'X';
      else
        str << ' ';
    }
    str << '\n';
  }
  return str.str();
}

std::string InputColumnStr(const ImageF& pathImage, size_t x) {
  std::ostringstream str;
  for (size_t y = 0; y != pathImage.Height(); ++y) {
    if (pathImage[y * pathImage.Width() + x] == 10.0f)
      str << ' ';
    else
      str << 'X';
  }
  return str.str();
}

std::string InputRowStr(const ImageF& pathImage, size_t y) {
  std::ostringstream str;
  for (size_t x = 0; x != pathImage.Width(); ++x) {
    if (pathImage[y * pathImage.Width() + x] == 10.0f)
      str << ' ';
    else
      str << 'X';
  }
  return str.str();
}

BOOST_AUTO_TEST_CASE(vertical) {
  ImageF image = MakeImage(10,
                           "X         "
                           " X        "
                           "  X       "
                           "   XXX    "
                           "     X    "
                           "         X"
                           "   X      "
                           "    XXXX  "
                           "        X "
                           "      XX  ");
  ImageF output(image.Width(), image.Height(), 0.0f);
  const DijkstraSplitter splitter(image.Width(), image.Height());
  splitter.DivideVertically(image.data(), output.data(), 0, image.Width());

  BOOST_CHECK_EQUAL(PathStr(output),
                    "X         \n"
                    " X        \n"
                    "  X       \n"
                    "   XX     \n"
                    "     X    \n"
                    "    X     \n"
                    "   X      \n"
                    "    XXXX  \n"
                    "        X \n"
                    "       X  \n");
}

BOOST_AUTO_TEST_CASE(vertical_constrained) {
  ImageF input = MakeImage(10,
                           " X  X     "
                           " X        "
                           "  X       "
                           "   XXX    "
                           "     X    "
                           "XX       X"
                           "  XX      "
                           "    XXXX  "
                           "        X "
                           "      XX  ");
  ImageF output(input);
  const DijkstraSplitter splitter(input.Width(), input.Height());
  splitter.DivideVertically(input.data(), output.data(), 2, 8);

  BOOST_CHECK_EQUAL(PathStr(output),
                    "XX  X   XX\n"
                    "XX X    XX\n"
                    "XXX     XX\n"
                    "XX XX   XX\n"
                    "XX   X  XX\n"
                    "XX  X   XX\n"
                    "XX X    XX\n"
                    "XX  X   XX\n"
                    "XX   X  XX\n"
                    "XX    X XX\n");
  // The input shouldn't have changed for the columns that were excluded:
  BOOST_CHECK_EQUAL(InputColumnStr(output, 0), "     X    ");
  BOOST_CHECK_EQUAL(InputColumnStr(output, 1), "XX   X    ");
  BOOST_CHECK_EQUAL(InputColumnStr(output, 8), "        X ");
  BOOST_CHECK_EQUAL(InputColumnStr(output, 9), "     X    ");
}

BOOST_AUTO_TEST_CASE(horizontal) {
  ImageF input = MakeImage(10,
                           "    X     "
                           "          "
                           "  X       "
                           "   XXXXXX "
                           "     X    "
                           " X   X   X"
                           " X    X   "
                           " X     X  "
                           " X      X "
                           "X     XX X");
  ImageF output(input.Width(), input.Height(), 0.0f);
  const DijkstraSplitter splitter(input.Width(), input.Height());
  splitter.DivideHorizontally(input.data(), output.data(), 0, input.Width());

  BOOST_CHECK_EQUAL(PathStr(output),
                    "          \n"
                    "          \n"
                    "          \n"
                    "   XX     \n"
                    "  X  X    \n"
                    " X   X    \n"
                    " X    X   \n"
                    " X     X  \n"
                    " X      X \n"
                    "X        X\n");
}

BOOST_AUTO_TEST_CASE(horizontal_constrained) {
  ImageF image = MakeImage(10,
                           "  XXX     "
                           " XXXXXX   "
                           " X     XXX"
                           "X   XXX   "
                           "   XX     "
                           "X        X"
                           "XX        "
                           "  X      X"
                           "   XXXXXX "
                           "    XXXX  ");
  ImageF output(image);
  const DijkstraSplitter splitter(image.Width(), image.Height());
  splitter.DivideHorizontally(image.data(), output.data(), 2, 8);

  BOOST_CHECK_EQUAL(PathStr(output),
                    "XXXXXXXXXX\n"
                    "XXXXXXXXXX\n"
                    " X     XXX\n"
                    "X X XXX   \n"
                    "   X      \n"
                    "          \n"
                    "          \n"
                    "          \n"
                    "XXXXXXXXXX\n"
                    "XXXXXXXXXX\n");
  // The input shouldn't have changed for the rows that were excluded:
  BOOST_CHECK_EQUAL(InputRowStr(output, 0), "  XXX     ");
  BOOST_CHECK_EQUAL(InputRowStr(output, 1), " XXXXXX   ");
  BOOST_CHECK_EQUAL(InputRowStr(output, 8), "   XXXXXX ");
  BOOST_CHECK_EQUAL(InputRowStr(output, 9), "    XXXX  ");
}

BOOST_AUTO_TEST_CASE(flood_vertical_area) {
  const size_t width = 9, height = 9;
  ImageF image = MakeImage(width,
                           "   X     "
                           "    X    "
                           "    X    "
                           "   X     "
                           "  X      "
                           "   XXX   "
                           "      X  "
                           "      X  "
                           "      X  ");
  DijkstraSplitter splitter(width, height);
  ImageF scratch(image), dividingLines(width, height, 0.0f);
  splitter.AddVerticalDivider(image.data(), scratch.data(),
                              dividingLines.data(), 2, 7);

  BOOST_CHECK_EQUAL(PathStr(dividingLines),
                    "   X     \n"
                    "    X    \n"
                    "    X    \n"
                    "   X     \n"
                    "  X      \n"
                    "   XXX   \n"
                    "      X  \n"
                    "      X  \n"
                    "      X  \n");

  aocommon::UVector<bool> mask(width * height);
  size_t subX, subWidth;
  splitter.FloodVerticalArea(dividingLines.data(), 1, mask.data(), subX,
                             subWidth);
  BOOST_CHECK_EQUAL(subX, 0u);
  BOOST_CHECK_EQUAL(subWidth, 6u);
  BOOST_CHECK_EQUAL(PathStr(mask, width),
                    "XXX      \n"
                    "XXXX     \n"
                    "XXXX     \n"
                    "XXX      \n"
                    "XX       \n"
                    "XXX      \n"
                    "XXXXXX   \n"
                    "XXXXXX   \n"
                    "XXXXXX   \n");

  splitter.FloodVerticalArea(dividingLines.data(), 7, mask.data(), subX,
                             subWidth);
  BOOST_CHECK_EQUAL(subX, 2u);
  BOOST_CHECK_EQUAL(subWidth, 7u);
  BOOST_CHECK_EQUAL(PathStr(mask, width),
                    "   XXXXXX\n"
                    "    XXXXX\n"
                    "    XXXXX\n"
                    "   XXXXXX\n"
                    "  XXXXXXX\n"
                    "   XXXXXX\n"
                    "      XXX\n"
                    "      XXX\n"
                    "      XXX\n");
}

BOOST_AUTO_TEST_CASE(flood_horizontal_area) {
  const size_t width = 9, height = 9;
  ImageF image = MakeImage(width,
                           "         "
                           "         "
                           "  XX    X"
                           " X  X  X "
                           " X   X X "
                           " X   X X "
                           "X     X  "
                           "         "
                           "         ");
  DijkstraSplitter splitter(width, height);
  ImageF scratch(image), dividingLines(width, height, 0.0f);
  splitter.AddHorizontalDivider(image.data(), scratch.data(),
                                dividingLines.data(), 2, 7);

  BOOST_CHECK_EQUAL(PathStr(dividingLines),
                    "         \n"
                    "         \n"
                    "  XX    X\n"
                    " X  X  X \n"
                    " X   X X \n"
                    " X   X X \n"
                    "X     X  \n"
                    "         \n"
                    "         \n");

  aocommon::UVector<bool> mask(width * height);
  size_t subY, subHeight;
  splitter.FloodHorizontalArea(dividingLines.data(), 1, mask.data(), subY,
                               subHeight);
  BOOST_CHECK_EQUAL(subY, 0u);
  BOOST_CHECK_EQUAL(subHeight, 6u);
  BOOST_CHECK_EQUAL(PathStr(mask, width),
                    "XXXXXXXXX\n"
                    "XXXXXXXXX\n"
                    "XX  XXXX \n"
                    "X    XX  \n"
                    "X     X  \n"
                    "X     X  \n"
                    "         \n"
                    "         \n"
                    "         \n");

  splitter.FloodHorizontalArea(dividingLines.data(), 7, mask.data(), subY,
                               subHeight);
  BOOST_CHECK_EQUAL(subY, 2u);
  BOOST_CHECK_EQUAL(subHeight, 7u);
  BOOST_CHECK_EQUAL(PathStr(mask, width),
                    "         \n"
                    "         \n"
                    "  XX    X\n"
                    " XXXX  XX\n"
                    " XXXXX XX\n"
                    " XXXXX XX\n"
                    "XXXXXXXXX\n"
                    "XXXXXXXXX\n"
                    "XXXXXXXXX\n");
}

BOOST_AUTO_TEST_CASE(get_bounding_mask) {
  const size_t width = 9, height = 9;
  ImageF image = MakeImage(width,
                           "    X    "
                           "    X    "
                           "    X    "
                           "    X    "
                           "XXXXXXXXX"
                           "    X    "
                           "    X    "
                           "    X    "
                           "    X    ");
  DijkstraSplitter splitter(width, height);
  ImageF dividingLines(width, height, 0.0f);
  splitter.DivideVertically(image.data(), dividingLines.data(), 3, 6);

  aocommon::UVector<bool> mask(width * height);
  size_t subXL, subXR, subY, subWidthL, subWidthR, subHeight;

  splitter.FloodVerticalArea(dividingLines.data(), 1, mask.data(), subXL,
                             subWidthL);
  BOOST_CHECK_EQUAL(subXL, 0u);
  BOOST_CHECK_EQUAL(subWidthL, 4u);
  aocommon::UVector<bool> maskL(subWidthL * height);
  Image::TrimBox(maskL.data(), subXL, 0, subWidthL, height, mask.data(), width,
                 height);
  BOOST_CHECK_EQUAL(PathStr(maskL, subWidthL),
                    "XXXX\n"
                    "XXXX\n"
                    "XXXX\n"
                    "XXXX\n"
                    "XXX \n"
                    "XXXX\n"
                    "XXXX\n"
                    "XXXX\n"
                    "XXXX\n");

  splitter.FloodVerticalArea(dividingLines.data(), 7, mask.data(), subXR,
                             subWidthR);
  BOOST_CHECK_EQUAL(subXR, 3u);
  BOOST_CHECK_EQUAL(subWidthR, 6u);
  aocommon::UVector<bool> maskR(subWidthR * height);
  Image::TrimBox(maskR.data(), subXR, 0, subWidthR, height, mask.data(), width,
                 height);
  BOOST_CHECK_EQUAL(PathStr(maskR, subWidthR),
                    " XXXXX\n"
                    " XXXXX\n"
                    " XXXXX\n"
                    " XXXXX\n"
                    "XXXXXX\n"
                    " XXXXX\n"
                    " XXXXX\n"
                    " XXXXX\n"
                    " XXXXX\n");

  dividingLines = 0.0f;
  splitter.DivideHorizontally(image.data(), dividingLines.data(), 3, 6);

  splitter.FloodHorizontalArea(dividingLines.data(), 1, mask.data(), subY,
                               subHeight);
  BOOST_CHECK_EQUAL(subY, 0u);
  BOOST_CHECK_EQUAL(subHeight, 4u);
  BOOST_CHECK_EQUAL(PathStr(mask, width),
                    "XXXXXXXXX\n"
                    "XXXXXXXXX\n"
                    "XXXXXXXXX\n"
                    "XXXX XXXX\n"
                    "         \n"
                    "         \n"
                    "         \n"
                    "         \n"
                    "         \n");

  aocommon::UVector<bool> output(width * height, false);
  size_t maskX, maskY, maskWidth, maskHeight;

  splitter.GetBoundingMask(maskL.data(), subXL, subWidthL, mask.data(),
                           output.data(), maskX, maskY, maskWidth, maskHeight);
  BOOST_CHECK_EQUAL(maskX, 0u);
  BOOST_CHECK_EQUAL(maskY, 0u);
  BOOST_CHECK_EQUAL(maskWidth, 4u);
  BOOST_CHECK_EQUAL(maskHeight, 4u);
  BOOST_CHECK_EQUAL(PathStr(output, width),
                    "XXXX     \n"
                    "XXXX     \n"
                    "XXXX     \n"
                    "XXXX     \n"
                    "         \n"
                    "         \n"
                    "         \n"
                    "         \n"
                    "         \n");

  output.assign(width * height, false);
  splitter.GetBoundingMask(maskR.data(), subXR, subWidthR, mask.data(),
                           output.data(), maskX, maskY, maskWidth, maskHeight);
  BOOST_CHECK_EQUAL(maskX, 4u);
  BOOST_CHECK_EQUAL(maskY, 0u);
  BOOST_CHECK_EQUAL(maskWidth, 5u);
  BOOST_CHECK_EQUAL(maskHeight, 4u);
  BOOST_CHECK_EQUAL(PathStr(output, width),
                    "    XXXXX\n"
                    "    XXXXX\n"
                    "    XXXXX\n"
                    "     XXXX\n"
                    "         \n"
                    "         \n"
                    "         \n"
                    "         \n"
                    "         \n");

  splitter.FloodHorizontalArea(dividingLines.data(), 7, mask.data(), subY,
                               subHeight);
  BOOST_CHECK_EQUAL(subY, 3u);
  BOOST_CHECK_EQUAL(subHeight, 6u);

  output.assign(width * height, false);
  splitter.GetBoundingMask(maskL.data(), subXL, subWidthL, mask.data(),
                           output.data(), maskX, maskY, maskWidth, maskHeight);
  BOOST_CHECK_EQUAL(maskX, 0u);
  BOOST_CHECK_EQUAL(maskY, 4u);
  BOOST_CHECK_EQUAL(maskWidth, 4u);
  BOOST_CHECK_EQUAL(maskHeight, 5u);
  BOOST_CHECK_EQUAL(PathStr(output, width),
                    "         \n"
                    "         \n"
                    "         \n"
                    "         \n"
                    "XXX      \n"
                    "XXXX     \n"
                    "XXXX     \n"
                    "XXXX     \n"
                    "XXXX     \n");

  output.assign(width * height, false);
  splitter.GetBoundingMask(maskR.data(), subXR, subWidthR, mask.data(),
                           output.data(), maskX, maskY, maskWidth, maskHeight);
  BOOST_CHECK_EQUAL(maskX, 3u);
  BOOST_CHECK_EQUAL(maskY, 3u);
  BOOST_CHECK_EQUAL(maskWidth, 6u);
  BOOST_CHECK_EQUAL(maskHeight, 6u);
  BOOST_CHECK_EQUAL(PathStr(output, width),
                    "         \n"
                    "         \n"
                    "         \n"
                    "    X    \n"
                    "   XXXXXX\n"
                    "    XXXXX\n"
                    "    XXXXX\n"
                    "    XXXXX\n"
                    "    XXXXX\n");
}

BOOST_AUTO_TEST_CASE(get_bounding_mask_on_noise) {
  const size_t width = 80, height = 80;
  ImageF image(width, height);
  std::mt19937 rnd;
  std::normal_distribution<float> gaus(0.0f, 1.0f);

  for (size_t repeat = 0; repeat != 1000; ++repeat) {
    for (size_t i = 0; i != width * height; ++i) image[i] = gaus(rnd);

    DijkstraSplitter splitter(width, height);
    ImageF dividingLinesV(width, height, 0.0f),
        dividingLinesH(width, height, 0.0f);
    splitter.DivideVertically(image.data(), dividingLinesV.data(), width / 4,
                              width * 3 / 4);
    splitter.DivideHorizontally(image.data(), dividingLinesH.data(), height / 4,
                                height * 3 / 4);

    aocommon::UVector<bool> maskL(width * height), maskR(width * height),
        maskT(width * height), maskB(width * height),
        mask1(width * height, false), mask2(width * height, false),
        mask3(width * height, false), mask4(width * height, false);
    size_t subX, subY, subWidth, subHeight;
    splitter.FloodVerticalArea(dividingLinesV.data(), width / 8, maskL.data(),
                               subX, subWidth);
    splitter.FloodVerticalArea(dividingLinesV.data(), width * 7 / 8,
                               maskR.data(), subX, subWidth);
    splitter.FloodHorizontalArea(dividingLinesH.data(), width / 8, maskT.data(),
                                 subY, subHeight);
    splitter.FloodHorizontalArea(dividingLinesH.data(), width * 7 / 8,
                                 maskB.data(), subY, subHeight);
    splitter.GetBoundingMask(maskL.data(), 0, width, maskT.data(), mask1.data(),
                             subX, subY, subWidth, subHeight);
    splitter.GetBoundingMask(maskR.data(), 0, width, maskT.data(), mask2.data(),
                             subX, subY, subWidth, subHeight);
    splitter.GetBoundingMask(maskL.data(), 0, width, maskB.data(), mask3.data(),
                             subX, subY, subWidth, subHeight);
    splitter.GetBoundingMask(maskR.data(), 0, width, maskB.data(), mask4.data(),
                             subX, subY, subWidth, subHeight);
    for (size_t i = 0; i != mask1.size(); ++i) {
      size_t n = 0;
      if (mask1[n]) ++n;
      if (mask2[n]) ++n;
      if (mask3[n]) ++n;
      if (mask4[n]) ++n;
      BOOST_CHECK_EQUAL(n, 1u);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
