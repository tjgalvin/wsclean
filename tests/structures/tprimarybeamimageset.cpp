
#include "../../structures/primarybeamimageset.h"

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(primary_beam_image_set)

BOOST_AUTO_TEST_CASE(apply_stokes_i) {
  constexpr size_t kWidth = 2;
  constexpr size_t kHeight = 2;
  PrimaryBeamImageSet beams(kWidth, kHeight);
  beams.SetToZero();

  // a complex Hermitian Mueller matrix consists of 16 values:
  BOOST_REQUIRE_EQUAL(beams.NImages(), 16);
  for (size_t i = 0; i != beams.NImages(); ++i) beams[i][0] = i + 3;

  beams[0][1] = 0.01;
  beams[15][1] = 0.01;

  float stokes_i[kWidth * kHeight] = {1.0, 1.0, 1.0, 1.0};
  beams.ApplyStokesI(stokes_i, 0.1);
  const float expected = 2.0 / (beams[0][0] + beams[15][0]);
  BOOST_CHECK_CLOSE_FRACTION(stokes_i[0], expected, 1e-6);

  // According to the beam, the other pixels do not have enough
  // sensitivity and should be set to NaN:
  for (size_t i = 1; i != 4; ++i) BOOST_CHECK(!std::isfinite(stokes_i[i]));
}

BOOST_AUTO_TEST_SUITE_END()
