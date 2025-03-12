#include "../../model/bbsmodel.h"

#include <boost/test/unit_test.hpp>

namespace wsclean {

namespace {
const std::string header =
    "FORMAT = Name, Type, Patch, Ra, Dec, I, Q, U, V, MajorAxis, MinorAxis, "
    "Orientation, ReferenceFrequency='147500000.0', SpectralIndex='[]'\n";
const std::string comments = R"(# LSMTool history:
# 2023-03-23 16:35:31: LOAD (from file '/dir/initial_skymodel.txt')
# 2023-03-23 16:39:14: SETPATCHPOSITIONS (method = 'mid') 
)";
const std::string patches = R"( , , Patch_176, 15:41:28.2192, 34.11.03.3360
 , , Patch_21, 15:31:25.2240, 35.33.40.6080
)";
const std::string gaussian_source =
    "J154128.2+341103, GAUSSIAN, Patch_176, 15:41:28.2192, 34.11.03.336, "
    "4.8242, 0, 0, 0, 65.4, 8.8, 27.6,  , [-0.73]\n";
}  // namespace

BOOST_AUTO_TEST_SUITE(bbs_model)

BOOST_AUTO_TEST_CASE(read_gaussian_source) {
  std::istringstream input(header + comments + patches + gaussian_source);
  const Model model = BBSModel::Read(input);
  BOOST_REQUIRE_EQUAL(model.SourceCount(), 1);
  const ModelSource source = model.Source(0);
  BOOST_CHECK_EQUAL(source.Name(), "J154128.2+341103");
  BOOST_REQUIRE_EQUAL(source.ComponentCount(), 1);
  const ModelComponent component = source.Component(0);
  BOOST_CHECK(component.Type() == ModelComponent::GaussianSource);
  BOOST_CHECK_CLOSE_FRACTION(component.PosRA(), 4.10793922345, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(component.PosDec(), 0.59662788936, 1e-6);
  BOOST_REQUIRE(component.HasPowerLawSED());
  const PowerLawSED& sed = static_cast<const PowerLawSED&>(component.SED());
  BOOST_REQUIRE_EQUAL(sed.NTerms(), 2);
  double reference_frequency = 0.0;
  double brightness[4] = {0.0, 0.0, 0.0, 0.0};
  std::vector<double> terms;
  sed.GetData(reference_frequency, brightness, terms);
  BOOST_CHECK_CLOSE_FRACTION(brightness[0], 4.8242, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(brightness[1], 0.0, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(brightness[2], 0.0, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(brightness[3], 0.0, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(terms[0], -0.73, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(reference_frequency, 147500000.0, 1e-6);
}

BOOST_AUTO_TEST_CASE(check_last_keyword) {
  const std::string header =
      "FORMAT = Name, Patch, Ra, Dec, I, Q, U, V, "
      "ReferenceFrequency='147500000.0', SpectralIndex='[]', Type\n";
  const std::string point_source =
      "J154128.2+341103, , 15:41:28.2192, 34.11.03.336, "
      "4.8242, 0, 0, 0, , [-0.73], POINT\n";
  const std::string bad_source =
      "J154128.2+341103, , 15:41:28.2192, 34.11.03.336, "
      "4.8242, 0, 0, 0, , [-0.73], NOT_A_TYPE\n";

  std::istringstream input1(header + point_source);
  const Model model = BBSModel::Read(input1);
  BOOST_REQUIRE_EQUAL(model.SourceCount(), 1);

  std::istringstream input2(header + bad_source);
  BOOST_CHECK_THROW(BBSModel::Read(input1), std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace wsclean
