
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include "../../gridding/visibilitymodifier.h"

namespace wsclean {

BOOST_AUTO_TEST_SUITE(visibility_modifier)

namespace {
constexpr size_t kNChannels = 1;
constexpr size_t kDataDescId = 0;
constexpr size_t kNStations = 2;
constexpr std::complex<float> kAntenna1BeamGainX = {0.5, 0.0};
constexpr std::complex<float> kAntenna2BeamGainX = {3.0, 4.0};
constexpr std::complex<float> kAntenna1ParmGainX = {0.25, -0.25};
constexpr std::complex<float> kAntenna2ParmGainX = {5.0, -1.0};
constexpr std::complex<float> kAntenna1BeamGainY = {1.0, 0.0};
constexpr std::complex<float> kAntenna2BeamGainY = {1.0, -1.0};
constexpr std::complex<float> kAntenna1ParmGainY = {3.0, -0.5};
constexpr std::complex<float> kAntenna2ParmGainY = {-1.0, 7.0};
constexpr size_t kAntenna1 = 0;
constexpr size_t kAntenna2 = 1;
constexpr size_t kMsIndex = 0;
constexpr size_t kTimeOffset = 0;
constexpr size_t kNPolarizations = 2;
constexpr float kWeights[kNPolarizations * kNChannels] = {1.0, 1.0};
constexpr float kImageWeights[kNChannels] = {1.0};
const std::vector<aocommon::ChannelInfo> kChannels(kNChannels, aocommon::ChannelInfo(150e6, 1e6));
const aocommon::BandData kBand(kChannels, 150e6);
const aocommon::MultiBandData kBands = []() -> aocommon::MultiBandData {
  aocommon::MultiBandData bands;
  bands.SetBand(kDataDescId, kBand);
  return bands;
}();
}  // namespace

template <size_t NTerms>
struct ModifierFixture {
 public:
  static_assert(NTerms == 2 || NTerms == 4);
  ModifierFixture() {
    const std::vector<aocommon::MC2x2F> beam_response = {
        aocommon::MC2x2F(kAntenna1BeamGainX, 0.0, 0.0,
                         kAntenna1BeamGainY),  // antenna 1
        aocommon::MC2x2F(kAntenna2BeamGainX, 0.0, 0.0,
                         kAntenna2BeamGainY)  // antenna 2
    };
    const std::vector<std::complex<float>> parm_response =
     NTerms == 2 ? std::vector<std::complex<float>>{
        kAntenna1ParmGainX, kAntenna1ParmGainY, // antenna 1
        kAntenna2ParmGainX, kAntenna2ParmGainY // antenna 2
    } : std::vector<std::complex<float>>{
        kAntenna1ParmGainX, 0.0, 0.0, kAntenna1ParmGainY, // antenna 1
        kAntenna2ParmGainX, 0.0, 0.0, kAntenna2ParmGainY // antenna 2
    };
    gain_types.push_back(NTerms == 2 ? schaapcommon::h5parm::GainType::kDiagonalComplex
                                     : schaapcommon::h5parm::GainType::kFullJones);
    modifier.SetH5Parm(dummy_h5parms, dummy_solutions, dummy_solutions, gain_types);
    modifier.InitializeMockResponse(kNStations, kBands, beam_response, parm_response);
  }

  void CheckClose(std::complex<float> a, std::complex<float> b) {
    BOOST_CHECK_CLOSE_FRACTION(a.real(), b.real(), 1e-5);
    BOOST_CHECK_CLOSE_FRACTION(a.imag(), b.imag(), 1e-5);
  }

#ifdef HAVE_EVERYBEAM
  void TestApplyConjugatedDual() {
    constexpr float kXx = -1.0;
    constexpr float kYy = 2.0;
    std::complex<float> data[] = {kXx, kYy};
    modifier.ApplyConjugatedDualForRow<ModifierBehaviour::kApplyAndSum, GainMode::k2VisDiagonal,
                                       NTerms, false>(data, kWeights, kImageWeights, kNChannels,
                                                      kDataDescId, kNStations, kAntenna1, kAntenna2,
                                                      kMsIndex);
    const float a1_parm_norm_x = std::norm(kAntenna1ParmGainX);
    const float a2_parm_norm_x = std::norm(kAntenna2ParmGainX);
    const float a1_beam_norm_x = std::norm(kAntenna1BeamGainX);
    const float a2_beam_norm_x = std::norm(kAntenna2BeamGainX);
    const float x_correction_reference =
        a1_parm_norm_x * a2_parm_norm_x * a1_beam_norm_x * a2_beam_norm_x;
    const float a1_parm_norm_y = std::norm(kAntenna1ParmGainY);
    const float a2_parm_norm_y = std::norm(kAntenna2ParmGainY);
    const float a1_beam_norm_y = std::norm(kAntenna1BeamGainY);
    const float a2_beam_norm_y = std::norm(kAntenna2BeamGainY);
    const float y_correction_reference =
        a1_parm_norm_y * a2_parm_norm_y * a1_beam_norm_y * a2_beam_norm_y;
    const aocommon::HMC4x4 correction = modifier.TotalCorrectionSum().GetMatrixValue();
    BOOST_CHECK_CLOSE_FRACTION(correction[0], x_correction_reference, 1e-5);
    BOOST_CHECK_CLOSE_FRACTION(correction[15], y_correction_reference, 1e-5);
    const float x_b_reference = a1_beam_norm_x * a2_beam_norm_x;
    const float y_b_reference = a1_beam_norm_y * a2_beam_norm_y;
    const aocommon::HMC4x4 beam_correction = modifier.BeamCorrectionSum().GetMatrixValue();
    BOOST_CHECK_CLOSE_FRACTION(beam_correction[0], x_b_reference, 1e-5);
    BOOST_CHECK_CLOSE_FRACTION(beam_correction[15], y_b_reference, 1e-5);
    const std::complex<float> backward_x = std::conj(kAntenna1BeamGainX) * kAntenna2BeamGainX *
                                           std::conj(kAntenna1ParmGainX) * kAntenna2ParmGainX;
    const std::complex<float> backward_y = std::conj(kAntenna1BeamGainY) * kAntenna2BeamGainY *
                                           std::conj(kAntenna1ParmGainY) * kAntenna2ParmGainY;
    CheckClose(data[0], kXx * backward_x);
    CheckClose(data[1], kYy * backward_y);
  }
#endif  // HAVE_EVERYBEAM

  std::vector<schaapcommon::h5parm::H5Parm> dummy_h5parms;
  std::vector<schaapcommon::h5parm::SolTab*> dummy_solutions;
  std::vector<schaapcommon::h5parm::GainType> gain_types;
  VisibilityModifier modifier;
};

BOOST_FIXTURE_TEST_CASE(get_h5parm, ModifierFixture<2>) {
  const auto [parm1, parm2] = modifier.GetParmResponseForRow(kMsIndex, kDataDescId, kNStations,
                                                             kAntenna1, kAntenna2, kTimeOffset);
  CheckClose(parm1[0], kAntenna1ParmGainX);
  CheckClose(parm1[1], kAntenna1ParmGainY);
  CheckClose(parm2[0], kAntenna2ParmGainX);
  CheckClose(parm2[1], kAntenna2ParmGainY);
}

BOOST_FIXTURE_TEST_CASE(apply_h5parm, ModifierFixture<2>) {
  constexpr float kXx = 7.0;
  constexpr float kYy = 8.0;
  std::complex<float> data[] = {kXx, kYy};
  modifier.ApplyParmResponseForRow<GainMode::k2VisDiagonal>(data, kMsIndex, kNChannels, kDataDescId,
                                                            kNStations, kAntenna1, kAntenna2);
  CheckClose(data[0], kXx * kAntenna1ParmGainX * std::conj(kAntenna2ParmGainX));
  CheckClose(data[1], kYy * kAntenna1ParmGainY * std::conj(kAntenna2ParmGainY));
}

BOOST_FIXTURE_TEST_CASE(apply_conjugate_h5parm, ModifierFixture<2>) {
  constexpr float kXx = -1.0;
  constexpr float kYy = 2.0;
  std::complex<float> data[] = {kXx, kYy};
  modifier.ApplyConjugatedParmResponseForRow<ModifierBehaviour::kApplyAndSum,
                                             GainMode::k2VisDiagonal, 2, false>(
      data, kWeights, kImageWeights, kMsIndex, kNChannels, kDataDescId, kNStations, kAntenna1,
      kAntenna2);
  const float reference_x_correction =
      std::norm(kAntenna1ParmGainX) * std::norm(kAntenna2ParmGainX);
  const float reference_y_correction =
      std::norm(kAntenna1ParmGainY) * std::norm(kAntenna2ParmGainY);
  const float reference_xy_correction =
      0.5 * (std::norm(kAntenna1ParmGainX * std::conj(kAntenna2ParmGainY)) +
             std::norm(kAntenna2ParmGainX * std::conj(kAntenna1ParmGainY)));
  const aocommon::HMC4x4 correction = modifier.TotalCorrectionSum().GetMatrixValue();
  BOOST_CHECK_CLOSE_FRACTION(correction[0], reference_x_correction, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(correction[5], reference_xy_correction, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(correction[10], reference_xy_correction, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(correction[15], reference_y_correction, 1e-5);
  CheckClose(data[0], kXx * std::conj(kAntenna1ParmGainX) * kAntenna2ParmGainX);
  CheckClose(data[1], kYy * std::conj(kAntenna1ParmGainY) * kAntenna2ParmGainY);
}

#ifdef HAVE_EVERYBEAM

BOOST_FIXTURE_TEST_CASE(apply_beam, ModifierFixture<2>) {
  constexpr float kXx = 3.14;
  constexpr float kYy = -1.0;
  std::complex<float> data[] = {kXx, kYy};
  modifier.ApplyBeamResponse<GainMode::k2VisDiagonal>(data, kBand, kDataDescId, kAntenna1,
                                                      kAntenna2);
  CheckClose(data[0], kXx * kAntenna1BeamGainX * std::conj(kAntenna2BeamGainX));
  CheckClose(data[1], kYy * kAntenna1BeamGainY * std::conj(kAntenna2BeamGainY));
}

BOOST_AUTO_TEST_CASE(get_beam_response_for_row) {
  const std::vector<aocommon::MC2x2F> beam_response = {
      aocommon::MC2x2F(5.0, 0.0, 0.0, 6.0),  // antenna 1
      aocommon::MC2x2F(7.0, 0.0, 0.0, 8.0)   // antenna 2
  };
  const std::vector<std::complex<float>> gain_response = {
      1.0, 0.0, 0.0, 1.0,  // antenna 1
      1.0, 0.0, 0.0, 1.0   // antenna 2
  };
  VisibilityModifier modifier;
  modifier.InitializeMockResponse(kNStations, kBands, beam_response, gain_response);
  auto [gains1, gains2] = modifier.GetBeamResponseForRow(kDataDescId, kAntenna1, kAntenna2);
  for (size_t ch = 0; ch != kNChannels; ++ch) {
    BOOST_CHECK_CLOSE_FRACTION(gains1[ch].Get(0).real(), 5.0, 1e-5);
    BOOST_CHECK_CLOSE_FRACTION(gains1[ch].Get(3).real(), 6.0, 1e-5);
    BOOST_CHECK_CLOSE_FRACTION(gains2[ch].Get(0).real(), 7.0, 1e-5);
    BOOST_CHECK_CLOSE_FRACTION(gains2[ch].Get(3).real(), 8.0, 1e-5);
  }
}

BOOST_AUTO_TEST_CASE(get_beam_response_for_visibility) {
  const std::vector<aocommon::MC2x2F> beam_response = {
      aocommon::MC2x2F(5.0, 0.0, 0.0, 6.0),  // antenna 1
      aocommon::MC2x2F(7.0, 0.0, 0.0, 8.0)   // antenna 2
  };
  const std::vector<std::complex<float>> gain_response = {
      1.0, 0.0, 0.0, 1.0,  // antenna 1
      1.0, 0.0, 0.0, 1.0   // antenna 2
  };
  VisibilityModifier modifier;
  modifier.InitializeMockResponse(kNStations, kBands, beam_response, gain_response);
  const BeamResponseMap& map = modifier.GetCachedBeamResponse();
  constexpr size_t kChannelIndex = 0;
  const auto [gains1, gains2] =
      modifier.GetBeamResponseForVisibility(map, kDataDescId, kAntenna1, kAntenna2, kChannelIndex);
  BOOST_CHECK_CLOSE_FRACTION(gains1.Get(0).real(), 5.0, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(gains1.Get(3).real(), 6.0, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(gains2.Get(0).real(), 7.0, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(gains2.Get(3).real(), 8.0, 1e-5);
}

BOOST_AUTO_TEST_CASE(get_solution_response_for_row) {
  const std::vector<aocommon::MC2x2F> beam_response = {
      aocommon::MC2x2F(1.0, 0.0, 0.0, 1.0),  // antenna 1
      aocommon::MC2x2F(1.0, 0.0, 0.0, 1.0)   // antenna 2
  };
  const std::vector<std::complex<float>> gain_response = {
      2.0, 0.0, 0.0, 3.0,  // antenna 1
      4.0, 0.0, 0.0, 5.0   // antenna 2
  };
  VisibilityModifier modifier;
  std::vector<schaapcommon::h5parm::GainType> gain_types{
      schaapcommon::h5parm::GainType::kFullJones};
  modifier.SetH5Parm({}, {}, {}, gain_types);
  modifier.InitializeMockResponse(kNStations, kBands, beam_response, gain_response);
  constexpr size_t kNParms = 4;
  auto [gains1, gains2] = modifier.GetParmResponseForRow(kMsIndex, kDataDescId, kNStations,
                                                         kAntenna1, kAntenna2, kTimeOffset);
  for (size_t ch = 0; ch != kNChannels; ++ch) {
    BOOST_CHECK_CLOSE_FRACTION(gains1[ch * kNParms + 0].real(), 2.0, 1e-5);
    BOOST_CHECK_CLOSE_FRACTION(gains1[ch * kNParms + 3].real(), 3.0, 1e-5);
    BOOST_CHECK_CLOSE_FRACTION(gains2[ch * kNParms + 0].real(), 4.0, 1e-5);
    BOOST_CHECK_CLOSE_FRACTION(gains2[ch * kNParms + 3].real(), 5.0, 1e-5);
  }
}

BOOST_AUTO_TEST_CASE(get_solution_response_for_visibility) {
  const std::vector<aocommon::MC2x2F> beam_response = {
      aocommon::MC2x2F(1.0, 0.0, 0.0, 1.0),  // antenna 1
      aocommon::MC2x2F(1.0, 0.0, 0.0, 1.0)   // antenna 2
  };
  const std::vector<std::complex<float>> gain_response = {
      2.0, 0.0, 0.0, 3.0,  // antenna 1
      4.0, 0.0, 0.0, 5.0   // antenna 2
  };
  VisibilityModifier modifier;
  std::vector<schaapcommon::h5parm::GainType> gain_types{
      schaapcommon::h5parm::GainType::kFullJones};
  modifier.SetH5Parm({}, {}, {}, gain_types);
  modifier.InitializeMockResponse(kNStations, kBands, beam_response, gain_response);
  const SolutionsResponseMap& map = modifier.GetCachedParmResponse(kMsIndex);
  constexpr size_t kChannelIndex = 0;
  constexpr size_t kNParms = 4;
  const auto [gains1, gains2] = modifier.GetParmResponseForVisibility<kNParms>(
      map[kDataDescId].data(), kDataDescId, kNStations, kAntenna1, kAntenna2, kChannelIndex,
      kTimeOffset);
  BOOST_CHECK_CLOSE_FRACTION(gains1[0].real(), 2.0, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(gains1[3].real(), 3.0, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(gains2[0].real(), 4.0, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(gains2[3].real(), 5.0, 1e-5);
}

BOOST_AUTO_TEST_CASE(apply_unit_beam) {
  constexpr float kXx = 1.0;
  constexpr float kYy = 1.0;
  const std::vector<aocommon::MC2x2F> beam_response = {
      aocommon::MC2x2F(1.0, 0.0, 0.0, 1.0),  // antenna 1
      aocommon::MC2x2F(1.0, 0.0, 0.0, 1.0)   // antenna 2
  };
  const std::vector<std::complex<float>> gain_response = {
      1.0, 0.0, 0.0, 1.0,  // antenna 1
      1.0, 0.0, 0.0, 1.0   // antenna 2
  };
  VisibilityModifier modifier;
  modifier.InitializeMockResponse(kNStations, kBands, beam_response, gain_response);
  std::complex<float> data[] = {kXx, kYy};
  modifier.ApplyConjugatedBeamResponseForRow<ModifierBehaviour::kApplyAndSum,
                                             GainMode::k2VisDiagonal, false>(
      data, kWeights, kImageWeights, kDataDescId, kBand, kAntenna1, kAntenna2);
  const AverageCorrection correction = modifier.TotalCorrectionSum();
  BOOST_CHECK_CLOSE_FRACTION(correction.GetMatrixValue()[0], 1.0, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(correction.GetMatrixValue()[15], 1.0, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(correction.GetStokesIValue(), 1.0, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(data[0].real(), 1.0, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(data[1].real(), 1.0, 1e-5);
}

BOOST_FIXTURE_TEST_CASE(apply_conjugate_beam, ModifierFixture<2>) {
  constexpr float kXx = 3.14;
  constexpr float kYy = -1.0;
  std::complex<float> data[] = {kXx, kYy};
  modifier.ApplyConjugatedBeamResponseForRow<ModifierBehaviour::kApplyAndSum,
                                             GainMode::k2VisDiagonal, false>(
      data, kWeights, kImageWeights, kDataDescId, kBand, kAntenna1, kAntenna2);
  const float reference_x_correction =
      std::norm(kAntenna1BeamGainX) * std::norm(kAntenna2BeamGainX);
  const float reference_y_correction =
      std::norm(kAntenna1BeamGainY) * std::norm(kAntenna2BeamGainY);
  const aocommon::HMC4x4 correction = modifier.TotalCorrectionSum().GetMatrixValue();
  BOOST_CHECK_CLOSE_FRACTION(correction[0], reference_x_correction, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(correction[15], reference_y_correction, 1e-5);
  CheckClose(data[0], kXx * std::conj(kAntenna1BeamGainX) * kAntenna2BeamGainX);
  CheckClose(data[1], kYy * std::conj(kAntenna1BeamGainY) * kAntenna2BeamGainY);
}

BOOST_FIXTURE_TEST_CASE(apply_dual, ModifierFixture<2>) {
  constexpr float kXx = 1e3;
  constexpr float kYy = 1e-3;
  std::complex<float> data[] = {kXx, kYy};
  modifier.ApplyBeamResponse<GainMode::k2VisDiagonal>(data, kBand, kDataDescId, kAntenna1,
                                                      kAntenna2);
  modifier.ApplyParmResponseForRow<GainMode::k2VisDiagonal>(data, 0, kNChannels, kDataDescId,
                                                            kNStations, kAntenna1, kAntenna2);
  const std::complex<float> forward_x = kAntenna1BeamGainX * std::conj(kAntenna2BeamGainX) *
                                        kAntenna1ParmGainX * std::conj(kAntenna2ParmGainX);
  const std::complex<float> forward_y = kAntenna1BeamGainY * std::conj(kAntenna2BeamGainY) *
                                        kAntenna1ParmGainY * std::conj(kAntenna2ParmGainY);
  CheckClose(data[0], kXx * forward_x);
  CheckClose(data[1], kYy * forward_y);
}

BOOST_FIXTURE_TEST_CASE(apply_conjugated_dual_2_terms, ModifierFixture<2>) {
  TestApplyConjugatedDual();
}

BOOST_FIXTURE_TEST_CASE(apply_conjugated_dual_4_terms, ModifierFixture<4>) {
  TestApplyConjugatedDual();
}

BOOST_FIXTURE_TEST_CASE(apply_conjugated_dual_forward, ModifierFixture<2>) {
  constexpr float kXx = -1.0;
  constexpr float kYy = 2.0;
  std::complex<float> data[] = {kXx, kYy};
  modifier
      .ApplyConjugatedDualForRow<ModifierBehaviour::kApplyAndSum, GainMode::k2VisDiagonal, 2, true>(
          data, kWeights, kImageWeights, kNChannels, kDataDescId, kNStations, kAntenna1, kAntenna2,
          kMsIndex);
  const float a1_parm_norm_x = std::norm(kAntenna1ParmGainX);
  const float a2_parm_norm_x = std::norm(kAntenna2ParmGainX);
  const float a1_beam_norm_x = std::norm(kAntenna1BeamGainX);
  const float a2_beam_norm_x = std::norm(kAntenna2BeamGainX);
  const float x_correction_reference =
      a1_parm_norm_x * a2_parm_norm_x * a1_beam_norm_x * a2_beam_norm_x;
  const float a1_parm_norm_y = std::norm(kAntenna1ParmGainY);
  const float a2_parm_norm_y = std::norm(kAntenna2ParmGainY);
  const float a1_beam_norm_y = std::norm(kAntenna1BeamGainY);
  const float a2_beam_norm_y = std::norm(kAntenna2BeamGainY);
  const float y_correction_reference =
      a1_parm_norm_y * a2_parm_norm_y * a1_beam_norm_y * a2_beam_norm_y;
  const aocommon::HMC4x4 correction = modifier.TotalCorrectionSum().GetMatrixValue();
  BOOST_CHECK_CLOSE_FRACTION(correction[0], x_correction_reference, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(correction[15], y_correction_reference, 1e-5);
  const float x_b_reference = a1_beam_norm_x * a2_beam_norm_x;
  const float y_b_reference = a1_beam_norm_y * a2_beam_norm_y;
  const aocommon::HMC4x4 beam_correction = modifier.BeamCorrectionSum().GetMatrixValue();
  BOOST_CHECK_CLOSE_FRACTION(beam_correction[0], x_b_reference, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(beam_correction[15], y_b_reference, 1e-5);
  const std::complex<float> backward_x = std::conj(kAntenna1BeamGainX) * kAntenna2BeamGainX *
                                         std::conj(kAntenna1ParmGainX) * kAntenna2ParmGainX;
  const std::complex<float> backward_y = std::conj(kAntenna1BeamGainY) * kAntenna2BeamGainY *
                                         std::conj(kAntenna1ParmGainY) * kAntenna2ParmGainY;
  const std::complex<float> forward_x = kAntenna1BeamGainX * std::conj(kAntenna2BeamGainX) *
                                        kAntenna1ParmGainX * std::conj(kAntenna2ParmGainX);
  const std::complex<float> forward_y = kAntenna1BeamGainY * std::conj(kAntenna2BeamGainY) *
                                        kAntenna1ParmGainY * std::conj(kAntenna2ParmGainY);
  CheckClose(data[0], kXx * backward_x * forward_x);
  CheckClose(data[1], kYy * backward_y * forward_y);
}

#endif  // HAVE_EVERYBEAM

BOOST_AUTO_TEST_SUITE_END()

}  // namespace wsclean
