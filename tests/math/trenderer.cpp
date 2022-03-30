#include <boost/test/unit_test.hpp>

#include "../../math/renderer.h"

#include "../../model/model.h"
#include "../../model/powerlawsed.h"

#include <aocommon/image.h>

#include <schaapcommon/fft/restoreimage.h>
#include <schaapcommon/fitters/gaussianfitter.h>

namespace {
constexpr size_t kThreadCount = 1;
constexpr size_t kWidth = 64;
constexpr size_t kHeight = 64;
constexpr long double kPixelSize = 1 /*amin*/ * (M_PI / 180.0 / 60.0);

struct RendererFixture {
  RendererFixture() : restored(kWidth, kHeight, 0.0) {
    imageSettings.ra = 0.0;
    imageSettings.dec = 0.0;
    imageSettings.pixel_scale_l = kPixelSize;
    imageSettings.pixel_scale_m = kPixelSize;
    imageSettings.l_shift = 0.0;
    imageSettings.m_shift = 0.0;
  }

  aocommon::Image restored;
  renderer::ImageCoordinateSettings imageSettings;
};

}  // namespace

BOOST_AUTO_TEST_SUITE(render_and_fit)

BOOST_FIXTURE_TEST_CASE(fit_with_bad_initial_value, RendererFixture) {
  PowerLawSED sed(150.0e6, 1.0);
  ModelComponent component;
  component.SetPosDec(0.0);
  component.SetPosRA(0.0);
  component.SetSED(sed);
  ModelSource source;
  source.AddComponent(component);
  Model model;
  model.AddSource(source);

  const long double beamMaj = 4.0L * kPixelSize;
  const long double beamMin = 4.0L * kPixelSize;
  const long double beamPA = 0.0;
  const long double estimatedBeamPx = 1.0;  // this is on purpose way off

  renderer::RestoreWithEllipticalBeam(
      restored, imageSettings, model, beamMaj, beamMin, beamPA, 100e6, 200e6,
      aocommon::Polarization::StokesI, kThreadCount);

  schaapcommon::fitters::GaussianFitter fitter;
  double fitMajor;
  double fitMinor;
  double fitPA;
  fitter.Fit2DGaussianCentred(restored.Data(), restored.Width(),
                              restored.Height(), estimatedBeamPx, fitMajor,
                              fitMinor, fitPA, 10.0, false);

  BOOST_CHECK_CLOSE_FRACTION(fitMajor, 4.0, 1e-4);
  BOOST_CHECK_CLOSE_FRACTION(fitMinor, 4.0, 1e-4);
}

BOOST_FIXTURE_TEST_CASE(fit_circular, RendererFixture) {
  PowerLawSED sed(150.0e6, 1.0);
  ModelComponent component;
  component.SetPosDec(0.0);
  component.SetPosRA(0.0);
  component.SetSED(sed);
  ModelSource source;
  source.AddComponent(component);
  Model model;
  model.AddSource(source);

  const long double beamMaj = 4.0L * kPixelSize;
  const long double beamMin = 4.0L * kPixelSize;
  const long double beamPA = 0.0;
  const long double estimatedBeamPx = 1.0;  // this is on purpose way off

  renderer::RestoreWithEllipticalBeam(
      restored, imageSettings, model, beamMaj, beamMin, beamPA, 100e6, 200e6,
      aocommon::Polarization::StokesI, kThreadCount);

  schaapcommon::fitters::GaussianFitter fitter;
  double fitMajor = estimatedBeamPx;
  fitter.Fit2DCircularGaussianCentred(restored.Data(), restored.Width(),
                                      restored.Height(), fitMajor);

  BOOST_CHECK_CLOSE_FRACTION(fitMajor, 4.0, 1e-4);
}

BOOST_FIXTURE_TEST_CASE(fit_small_beam, RendererFixture) {
  PowerLawSED sed(150.0e6, 1.0);
  ModelComponent component;
  component.SetPosDec(0.0);
  component.SetPosRA(0.0);
  component.SetSED(sed);
  ModelSource source;
  source.AddComponent(component);
  Model model;
  model.AddSource(source);

  const long double beamMaj = 4.0L * kPixelSize;
  const long double beamMin = 0.5L * kPixelSize;
  const long double beamPA = 0.0;
  const long double estimatedBeamPx = 1.0;  // this is on purpose way off

  renderer::RestoreWithEllipticalBeam(
      restored, imageSettings, model, beamMaj, beamMin, beamPA, 100e6, 200e6,
      aocommon::Polarization::StokesI, kThreadCount);

  schaapcommon::fitters::GaussianFitter fitter;
  double fitMajor = estimatedBeamPx;
  double fitMinor = estimatedBeamPx;
  double fitPA = 0.0;
  fitter.Fit2DGaussianCentred(restored.Data(), restored.Width(),
                              restored.Height(), estimatedBeamPx, fitMajor,
                              fitMinor, fitPA, 10.0, false);

  BOOST_CHECK_CLOSE_FRACTION(fitMinor, 0.5, 1e-4);
}

BOOST_AUTO_TEST_SUITE_END()
