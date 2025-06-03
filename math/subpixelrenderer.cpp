#include "subpixelrenderer.h"

#include <cmath>
#include <vector>

#include <aocommon/lane.h>
#include <aocommon/logger.h>
#include <aocommon/uvector.h>
#include <aocommon/threadpool.h>

#include <schaapcommon/math/drawgaussian.h>
#include <schaapcommon/fitters/polynomialfitter.h>

#include "../model/bbsmodel.h"
#include "../model/modelsource.h"
#include "../structures/resources.h"
#include "fourierdomainrenderer.h"

using aocommon::Image;

namespace wsclean::math {
namespace {

struct RenderingInfo {
  double central_frequency;
  double start_frequency;
  double end_frequency;
  const aocommon::CoordinateSystem& coordinate_system;
  std::vector<std::vector<Image>> images;
};

aocommon::UVector<float> MakeSincKernel(double value, size_t size) {
  aocommon::UVector<float> kernel(size + 1);
  const int mid = kernel.size() / 2;
  const double fraction = value - std::floor(value);
  for (size_t i = 0; i != kernel.size(); ++i) {
    const double value = (int(i) - mid - fraction) * M_PI;
    kernel[i] = (value == 0) ? 1.0 : std::sin(value) / value;
  }
  return kernel;
}

void MakeWindowedKernel(double value, size_t n,
                        aocommon::UVector<float>& kernel) {
  const int midH = n / 2;
  const float fraction = value - std::floor(value);
  for (size_t i = 0; i != n; ++i) {
    const float xi = (int(i) - midH - fraction) * M_PI;
    // The Hann window is cos^2 (xi/n).
    const float hann_term = std::cos(xi / n);
    kernel[i] = (xi == 0.0) ? 1.0 : hann_term * hann_term * std::sin(xi) / xi;
  }
}

// This function reads the spectral terms (including Stokes I flux) from the
// component model data, and, if required, the terms are converted to those of
// the desired (ordinary) polynomial.
std::vector<float> GetSpectralTerms(const ModelComponent& component,
                                    const RenderingInfo& settings) {
  if (!component.HasPowerLawSED())
    throw std::runtime_error(
        "Spectral term imaging requires the model to specify functional terms");
  const PowerLawSED& sed = static_cast<const PowerLawSED&>(component.SED());
  double component_reference_frequency = 0.0;
  double component_brightness[4] = {0.0, 0.0, 0.0, 0.0};
  std::vector<double> component_terms;
  sed.GetData(component_reference_frequency, component_brightness,
              component_terms);

  // Combine Stokes I flux with the higher-order spectral terms.
  const float stokes_i = component_brightness[0];
  std::vector<float> spectral_terms = {stokes_i};
  spectral_terms.insert(spectral_terms.end(), component_terms.begin(),
                        component_terms.end());

  if (sed.IsLogarithmic()) {
    std::vector<float> polynomial_terms(spectral_terms.size(), 0.0);
    schaapcommon::fitters::PowerLawToPolynomialCoefficients(
        polynomial_terms, spectral_terms, component_reference_frequency,
        settings.central_frequency, settings.start_frequency,
        settings.end_frequency);

    spectral_terms = polynomial_terms;
  } else if (component_reference_frequency != settings.central_frequency) {
    schaapcommon::fitters::ShiftPolynomialReferenceFrequency(
        spectral_terms, component_reference_frequency,
        settings.central_frequency);
  }

  return spectral_terms;
}

// This function waits for sources to be placed in the source lane, and renders
// these to the image. It runs until the lane receives a write_end() call.
void RenderSourceTasks(aocommon::Lane<ModelSource>& source_lane,
                       std::vector<Image>& images, size_t window_size,
                       size_t n_terms, const RenderingInfo& settings) {
  const double pixel_size = settings.coordinate_system.dl;
  FourierDomainRenderer fd_renderer(window_size, pixel_size);
  SubPixelRenderer renderer(window_size);
  ModelSource source;
  while (source_lane.read(source)) {
    for (const ModelComponent& component : source) {
      std::vector<float> spectral_terms = GetSpectralTerms(component, settings);
      if (n_terms < spectral_terms.size()) {
        aocommon::Logger::Warn << "Consider increasing the number of spectral "
                                  "terms set with -draw-spectral-terms (" +
                                      std::to_string(spectral_terms.size()) +
                                      " terms available, but only " +
                                      std::to_string(n_terms) + " requested)";
      }
      for (size_t image_index = 0; image_index < n_terms; ++image_index) {
        const float term = spectral_terms[image_index];
        Image& image = images[image_index];
        if (component.Type() == ModelComponent::GaussianSource) {
          bool is_small_gaussian = component.MinorAxis() < 50.0 * pixel_size;
          if (is_small_gaussian) {
            double l, m;
            aocommon::ImageCoordinates::RaDecToLM<double>(
                component.PosRA(), component.PosDec(),
                settings.coordinate_system.ra, settings.coordinate_system.dec,
                l, m);
            l += settings.coordinate_system.l_shift;
            m += settings.coordinate_system.m_shift;
            float x, y;
            aocommon::ImageCoordinates::LMToXYfloat<float>(
                l, m, settings.coordinate_system.dl,
                settings.coordinate_system.dm, settings.coordinate_system.width,
                settings.coordinate_system.height, x, y);

            fd_renderer.RenderModelComponent(
                image.Data(), component, settings.coordinate_system.width,
                settings.coordinate_system.height, term, x, y);
          } else {
            const aocommon::CoordinateSystem& cs = settings.coordinate_system;
            const schaapcommon::math::Ellipse shape(component.MajorAxis(),
                                                    component.MinorAxis(),
                                                    component.PositionAngle());
            DrawGaussianToLm(image.Data(), cs.width, cs.height, cs.ra, cs.dec,
                             cs.dl, cs.dm, cs.l_shift, cs.m_shift,
                             component.PosRA(), component.PosDec(), shape,
                             term);
          }
        } else {
          double l, m;
          aocommon::ImageCoordinates::RaDecToLM<double>(
              component.PosRA(), component.PosDec(),
              settings.coordinate_system.ra, settings.coordinate_system.dec, l,
              m);
          l += settings.coordinate_system.l_shift;
          m += settings.coordinate_system.m_shift;
          float x, y;
          aocommon::ImageCoordinates::LMToXYfloat<float>(
              l, m, settings.coordinate_system.dl,
              settings.coordinate_system.dm, settings.coordinate_system.width,
              settings.coordinate_system.height, x, y);
          if (window_size)
            renderer.RenderWindowedSource(
                image.Data(), settings.coordinate_system.width,
                settings.coordinate_system.height, term, x, y);
          else
            renderer.RenderSource(
                image.Data(), settings.coordinate_system.width,
                settings.coordinate_system.height, term, x, y);
        }
      }
    }
  }
}

}  // namespace

void SubPixelRenderer::RenderSource(float* image, size_t width, size_t height,
                                    float brightness, double x, double y) {
  const aocommon::UVector<float> x_sinc = MakeSincKernel(x, width);
  const aocommon::UVector<float> y_sinc = MakeSincKernel(y, height);

  const int mid_x = x_sinc.size() / 2;
  const int mid_y = y_sinc.size() / 2;
  const int x_offset = std::floor(x) - mid_x;
  const int y_offset = std::floor(y) - mid_y;
  const size_t start_x = std::max<int>(x_offset, 0);
  const size_t start_y = std::max<int>(y_offset, 0);
  const size_t end_x = std::min<size_t>(x_offset + x_sinc.size(), width);
  const size_t end_y = std::min<size_t>(y_offset + y_sinc.size(), height);

  for (size_t yi = start_y; yi != end_y; ++yi) {
    float* ptr = &image[yi * width];
    const double y_value = brightness * y_sinc[yi - y_offset];
    for (size_t xi = start_x; xi != end_x; ++xi) {
      ptr[xi] += y_value * x_sinc[xi - x_offset];
    }
  }
}

void SubPixelRenderer::RenderWindowedSource(float* image, size_t width,
                                            size_t height, float flux, float x,
                                            float y) {
  const size_t n = x_kernel_.size();

  MakeWindowedKernel(x, n, x_kernel_);
  MakeWindowedKernel(y, n, y_kernel_);

  const int mid_x = n / 2;
  const int mid_y = n / 2;
  const int x_offset = std::floor(x) - mid_x;
  const int y_offset = std::floor(y) - mid_y;
  const size_t start_x = std::max<int>(x_offset, 0);
  const size_t start_y = std::max<int>(y_offset, 0);
  const size_t end_x =
      std::max<int>(std::min<int>(x_offset + int(n), int(width)), int(start_x));
  const size_t end_y = std::max<int>(
      std::min<int>(y_offset + int(n), int(height)), int(start_y));
  for (size_t yi = start_y; yi != end_y; ++yi) {
    float* ptr = &image[yi * width];
    const float y_corrected_flux = flux * y_kernel_[yi - y_offset];
    for (size_t xi = start_x; xi != end_x; ++xi) {
      ptr[xi] += y_corrected_flux * x_kernel_[xi - x_offset];
    }
  }
}

std::vector<aocommon::Image> RenderSubPixelModel(
    const std::string& model_filename,
    const aocommon::CoordinateSystem& coordinate_system, double frequency,
    double bandwidth, size_t window_size, size_t n_terms, double mem_fraction,
    double mem_limit) {
  aocommon::Logger::Info << "Rendering sources...\n";
  aocommon::ThreadPool& pool = aocommon::ThreadPool::GetInstance();
  // Each thread will get their own list of images, to prevent having to
  // synchronize. These lists consist of the Stokes I image and images for each
  // of the higher order spectral terms, as requested.
  std::vector<std::vector<Image>> images;

  // Limit the number of threads based on the available memory.
  const size_t thread_limit =
      std::max<size_t>(1, GetAvailableMemory(0.5 * mem_fraction, mem_limit) /
                              (n_terms * coordinate_system.width *
                               coordinate_system.height * sizeof(float)));
  const size_t n_threads = std::min(pool.NThreads(), thread_limit);
  if (n_threads < pool.NThreads()) {
    aocommon::Logger::Warn << "Rendering sky model with fewer threads ("
                           << n_threads << ") than requested ("
                           << pool.NThreads() << ").";
  }

  // Thread pool is only used for renderering, hence, it's safe to resize.
  pool.SetNThreads(n_threads);

  images.resize(n_threads);
  const RenderingInfo settings(frequency, frequency - bandwidth * 0.5,
                               frequency + bandwidth * 0.5, coordinate_system,
                               images);
  aocommon::Lane<ModelSource> source_lane(n_threads);
  for (size_t i = 0; i != n_threads; ++i) {
    for (size_t t = 0; t != n_terms; ++t) {
      images[i].emplace_back(coordinate_system.width, coordinate_system.height,
                             0.0f);
    }
  }

  pool.StartParallelExecution([&](size_t thread_index) {
    RenderSourceTasks(source_lane, images[thread_index], window_size, n_terms,
                      settings);
  });

  // We use the "streaming" reading function from BBSModel. This allows using
  // very large files without keeping it in memory all the time. The streaming
  // read will, after reading a source, call the provided function. That
  // function will place the source in the queue, which will be rendered by the
  // threads running in the thread pool.
  auto process_function = [&source_lane](const ModelSource& s) {
    source_lane.write(s);
  };

  std::exception_ptr exception;
  try {
    BBSModel::Read(model_filename, process_function);
  } catch (...) {
    // In case of an exception (e.g. file not found), the threads still need to
    // receive a signal to end.
    exception = std::current_exception();
  }

  aocommon::Logger::Info << "Finishing...\n";
  source_lane.write_end();
  pool.FinishParallelExecution();

  if (exception) rethrow_exception(exception);

  // Add all images together
  for (size_t thread_index = 1; thread_index != images.size(); ++thread_index) {
    for (size_t image_index = 0; image_index < n_terms; ++image_index) {
      images[0][image_index] += images[thread_index][image_index];
    }
  }

  std::vector<aocommon::Image> result(std::move(images.front()));
  return result;
}

}  // namespace wsclean::math
