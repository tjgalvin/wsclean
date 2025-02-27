#include "subpixelrenderer.h"

#include <cmath>
#include <vector>

#include <aocommon/lane.h>
#include <aocommon/logger.h>
#include <aocommon/uvector.h>
#include <aocommon/threadpool.h>

#include <schaapcommon/math/drawgaussian.h>

#include "../model/bbsmodel.h"
#include "../model/modelsource.h"

using aocommon::Image;

namespace wsclean::math {
namespace {

struct RenderingInfo {
  double start_frequency;
  double end_frequency;
  const aocommon::CoordinateSystem& coordinate_system;
  std::vector<aocommon::Image> images;
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

// This function waits for sources to be placed in the source lane, and renders
// these to the image. It runs until the lane receives a write_end() call.
void RenderSourceTasks(aocommon::Lane<ModelSource>& source_lane, Image& image,
                       size_t window_size, const RenderingInfo& settings) {
  SubPixelRenderer renderer(window_size);
  ModelSource source;
  while (source_lane.read(source)) {
    for (const ModelComponent& comp : source) {
      const float flux = comp.SED().IntegratedFlux(
          settings.start_frequency, settings.end_frequency,
          aocommon::Polarization::StokesI);
      if (!std::isfinite(flux)) {
        std::cout << "Evaluating the spectrum for source " + source.Name()
                  << " resulted in a non-finite value.\n";
        throw std::runtime_error("Evaluating the spectrum for source " +
                                 source.Name() +
                                 " resulted in a non-finite value");
      }
      if (comp.Type() != ModelComponent::PointSource) {
        // TODO: also (small) Gaussian sources should be sinc-convolved
        const aocommon::CoordinateSystem& cs = settings.coordinate_system;
        const schaapcommon::math::Ellipse shape(
            comp.MajorAxis(), comp.MinorAxis(), comp.PositionAngle());
        DrawGaussianToLm(image.Data(), cs.width, cs.height, cs.ra, cs.dec,
                         cs.dl, cs.dm, cs.l_shift, cs.m_shift, comp.PosRA(),
                         comp.PosDec(), shape, flux);
      } else {
        double l, m;
        float x, y;
        aocommon::ImageCoordinates::RaDecToLM<double>(
            comp.PosRA(), comp.PosDec(), settings.coordinate_system.ra,
            settings.coordinate_system.dec, l, m);
        l += settings.coordinate_system.l_shift;
        m += settings.coordinate_system.m_shift;
        aocommon::ImageCoordinates::LMToXYfloat<float>(
            l, m, settings.coordinate_system.dl, settings.coordinate_system.dm,
            settings.coordinate_system.width, settings.coordinate_system.height,
            x, y);
        if (window_size)
          renderer.RenderWindowedSource(
              image.Data(), settings.coordinate_system.width,
              settings.coordinate_system.height, flux, x, y);
        else
          renderer.RenderSource(image.Data(), settings.coordinate_system.width,
                                settings.coordinate_system.height, flux, x, y);
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

aocommon::Image RenderSubPixelModel(
    const std::string& model_filename,
    const aocommon::CoordinateSystem& coordinate_system, double frequency,
    double bandwidth, size_t window_size) {
  aocommon::Logger::Info << "Rendering sources...\n";
  aocommon::ThreadPool& pool = aocommon::ThreadPool::GetInstance();
  // Each thread will get their own image, to prevent having to synchronize.
  std::vector<Image> images;
  const RenderingInfo settings(frequency - bandwidth * 0.5,
                               frequency + bandwidth * 0.5, coordinate_system,
                               images);
  aocommon::Lane<ModelSource> source_lane(pool.NThreads());
  for (size_t i = 0; i != pool.NThreads(); ++i)
    images.emplace_back(coordinate_system.width, coordinate_system.height,
                        0.0f);

  pool.StartParallelExecution([&](size_t thread_index) {
    RenderSourceTasks(source_lane, images[thread_index], window_size, settings);
  });

  // We use the "streaming" reading function from BBSModel. This allows using
  // very large files without keeping it in memory all the time. The streaming
  // read will, after reading a source, call the provided function. That
  // function will place the source in the queue, which will be rendered by the
  // threads running in the thread pool.
  auto process_function = [&source_lane](const ModelSource& s) {
    source_lane.write(s);
  };
  BBSModel::Read(model_filename, process_function);

  aocommon::Logger::Info << "Finishing...\n";
  source_lane.write_end();
  pool.FinishParallelExecution();

  // Add all images together
  for (size_t image_index = 1; image_index != images.size(); ++image_index) {
    images[0] += images[image_index];
  }

  aocommon::Image result(std::move(images.front()));
  return result;
}

}  // namespace wsclean::math
