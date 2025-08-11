#include "directmsgridder.h"

#include "msgriddermanager.h"

#include "../main/progressbar.h"
#include "../msproviders/msprovider.h"
#include "../msproviders/msreaders/msreader.h"

#include <thread>
#include <vector>

#include <aocommon/threadpool.h>

namespace wsclean {

template <typename num_t>
DirectMSGridder<num_t>::DirectMSGridder(
    const Settings& settings, const Resources& resources,
    MsProviderCollection& ms_provider_collection)
    : MsGridder(settings, ms_provider_collection), _resources(resources) {}

template <typename num_t>
void DirectMSGridder<num_t>::StartInversion() {
  _sqrtLMTable = GetSqrtLMLookupTable();
  ResetVisibilityCounters();

  progress_bar_ =
      std::make_unique<ProgressBar>("Performing direct Fourier transform");

  _inversionLane.resize(_resources.NCpus() * 1024);

  for (size_t t = 0; t != _resources.NCpus(); ++t) {
    _layers.emplace_back(aocommon::ImageBase<num_t>(TrimWidth(), TrimHeight(),
                                                    static_cast<num_t>(0.0)));
  }

  aocommon::ThreadPool& thread_pool = aocommon::ThreadPool::GetInstance();
  thread_pool.SetNThreads(_resources.NCpus() + 1);
  thread_pool.StartParallelExecution([&](size_t thread_index) {
    const size_t layer = thread_index - 1;
    InversionSample sample;
    while (_inversionLane.read(sample)) {
      gridSample(sample, layer);
    }
  });
}

template <typename num_t>
size_t DirectMSGridder<num_t>::GridMeasurementSet(
    const MsProviderCollection::MsData& ms_data) {
  InvertMeasurementSet(ms_data, ms_data.internal_ms_index);
  return 0;
}

template <typename num_t>
void DirectMSGridder<num_t>::FinishInversion() {
  _inversionLane.write_end();
  aocommon::ThreadPool::GetInstance().FinishParallelExecution();

  aocommon::ImageBase<num_t> scratch = std::move(_layers.back());
  _layers.pop_back();

  for (const aocommon::ImageBase<num_t>& layer : _layers) {
    scratch += layer;
  }

  _layers.clear();
  _sqrtLMTable.Reset();

  const size_t width = TrimWidth();
  const size_t height = TrimHeight();

  // Wrap the image correctly and normalize it
  _image = aocommon::Image(TrimWidth(), TrimHeight());
  const double weight_factor = 1.0 / ImageWeight();
  for (size_t y = 0; y != height; ++y) {
    size_t ySrc = (height - y) + height / 2;
    if (ySrc >= height) ySrc -= height;

    for (size_t x = 0; x != width; ++x) {
      size_t xSrc = x + width / 2;
      if (xSrc >= width) xSrc -= width;

      _image[x + y * width] = scratch[xSrc + ySrc * width] * weight_factor;
    }
  }

  progress_bar_.reset();
}

template <typename num_t>
inline void DirectMSGridder<num_t>::gridSample(const InversionSample& sample,
                                               size_t layerIndex) {
  // Contribution of one visibility:
  //  I(l, m) = V(u, v, w) exp (2 pi i (ul + vm + w (sqrt(1 - l^2 - m^2) - 1)))
  //   Since every visibility has a conjugate visibility for (-u, -v, -w), we
  //   can simultaneously add:
  // Ic(l, m) = V^*(u, v, w) exp (-2 pi i (ul + vm + w (sqrt(1 - l^2 - m^2) -
  // 1)))
  //   Adding those together gives one real value:
  //     I+Ic = real(V) 2 cos (2 pi (ul + vm + w (sqrt(1 - l^2 - m^2) - 1))) -
  //            imag(V) 2 sin (2 pi (ul + vm + w (sqrt(1 - l^2 - m^2) - 1)))
  aocommon::ImageBase<num_t>& layer = _layers[layerIndex];
  const std::complex<num_t> val = sample.sample;
  const size_t width = TrimWidth();
  const size_t height = TrimHeight();
  constexpr num_t minTwoPi = num_t(-2.0 * M_PI);
  const num_t u = sample.uInLambda;
  const num_t v = sample.vInLambda;
  const num_t w = sample.wInLambda;

  for (size_t y = 0; y != height; ++y) {
    const size_t yIndex = y * width;

    size_t ySrc = (height - y) + height / 2;
    if (ySrc >= height) ySrc -= height;
    const num_t m =
        num_t(((num_t)ySrc - (height / 2)) * PixelSizeY() + MShift());

    for (size_t x = 0; x != width; ++x) {
      size_t xSrc = x + width / 2;
      if (xSrc >= width) xSrc -= width;
      const num_t l =
          num_t(((width / 2) - (num_t)xSrc) * PixelSizeX() + LShift());

      const size_t index = yIndex + x;
      const num_t angle = minTwoPi * (u * l + v * m + w * _sqrtLMTable[index]);
      layer[index] +=
          val.real() * std::cos(angle) - val.imag() * std::sin(angle);
    }
  }
}

template <typename num_t>
void DirectMSGridder<num_t>::InvertMeasurementSet(
    const MsProviderCollection::MsData& ms_data, size_t ms_index) {
  const size_t n_vis_polarizations = ms_data.ms_provider->NPolarizations();
  const aocommon::MultiBandData& selected_bands(
      ms_data.ms_provider->SelectedBands());

  const size_t max_data_size =
      selected_bands.MaxBandChannels() * n_vis_polarizations;
  const size_t n_parms = NumValuesPerSolution();
  aocommon::UVector<std::complex<float>> model_buffer(max_data_size);
  aocommon::UVector<float> weight_buffer(max_data_size);
  aocommon::UVector<bool> selection_buffer(max_data_size, true);

  InversionRow row_data;
  aocommon::UVector<std::complex<float>> row_visibilities(max_data_size);
  row_data.data = row_visibilities.data();

  size_t rowIndex = 0;
  std::unique_ptr<MSReader> ms_reader = ms_data.ms_provider->MakeReader();
  while (ms_reader->CurrentRowAvailable()) {
    // This is temporarily disabled during BDA changes. Once MsProvider provides
    // a way to get the nr of selected rows, this can be enabled.
    // progress_bar_->SetProgress(ms_index * id_to_ms_row.size() + rowIndex,
    //                            GetMsCount() * id_to_ms_row.size());

    MSProvider::MetaData metadata;
    ms_reader->ReadMeta(metadata);
    row_data.uvw[0] = metadata.u_in_m;
    row_data.uvw[1] = metadata.v_in_m;
    row_data.uvw[2] = metadata.w_in_m;

    if (n_parms == 2) {
      GetCollapsedVisibilities<2>(*ms_reader, ms_data.antenna_names.size(),
                                  row_data, weight_buffer.data(),
                                  model_buffer.data(), selection_buffer.data(),
                                  metadata);
    } else {
      GetCollapsedVisibilities<4>(*ms_reader, ms_data.antenna_names.size(),
                                  row_data, weight_buffer.data(),
                                  model_buffer.data(), selection_buffer.data(),
                                  metadata);
    }
    const aocommon::BandData& band = selected_bands[metadata.data_desc_id];
    InversionSample sample;
    for (size_t channel = 0; channel != band.ChannelCount(); ++channel) {
      const double wavelength = band.ChannelWavelength(channel);
      sample.uInLambda = row_data.uvw[0] / wavelength;
      sample.vInLambda = row_data.uvw[1] / wavelength;
      sample.wInLambda = row_data.uvw[2] / wavelength;
      sample.sample = row_data.data[channel];
      _inversionLane.write(sample);
    }

    ms_reader->NextInputRow();
    ++rowIndex;
  }
}

template <typename num_t>
void DirectMSGridder<num_t>::StartPredict(
    std::vector<aocommon::Image>&& /*image*/) {
  throw std::runtime_error(
      "Prediction not yet implemented for direct FT gridding");
}

template <typename num_t>
size_t DirectMSGridder<num_t>::PredictMeasurementSet(
    const MsProviderCollection::MsData& /*ms_data*/) {
  throw std::runtime_error(
      "Prediction not yet implemented for direct FT gridding");
  return 0;
}

template <typename num_t>
void DirectMSGridder<num_t>::FinishPredict() {
  throw std::runtime_error(
      "Prediction not yet implemented for direct FT gridding");
}

template <typename num_t>
aocommon::ImageBase<num_t> DirectMSGridder<num_t>::GetSqrtLMLookupTable()
    const {
  const size_t width = TrimWidth();
  const size_t height = TrimHeight();
  aocommon::ImageBase<num_t> sqrtLMTable(ImageWidth(), ImageHeight());
  num_t* iter = sqrtLMTable.Data();
  for (size_t y = 0; y != height; ++y) {
    size_t ySrc = (height - y) + height / 2;
    if (ySrc >= height) ySrc -= height;
    num_t m = num_t(((num_t)ySrc - (height / 2)) * PixelSizeY() + MShift());

    for (size_t x = 0; x != width; ++x) {
      size_t xSrc = x + width / 2;
      if (xSrc >= width) xSrc -= width;
      num_t l = num_t(((width / 2) - (num_t)xSrc) * PixelSizeX() + LShift());

      if (l * l + m * m < 1.0)
        *iter = std::sqrt(1.0 - l * l - m * m) - 1.0;
      else
        *iter = 0.0;
      ++iter;
    }
  }
  return sqrtLMTable;
}

template class DirectMSGridder<float>;
template class DirectMSGridder<double>;
template class DirectMSGridder<long double>;

}  // namespace wsclean
