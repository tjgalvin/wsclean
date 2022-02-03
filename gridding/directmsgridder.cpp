#include "directmsgridder.h"

#include "../main/progressbar.h"

#include "../msproviders/msprovider.h"
#include "../msproviders/msreaders/msreader.h"

#include <thread>
#include <vector>

template <typename num_t>
DirectMSGridder<num_t>::DirectMSGridder(const Settings& settings)
    : MSGridderBase(settings), _nThreads(_settings.threadCount) {}

template <typename num_t>
void DirectMSGridder<num_t>::Invert() {
  initializeSqrtLMLookupTable();
  const size_t width = TrimWidth(), height = TrimHeight();

  std::vector<MSData> msDataVector;
  initializeMSDataVector(msDataVector);
  resetVisibilityCounters();

  ProgressBar progress("Performing direct Fourier transform");

  _inversionLane.resize(_nThreads * 1024);

  std::vector<std::thread> threads;
  threads.reserve(_nThreads);
  for (size_t t = 0; t != _nThreads; ++t) {
    _layers.emplace_back(allocate());
    std::fill(_layers[t], _layers[t] + width * height, num_t(0.0));
    threads.emplace_back(
        [&](size_t threadIndex) { inversionWorker(threadIndex); }, t);
  }

  for (size_t i = 0; i != MeasurementSetCount(); ++i) {
    MSData& msData = msDataVector[i];

    if (Polarization() == aocommon::Polarization::XX) {
      invertMeasurementSet<DDGainMatrix::kXX>(msData, progress, i);
    } else if (Polarization() == aocommon::Polarization::YY) {
      invertMeasurementSet<DDGainMatrix::kYY>(msData, progress, i);
    } else {
      invertMeasurementSet<DDGainMatrix::kTrace>(msData, progress, i);
    }
  }

  _inversionLane.write_end();
  for (std::thread& t : threads) t.join();
  threads.clear();

  num_t* scratch;
  scratch = std::move(_layers.back());
  _layers.pop_back();

  for (const num_t* layer : _layers) {
    for (size_t i = 0; i != width * height; ++i) scratch[i] += layer[i];
  }

  // Wrap the image correctly and normalize it
  _image = aocommon::Image(TrimWidth(), TrimHeight());
  double wFactor = 1.0 / totalWeight();
  for (size_t y = 0; y != height; ++y) {
    size_t ySrc = (height - y) + height / 2;
    if (ySrc >= height) ySrc -= height;

    for (size_t x = 0; x != width; ++x) {
      size_t xSrc = x + width / 2;
      if (xSrc >= width) xSrc -= width;

      _image[x + y * width] = scratch[xSrc + ySrc * width] * wFactor;
    }
  }

  freeImg(scratch);
  for (num_t* layer : _layers) freeImg(layer);
  _layers.clear();
  freeImg(_sqrtLMTable);
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
  num_t* layer = _layers[layerIndex];
  const std::complex<num_t> val = sample.sample;
  const size_t width = TrimWidth(), height = TrimHeight();
  const num_t minTwoPi = num_t(-2.0 * M_PI), u = sample.uInLambda,
              v = sample.vInLambda, w = sample.wInLambda;

  for (size_t y = 0; y != height; ++y) {
    size_t yIndex = y * height;

    size_t ySrc = (height - y) + height / 2;
    if (ySrc >= height) ySrc -= height;
    const num_t m =
        num_t(((num_t)ySrc - (height / 2)) * PixelSizeY() + PhaseCentreDM());

    for (size_t x = 0; x != width; ++x) {
      size_t xSrc = x + width / 2;
      if (xSrc >= width) xSrc -= width;
      const num_t l =
          num_t(((width / 2) - (num_t)xSrc) * PixelSizeX() + PhaseCentreDL());

      size_t index = yIndex + x;
      num_t angle = minTwoPi * (u * l + v * m + w * _sqrtLMTable[index]);
      layer[index] +=
          val.real() * std::cos(angle) - val.imag() * std::sin(angle);
    }
  }
}

template <typename num_t>
void DirectMSGridder<num_t>::inversionWorker(size_t layer) {
  InversionSample sample;
  while (_inversionLane.read(sample)) {
    gridSample(sample, layer);
  }
}

template <typename num_t>
template <DDGainMatrix GainEntry>
void DirectMSGridder<num_t>::invertMeasurementSet(
    const MSGridderBase::MSData& msData, ProgressBar& progress,
    size_t msIndex) {
  StartMeasurementSet(msData, false);
  const aocommon::BandData selectedBand(msData.SelectedBand());
  aocommon::UVector<std::complex<float>> modelBuffer(
      selectedBand.ChannelCount());
  aocommon::UVector<float> weightBuffer(selectedBand.ChannelCount());
  aocommon::UVector<bool> isSelected(selectedBand.ChannelCount(), true);

  InversionRow newItem;
  aocommon::UVector<std::complex<float>> newItemData(
      selectedBand.ChannelCount());
  newItem.data = newItemData.data();

  std::vector<size_t> idToMSRow;
  msData.msProvider->MakeIdToMSRowMapping(idToMSRow);
  size_t rowIndex = 0;
  std::unique_ptr<MSReader> msReader = msData.msProvider->MakeReader();
  while (msReader->CurrentRowAvailable()) {
    progress.SetProgress(msIndex * idToMSRow.size() + rowIndex,
                         MeasurementSetCount() * idToMSRow.size());

    msReader->ReadMeta(newItem.uvw[0], newItem.uvw[1], newItem.uvw[2]);

    readAndWeightVisibilities<1, GainEntry>(
        *msReader, msData.antennaNames, newItem, selectedBand,
        weightBuffer.data(), modelBuffer.data(), isSelected.data());
    InversionSample sample;
    for (size_t ch = 0; ch != selectedBand.ChannelCount(); ++ch) {
      const double wl = selectedBand.ChannelWavelength(ch);
      sample.uInLambda = newItem.uvw[0] / wl;
      sample.vInLambda = newItem.uvw[1] / wl;
      sample.wInLambda = newItem.uvw[2] / wl;
      sample.sample = newItem.data[ch];
      _inversionLane.write(sample);
    }

    msReader->NextInputRow();
    ++rowIndex;
  }
}

template <typename num_t>
void DirectMSGridder<num_t>::Predict(std::vector<aocommon::Image>&& /*image*/) {
  throw std::runtime_error(
      "Prediction not yet implemented for direct FT gridding");
}

template <typename num_t>
void DirectMSGridder<num_t>::initializeSqrtLMLookupTable() {
  const size_t width = TrimWidth(), height = TrimHeight();
  _sqrtLMTable = allocate();
  num_t* iter = _sqrtLMTable;
  for (size_t y = 0; y != height; ++y) {
    size_t ySrc = (height - y) + height / 2;
    if (ySrc >= height) ySrc -= height;
    num_t m =
        num_t(((num_t)ySrc - (height / 2)) * PixelSizeY() + PhaseCentreDM());

    for (size_t x = 0; x != width; ++x) {
      size_t xSrc = x + width / 2;
      if (xSrc >= width) xSrc -= width;
      num_t l =
          num_t(((width / 2) - (num_t)xSrc) * PixelSizeX() + PhaseCentreDL());

      if (l * l + m * m < 1.0)
        *iter = std::sqrt(1.0 - l * l - m * m) - 1.0;
      else
        *iter = 0.0;
      ++iter;
    }
  }
}

template class DirectMSGridder<float>;
template class DirectMSGridder<double>;
template class DirectMSGridder<long double>;

template void DirectMSGridder<float>::invertMeasurementSet<DDGainMatrix::kXX>(
    const MSData& msData, class ProgressBar& progress, size_t msIndex);
template void DirectMSGridder<float>::invertMeasurementSet<DDGainMatrix::kYY>(
    const MSData& msData, class ProgressBar& progress, size_t msIndex);
template void
DirectMSGridder<float>::invertMeasurementSet<DDGainMatrix::kTrace>(
    const MSData& msData, class ProgressBar& progress, size_t msIndex);

template void DirectMSGridder<double>::invertMeasurementSet<DDGainMatrix::kXX>(
    const MSData& msData, class ProgressBar& progress, size_t msIndex);
template void DirectMSGridder<double>::invertMeasurementSet<DDGainMatrix::kYY>(
    const MSData& msData, class ProgressBar& progress, size_t msIndex);
template void
DirectMSGridder<double>::invertMeasurementSet<DDGainMatrix::kTrace>(
    const MSData& msData, class ProgressBar& progress, size_t msIndex);

template void
DirectMSGridder<long double>::invertMeasurementSet<DDGainMatrix::kXX>(
    const MSData& msData, class ProgressBar& progress, size_t msIndex);
template void
DirectMSGridder<long double>::invertMeasurementSet<DDGainMatrix::kYY>(
    const MSData& msData, class ProgressBar& progress, size_t msIndex);
template void
DirectMSGridder<long double>::invertMeasurementSet<DDGainMatrix::kTrace>(
    const MSData& msData, class ProgressBar& progress, size_t msIndex);
