#include "msprovidercollection.h"

#include <vector>

#include <aocommon/logger.h>

#include <schaapcommon/reordering/msselection.h>

#include "msgridder.h"
#include "msgriddermanager.h"

using aocommon::Logger;

namespace wsclean {

namespace {
template <size_t NPolInMSProvider>
inline void CalculateMsLimits(MsProviderCollection::MsData& ms_data, double u,
                              double v, double w, double baseline_in_m,
                              double wavelength, double pixel_size_x,
                              double pixel_size_y, size_t image_width,
                              size_t image_height,
                              const ImageWeights* image_weights) {
  const double half_width = 0.5 * image_width;
  const double half_height = 0.5 * image_height;
  const double x = u * pixel_size_x * image_width;
  const double y = u * pixel_size_y * image_height;
  const double imaging_weight = image_weights->GetWeight(u, v);
  if (imaging_weight != 0.0) {
    if (std::floor(x) > -half_width && std::ceil(x) < half_width &&
        std::floor(y) > -half_height && std::ceil(y) < half_height) {
      ms_data.max_w = std::max(ms_data.max_w, std::fabs(w));
      ms_data.min_w = std::min(ms_data.min_w, std::fabs(w));
      ms_data.max_baseline_uvw =
          std::max(ms_data.max_baseline_uvw, baseline_in_m / wavelength);
      ms_data.max_baseline_in_m =
          std::max(ms_data.max_baseline_in_m, baseline_in_m);
    }
  }
}

std::vector<double> SelectH5parmTimes(MSProvider& ms_provider) {
  std::unique_ptr<MSReader> ms_reader = ms_provider.MakeReader();
  std::vector<double> unique_times;
  while (ms_reader->CurrentRowAvailable()) {
    MSProvider::MetaData meta_data;
    ms_reader->ReadMeta(meta_data);
    if (!unique_times.empty() && meta_data.time < unique_times.back()) {
      throw std::runtime_error(
          "The measurement set is not sorted in time. To apply h5parm "
          "solutions, this is currently required.");
    }
    if (unique_times.empty() || meta_data.time != unique_times.back()) {
      unique_times.emplace_back(meta_data.time);
    }
    ms_reader->NextInputRow();
  }
  return unique_times;
}
}  // namespace

void MsProviderCollection::InitializeMSDataVector(
    const std::vector<MsGridder*>& gridders, double w_limit,
    bool has_solution_data) {
  assert(Count() != 0);

  bool has_cache = false;
  for (MsGridder* facet_gridder : gridders) {
    has_cache = facet_gridder->HasMetaDataCache();
    if (!has_cache) facet_gridder->AllocateMetaDataCache(Count());

    facet_gridder->ResetVisibilityModifierCache(Count());
  }

  ms_limits_.max_baseline = 0.0;
  ms_limits_.max_w = 0.0;
  ms_limits_.min_w = std::numeric_limits<double>::max();
  for (size_t i = 0; i != Count(); ++i) {
    MsData& ms_data = ms_data_vector_[i];
    ms_data.internal_ms_index = i;
    ms_data.original_ms_index = Index(i);
    InitializeMeasurementSet(ms_data, gridders, has_cache, has_solution_data);

    ms_limits_.Calculate(ms_data.ms_provider->SelectedBands(),
                         ms_data.ms_provider->StartTime());
    ms_limits_.max_baseline =
        std::max(ms_limits_.max_baseline, ms_data.max_baseline_uvw);
    ms_limits_.max_w = std::max(ms_limits_.max_w, ms_data.max_w);
    ms_limits_.min_w = std::min(ms_limits_.min_w, ms_data.min_w);
  }

  ms_limits_.Validate();

  if (w_limit != 0.0) {
    ms_limits_.max_w *= (1.0 - w_limit);
    if (ms_limits_.max_w < ms_limits_.min_w)
      ms_limits_.max_w = ms_limits_.min_w;
  }

  for (MsGridder* facet_gridder : gridders) {
    facet_gridder->SetMaxW(ms_limits_.max_w);
    facet_gridder->SetMinW(ms_limits_.min_w);
    facet_gridder->SetMaxBaseline(ms_limits_.max_baseline);
  }
}

void MsProviderCollection::InitializeMS() { ms_data_vector_.resize(Count()); }

void MsProviderCollection::InitializeMeasurementSet(
    MsData& ms_data, const std::vector<MsGridder*>& gridders, bool is_cached,
    bool has_solution_data) {
  MSProvider& ms_provider = MeasurementSet(ms_data.internal_ms_index);

  ms_data.ms_provider = &ms_provider;
  ms_data.antenna_names = ms_provider.GetAntennaNames();

  // wlimits will vary across facets in a facet group, however these limits are
  // "estimates".
  // They should not be too small but can be too large (though this can have
  // performance implications.
  // As a code simplification and performance improvement we select the smallest
  // width and height out of all facets in a facet group to calculate the
  // wlimits rather than doing it individually for each one.
  size_t min_image_width = gridders[0]->ImageWidth();
  size_t min_image_height = gridders[0]->ImageHeight();
  for (const MsGridder* gridder : gridders) {
    min_image_width = std::min(min_image_width, gridder->ImageWidth());
    min_image_height = std::min(min_image_height, gridder->ImageHeight());
  }

  MetaDataCache::Entry& cache_entry =
      gridders[0]->GetMetaDataCacheItem(ms_data.internal_ms_index);

  if (is_cached) {
    ms_data.max_w = cache_entry.max_w;
    ms_data.max_w_with_flags = cache_entry.max_w_with_flags;
    ms_data.min_w = cache_entry.min_w;
    ms_data.max_baseline_uvw = cache_entry.max_baseline_uvw;
    ms_data.max_baseline_in_m = cache_entry.max_baseline_in_m;
    ms_data.integration_time = cache_entry.integration_time;
  } else {
    if (ms_provider.NPolarizations() == 4)
      CalculateMsLimits<4>(ms_data, gridders[0]->PixelSizeX(),
                           gridders[0]->PixelSizeY(), min_image_width,
                           min_image_height, gridders[0]->GetImageWeights());
    else if (ms_provider.NPolarizations() == 2)
      CalculateMsLimits<2>(ms_data, gridders[0]->PixelSizeX(),
                           gridders[0]->PixelSizeY(), min_image_width,
                           min_image_height, gridders[0]->GetImageWeights());
    else
      CalculateMsLimits<1>(ms_data, gridders[0]->PixelSizeX(),
                           gridders[0]->PixelSizeY(), min_image_width,
                           min_image_height, gridders[0]->GetImageWeights());
    cache_entry.max_w = ms_data.max_w;
    cache_entry.max_w_with_flags = ms_data.max_w_with_flags;
    cache_entry.min_w = ms_data.min_w;
    cache_entry.max_baseline_uvw = ms_data.max_baseline_uvw;
    cache_entry.max_baseline_in_m = ms_data.max_baseline_in_m;
    cache_entry.integration_time = ms_data.integration_time;
  }

  if (has_solution_data) {
    ms_data.unique_times =
        std::make_shared<std::vector<double>>(SelectH5parmTimes(ms_provider));
    for (MsGridder* gridder : gridders) {
      gridder->GetVisibilityModifier().SetMSTimes(ms_data.original_ms_index,
                                                  ms_data.unique_times);
    }
  }
}

template <size_t NPolInMSProvider>
void MsProviderCollection::CalculateMsLimits(
    MsData& ms_data, double pixel_size_x, double pixel_size_y,
    size_t image_width, size_t image_height,
    const ImageWeights* image_weights) {
  Logger::Info << "Determining min and max w & theoretical beam size... ";
  Logger::Info.Flush();
  ms_data.max_w = 0.0;
  ms_data.max_w_with_flags = 0.0;
  ms_data.min_w = 1e100;
  ms_data.max_baseline_uvw = 0.0;
  ms_data.max_baseline_in_m = 0.0;
  const aocommon::MultiBandData& selected_bands =
      ms_data.ms_provider->SelectedBands();
  std::vector<float> weightArray(selected_bands.MaxBandChannels() *
                                 NPolInMSProvider);
  double curTimestep = -1, firstTime = -1, lastTime = -1;
  size_t nTimesteps = 0;
  std::unique_ptr<MSReader> msReader = ms_data.ms_provider->MakeReader();
  const double shortest_wavelength = selected_bands.ShortestWavelength();
  const double longest_wavelength = selected_bands.LongestWavelength();
  while (msReader->CurrentRowAvailable()) {
    MSProvider::MetaData metadata;
    msReader->ReadMeta(metadata);

    if (curTimestep != metadata.time) {
      curTimestep = metadata.time;
      ++nTimesteps;
      if (firstTime == -1) firstTime = curTimestep;
      lastTime = curTimestep;
    }

    const aocommon::BandData& band = selected_bands[metadata.data_desc_id];
    const double wHi = std::fabs(metadata.w_in_m / shortest_wavelength);
    const double wLo = std::fabs(metadata.w_in_m / longest_wavelength);
    const double baselineInM = std::sqrt(metadata.u_in_m * metadata.u_in_m +
                                         metadata.v_in_m * metadata.v_in_m +
                                         metadata.w_in_m * metadata.w_in_m);
    if (wHi > ms_data.max_w || wLo < ms_data.min_w ||
        baselineInM / band.SmallestWavelength() > ms_data.max_baseline_uvw) {
      msReader->ReadWeights(weightArray.data());
      const float* weightPtr = weightArray.data();

      for (size_t ch = 0; ch != band.ChannelCount(); ++ch) {
        const double wavelength = band.ChannelWavelength(ch);
        double wInL = metadata.w_in_m / wavelength;
        ms_data.max_w_with_flags =
            std::max(ms_data.max_w_with_flags, fabs(wInL));
        if (*weightPtr != 0.0) {
          double uInL = metadata.u_in_m / wavelength;
          double vInL = metadata.v_in_m / wavelength;
          wsclean::CalculateMsLimits<NPolInMSProvider>(
              ms_data, uInL, vInL, wInL, baselineInM, wavelength, pixel_size_x,
              pixel_size_y, image_width, image_height, image_weights);
        }
        weightPtr += NPolInMSProvider;
      }
    }
    msReader->NextInputRow();
  }

  if (ms_data.min_w == 1e100) {
    ms_data.min_w = 0.0;
    ms_data.max_w_with_flags = 0.0;
    ms_data.max_w = 0.0;
  }

  Logger::Info << "DONE (w=[" << ms_data.min_w << ":" << ms_data.max_w
               << "] lambdas, maxuvw=" << ms_data.max_baseline_uvw
               << " lambda)\n";
  if (ms_data.max_w_with_flags != ms_data.max_w) {
    Logger::Debug << "Discarded data has higher w value of "
                  << ms_data.max_w_with_flags << " lambda.\n";
  }

  if (lastTime == firstTime || nTimesteps < 2)
    ms_data.integration_time = 1;
  else
    ms_data.integration_time = (lastTime - firstTime) / (nTimesteps - 1);
}

template void MsProviderCollection::CalculateMsLimits<1>(
    MsData& ms_data, double pixel_size_x, double pixel_size_y,
    size_t image_width, size_t image_height, const ImageWeights* image_weights);
template void MsProviderCollection::CalculateMsLimits<2>(
    MsData& ms_data, double pixel_size_x, double pixel_size_y,
    size_t image_width, size_t image_height, const ImageWeights* image_weights);
template void MsProviderCollection::CalculateMsLimits<4>(
    MsData& ms_data, double pixel_size_x, double pixel_size_y,
    size_t image_width, size_t image_height, const ImageWeights* image_weights);

MSProvider& MsProviderCollection::MeasurementSet(
    size_t internal_ms_index) const {
  return *std::get<0>(ms_provider_collection_[internal_ms_index]);
}
const MSSelection& MsProviderCollection::Selection(
    size_t internal_ms_index) const {
  return std::get<1>(ms_provider_collection_[internal_ms_index]);
}
/**
 * Maps an internal ms index to its original (command line) index.
 */
size_t MsProviderCollection::Index(size_t internal_ms_index) const {
  return std::get<2>(ms_provider_collection_[internal_ms_index]);
}

size_t MsProviderCollection::Count() const {
  return ms_provider_collection_.size();
}

/**
 * @param index the original command line index of this measurement set. This
 * allows associating the measurement set with the right h5parm solution file.
 */
void MsProviderCollection::Add(std::unique_ptr<MSProvider> ms_provider,
                               const MSSelection& selection, size_t index) {
  ms_provider_collection_.emplace_back(std::move(ms_provider), selection,
                                       index);
}

}  // namespace wsclean
