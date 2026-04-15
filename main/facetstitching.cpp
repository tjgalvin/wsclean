#include "facetstitching.h"

#include "../math/imageoperations.h"
#include "../math/tophatconvolution.h"
#include "../structures/resources.h"

#include <aocommon/logger.h>
#include <aocommon/recursivefor.h>
#include <aocommon/staticfor.h>

#include <schaapcommon/facets/facetimage.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <vector>

using aocommon::Image;
using aocommon::Logger;
using aocommon::Polarization;
using aocommon::PolarizationEnum;

namespace wsclean {

namespace {

constexpr int64_t kBytesPerGiB = 1024ll * 1024ll * 1024ll;

struct StitchSqGroupJob {
  size_t squared_group_index = 0;
  ImagingTable squared_group;
  StitchImageType type = StitchImageType::kOther;
  bool apply_correction = false;
  bool apply_matrix = false;
  size_t estimated_job_bytes = 0;
  size_t required_memory_slots = 1;
};

struct StitchSqGroupWorkspace {
  StitchSqGroupWorkspace(size_t image_width, size_t image_height)
      : full_image(image_width, image_height), facet_image(image_width, image_height, 1) {}

  Image full_image;
  std::unique_ptr<Image> weight_image;
  schaapcommon::facets::FacetImage facet_image;
};

struct StitchExecutionPlan {
  std::vector<StitchSqGroupJob> jobs;
  size_t cpu_parallelism = 1;
  size_t image_slot_bytes = sizeof(float);
  size_t total_memory_slots = 1;
};

std::map<PolarizationEnum, Image> LoadAndApplyFacetCorrection(
    const ImagingTable& polarized_facet_entries, CachedImageSet& image_cache,
    const Settings& settings, const std::vector<OutputChannelInfo>& info_per_channel) {
  std::map<PolarizationEnum, Image> corrected_images;
  std::vector<Image> facet_images;
  facet_images.reserve(polarized_facet_entries.EntryCount());

  for (const ImagingTableEntry& facet_entry : polarized_facet_entries) {
    const schaapcommon::facets::BoundingBox& box = facet_entry.facet->GetTrimmedBoundingBox();
    Image& facet_image = facet_images.emplace_back(box.Width(), box.Height());
    image_cache.LoadFacet(facet_image.Data(), facet_entry.polarization,
                          facet_entry.outputChannelIndex, facet_entry.facetIndex, facet_entry.facet,
                          false);
  }

  const size_t channel_index = polarized_facet_entries.begin()->outputChannelIndex;
  const size_t facet_index = polarized_facet_entries.begin()->facetIndex;
  const double stokes_i_sqrt = std::sqrt(
      info_per_channel[channel_index].averageFacetCorrection[facet_index].GetStokesIValue());
  aocommon::HMC4x4 matrix =
      info_per_channel[channel_index].averageFacetCorrection[facet_index].GetMatrixValue();
  if (matrix.Invert()) {
    if (Polarization::HasFullStokesPolarization(settings.polarizations)) {
      std::array<Image*, 4> images{&facet_images[0], &facet_images[1], &facet_images[2],
                                   &facet_images[3]};
      math::CorrectImagesForMuellerMatrix(matrix * stokes_i_sqrt, images);
    } else if (settings.polarizations == std::set{PolarizationEnum::XX, PolarizationEnum::YY}) {
      std::array<Image*, 2> images{&facet_images[0], &facet_images[1]};
      math::CorrectDualImagesForMuellerMatrix(matrix * stokes_i_sqrt, images);
    } else {
      Logger::Warn << "Warning: Performing facet corrections with incomplete "
                      "Stokes images. This may not\n"
                      "correct for polarization leakage well.\n";

      for (PolarizationEnum p : settings.polarizations) {
        if (!Polarization::IsStokes(p)) {
          throw std::runtime_error(
              "A polarized facet correction was requested that is not "
              "implemented: facet corrections can only be applied on Stokes "
              "images or on the combination of XX and YY.");
        }
      }

      std::array<Image*, 4> images{nullptr, nullptr, nullptr, nullptr};
      auto facet_image_iterator = facet_images.begin();
      for (PolarizationEnum p : settings.polarizations) {
        const size_t pol_index = Polarization::StokesToIndex(p);
        images[pol_index] = &*facet_image_iterator;
        ++facet_image_iterator;
      }

      Image zero_images[4];
      for (size_t i = 0; i != images.size(); ++i) {
        if (!images[i]) {
          zero_images[i] = Image(facet_images.front().Width(), facet_images.front().Height(), 0.0f);
          images[i] = &zero_images[i];
        }
      }
      math::CorrectImagesForMuellerMatrix(matrix * stokes_i_sqrt, images);
    }
  } else {
    Logger::Warn << "Warning: Facet correction matrix is singular: can't be "
                    "inverted.\n";
    for (Image& image : facet_images) {
      image = std::numeric_limits<float>::quiet_NaN();
    }
  }

  for (size_t polarization_index = 0; polarization_index != facet_images.size();
       ++polarization_index) {
    corrected_images.emplace(polarized_facet_entries[polarization_index].polarization,
                             std::move(facet_images[polarization_index]));
  }
  return corrected_images;
}

template <typename ProcessFullPsfFn, typename CreateWriterFn, typename CreateWriterForPolFn>
void StitchSingleGroup(const ImagingTable::Group& facetGroup, size_t imageIndex,
                       CachedImageSet& imageCache, StitchImageType type, Image& fullImage,
                       std::unique_ptr<Image>& weight_image,
                       schaapcommon::facets::FacetImage& facetImage, bool apply_scalar,
                       const std::map<size_t, ImagingTable>& facets_by_index,
                       const Settings& settings,
                       const std::vector<OutputChannelInfo>& info_per_channel,
                       ProcessFullPsfFn&& process_full_psf, CreateWriterFn&& create_writer,
                       CreateWriterForPolFn&& create_writer_for_pol) {
  const size_t feather_size =
      (type == StitchImageType::kFacetPbModel) ? 0 : settings.GetFeatherSize();
  const bool isImaginary = (imageIndex == 1);
  const bool apply_matrix = !facets_by_index.empty() && !isImaginary;
  fullImage = 0.0f;
  if (feather_size != 0) {
    if (weight_image)
      *weight_image = 0.0f;
    else
      weight_image =
          std::make_unique<Image>(settings.trimmedImageWidth, settings.trimmedImageHeight, 0.0f);
  }
  for (const std::shared_ptr<ImagingTableEntry>& facetEntry : facetGroup) {
    facetImage.SetFacet(*facetEntry->facet, true);
    if (apply_matrix) {
      const auto facet_iter = facets_by_index.find(facetEntry->facetIndex);
      if (facet_iter == facets_by_index.end()) {
        throw std::runtime_error("Missing facet correction input for stitch task");
      }
      std::map<PolarizationEnum, Image> corrected_images =
          LoadAndApplyFacetCorrection(facet_iter->second, imageCache, settings, info_per_channel);
      const auto image_iter = corrected_images.find(facetEntry->polarization);
      if (image_iter == corrected_images.end()) {
        throw std::runtime_error("Missing corrected polarization image for stitch task");
      }
      std::copy(image_iter->second.begin(), image_iter->second.end(), facetImage.Data(0));
    } else {
      imageCache.LoadFacet(facetImage.Data(0), facetEntry->polarization,
                           facetEntry->outputChannelIndex, facetEntry->facetIndex,
                           facetEntry->facet, isImaginary);
    }

    if (apply_scalar) {
      const size_t channelIndex = facetEntry->outputChannelIndex;
      const double factor = info_per_channel[channelIndex]
                                .averageFacetCorrection[facetEntry->facetIndex]
                                .GetStokesIValue();
      facetImage *= 1.0 / std::sqrt(factor);
    }

    if (feather_size != 0) {
      aocommon::Image mask = facetImage.MakeMask();
      const schaapcommon::facets::Facet& facet = *facetEntry->facet;
      const bool needs_padding = facet.GetConvolutionBox() != facet.GetTrimmedBoundingBox();
      size_t left = 0;
      size_t top = 0;
      if (needs_padding) {
        const auto clipped_subtract = [](size_t a, size_t b) { return a > b ? a - b : 0; };
        left = clipped_subtract(facet.GetTrimmedBoundingBox().Min().x,
                                facet.GetConvolutionBox().Min().x);
        top = clipped_subtract(facet.GetTrimmedBoundingBox().Min().y,
                               facet.GetConvolutionBox().Min().y);
        const size_t right = clipped_subtract(facet.GetConvolutionBox().Max().x,
                                              facet.GetTrimmedBoundingBox().Max().x);
        const size_t bottom = clipped_subtract(facet.GetConvolutionBox().Max().y,
                                               facet.GetTrimmedBoundingBox().Max().y);
        mask = mask.Pad(left, top, right, bottom);
      }
      tophat_convolution::Convolve(mask, feather_size);
      if (needs_padding) {
        mask = mask.TrimBox(left, top, facet.GetTrimmedBoundingBox().Width(),
                            facet.GetTrimmedBoundingBox().Height());
      }
      facetImage.AddWithMask(fullImage.Data(), weight_image->Data(), mask, 0);
    } else {
      facetImage.AddToImage({fullImage.Data()});
    }
  }
  if (weight_image) {
    for (size_t i = 0; i != fullImage.Size(); ++i) {
      if ((*weight_image)[i] != 0.0)
        fullImage[i] /= (*weight_image)[i];
      else
        fullImage[i] = 0.0;
    }
  }
  if (type == StitchImageType::kDirty) {
    if (settings.isDirtySaved) {
      WSCFitsWriter writer(create_writer(*facetGroup.front(), isImaginary, false));
      Logger::Info << "Writing dirty image...\n";
      writer.WriteImage("dirty.fits", fullImage);
    }
  }
  if (type == StitchImageType::kFacetPbModel) {
    WSCFitsWriter writer(create_writer(*facetGroup.front(), isImaginary, false));
    Logger::Info << "Writing facet-corrected primary-beam model image...\n";
    writer.WriteImage("model-fpb.fits", fullImage);
  }

  if (type == StitchImageType::kPSF) {
    const ImagingTableEntry& entry = *facetGroup.front();
    process_full_psf(fullImage, entry);
  }

  if (type != StitchImageType::kFacetPbModel) {
    const size_t channelIndex = facetGroup.front()->outputChannelIndex;
    const PolarizationEnum polarization = facetGroup.front()->polarization;
    WSCFitsWriter cache_writer =
        create_writer_for_pol(*facetGroup.front(), polarization, false, false);
    imageCache.Store(fullImage, cache_writer, polarization, channelIndex, isImaginary);
  }
}

}  // namespace

void StitchFacets(const ImagingTable& table, CachedImageSet& image_cache, StitchImageType type,
                  const StitchFacetsConfig& config) {
  const Settings& settings = *config.settings;
  const std::vector<OutputChannelInfo>& info_per_channel = *config.info_per_channel;

  if (config.facet_count == 0) return;

  Logger::Info << "Stitching facets onto full image...\n";

  if (type == StitchImageType::kDirty) {
    for (const ImagingTable::Group& facet_group : table.FacetGroups()) {
      config.initialize_model_images(*facet_group.front(), facet_group.front()->polarization,
                                     table.MaxFacetGroupIndex());
    }
  }

  StitchExecutionPlan execution_plan;
  const size_t max_parallel_jobs = std::max<size_t>(1, table.SquaredGroupCount());
  execution_plan.cpu_parallelism =
      std::clamp<size_t>(settings.threadCount, size_t{1}, max_parallel_jobs);
  execution_plan.image_slot_bytes = std::max<size_t>(
      sizeof(float), settings.trimmedImageWidth * settings.trimmedImageHeight * sizeof(float));

  const int64_t available_memory = GetAvailableMemory(settings.memFraction, settings.absMemLimit);
  const int64_t reserve_memory = std::max<int64_t>(kBytesPerGiB, available_memory / 5);
  const int64_t usable_stitch_memory = std::max<int64_t>(
      static_cast<int64_t>(execution_plan.image_slot_bytes), available_memory - reserve_memory);
  execution_plan.total_memory_slots = std::max<size_t>(
      1, static_cast<size_t>(usable_stitch_memory / execution_plan.image_slot_bytes));
  const auto process_full_psf = [&](Image& image, const ImagingTableEntry& entry) {
    config.process_full_psf(image, entry);
  };
  const auto create_writer = [&](const ImagingTableEntry& entry, bool is_imaginary, bool is_model) {
    return config.create_writer(entry, is_imaginary, is_model);
  };
  const auto create_writer_for_pol = [&](const ImagingTableEntry& entry, PolarizationEnum pol,
                                         bool is_imaginary, bool is_model) {
    return config.create_writer_for_pol(entry, pol, is_imaginary, is_model);
  };

  execution_plan.jobs.reserve(table.SquaredGroupCount());
  for (size_t sq_group = 0; sq_group != table.SquaredGroupCount(); ++sq_group) {
    const ImagingTable squared_group = table.GetSquaredGroup(sq_group);
    const bool apply_correction = type != StitchImageType::kPSF &&
                                  type != StitchImageType::kFacetPbModel &&
                                  settings.UseFacetCorrections();
    const bool apply_matrix = apply_correction && settings.polarizations.size() != 1;
    const detail::StitchMemoryModel memory_model{
        settings.trimmedImageWidth,
        settings.trimmedImageHeight,
        settings.polarizations.size(),
        type,
        apply_matrix,
        settings.GetFeatherSize() != 0 && type != StitchImageType::kFacetPbModel,
        type == StitchImageType::kPSF && settings.isUVImageSaved};
    const size_t estimated_job_bytes = detail::EstimateStitchJobBytes(memory_model);
    const size_t required_memory_slots =
        std::min(execution_plan.total_memory_slots,
                 detail::CalculateRequiredStitchMemorySlots(estimated_job_bytes,
                                                            execution_plan.image_slot_bytes));
    execution_plan.jobs.emplace_back(StitchSqGroupJob{sq_group, squared_group, type,
                                                      apply_correction, apply_matrix,
                                                      estimated_job_bytes, required_memory_slots});
  }

  detail::MemoryGate memory_gate(execution_plan.total_memory_slots);
  std::optional<aocommon::RecursiveFor> recursive_for;
  if (execution_plan.cpu_parallelism > 1 && aocommon::RecursiveFor::GetInstance() == nullptr) {
    recursive_for.emplace();
  }
  aocommon::StaticFor<size_t> loop;
  loop.ConstrainedRun(
      0, execution_plan.jobs.size(), execution_plan.cpu_parallelism, [&](size_t start, size_t end) {
        for (size_t job_index = start; job_index != end; ++job_index) {
          const StitchSqGroupJob& job = execution_plan.jobs[job_index];
          [[maybe_unused]] detail::MemoryReservation reservation =
              memory_gate.Acquire(job.required_memory_slots);

          StitchSqGroupWorkspace workspace(settings.trimmedImageWidth, settings.trimmedImageHeight);
          std::map<size_t, ImagingTable> facets_by_index;
          if (job.apply_matrix) {
            for (size_t facet_index = 0; facet_index != job.squared_group.FacetCount();
                 ++facet_index) {
              const ImagingTable polarized_entries = job.squared_group.GetFacet(facet_index);
              facets_by_index.emplace(polarized_entries.Front().facetIndex, polarized_entries);
            }
          }

          for (const ImagingTable::Group& facet_group : job.squared_group.FacetGroups()) {
            if (job.type == StitchImageType::kPSF &&
                facet_group.front()->polarization != *settings.polarizations.begin()) {
              continue;
            }
            const size_t image_count = facet_group.front()->imageCount;
            for (size_t image_index = 0; image_index != image_count; ++image_index) {
              StitchSingleGroup(facet_group, image_index, image_cache, job.type,
                                workspace.full_image, workspace.weight_image, workspace.facet_image,
                                job.apply_correction && !job.apply_matrix, facets_by_index,
                                settings, info_per_channel, process_full_psf, create_writer,
                                create_writer_for_pol);
            }
          }
        }
      });
}

}  // namespace wsclean
