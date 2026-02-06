#include "idgmsgridder.h"

#include "../msproviders/msreaders/timestepbufferreader.h"

#include <cmath>
#include <thread>

#include <idg-api.h>

#include <aocommon/coordinatesystem.h>
#include <aocommon/fits/fitsreader.h>
#include <aocommon/logger.h>

#include "../msproviders/msprovider.h"
#include "../msproviders/timestepbuffer.h"

#include "../io/findmwacoefffile.h"
#include "../io/imagefilename.h"
#include "../io/parsetreader.h"

#include "../structures/imagingtable.h"

#include "../main/settings.h"

#include "averagebeam.h"
#include "idgcommon.h"
#include "idgconfiguration.h"

#ifdef HAVE_EVERYBEAM
#include <EveryBeam/aterms/atermconfig.h>
#include <EveryBeam/options.h>
#include <EveryBeam/load.h>

using everybeam::ATermSettings;
using everybeam::aterms::ATermBase;
using everybeam::aterms::ATermBeam;
using everybeam::aterms::ATermConfig;
#endif  // HAVE_EVERYBEAM

using aocommon::CoordinateSystem;
using aocommon::Image;
using aocommon::Logger;

namespace wsclean {

namespace {
constexpr const size_t kGridderIndex = 0;
}

IdgMsGridder::IdgMsGridder(const Settings& settings, const Resources& resources,
                           MsProviderCollection& ms_provider_collection)
    : MsGridder(settings, ms_provider_collection),
      _averageBeam(nullptr),
      _outputProvider(nullptr),
      _proxyType(idg::api::Type::CPU_OPTIMIZED),
      _buffersize(0),
      _resources(resources) {
  IdgConfiguration::Read(_proxyType, _buffersize, _options);
  _proxyType = GetIdgType(GetSettings());
  _bufferset = std::unique_ptr<idg::api::BufferSet>(
      idg::api::BufferSet::create(_proxyType));
  if (settings.gridWithBeam || !settings.atermConfigFilename.empty())
    _options["a_term_kernel_size"] = float(GetSettings().atermKernelSize);
  _options["max_threads"] = int(resources.NCpus());
  if (settings.gridMode == GriddingKernelMode::BlackmanHarris)
    _options["taper"] = std::string("blackman-harris");
}

size_t IdgMsGridder::GetNInversionPasses() const {
  if (GetPsfMode() == PsfMode::kNone && _averageBeam->Empty()) {
    // Two passes are used for this case: i) average beam creation; ii)
    // gridding.
    return 2;
  }
  return 1;
}

void IdgMsGridder::StartInversion() {
  const size_t untrimmed_width = ImageWidth();

  assert(TrimWidth() == TrimHeight());
  assert(untrimmed_width == ImageHeight());

  _options["padded_size"] = untrimmed_width;

  _options["stokes_I_only"] =
      (Polarization() == aocommon::Polarization::StokesI);

  if (!_averageBeam) _averageBeam.reset(new AverageBeam());

  double max_w = 0;
  for (size_t i = 0; i != GetMsCount(); ++i) {
    max_w = std::max(max_w, GetMsData(i).max_w_with_flags);
  }

  const double shift_l = LShift();
  const double shift_m = MShift();
  const double shift_p =
      std::sqrt(1.0 - shift_l * shift_l - shift_m * shift_m) - 1.0;
  _bufferset->init(TrimWidth(), ActualPixelSizeX(), max_w + 1.0, shift_l,
                   shift_m, shift_p, _options);
  Logger::Debug << "IDG subgrid size: " << _bufferset->get_subgridsize()
                << '\n';

  if (GetPsfMode() != PsfMode::kNone) {
    // Computing the PSF
    // For the PSF the aterm is not applied
    _bufferset->set_apply_aterm(false);
    _bufferset->unset_matrix_inverse_beam();
  } else {
    // Compute a dirty/residual image
    // with application of the a term
    _bufferset->set_apply_aterm(true);

    // Because compensation for the average beam happens at subgrid level
    // it needs to be known in advance.
    // If it is not in the cache it needs to be computed first, which will occur
    // inside StartInversionPass().
    if (!_averageBeam->Empty()) {
      // Set avg beam from cache
      Logger::Debug << "Using average beam from cache.\n";
      _bufferset->set_scalar_beam(_averageBeam->ScalarBeam());
      _bufferset->set_matrix_inverse_beam(_averageBeam->MatrixInverseBeam());
    }
  }
}

void IdgMsGridder::StartInversionPass(size_t pass_index) {
  if (pass_index == 0 && GetNInversionPasses() == 2) {
    // Compute avg beam
    Logger::Debug << "Computing average beam.\n";
    _bufferset->init_compute_avg_beam(idg::api::compute_flags::compute_only);
  } else {
    ResetVisibilityCounters();
  }
}

void IdgMsGridder::FinishInversionPass(size_t pass_index) {
  if (pass_index == 0 && GetNInversionPasses() == 2) {
    _bufferset->finalize_compute_avg_beam();
    Logger::Debug << "Finished computing average beam.\n";
    _averageBeam->SetScalarBeam(_bufferset->get_scalar_beam(), TrimWidth(),
                                TrimHeight());
    _averageBeam->SetMatrixInverseBeam(_bufferset->get_matrix_inverse_beam(),
                                       _bufferset->get_subgridsize(),
                                       _bufferset->get_subgridsize());
  }
}

void IdgMsGridder::FinishInversion() {
  const size_t n_image_polarizations = _options["stokes_I_only"] ? 1 : 4;

  // GridMeasurementSet calls have added the gridding result to _image member
  _image.assign(n_image_polarizations * TrimWidth() * TrimHeight(), 0.0);
  _bufferset->get_image(_image.data());
  if (GetPsfMode() != PsfMode::kNone) {
    Logger::Debug << "Total weight: " << ImageWeight() << '\n';
  }
  // result is now in _image member
  // Can be accessed by subsequent calls to ResultImages()
}

size_t IdgMsGridder::GridMeasurementSet(
    const MsProviderCollection::MsData& ms_data) {
  aocommon::UVector<std::complex<float>> aterm_buffer;
#ifdef HAVE_EVERYBEAM
  std::unique_ptr<ATermBase> aterm_maker;
  if (!prepareForMeasurementSet(ms_data, aterm_maker, aterm_buffer,
                                idg::api::BufferSetType::gridding))
    return 0;
#else
  if (!prepareForMeasurementSet(ms_data, aterm_buffer,
                                idg::api::BufferSetType::gridding))
    return 0;
#endif

  // TODO For now we do not allow multiple bands in one msprovider
  if (!ms_data.ms_provider->IsRegular()) {
    throw std::runtime_error(
        "IDG implementation does not support irregular data yet");
  }
  // Regular data should always have one band...
  assert(_selectedBands.BandCount() == 1);
  const aocommon::BandData& selected_band = *_selectedBands.begin();

  const size_t n_vis_polarizations = ms_data.ms_provider->NPolarizations();
  constexpr size_t n_idg_polarizations = 4;
  const size_t data_size = selected_band.ChannelCount() * n_idg_polarizations;
  aocommon::UVector<float> weight_buffer(data_size);
  aocommon::UVector<std::complex<float>> model_buffer(data_size);
  aocommon::UVector<bool> selection_buffer(selected_band.ChannelCount(), true);

  // The gridder doesn't need to know the absolute time index; this value
  // indexes relatively to where we start in the measurement set, and only
  // increases when the time changes.
  int time_index = -1;
  double current_time = -1.0;
  aocommon::UVector<double> uvws(ms_data.ms_provider->NAntennas() * 3, 0.0);

  TimestepBuffer timestep_buffer(ms_data.ms_provider, DoSubtractModel());
  std::unique_ptr<MSReader> ms_reader = timestep_buffer.MakeReader();
  TimestepBufferReader& timestep_reader =
      static_cast<TimestepBufferReader&>(*ms_reader);
  aocommon::UVector<std::complex<float>> row_visibilities(data_size);
  IDGInversionRow row_data;
  row_data.data = row_visibilities.data();

  while (ms_reader->CurrentRowAvailable()) {
    MSProvider::MetaData metadata;
    timestep_reader.ReadMeta(metadata);

    if (current_time != metadata.time) {
      current_time = metadata.time;
      time_index++;
#ifdef HAVE_EVERYBEAM
      if (aterm_maker) {
        timestep_reader.GetUVWsForTimestep(uvws);
        if (aterm_maker->Calculate(aterm_buffer.data(), current_time,
                                   selected_band.CentreFrequency(),
                                   metadata.field_id, uvws.data())) {
          _bufferset->get_gridder(kGridderIndex)
              ->set_aterm(time_index, aterm_buffer.data());
          Logger::Debug << "Calculated a-terms for timestep " << time_index
                        << "\n";
        }
      }
#endif
    }

    row_data.uvw[0] = metadata.u_in_m;
    row_data.uvw[1] = metadata.v_in_m;
    row_data.uvw[2] = metadata.w_in_m;

    row_data.antenna1 = metadata.antenna1;
    row_data.antenna2 = metadata.antenna2;
    row_data.time_index = time_index;

    const size_t n_parms = NumValuesPerSolution();
    if (n_vis_polarizations == 1) {
      if (n_parms == 2) {
        GetInstrumentalVisibilities<1, 2>(
            *ms_reader, ms_data.antenna_names.size(), row_data,
            weight_buffer.data(), model_buffer.data(), selection_buffer.data(),
            metadata);
      } else {
        GetInstrumentalVisibilities<1, 4>(
            *ms_reader, ms_data.antenna_names.size(), row_data,
            weight_buffer.data(), model_buffer.data(), selection_buffer.data(),
            metadata);
      }
      // The data is placed in the first quarter of the buffers: reverse copy it
      // and expand it to 4 polarizations. TODO at a later time, IDG should
      // be able to directly accept 1 polarization instead of 4.
      size_t source_index = row_visibilities.size() / 4;
      for (size_t i = row_visibilities.size(); i != 0; i -= 4) {
        row_visibilities[i - 1] = row_visibilities[source_index - 1];
        row_visibilities[i - 2] = 0.0;
        row_visibilities[i - 3] = 0.0;
        row_visibilities[i - 4] = row_visibilities[source_index - 1];
        weight_buffer[i - 1] = weight_buffer[source_index - 1];
        weight_buffer[i - 2] = weight_buffer[source_index - 1];
        weight_buffer[i - 3] = weight_buffer[source_index - 1];
        weight_buffer[i - 4] = weight_buffer[source_index - 1];
        source_index--;
      }
    } else if (n_vis_polarizations == 2) {
      if (n_parms == 2) {
        GetInstrumentalVisibilities<2, 2>(
            *ms_reader, ms_data.antenna_names.size(), row_data,
            weight_buffer.data(), model_buffer.data(), selection_buffer.data(),
            metadata);
      } else {
        GetInstrumentalVisibilities<2, 4>(
            *ms_reader, ms_data.antenna_names.size(), row_data,
            weight_buffer.data(), model_buffer.data(), selection_buffer.data(),
            metadata);
      }
      // The data is placed in the first half of the buffers: reverse copy it
      // and expand it to 4 polarizations. TODO at a later time, IDG should
      // be able to directly accept 2 pols instead of 4.
      size_t source_index = row_visibilities.size() / 2;
      for (size_t i = row_visibilities.size(); i != 0; i -= 4) {
        row_visibilities[i - 1] = row_visibilities[source_index - 1];
        row_visibilities[i - 2] = 0.0;
        row_visibilities[i - 3] = 0.0;
        row_visibilities[i - 4] = row_visibilities[source_index - 2];
        weight_buffer[i - 1] = weight_buffer[source_index - 1];
        weight_buffer[i - 2] = weight_buffer[source_index - 1];
        weight_buffer[i - 3] = weight_buffer[source_index - 2];
        weight_buffer[i - 4] = weight_buffer[source_index - 2];
        source_index -= 2;
      }
    } else {
      assert(n_vis_polarizations == 4);
      if (n_parms == 2) {
        GetInstrumentalVisibilities<4, 2>(
            *ms_reader, ms_data.antenna_names.size(), row_data,
            weight_buffer.data(), model_buffer.data(), selection_buffer.data(),
            metadata);
      } else {
        GetInstrumentalVisibilities<4, 4>(
            *ms_reader, ms_data.antenna_names.size(), row_data,
            weight_buffer.data(), model_buffer.data(), selection_buffer.data(),
            metadata);
      }
    }

    row_data.uvw[1] = -metadata.v_in_m;  // DEBUG vdtol, flip axis
    row_data.uvw[2] = -metadata.w_in_m;  //

    _bufferset->get_gridder(kGridderIndex)
        ->grid_visibilities(time_index, metadata.antenna1, metadata.antenna2,
                            row_data.uvw, row_data.data, weight_buffer.data());

    ms_reader->NextInputRow();
  }
  _bufferset->finished();

  return 0;
}

void IdgMsGridder::StartPredict(std::vector<Image>&& images) {
  if (images.size() == 2)
    throw std::runtime_error("IDG gridder cannot make complex images");
  const size_t untrimmed_width = ImageWidth();
  const size_t width = TrimWidth();
  const size_t height = TrimHeight();

  assert(width == height);
  assert(untrimmed_width == ImageHeight());

  _options["padded_size"] = untrimmed_width;

  const bool stokes_I_only =
      (Polarization() == aocommon::Polarization::StokesI);
  _options["stokes_I_only"] = stokes_I_only;
  const size_t n_image_polarizations = stokes_I_only ? 1 : 4;

  _image.assign(n_image_polarizations * width * height, 0.0);
  if (!_averageBeam) {
    Logger::Debug << "No average beam in cache, creating an empty one.\n";
    _averageBeam.reset(new AverageBeam());
  }

  assert(images.size() == n_image_polarizations);
  if (Polarization() == aocommon::Polarization::FullStokes) {
    for (size_t polarization_index = 0;
         polarization_index != n_image_polarizations; ++polarization_index) {
      std::copy_n(images[polarization_index].Data(), width * height,
                  _image.data() + polarization_index * width * height);
    }
  } else {
    const size_t stokes_index =
        aocommon::Polarization::StokesToIndex(Polarization());
    std::copy_n(images[0].Data(), width * height,
                _image.data() + stokes_index * width * height);
  }

  bool do_scale = false;
  if (!_averageBeam->Empty()) {
    // Set avg beam from cache
    Logger::Debug << "Average beam is already in cache.\n";
    _bufferset->set_scalar_beam(_averageBeam->ScalarBeam());
    do_scale = true;
  }

  double max_w = 0;
  for (size_t i = 0; i != GetMsCount(); ++i) {
    max_w = std::max(max_w, GetMsData(i).max_w_with_flags);
  }

  const double shift_l = LShift();
  const double shift_m = MShift();
  const double shift_p =
      std::sqrt(1.0 - shift_l * shift_l - shift_m * shift_m) - 1.0;
  _bufferset->init(width, ActualPixelSizeX(), max_w + 1.0, shift_l, shift_m,
                   shift_p, _options);
  _bufferset->set_image(_image.data(), do_scale);
}

void IdgMsGridder::FinishPredict() {}

size_t IdgMsGridder::PredictMeasurementSet(
    const MsProviderCollection::MsData& ms_data) {
  aocommon::UVector<std::complex<float>> aterm_buffer;
#ifdef HAVE_EVERYBEAM
  std::unique_ptr<ATermBase> aterm_maker;
  if (!prepareForMeasurementSet(ms_data, aterm_maker, aterm_buffer,
                                idg::api::BufferSetType::degridding))
    return 0;
#else
  if (!prepareForMeasurementSet(ms_data, aterm_buffer,
                                idg::api::BufferSetType::degridding))
    return 0;
#endif

  ms_data.ms_provider->ReopenRW();

  _outputProvider = ms_data.ms_provider;

  // TODO For now we do not allow multiple bands in one msprovider
  if (!ms_data.ms_provider->IsRegular())
    throw std::runtime_error(
        "IDG implementation does not support irregular data yet");
  // Regular data should always have one band...
  assert(_selectedBands.BandCount() == 1);
  const aocommon::BandData& selected_band = *_selectedBands.begin();

  constexpr size_t n_idg_polarizations = 4;
  aocommon::UVector<std::complex<float>> buffer(selected_band.ChannelCount() *
                                                n_idg_polarizations);

  int time_index = -1;
  double current_time = -1.0;
  aocommon::UVector<double> uvws(ms_data.ms_provider->NAntennas() * 3, 0.0);

  TimestepBuffer timestep_buffer(ms_data.ms_provider, false);
  timestep_buffer.ResetWritePosition();
  for (std::unique_ptr<MSReader> ms_reader = timestep_buffer.MakeReader();
       ms_reader->CurrentRowAvailable(); ms_reader->NextInputRow()) {
    TimestepBufferReader& timestep_reader =
        static_cast<TimestepBufferReader&>(*ms_reader);

    MSProvider::MetaData metadata;
    timestep_reader.ReadMeta(metadata);

    const size_t provRowId = timestep_reader.RowId();
    if (current_time != metadata.time) {
      current_time = metadata.time;
      time_index++;

#ifdef HAVE_EVERYBEAM
      if (aterm_maker) {
        timestep_reader.GetUVWsForTimestep(uvws);
        if (aterm_maker->Calculate(aterm_buffer.data(), current_time,
                                   selected_band.CentreFrequency(),
                                   metadata.field_id, uvws.data())) {
          _bufferset->get_degridder(kGridderIndex)
              ->set_aterm(time_index, aterm_buffer.data());
          Logger::Debug << "Calculated new a-terms for timestep " << time_index
                        << "\n";
        }
      }
#endif
    }

    IDGPredictionRow row;
    row.uvw[0] = metadata.u_in_m;
    row.uvw[1] = -metadata.v_in_m;
    row.uvw[2] = -metadata.w_in_m;
    row.antenna1 = metadata.antenna1;
    row.antenna2 = metadata.antenna2;
    row.time_index = time_index;
    row.row_id = provRowId;
    predictRow(row, ms_data.antenna_names);
  }

  computePredictionBuffer(ms_data.antenna_names);
  return 0;
}

void IdgMsGridder::predictRow(IDGPredictionRow& row,
                              const std::vector<std::string>& antenna_names) {
  while (_bufferset->get_degridder(kGridderIndex)
             ->request_visibilities(row.row_id, row.time_index, row.antenna1,
                                    row.antenna2, row.uvw)) {
    computePredictionBuffer(antenna_names);
  }
}

void IdgMsGridder::computePredictionBuffer(
    const std::vector<std::string>& antenna_names) {
  auto available_row_ids = _bufferset->get_degridder(kGridderIndex)->compute();
  Logger::Debug << "Computed " << available_row_ids.size() << " rows.\n";
  const size_t n_vis_polarizations = _outputProvider->NPolarizations();
  for (std::pair<long unsigned, std::complex<float>*>& row :
       available_row_ids) {
    MSProvider::MetaData metadata;
    ReadPredictMetaData(metadata);
    const aocommon::BandData& band = _selectedBands[metadata.data_desc_id];
    if (n_vis_polarizations == 1) {
      // Place Stokes I in the first quarter of the array
      for (size_t i = 0; i != band.ChannelCount(); ++i) {
        row.second[i] = (row.second[i * 4] + row.second[i * 4 + 3]) / 2.0f;
      }

    } else if (n_vis_polarizations == 2) {
      // Remove the XY/YX pols from the data and place the result in the first
      // half of the array
      for (size_t i = 0; i != band.ChannelCount(); ++i) {
        row.second[i * 2] = row.second[i * 4];
        row.second[i * 2 + 1] = row.second[i * 4 + 3];
      }
    } else {
      assert(n_vis_polarizations == 4);
    }
    const double* uvw = nullptr;
    WriteInstrumentalVisibilities(*_outputProvider, antenna_names.size(),
                                  metadata.data_desc_id, row.second, uvw,
                                  metadata.field_id, metadata.antenna1,
                                  metadata.antenna2, metadata.time);
  }
  _bufferset->get_degridder(kGridderIndex)->finished_reading();
}

std::vector<Image> IdgMsGridder::ResultImages() {
  const size_t width = TrimWidth();
  const size_t height = TrimHeight();
  std::vector<Image> images;
  if (Polarization() == aocommon::Polarization::FullStokes) {
    images.reserve(4);
    for (size_t polIndex = 0; polIndex != 4; ++polIndex) {
      images.emplace_back(width, height);
      std::copy_n(_image.data() + polIndex * width * height, width * height,
                  images[polIndex].Data());
    }
  } else {
    size_t polIndex = aocommon::Polarization::StokesToIndex(Polarization());
    images.emplace_back(width, height);
    std::copy_n(_image.data() + polIndex * width * height, width * height,
                images[0].Data());
  }
  return images;
}

void IdgMsGridder::SetAverageBeam(std::unique_ptr<AverageBeam> average_beam) {
  _averageBeam = std::move(average_beam);
}

std::unique_ptr<AverageBeam> IdgMsGridder::ReleaseAverageBeam() {
  return std::move(_averageBeam);
}

void IdgMsGridder::SaveBeamImage(const ImagingTableEntry& entry,
                                 ImageFilename& filename,
                                 const Settings& settings, double ra,
                                 double dec, double pdl, double pdm,
                                 const AverageBeam& average_beam) {
  aocommon::FitsWriter writer;
  writer.SetImageDimensions(settings.trimmedImageWidth,
                            settings.trimmedImageHeight, ra, dec,
                            settings.pixelScaleX, settings.pixelScaleY);
  writer.SetPhaseCentreShift(pdl, pdm);
  ImageFilename polName(filename);
  polName.SetPolarization(aocommon::Polarization::StokesI);
  writer.SetPolarization(aocommon::Polarization::StokesI);
  writer.SetFrequency(entry.CentralFrequency(),
                      entry.bandEndFrequency - entry.bandStartFrequency);
  writer.Write(polName.GetBeamPrefix(settings) + ".fits",
               average_beam.ScalarBeam()->data());
}

void IdgMsGridder::SavePBCorrectedImages(aocommon::FitsWriter& writer,
                                         const ImageFilename& filename,
                                         const std::string& filenameKind,
                                         const Settings& settings) {
  ImageFilename beamName(filename);
  beamName.SetPolarization(aocommon::Polarization::StokesI);
  aocommon::FitsReader reader(beamName.GetBeamPrefix(settings) + ".fits");

  Image beam(reader.ImageWidth(), reader.ImageHeight());
  reader.Read(beam.Data());

  Image image;
  for (size_t polIndex = 0; polIndex != 4; ++polIndex) {
    aocommon::PolarizationEnum pol =
        aocommon::Polarization::IndexToStokes(polIndex);
    ImageFilename name(filename);
    name.SetPolarization(pol);
    aocommon::FitsReader reader(name.GetPrefix(settings) + "-" + filenameKind +
                                ".fits");
    if (image.Empty()) image = Image(reader.ImageWidth(), reader.ImageHeight());
    reader.Read(image.Data());

    for (size_t i = 0; i != reader.ImageWidth() * reader.ImageHeight(); ++i) {
      if (beam[i] > 1e-6)
        image[i] /= beam[i];
      else
        image[i] = std::numeric_limits<double>::quiet_NaN();
    }

    writer.SetPolarization(pol);
    writer.Write(name.GetPrefix(settings) + "-" + filenameKind + "-pb.fits",
                 image.Data());
  }
}

#ifdef HAVE_EVERYBEAM
bool IdgMsGridder::prepareForMeasurementSet(
    const MsProviderCollection::MsData& ms_data,
    std::unique_ptr<ATermBase>& aTermMaker,
    aocommon::UVector<std::complex<float>>& aTermBuffer,
    idg::api::BufferSetType bufferSetType) {
#else
bool IdgMsGridder::prepareForMeasurementSet(
    const MsProviderCollection::MsData& ms_data,
    aocommon::UVector<std::complex<float>>& aTermBuffer,
    idg::api::BufferSetType bufferSetType) {
#endif  // HAVE_EVERYBEAM
  const float max_baseline = ms_data.max_baseline_in_m;
  // Skip this ms if there is no data in it
  if (!max_baseline) return false;

  _selectedBands = ms_data.ms_provider->SelectedBands();

  // TODO For now we do not allow multiple bands in one msprovider
  if (!ms_data.ms_provider->IsRegular())
    throw std::runtime_error(
        "IDG implementation does not support irregular data yet");
  // Regular data should always have one band...
  assert(_selectedBands.BandCount() == 1);
  const aocommon::BandData& selected_band = *_selectedBands.begin();

  // TODO for now we map the ms antennas directly to the gridder's antenna,
  // including non-selected antennas. Later this can be made more efficient.
  const size_t nStations = ms_data.ms_provider->NAntennas();

  std::vector<std::vector<double>> bands;
  bands.emplace_back(selected_band.begin(), selected_band.end());
  const size_t nChannels = selected_band.ChannelCount();

  // Only one-third of the mem is allocated to the buffers, so that memory
  // remains available for the images and other things done by IDG.
  // Never use more than 16 GB
  const size_t memSize = std::min<uint64_t>(16ul * 1024ul * 1024ul * 1024ul,
                                            _resources.Memory() / 3);
  uint64_t memPerTimestep =
      idg::api::BufferSet::get_memory_per_timestep(nStations, nChannels);
#ifdef HAVE_EVERYBEAM
  aTermMaker = getATermMaker(ms_data);
  if (aTermMaker) {
    const size_t subgridsize = _bufferset->get_subgridsize();
    aTermBuffer.resize(subgridsize * subgridsize * 4 * nStations);
    // When a-terms are used, they will also take memory. Here we calculate
    // their approx contribution.
    double avgUpdate = aTermMaker->AverageUpdateTime();
    Logger::Debug << "A-terms change on average every " << avgUpdate
                  << " s, once every " << (avgUpdate / ms_data.integration_time)
                  << " timesteps.\n";
    const uint64_t atermMemPerTimestep =
        subgridsize * subgridsize * nStations *  // size of grid x nr of grids
        (4 * 8) *  // 4 pol, 8 bytes per complex value
        (ms_data.integration_time /
         avgUpdate);  // Average number of aterms per timestep
    Logger::Debug << "A-terms increase mem per timestep from " << memPerTimestep
                  << " bytes to " << (memPerTimestep + atermMemPerTimestep)
                  << " bytes.\n";
    memPerTimestep += atermMemPerTimestep;
  }
#else
  if (!GetSettings().atermConfigFilename.empty() ||
      GetSettings().gridWithBeam) {
    throw std::runtime_error(
        "ATerm correction requested, but the software has been compiled "
        "without EveryBeam. Recompile your software and make sure that "
        "cmake finds the EveryBeam library.");
  }
#endif  // HAVE_EVERYBEAM

  // IDG can allocate two visibility buffers: (for parallel processing)
  memPerTimestep *= 2;

  _buffersize = std::max<size_t>(1, memSize / memPerTimestep);

  Logger::Debug << "Allocatable timesteps (" << nStations << " stations, "
                << nChannels << " channels, " << memSize / (1024 * 1024 * 1024)
                << " GB mem): " << _buffersize << '\n';
  _bufferset->init_buffers(_buffersize, bands, nStations, max_baseline,
                           _options, bufferSetType);

  return true;
}

#ifdef HAVE_EVERYBEAM
std::unique_ptr<class ATermBase> IdgMsGridder::getATermMaker(
    const MsProviderCollection::MsData& ms_data) {
  const size_t nr_stations = ms_data.ms_provider->NAntennas();
  if (!GetSettings().atermConfigFilename.empty() ||
      GetSettings().gridWithBeam) {
    // IDG uses a flipped coordinate system which is moved by half a pixel:
    CoordinateSystem system;
    system.width = _bufferset->get_subgridsize();
    system.height = system.width;
    system.ra = PhaseCentreRA();
    system.dec = PhaseCentreDec();
    system.dl = -_bufferset->get_subgrid_pixelsize();
    system.dm = -_bufferset->get_subgrid_pixelsize();
    system.l_shift = LShift() - 0.5 * system.dl;
    system.m_shift = MShift() + 0.5 * system.dm;

    everybeam::ATermSettings aterm_settings;
    aterm_settings.save_aterms_prefix = GetSettings().prefixName;
    aterm_settings.data_column_name = GetSettings().dataColumnName;
    aterm_settings.filenames = GetSettings().filenames;
    aterm_settings.aterm_update_interval = GetSettings().beamAtermUpdateTime;
    aterm_settings.padded_image_width = GetSettings().paddedImageWidth;
    aterm_settings.trimmed_image_width = GetSettings().trimmedImageWidth;
    aterm_settings.max_support = GetSettings().atermKernelSize;

    std::optional<SynchronizedMS> optional_ms =
        ms_data.ms_provider->MsIfAvailable();
    if (!optional_ms) {
      throw std::runtime_error(
          "Gridding with beam or a-term configuration file is currently not "
          "supported without an on-disk measurement set");
    }
    if (!GetSettings().atermConfigFilename.empty()) {
      ParsetReader parset_aterms(GetSettings().atermConfigFilename);
      // If MWA MS and "beam" in aterms, get the path to the coefficient file
      if (everybeam::GetTelescopeType(**optional_ms) ==
          everybeam::TelescopeType::kMWATelescope) {
        std::vector<std::string> aterms = parset_aterms.GetStringList("aterms");
        if (std::find(aterms.begin(), aterms.end(), "beam") != aterms.end()) {
          aterm_settings.coeff_path =
              wsclean::mwa::FindCoeffFile(GetSettings().mwaPath);
        }
      }

      std::unique_ptr<ATermConfig> config(
          new ATermConfig(nr_stations, system, aterm_settings));
      config->SetSaveATerms(GetSettings().saveATerms, GetSettings().prefixName);
      config->Read(**optional_ms, parset_aterms, optional_ms->Filename());
      return config;
    } else {
      // If MWA MS, get the path to the coefficient file
      if (everybeam::GetTelescopeType(**optional_ms) ==
          everybeam::TelescopeType::kMWATelescope) {
        aterm_settings.coeff_path =
            wsclean::mwa::FindCoeffFile(GetSettings().mwaPath);
      }
      bool frequencyInterpolation = true;
      bool useChannelFrequency = true;
      std::string elementResponseModel = !GetSettings().beamModel.empty()
                                             ? GetSettings().beamModel
                                             : "DEFAULT";

      std::unique_ptr<ATermBeam> beam = ATermConfig::GetATermBeam(
          **optional_ms, system, aterm_settings, frequencyInterpolation,
          GetSettings().beamNormalisationMode, useChannelFrequency,
          elementResponseModel, GetSettings().beamMode);
      beam->SetSaveATerms(GetSettings().saveATerms, GetSettings().prefixName);
      beam->SetUpdateInterval(GetSettings().beamAtermUpdateTime);
      return beam;
    }
  } else {
    return std::unique_ptr<ATermBase>();
  }
}

#endif  // HAVE_EVERYBEAM

}  // namespace wsclean
