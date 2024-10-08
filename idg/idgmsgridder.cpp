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

namespace {
constexpr const size_t kGridderIndex = 0;
}

IdgMsGridder::IdgMsGridder(const Settings& settings, const Resources& resources)
    : MSGridderBase(settings),
      _averageBeam(nullptr),
      _outputProvider(nullptr),
      _proxyType(idg::api::Type::CPU_OPTIMIZED),
      _buffersize(0),
      _resources(resources) {
  IdgConfiguration::Read(_proxyType, _buffersize, _options);
  setIdgType();
  _bufferset = std::unique_ptr<idg::api::BufferSet>(
      idg::api::BufferSet::create(_proxyType));
  if (settings.gridWithBeam || !settings.atermConfigFilename.empty())
    _options["a_term_kernel_size"] = float(_settings.atermKernelSize);
  _options["max_threads"] = int(resources.NCpus());
  if (settings.gridMode == GriddingKernelMode::BlackmanHarris)
    _options["taper"] = std::string("blackman-harris");
}

IdgMsGridder::~IdgMsGridder() {
  Logger::Info << "Gridding: " << _griddingWatch.ToString()
               << ", degridding: " << _degriddingWatch.ToString() << '\n';
}

void IdgMsGridder::Invert() {
  const size_t untrimmedWidth = ImageWidth();
  const size_t width = TrimWidth(), height = TrimHeight();

  assert(width == height);
  assert(untrimmedWidth == ImageHeight());

  _options["padded_size"] = untrimmedWidth;

  const bool stokes_I_only =
      (Polarization() == aocommon::Polarization::StokesI);
  _options["stokes_I_only"] = stokes_I_only;
  const size_t n_image_polarizations = stokes_I_only ? 1 : 4;

  if (!_averageBeam) _averageBeam.reset(new AverageBeam());

  std::vector<MSData> msDataVector;
  initializeMSDataVector(msDataVector);

  double max_w = 0;
  for (const MSData& msData : msDataVector) {
    max_w = std::max(max_w, msData.maxWWithFlags);
  }

  const double shiftl = LShift();
  const double shiftm = MShift();
  const double shiftp =
      std::sqrt(1.0 - shiftl * shiftl - shiftm * shiftm) - 1.0;
  _bufferset->init(width, ActualPixelSizeX(), max_w + 1.0, shiftl, shiftm,
                   shiftp, _options);
  Logger::Debug << "IDG subgrid size: " << _bufferset->get_subgridsize()
                << '\n';

  if (GetPsfMode() != PsfMode::kNone) {
    // Computing the PSF
    // For the PSF the aterm is not applied
    _bufferset->set_apply_aterm(false);
    _bufferset->unset_matrix_inverse_beam();
    resetVisibilityCounters();
    for (const MSData& msData : msDataVector) {
      // Adds the gridding result to _image member
      gridMeasurementSet(msData);
    }
    _image.assign(n_image_polarizations * width * height, 0.0);
    _bufferset->get_image(_image.data());

    Logger::Debug << "Total weight: " << totalWeight() << '\n';
  } else {
    // Compute a dirty/residual image
    // with application of the a term
    _bufferset->set_apply_aterm(true);

    // Because compensation for the average beam happens at subgrid level
    // it needs to be known in advance.
    // If it is not in the cache it needs to be computed first
    if (!_averageBeam->Empty()) {
      // Set avg beam from cache
      Logger::Debug << "Using average beam from cache.\n";
      _bufferset->set_scalar_beam(_averageBeam->ScalarBeam());
      _bufferset->set_matrix_inverse_beam(_averageBeam->MatrixInverseBeam());
    } else {
      // Compute avg beam
      Logger::Debug << "Computing average beam.\n";
      _bufferset->init_compute_avg_beam(idg::api::compute_flags::compute_only);
      for (const MSData& msData : msDataVector) {
        gridMeasurementSet(msData);
      }
      _bufferset->finalize_compute_avg_beam();
      Logger::Debug << "Finished computing average beam.\n";
      _averageBeam->SetScalarBeam(_bufferset->get_scalar_beam(), width, height);
      _averageBeam->SetMatrixInverseBeam(_bufferset->get_matrix_inverse_beam(),
                                         _bufferset->get_subgridsize(),
                                         _bufferset->get_subgridsize());
    }

    resetVisibilityCounters();
    for (const MSData& msData : msDataVector) {
      // Adds the gridding result to _image member
      gridMeasurementSet(msData);
    }
    _image.assign(n_image_polarizations * width * height, 0.0);
    _bufferset->get_image(_image.data());
  }

  // result is now in _image member
  // Can be accessed by subsequent calls to ResultImages()
}

void IdgMsGridder::gridMeasurementSet(const MSGridderBase::MSData& msData) {
  aocommon::UVector<std::complex<float>> aTermBuffer;

#ifdef HAVE_EVERYBEAM
  std::unique_ptr<ATermBase> aTermMaker;
  if (!prepareForMeasurementSet(msData, aTermMaker, aTermBuffer,
                                idg::api::BufferSetType::gridding))
    return;
#else
  if (!prepareForMeasurementSet(msData, aTermBuffer,
                                idg::api::BufferSetType::gridding))
    return;
#endif

  StartMeasurementSet(msData, false);

  const size_t n_vis_polarizations = msData.msProvider->NPolarizations();
  constexpr size_t n_idg_polarizations = 4;
  aocommon::UVector<float> weightBuffer(_selectedBand.ChannelCount() *
                                        n_idg_polarizations);
  aocommon::UVector<std::complex<float>> modelBuffer(
      _selectedBand.ChannelCount() * n_idg_polarizations);
  aocommon::UVector<bool> isSelected(
      _selectedBand.ChannelCount() * n_idg_polarizations, true);
  aocommon::UVector<std::complex<float>> dataBuffer(
      _selectedBand.ChannelCount() * n_idg_polarizations);

  _griddingWatch.Start();

  // The gridder doesn't need to know the absolute time index; this value
  // indexes relatively to where we start in the measurement set, and only
  // increases when the time changes.
  int timeIndex = -1;
  double currentTime = -1.0;
  aocommon::UVector<double> uvws(msData.msProvider->NAntennas() * 3, 0.0);

  TimestepBuffer timestepBuffer(msData.msProvider, DoSubtractModel());
  for (std::unique_ptr<MSReader> msReader = timestepBuffer.MakeReader();
       msReader->CurrentRowAvailable(); msReader->NextInputRow()) {
    TimestepBufferReader& timestepReader =
        static_cast<TimestepBufferReader&>(*msReader);
    MSProvider::MetaData metaData;
    timestepReader.ReadMeta(metaData);

    if (currentTime != metaData.time) {
      currentTime = metaData.time;
      timeIndex++;
#ifdef HAVE_EVERYBEAM
      if (aTermMaker) {
        timestepReader.GetUVWsForTimestep(uvws);
        if (aTermMaker->Calculate(aTermBuffer.data(), currentTime,
                                  _selectedBand.CentreFrequency(),
                                  metaData.fieldId, uvws.data())) {
          _bufferset->get_gridder(kGridderIndex)
              ->set_aterm(timeIndex, aTermBuffer.data());
          Logger::Debug << "Calculated a-terms for timestep " << timeIndex
                        << "\n";
        }
      }
#endif
    }
    IDGInversionRow rowData;

    rowData.data = dataBuffer.data();
    rowData.uvw[0] = metaData.uInM;
    rowData.uvw[1] = metaData.vInM;
    rowData.uvw[2] = metaData.wInM;

    rowData.antenna1 = metaData.antenna1;
    rowData.antenna2 = metaData.antenna2;
    rowData.timeIndex = timeIndex;

    if (n_vis_polarizations == 1) {
      readAndWeightVisibilities<1>(*msReader, msData.antennaNames, rowData,
                                   _selectedBand, weightBuffer.data(),
                                   modelBuffer.data(), isSelected.data(),
                                   GainMode::kDiagonal);
      // The data is placed in the first quarter of the buffers: reverse copy it
      // and expand it to 4 polarizations. TODO at a later time, IDG should
      // be able to directly accept 1 polarization instead of 4.
      size_t source_index = dataBuffer.size() / 4;
      for (size_t i = dataBuffer.size(); i != 0; i -= 4) {
        dataBuffer[i - 1] = dataBuffer[source_index - 1];
        dataBuffer[i - 2] = 0.0;
        dataBuffer[i - 3] = 0.0;
        dataBuffer[i - 4] = dataBuffer[source_index - 1];
        weightBuffer[i - 1] = weightBuffer[source_index - 1];
        weightBuffer[i - 2] = weightBuffer[source_index - 1];
        weightBuffer[i - 3] = weightBuffer[source_index - 1];
        weightBuffer[i - 4] = weightBuffer[source_index - 1];
        source_index--;
      }
    } else if (n_vis_polarizations == 2) {
      readAndWeightVisibilities<2>(*msReader, msData.antennaNames, rowData,
                                   _selectedBand, weightBuffer.data(),
                                   modelBuffer.data(), isSelected.data(),
                                   GainMode::kFull);
      // The data is placed in the first half of the buffers: reverse copy it
      // and expand it to 4 polarizations. TODO at a later time, IDG should
      // be able to directly accept 2 pols instead of 4.
      size_t source_index = dataBuffer.size() / 2;
      for (size_t i = dataBuffer.size(); i != 0; i -= 4) {
        rowData.data[i - 1] = rowData.data[source_index - 1];
        rowData.data[i - 2] = 0.0;
        rowData.data[i - 3] = 0.0;
        rowData.data[i - 4] = rowData.data[source_index - 2];
        weightBuffer[i - 1] = weightBuffer[source_index - 1];
        weightBuffer[i - 2] = weightBuffer[source_index - 1];
        weightBuffer[i - 3] = weightBuffer[source_index - 2];
        weightBuffer[i - 4] = weightBuffer[source_index - 2];
        source_index -= 2;
      }
    } else {
      assert(n_vis_polarizations == 4);
      readAndWeightVisibilities<4>(*msReader, msData.antennaNames, rowData,
                                   _selectedBand, weightBuffer.data(),
                                   modelBuffer.data(), isSelected.data(),
                                   GainMode::kFull);
    }

    rowData.uvw[1] = -metaData.vInM;  // DEBUG vdtol, flip axis
    rowData.uvw[2] = -metaData.wInM;  //

    _bufferset->get_gridder(kGridderIndex)
        ->grid_visibilities(timeIndex, metaData.antenna1, metaData.antenna2,
                            rowData.uvw, rowData.data, weightBuffer.data());
  }
  _bufferset->finished();

  _griddingWatch.Pause();
}

void IdgMsGridder::Predict(std::vector<Image>&& images) {
  if (images.size() == 2)
    throw std::runtime_error("IDG gridder cannot make complex images");
  const size_t untrimmedWidth = ImageWidth();
  const size_t width = TrimWidth(), height = TrimHeight();

  assert(width == height);
  assert(untrimmedWidth == ImageHeight());

  _options["padded_size"] = untrimmedWidth;

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
    for (size_t polIndex = 0; polIndex != n_image_polarizations; ++polIndex) {
      std::copy_n(images[polIndex].Data(), width * height,
                  _image.data() + polIndex * width * height);
    }
  } else {
    const size_t stokesIndex =
        aocommon::Polarization::StokesToIndex(Polarization());
    std::copy_n(images[0].Data(), width * height,
                _image.data() + stokesIndex * width * height);
  }

  bool do_scale = false;
  if (!_averageBeam->Empty()) {
    // Set avg beam from cache
    Logger::Debug << "Average beam is already in cache.\n";
    _bufferset->set_scalar_beam(_averageBeam->ScalarBeam());
    do_scale = true;
  }

  std::vector<MSData> msDataVector;
  initializeMSDataVector(msDataVector);

  double max_w = 0;
  for (const MSData& msData : msDataVector) {
    max_w = std::max(max_w, msData.maxWWithFlags);
  }

  const double shiftl = LShift();
  const double shiftm = MShift();
  const double shiftp =
      std::sqrt(1.0 - shiftl * shiftl - shiftm * shiftm) - 1.0;
  _bufferset->init(width, ActualPixelSizeX(), max_w + 1.0, shiftl, shiftm,
                   shiftp, _options);
  _bufferset->set_image(_image.data(), do_scale);

  for (const MSData& msData : msDataVector) predictMeasurementSet(msData);
}

void IdgMsGridder::setIdgType() {
  switch (_settings.idgMode) {
    default:
      return;
    case Settings::IDG_CPU:
      _proxyType = idg::api::Type::CPU_OPTIMIZED;
      return;
    case Settings::IDG_GPU:
      _proxyType = idg::api::Type::CUDA_GENERIC;
      return;
    case Settings::IDG_HYBRID:
      _proxyType = idg::api::Type::HYBRID_CUDA_CPU_OPTIMIZED;
      return;
  }
}

void IdgMsGridder::predictMeasurementSet(const MSGridderBase::MSData& msData) {
  aocommon::UVector<std::complex<float>> aTermBuffer;
#ifdef HAVE_EVERYBEAM
  std::unique_ptr<ATermBase> aTermMaker;
  if (!prepareForMeasurementSet(msData, aTermMaker, aTermBuffer,
                                idg::api::BufferSetType::degridding))
    return;
#else
  if (!prepareForMeasurementSet(msData, aTermBuffer,
                                idg::api::BufferSetType::degridding))
    return;
#endif

  msData.msProvider->ReopenRW();

  _outputProvider = msData.msProvider;
  StartMeasurementSet(msData, true);

  constexpr size_t n_idg_polarizations = 4;
  aocommon::UVector<std::complex<float>> buffer(_selectedBand.ChannelCount() *
                                                n_idg_polarizations);
  _degriddingWatch.Start();

  int timeIndex = -1;
  double currentTime = -1.0;
  aocommon::UVector<double> uvws(msData.msProvider->NAntennas() * 3, 0.0);

  TimestepBuffer timestepBuffer(msData.msProvider, false);
  timestepBuffer.ResetWritePosition();
  for (std::unique_ptr<MSReader> msReader = timestepBuffer.MakeReader();
       msReader->CurrentRowAvailable(); msReader->NextInputRow()) {
    TimestepBufferReader& timestepReader =
        static_cast<TimestepBufferReader&>(*msReader);

    MSProvider::MetaData metaData;
    timestepReader.ReadMeta(metaData);

    const size_t provRowId = timestepReader.RowId();
    if (currentTime != metaData.time) {
      currentTime = metaData.time;
      timeIndex++;

#ifdef HAVE_EVERYBEAM
      if (aTermMaker) {
        timestepReader.GetUVWsForTimestep(uvws);
        if (aTermMaker->Calculate(aTermBuffer.data(), currentTime,
                                  _selectedBand.CentreFrequency(),
                                  metaData.fieldId, uvws.data())) {
          _bufferset->get_degridder(kGridderIndex)
              ->set_aterm(timeIndex, aTermBuffer.data());
          Logger::Debug << "Calculated new a-terms for timestep " << timeIndex
                        << "\n";
        }
      }
#endif
    }

    IDGPredictionRow row;
    row.uvw[0] = metaData.uInM;
    row.uvw[1] = -metaData.vInM;
    row.uvw[2] = -metaData.wInM;
    row.antenna1 = metaData.antenna1;
    row.antenna2 = metaData.antenna2;
    row.timeIndex = timeIndex;
    row.rowId = provRowId;
    predictRow(row, msData.antennaNames);
  }

  computePredictionBuffer(msData.antennaNames);
}

void IdgMsGridder::predictRow(IDGPredictionRow& row,
                              const std::vector<std::string>& antennaNames) {
  while (_bufferset->get_degridder(kGridderIndex)
             ->request_visibilities(row.rowId, row.timeIndex, row.antenna1,
                                    row.antenna2, row.uvw)) {
    computePredictionBuffer(antennaNames);
  }
}

void IdgMsGridder::computePredictionBuffer(
    const std::vector<std::string>& antennaNames) {
  auto available_row_ids = _bufferset->get_degridder(kGridderIndex)->compute();
  Logger::Debug << "Computed " << available_row_ids.size() << " rows.\n";
  const size_t n_vis_polarizations = _outputProvider->NPolarizations();
  for (std::pair<long unsigned, std::complex<float>*>& row :
       available_row_ids) {
    if (n_vis_polarizations == 1) {
      // Remove the XY/YX pols from the data and place the result in the first
      // quarter of the array
      for (size_t i = 0; i != _selectedBand.ChannelCount(); ++i) {
        row.second[i] = (row.second[i * 4] + row.second[i * 4 + 3]) / 2.0f;
      }
      writeVisibilities<1>(*_outputProvider, antennaNames, _selectedBand,
                           row.second, GainMode::kDiagonal);
    } else if (n_vis_polarizations == 2) {
      // Remove the XY/YX pols from the data and place the result in the first
      // half of the array
      for (size_t i = 0; i != _selectedBand.ChannelCount(); ++i) {
        row.second[i * 2] = row.second[i * 4];
        row.second[i * 2 + 1] = row.second[i * 4 + 3];
      }
      writeVisibilities<2>(*_outputProvider, antennaNames, _selectedBand,
                           row.second, GainMode::kFull);
    } else {
      assert(n_vis_polarizations == 4);
      writeVisibilities<4>(*_outputProvider, antennaNames, _selectedBand,
                           row.second, GainMode::kFull);
    }
  }
  _bufferset->get_degridder(kGridderIndex)->finished_reading();
  _degriddingWatch.Pause();
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
    const MSGridderBase::MSData& msData, std::unique_ptr<ATermBase>& aTermMaker,
    aocommon::UVector<std::complex<float>>& aTermBuffer,
    idg::api::BufferSetType bufferSetType) {
#else
bool IdgMsGridder::prepareForMeasurementSet(
    const MSGridderBase::MSData& msData,
    aocommon::UVector<std::complex<float>>& aTermBuffer,
    idg::api::BufferSetType bufferSetType) {
#endif  // HAVE_EVERYBEAM
  const float max_baseline = msData.maxBaselineInM;
  // Skip this ms if there is no data in it
  if (!max_baseline) return false;

  _selectedBand = msData.SelectedBand();

  // TODO for now we map the ms antennas directly to the gridder's antenna,
  // including non-selected antennas. Later this can be made more efficient.
  const size_t nStations = msData.msProvider->MS()->antenna().nrow();

  std::vector<std::vector<double>> bands;
  bands.emplace_back(_selectedBand.begin(), _selectedBand.end());
  const size_t nChannels = _selectedBand.ChannelCount();

  // Only one-third of the mem is allocated to the buffers, so that memory
  // remains available for the images and other things done by IDG.
  // Never use more than 16 GB
  const size_t memSize = std::min<uint64_t>(16ul * 1024ul * 1024ul * 1024ul,
                                            _resources.Memory() / 3);
  uint64_t memPerTimestep =
      idg::api::BufferSet::get_memory_per_timestep(nStations, nChannels);
#ifdef HAVE_EVERYBEAM
  aTermMaker = getATermMaker(msData);
  if (aTermMaker) {
    const size_t subgridsize = _bufferset->get_subgridsize();
    aTermBuffer.resize(subgridsize * subgridsize * 4 * nStations);
    // When a-terms are used, they will also take memory. Here we calculate
    // their approx contribution.
    double avgUpdate = aTermMaker->AverageUpdateTime();
    Logger::Debug << "A-terms change on average every " << avgUpdate
                  << " s, once every " << (avgUpdate / msData.integrationTime)
                  << " timesteps.\n";
    const uint64_t atermMemPerTimestep =
        subgridsize * subgridsize * nStations *  // size of grid x nr of grids
        (4 * 8) *  // 4 pol, 8 bytes per complex value
        (msData.integrationTime /
         avgUpdate);  // Average number of aterms per timestep
    Logger::Debug << "A-terms increase mem per timestep from " << memPerTimestep
                  << " bytes to " << (memPerTimestep + atermMemPerTimestep)
                  << " bytes.\n";
    memPerTimestep += atermMemPerTimestep;
  }
#else
  if (!_settings.atermConfigFilename.empty() || _settings.gridWithBeam) {
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
    const MSGridderBase::MSData& msData) {
  SynchronizedMS ms = msData.msProvider->MS();
  size_t nr_stations = ms->antenna().nrow();
  if (!_settings.atermConfigFilename.empty() || _settings.gridWithBeam) {
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
    aterm_settings.save_aterms_prefix = _settings.prefixName;
    aterm_settings.data_column_name = _settings.dataColumnName;
    aterm_settings.filenames = _settings.filenames;
    aterm_settings.aterm_update_interval = _settings.beamAtermUpdateTime;
    aterm_settings.padded_image_width = _settings.paddedImageWidth;
    aterm_settings.trimmed_image_width = _settings.trimmedImageWidth;
    aterm_settings.max_support = _settings.atermKernelSize;

    if (!_settings.atermConfigFilename.empty()) {
      ParsetReader parset_aterms(_settings.atermConfigFilename);
      // If MWA MS and "beam" in aterms, get the path to the coefficient file
      if (everybeam::GetTelescopeType(*ms) ==
          everybeam::TelescopeType::kMWATelescope) {
        std::vector<std::string> aterms = parset_aterms.GetStringList("aterms");
        if (std::find(aterms.begin(), aterms.end(), "beam") != aterms.end()) {
          aterm_settings.coeff_path =
              wsclean::mwa::FindCoeffFile(_settings.mwaPath);
        }
      }

      std::unique_ptr<ATermConfig> config(
          new ATermConfig(nr_stations, system, aterm_settings));
      config->SetSaveATerms(_settings.saveATerms, _settings.prefixName);
      config->Read(*ms, parset_aterms, ms.Filename());
      return config;
    } else {
      // If MWA MS, get the path to the coefficient file
      if (everybeam::GetTelescopeType(*ms) ==
          everybeam::TelescopeType::kMWATelescope) {
        aterm_settings.coeff_path =
            wsclean::mwa::FindCoeffFile(_settings.mwaPath);
      }
      bool frequencyInterpolation = true;
      bool useChannelFrequency = true;
      std::string elementResponseModel =
          !_settings.beamModel.empty() ? _settings.beamModel : "DEFAULT";

      std::unique_ptr<ATermBeam> beam = ATermConfig::GetATermBeam(
          *ms, system, aterm_settings, frequencyInterpolation,
          _settings.beamNormalisationMode, useChannelFrequency,
          elementResponseModel, _settings.beamMode);
      beam->SetSaveATerms(_settings.saveATerms, _settings.prefixName);
      beam->SetUpdateInterval(_settings.beamAtermUpdateTime);
      return beam;
    }
  } else {
    return std::unique_ptr<ATermBase>();
  }
}
#endif  // HAVE_EVERYBEAM
