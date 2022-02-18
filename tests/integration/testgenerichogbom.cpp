#include "fftwmanager.h"
#include <aocommon/fits/fitswriter.h>
#include <aocommon/uvector.h>
#include "wsclean/imagingtable.h"
#include "deconvolution/imageset.h"
#include "deconvolution/genericclean.h"
#include "stopwatch.h"

#include <aocommon/logger.h>
#include <aocommon/fits/fitsreader.h>

int main(int argc, char* argv[]) {
  if (argc <= 1) {
    std::cout << "Syntax: testgenerichogbom [-no-clark] <freqcount> <niter> "
                 "{<image> <psf>} [..]\n";
    return -1;
  }

  aocommon::Logger::SetVerbosity(aocommon::Logger::kVerboseVerbosity);

  bool useClark = true;
  int argi = 1;

  while (argi < argc && argv[argi][0] == '-') {
    std::string p(&argv[argi][1]);
    if (p == "no-clark")
      useClark = false;
    else
      throw std::runtime_error("Bad parameters");
    ++argi;
  }

  const size_t freqCount = atoi(argv[argi]);
  const size_t nIter = atoi(argv[argi + 1]);
  argi += 2;

  size_t width = 0, height = 0;
  ImagingTable table;
  for (size_t i = 0; i != freqCount; ++i) {
    ImagingTableEntry& e = table.AddEntry();
    e.index = i;
    e.squaredDeconvolutionIndex = i;
    e.outputChannelIndex = i;
    e.outputIntervalIndex = 0;
    e.imageCount = 1;
    e.joinedGroupIndex = 0;
    e.polarization = aocommon::Polarization::StokesI;
  }
  table.Update();

  Settings settings;
  settings.deconvolutionChannelCount = 0;
  settings.squaredJoins = false;
  settings.linkedPolarizations = std::set<aocommon::PolarizationEnum>();
  ImageSet dirtySet(&table, settings), modelSet(&table, settings);
  std::vector<aocommon::UVector<double>> psfs(freqCount);

  aocommon::FitsWriter writer;

  for (size_t i = 0; i != freqCount; ++i) {
    std::string imageName(argv[argi]), psfName(argv[argi + 1]);
    argi += 2;
    aocommon::FitsReader reader(imageName);
    if (i == 0) {
      writer = aocommon::FitsWriter(reader);
      width = reader.ImageWidth();
      height = reader.ImageHeight();
    }

    aocommon::UVector<double>& psf = psfs[i];
    if (!dirtySet.IsAllocated()) {
      dirtySet.AllocateImages(width, height);
      modelSet.AllocateImages(width, height);
    }

    psf.resize(width * height);

    std::cout << "Reading " << imageName << "...\n";
    reader.Read(dirtySet[i]);

    std::cout << "Reading " << psfName << "...\n";
    reader = aocommon::FitsReader(psfName);
    reader.Read(psf.data());
  }

  modelSet = 0.0;

  FFTWManager fftw;
  GenericClean clean(fftw, useClark);

  aocommon::UVector<const double*> psfVec(psfs.size());
  for (size_t i = 0; i != psfs.size(); ++i) psfVec[i] = psfs[i].data();

  Stopwatch watch(true);

  bool reachedMajorThreshold = false;
  clean.SetMaxNIter(nIter);
  clean.SetMGain(0.8);
  clean.SetCleanBorderRatio(0.0);
  clean.ExecuteMajorIteration(dirtySet, modelSet, psfVec, width, height,
                              reachedMajorThreshold);

  std::cout << "Time taken: " << watch.ToString()
            << ", iterations per second: " << double(nIter) / watch.Seconds()
            << '\n';

  std::cout << "Saving " << modelSet.size() << " output images...\n";
  for (size_t i = 0; i != modelSet.size(); ++i) {
    std::ostringstream sm, sr;
    sm << "model" << i << ".fits";
    writer.Write(sm.str(), modelSet[i]);
    sr << "residual" << i << ".fits";
    writer.Write(sr.str(), dirtySet[i]);
  }

  aocommon::UVector<double> image(width * height);
  dirtySet.GetLinearIntegrated(image.data());
  writer.Write("residual-mfs.fits", image.data());
  modelSet.GetLinearIntegrated(image.data());
  writer.Write("model-mfs.fits", image.data());

  return 0;
}
