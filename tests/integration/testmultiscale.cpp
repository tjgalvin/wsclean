#include "fftwmanager.h"
#include "fitsreader.h"
#include "fitswriter.h"
#include <aocommon/uvector.h>
#include "stopwatch.h"
#include "wsclean/imagingtable.h"
#include "deconvolution/imageset.h"
#include "multiscale/multiscalealgorithm.h"

int main(int argc, char* argv[]) {
  if (argc <= 1) {
    std::cout
        << "Syntax: testmultiscale <freqcount> <niter> {<image> <psf>} [..]\n";
  }
  const size_t freqCount = atoi(argv[1]);
  const size_t nIter = atoi(argv[2]);

  int argi = 3;

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

  WSCleanSettings settings;
  settings.deconvolutionChannelCount = 0;
  settings.squaredJoins = false;
  settings.linkedPolarizations = std::set<aocommon::PolarizationEnum>();
  ImageSet dirtySet(&table, settings), modelSet(&table, settings);
  std::vector<aocommon::UVector<double>> psfs(freqCount);

  FitsWriter writer;
  double beamSize = 0.0, pixelScaleX = 0.0, pixelScaleY = 0.0;

  for (size_t i = 0; i != freqCount; ++i) {
    std::string imageName(argv[argi]), psfName(argv[argi + 1]);
    argi += 2;
    FitsReader reader(imageName);
    if (i == 0) {
      writer = FitsWriter(reader);
      width = reader.ImageWidth();
      height = reader.ImageHeight();
      beamSize = reader.BeamMinorAxisRad();
      pixelScaleX = reader.PixelSizeX();
      pixelScaleY = reader.PixelSizeY();
      std::cout << "Beam size: " << beamSize / pixelScaleX << " pixels.\n";
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
    reader = FitsReader(psfName);
    reader.Read(psf.data());
  }

  double* model(modelSet[0]);

  // bool hasInitialModel = false;
  if (argc > argi) {
    std::string modelName(argv[4]);
    std::cout << "Reading " << modelName << "...\n";
    FitsReader reader(modelName);
    reader.Read(model);
    // hasInitialModel = true;
    for (size_t i = 1; i < freqCount; ++i) {
      for (size_t j = 0; j != width * height; ++j) modelSet[i][j] = model[j];
    }
  } else {
    std::cout << "Using empty model.\n";
    modelSet = 0.0;
  }

  FFTWManager fftwManager;
  MultiScaleAlgorithm multiscale(fftwManager, beamSize, pixelScaleX,
                                 pixelScaleY);

  /*if(hasInitialModel)
  {
          for(size_t i=0; i!=freqCount; ++i)
          {
                  aocommon::UVector<double> psfKernel(psfs[0].size());
                  FFTConvolver::PrepareKernel(psfKernel.data(), psfs[i].data(),
  width, height);

                  // Calculate: residual = dirty - model (x) psf
                  aocommon::UVector<double> tmp(modelSet[i],
  modelSet[i]+width*height); FFTConvolver::ConvolveSameSize(tmp.data(),
  psfKernel.data(), width, height);

                  // residual = residual - scratch
                  algorithm.Subtract(dirtySet[i], tmp);
          }
  }*/

  aocommon::UVector<const double*> psfVec(psfs.size());
  for (size_t i = 0; i != psfs.size(); ++i) psfVec[i] = psfs[i].data();

  Stopwatch watch(true);

  bool reachedMajorThreshold = false;
  multiscale.SetMaxNIter(nIter);
  multiscale.SetCleanBorderRatio(0.0);
  multiscale.ExecuteMajorIteration(dirtySet, modelSet, psfVec, width, height,
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
