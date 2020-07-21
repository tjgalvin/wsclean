#ifndef WSCLEAN_H
#define WSCLEAN_H

#include "../msselection.h"
#include <aocommon/polarization.h>
#include "../weightmode.h"
#include "../stopwatch.h"

#include "../deconvolution/deconvolution.h"

#include "../scheduling/griddingresult.h"

#include "cachedimageset.h"
#include "imagingtable.h"
#include "msgridderbase.h"
#include "observationinfo.h"
#include "outputchannelinfo.h"
#include "wscfitswriter.h"
#include "wscleansettings.h"

#include <set>

class WSClean {
 public:
  WSClean();
  ~WSClean();

  WSCleanSettings& Settings() { return _settings; }
  const WSCleanSettings& Settings() const { return _settings; }

  void SetCommandLine(const std::string& cmdLine) { _commandLine = cmdLine; }

  void RunClean();

  void RunPredict();

 private:
  void runIndependentGroup(ImagingTable& groupTable,
                           std::unique_ptr<class PrimaryBeam>& primaryBeam);
  void saveRestoredImagesForGroup(
      const ImagingTableEntry& tableEntry,
      std::unique_ptr<class PrimaryBeam>& primaryBeam) const;
  void predictGroup(const ImagingTable& imagingGroup);

  void runFirstInversion(ImagingTableEntry& entry,
                         std::unique_ptr<class PrimaryBeam>& primaryBeam);

  void performReordering(bool isPredictMode);

  std::shared_ptr<ImageWeights> initializeImageWeights(
      const ImagingTableEntry& entry,
      std::vector<std::unique_ptr<class MSDataDescription>>& msList);
  void initializeMFSImageWeights();
  void initializeMSList(
      const ImagingTableEntry& entry,
      std::vector<std::unique_ptr<MSDataDescription>>& msList);
  void storeAndCombineXYandYX(CachedImageSet& dest,
                              aocommon::PolarizationEnum polarization,
                              size_t joinedChannelIndex, bool isImaginary,
                              const double* image);
  bool selectChannels(MSSelection& selection, size_t msIndex, size_t bandIndex,
                      const ImagingTableEntry& entry);
  MSSelection selectInterval(MSSelection& fullSelection, size_t intervalIndex);
  void readEarlierModelImages(const ImagingTableEntry& entry);

  void makeImagingTable(size_t outputIntervalIndex);
  void makeImagingTableEntry(const std::vector<aocommon::ChannelInfo>& channels,
                             size_t outIntervalIndex, size_t outChannelIndex,
                             ImagingTableEntry& entry);
  void makeImagingTableEntryChannelSettings(
      const std::vector<aocommon::ChannelInfo>& channels,
      size_t outIntervalIndex, size_t outChannelIndex, size_t nOutChannels,
      ImagingTableEntry& entry);
  void addPolarizationsToImagingTable(size_t& joinedGroupIndex,
                                      size_t& squaredGroupIndex,
                                      size_t outChannelIndex,
                                      const ImagingTableEntry& templateEntry);
  std::unique_ptr<class ImageWeightCache> createWeightCache();

  void multiplyImage(double factor, double* image) const;
  void multiplyImage(double factor, Image& image) const {
    multiplyImage(factor, image.data());
  }

  GriddingResult loadExistingImage(ImagingTableEntry& entry, bool isPSF);
  void loadExistingPSF(ImagingTableEntry& entry);
  void loadExistingDirty(ImagingTableEntry& entry, bool updateBeamInfo);

  void imagePSF(ImagingTableEntry& entry);
  void imagePSFCallback(ImagingTableEntry& entry,
                        struct GriddingResult& result);

  void imageMain(ImagingTableEntry& entry, bool isFirstInversion,
                 bool updateBeamInfo);
  void imageMainCallback(ImagingTableEntry& entry,
                         struct GriddingResult& result, bool updateBeamInfo,
                         bool isInitialInversion);

  void predict(const ImagingTableEntry& entry);
  void predictCallback(const ImagingTableEntry& entry,
                       struct GriddingResult& result);

  // void makeMFSImage(const string& suffix, size_t intervalIndex,
  // aocommon::PolarizationEnum pol, bool isImaginary, bool isPSF = false); void
  // renderMFSImage(size_t intervalIndex, aocommon::PolarizationEnum pol, bool
  // isImaginary, bool isPBCorrected) const;
  void saveUVImage(const double* image, const ImagingTableEntry& entry,
                   bool isImaginary, const std::string& prefix) const;
  void writeFirstResidualImages(const ImagingTable& groupTable) const;
  void writeModelImages(const ImagingTable& groupTable) const;

  double minTheoreticalBeamSize(const ImagingTable& table) const;

  void makeBeam();

  WSCFitsWriter createWSCFitsWriter(const ImagingTableEntry& entry,
                                    bool isImaginary, bool isModel) const;

  WSCFitsWriter createWSCFitsWriter(const ImagingTableEntry& entry,
                                    aocommon::PolarizationEnum polarization,
                                    bool isImaginary, bool isModel) const;

  MSSelection _globalSelection;
  std::string _commandLine;

  WSCleanSettings _settings;

  std::vector<OutputChannelInfo> _infoPerChannel;
  OutputChannelInfo _infoForMFS;
  std::map<size_t, std::unique_ptr<MetaDataCache>> _msGridderMetaCache;

  std::unique_ptr<class GriddingTaskManager> _griddingTaskManager;
  std::unique_ptr<class ImageWeightCache> _imageWeightCache;
  Stopwatch _inversionWatch, _predictingWatch, _deconvolutionWatch;
  bool _isFirstInversion;
  size_t _majorIterationNr;
  CachedImageSet _psfImages, _modelImages, _residualImages;
  std::vector<PartitionedMS::Handle> _partitionedMSHandles;
  std::vector<MultiBandData> _msBands;
  Deconvolution _deconvolution;
  ImagingTable _imagingTable;
  ObservationInfo _observationInfo;
};

#endif
