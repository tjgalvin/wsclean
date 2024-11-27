#ifndef WSCLEAN_H
#define WSCLEAN_H

#include <optional>
#include <set>

#include <aocommon/image.h>
#include <aocommon/fits/fitsreader.h>
#include <aocommon/fits/fitswriter.h>
#include <aocommon/multibanddata.h>
#include <aocommon/polarization.h>

#include <schaapcommon/facets/facet.h>

#include <radler/radler.h>

#include "../scheduling/griddingresult.h"
#include "../scheduling/griddingtaskfactory.h"

#include "../io/cachedimageset.h"
#include "../io/wscfitswriter.h"

#include "../structures/imagingtable.h"
#include "../structures/msselection.h"
#include "../structures/observationinfo.h"
#include "../structures/outputchannelinfo.h"
#include "../structures/weightmode.h"

#include "../gridding/msgridder.h"

#include "../msproviders/reorderedmsprovider.h"

#include "imageweightinitializer.h"
#include "mshelper.h"
#include "stopwatch.h"
#include "settings.h"

namespace schaapcommon::facets {
class FacetImage;
}  // namespace schaapcommon::facets

namespace wsclean {

class ImageWeightCache;
class PrimaryBeam;

class WSClean {
 public:
  WSClean();
  ~WSClean();

  Settings& GetSettings() { return _settings; }
  const Settings& GetSettings() const { return _settings; }
  void ResetSettings() { _settings = Settings(); }
  void SetCommandLine(const std::string& cmdLine) { _commandLine = cmdLine; }

  void RunClean();

  /**
   * Entry point for performing a single prediction for an existing model image.
   *
   * In case of a facet-based prediction, the provided model images are assumed
   * to have the same size, so that the image size of the full image can be
   * inferred from the first entry in the _imagingTable in an early stage.
   */
  void RunPredict();

 private:
  void runIndependentGroup(ImagingTable& groupTable,
                           std::unique_ptr<PrimaryBeam>& primaryBeam);
  void saveRestoredImagesForGroup(
      const ImagingTable::Group& group,
      std::unique_ptr<PrimaryBeam>& primaryBeam) const;
  void predictGroup(const ImagingTable& groupTable);

  void runFirstInversions(ImagingTable& groupTable,
                          std::unique_ptr<PrimaryBeam>& primaryBeam,
                          bool requestPolarizationsAtOnce,
                          bool parallelizePolarizations);
  /**
   * @brief Run first inversion on all entries within a group.
   * @details A group should contain all facets of a single image.
   */
  void runFirstInversionGroup(ImagingTable::Group& facetGroup,
                              std::unique_ptr<PrimaryBeam>& primaryBeam);
  void runMajorIterations(ImagingTable& groupTable,
                          std::unique_ptr<PrimaryBeam>& primaryBeam,
                          bool requestPolarizationsAtOnce,
                          bool parallelizePolarizations);

  /**
   * Returns true when gridding is done with a-terms. This can either
   * be enabled by setting the gridWithBeam setting to true or by providing
   * an aterm config file. */
  bool griddingUsesATerms() const {
    return _settings.gridWithBeam || !_settings.atermConfigFilename.empty();
  }

  /**
   * True when the imaging uses any of the methods to apply a beam.
   * A beam can be applied through facetting (with solutions or beam),
   * through gridding with the beam using IDG or by correcting for the beam
   * in image space after imaging.
   */
  bool usesBeam() const {
    return _settings.applyPrimaryBeam || _settings.applyFacetBeam ||
           !_settings.facetSolutionFiles.empty() || griddingUsesATerms();
  }

  std::string PredictModelFileSuffix() const {
    return _settings.applyFacetBeam
               ? "-model-fpb.fits"
               : ((_settings.UseFacetCorrections() || griddingUsesATerms())
                      ? "-model-pb.fits"
                      : "-model.fits");
  }

  ObservationInfo getObservationInfo() const;
  std::pair<double, double> getLMShift() const;

  void resetModelColumns(const ImagingTable::Groups& facet_groups);
  void resetModelColumns(const ImagingTableEntry& entry);
  void storeAndCombineXYandYX(CachedImageSet& dest, size_t joinedChannelIndex,
                              const ImagingTableEntry& entry,
                              aocommon::PolarizationEnum polarization,
                              bool isImaginary, const aocommon::Image& image);
  MSSelection selectInterval(MSSelection& fullSelection, size_t intervalIndex);

  void makeImagingTable(size_t outputIntervalIndex);
  void makeImagingTableEntry(const std::vector<aocommon::ChannelInfo>& channels,
                             size_t outIntervalIndex, size_t outChannelIndex,
                             ImagingTableEntry& entry);
  void makeImagingTableEntryChannelSettings(
      const std::vector<aocommon::ChannelInfo>& channels,
      size_t outIntervalIndex, size_t outChannelIndex, size_t nOutChannels,
      ImagingTableEntry& entry);
  void addPolarizationsToImagingTable(ImagingTableEntry& templateEntry);
  void addFacetsToImagingTable(ImagingTableEntry& templateEntry,
                               const size_t facet_count);
  void updateFacetsInImagingTable(
      const std::vector<std::shared_ptr<schaapcommon::facets::Facet>>& facets,
      bool updateDdPsfs);
  std::unique_ptr<ImageWeightCache> createWeightCache();

  /**
   * Initializes full-size model images for the given entry. Depending on the
   * settings, this might load existing images from disk or initialize
   * them to zero.
   */
  void initializeModelImages(const ImagingTableEntry& entry,
                             aocommon::PolarizationEnum polarization,
                             size_t nFacetGroups);
  void readExistingModelImages(const ImagingTableEntry& entry,
                               aocommon::PolarizationEnum polarization,
                               size_t nFacetGroups);
  /**
   * Override the image settings given a FitsReader object.
   * The boolean return value indicates whether the gridder needs
   * to be reset.
   */
  bool overrideImageSettings(const aocommon::FitsReader& reader);
  void loadExistingImage(ImagingTableEntry& entry, bool isPSF);

  void ImagePsf(ImagingTable::Group&& facet_group);
  void ImagePsfCallback(ImagingTable::Group facet_group,
                        GriddingResult& result);

  void ImageMain(ImagingTable::Group& facet_group, bool is_first_inversion,
                 bool update_beam_info);
  void ImageMainCallback(ImagingTable::Group facet_group,
                         GriddingResult& result, bool update_beam_info,
                         bool is_first_inversion);

  void Predict(const ImagingTable::Group& facet_group);

  void saveUVImage(const aocommon::Image& image, const ImagingTableEntry& entry,
                   bool isImaginary, const std::string& prefix) const;

  void saveUVImage(const aocommon::Image& image, const ImagingTableEntry& entry,
                   const OutputChannelInfo& channel_info, bool isImaginary,
                   const std::string& prefix) const;

  void processFullPSF(aocommon::Image& image, const ImagingTableEntry& entry);

  void ApplyFacetCorrectionForSingleChannel(const ImagingTable& squared_group,
                                            CachedImageSet& image_cache);

  /**
   * @brief Stitch facets for all FacetGroups
   */
  void stitchFacets(const ImagingTable& table, CachedImageSet& image_cache,
                    bool write_dirty, bool is_psf, bool is_facet_pb_model);

  /**
   * Stitch facet for a single (Facet)Group
   * @param weight_image weight image pointer that should be either empty
   * or should be an image with the right size. This can be used to reuse
   * the same weight image over multiple calls and prevent re-allocation.
   */
  void stitchSingleGroup(const ImagingTable::Group& facetGroup,
                         size_t imageIndex, CachedImageSet& imageCache,
                         bool writeDirty, bool isPSF,
                         aocommon::Image& fullImage,
                         std::unique_ptr<aocommon::Image>& weight_image,
                         schaapcommon::facets::FacetImage& facetImage,
                         size_t maxFacetGroupIndex, bool apply_scalar,
                         bool is_facet_pb_model);
  /**
   * Partition model image into facets and save them into fits files
   */
  void partitionModelIntoFacets(const ImagingTable::Groups& facetGroups,
                                bool isPredictOnly);

  /**
   * Partition image into facets for a single (Facet)Group
   */
  void partitionSingleGroup(const ImagingTable::Group& facetGroup,
                            size_t imageIndex, CachedImageSet& imageCache,
                            const aocommon::Image& fullImage,
                            schaapcommon::facets::FacetImage& facetImage,
                            bool isPredictOnly);

  void writeFirstResidualImages(const ImagingTable& groupTable) const;
  void writeModelImages(const ImagingTable& groupTable) const;

  double minTheoreticalBeamSize(const ImagingTable& table) const;

  void makeBeam();

  WSCFitsWriter createWSCFitsWriter(const ImagingTableEntry& entry,
                                    const OutputChannelInfo& channel_info,
                                    bool isImaginary, bool isModel) const;

  WSCFitsWriter createWSCFitsWriter(const ImagingTableEntry& entry,
                                    bool isImaginary, bool isModel) const;

  WSCFitsWriter createWSCFitsWriter(const ImagingTableEntry& entry,
                                    aocommon::PolarizationEnum polarization,
                                    bool isImaginary, bool isModel) const;
  /**
   * @brief Apply the H5 solution to the (restored) image and save as -pb.fits
   * file. Method is only invoked in case no beam corrections are applied.
   */
  void correctImagesH5(aocommon::FitsWriter& writer,
                       const ImagingTable::Group& group,
                       const ImageFilename& imageName,
                       const std::string& filenameKind) const;

  void storeAverageBeam(const ImagingTableEntry& entry,
                        std::unique_ptr<AverageBeam>& averageBeam);

  /**
   * @brief Compute the total amount of MSProviders that will be generated.
   * This number is needed to initialize the writer locks in the prediction
   * tasks, which are set via a call to _griddingTaskManager->Start(). The
   * number of @p MSProviders is the acummulated number of bands per MS.
   *
   * @return size_t Number of MSProviders
   */
  size_t getMaxNrMSProviders() const {
    size_t msCount = 0;
    for (const auto& msBand : _msBands) msCount += msBand.DataDescCount();
    return msCount;
  }

  MSSelection _globalSelection;
  std::string _commandLine;

  Settings _settings;

  std::vector<OutputChannelInfo> _infoPerChannel;
  OutputChannelInfo _infoForMFS;

  std::unique_ptr<MsHelper> _msHelper;
  std::unique_ptr<ImageWeightInitializer> _image_weight_initializer;
  std::unique_ptr<GriddingTaskFactory> _griddingTaskFactory;
  std::unique_ptr<GriddingTaskManager> _griddingTaskManager;
  std::unique_ptr<ImageWeightCache> _imageWeightCache;
  Stopwatch _inversionWatch, _predictingWatch, _deconvolutionWatch;
  bool _isFirstInversionTask;  // Becomes false after the first inversion task.
  size_t _majorIterationNr;
  CachedImageSet _psfImages;
  CachedImageSet _modelImages;
  CachedImageSet _residualImages;
  CachedImageSet _scalarBeamImages;
  CachedImageSet _matrixBeamImages;
  std::vector<aocommon::MultiBandData> _msBands;
  // Radler object only needed in RunClean runs.
  std::optional<radler::Radler> _deconvolution;
  ImagingTable _imagingTable;
  ObservationInfo _observationInfo;
  std::size_t _facetCount;  // 0 means facets are not used.
  std::size_t _ddPsfCount;  // 0 means dd-psfs are not used.
  /// These contain the user-requested image shift values converted from ra,dec
  /// to l,m units
  /// @{
  double _l_shift;
  double _m_shift;
  /// @}

  double _lastStartTime;
};

}  // namespace wsclean

#endif
