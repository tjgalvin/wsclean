#ifndef PRIMARY_BEAM_H
#define PRIMARY_BEAM_H

#include <string>

#include "imagefilename.h"
#include "imagingtable.h"
#include "primarybeamimageset.h"

#include <aocommon/polarization.h>

#include "../aterms/telescope.h"

class PrimaryBeam {
 public:
  PrimaryBeam(const class WSCleanSettings& settings);
  ~PrimaryBeam();

  void SetPhaseCentre(double ra, double dec, double dl, double dm) {
    _phaseCentreRA = ra;
    _phaseCentreDec = dec;
    _phaseCentreDL = dl;
    _phaseCentreDM = dm;
  }

  void CorrectImages(const ImageFilename& imageName,
                     std::vector<double*>& images) {
    PrimaryBeamImageSet beamImages = load(imageName, _settings);
    if (_settings.polarizations.size() == 1 &&
        *_settings.polarizations.begin() == aocommon::Polarization::StokesI) {
      beamImages.ApplyStokesI(images[0]);
    } else if (_settings.polarizations.size() == 4 &&
               aocommon::Polarization::HasFullStokesPolarization(
                   _settings.polarizations)) {
      beamImages.ApplyFullStokes(images.data());
    }
  }

  PrimaryBeamImageSet Load(const ImageFilename& imageName) {
    return load(imageName, _settings);
  }

  void AddMS(std::unique_ptr<class MSDataDescription> description);

  void MakeBeamImages(const ImageFilename& imageName,
                      const ImagingTableEntry& entry,
                      std::shared_ptr<class ImageWeights> imageWeights);

  void CorrectImages(class FitsWriter& writer, const ImageFilename& imageName,
                     const std::string& filenameKind);

 private:
  const WSCleanSettings& _settings;
  std::vector<std::unique_ptr<class MSDataDescription>> _msList;
  double _phaseCentreRA, _phaseCentreDec, _phaseCentreDL, _phaseCentreDM;

  static PrimaryBeamImageSet load(const ImageFilename& imageName,
                                  const WSCleanSettings& settings);

  PrimaryBeamImageSet makeLOFARImage(
      const ImagingTableEntry& entry,
      std::shared_ptr<class ImageWeights> imageWeights);

  void makeMWAImage(PrimaryBeamImageSet& beamImages,
                    const ImagingTableEntry& entry);

  void makeATCAImage(PrimaryBeamImageSet& beamImages,
                     const ImagingTableEntry& entry);

  void makeVLAImage(PrimaryBeamImageSet& beamImages,
                    const ImagingTableEntry& entry);

  void makeFromVoltagePattern(PrimaryBeamImageSet& beamImages,
                              const ImagingTableEntry& entry,
                              class VoltagePattern& vp);
};

#endif
