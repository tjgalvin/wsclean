#ifndef MSPROVIDERS_MSREADERS_CONTIGUOUSMSREADER_H
#define MSPROVIDERS_MSREADERS_CONTIGUOUSMSREADER_H

#include "msreader.h"

class ContiguousMS;

class ContiguousMSReader final : public MSReader {
 public:
  explicit ContiguousMSReader(ContiguousMS* contiguousMS);

  virtual ~ContiguousMSReader(){};

  size_t RowId() const final override { return _currentRowId; }

  bool CurrentRowAvailable() final override;

  void NextInputRow() final override;

  void ReadMeta(double& u, double& v, double& w,
                size_t& dataDescId) final override;

  void ReadMeta(MSProvider::MetaData& metaData) final override;

  void ReadData(std::complex<float>* buffer) final override;

  void ReadModel(std::complex<float>* buffer) final override;

  void ReadWeights(std::complex<float>* buffer) final override;

  void ReadWeights(float* buffer) final override;

  void WriteImagingWeights(const float* buffer) final override;

 private:
  size_t _currentInputRow;
  size_t _currentInputTimestep;
  double _currentInputTime;
  size_t _currentRowId;

  bool _isDataRead, _isModelRead, _isWeightRead;
  std::unique_ptr<casacore::ArrayColumn<float>> _imagingWeightsColumn;

  void readData();

  void readWeights();

  void readModel();
};

#endif