#ifndef MSPROVIDERS_MSREADERS_CONTIGUOUSMSREADER_H
#define MSPROVIDERS_MSREADERS_CONTIGUOUSMSREADER_H

#include "msreader.h"

class ContiguousMS;

class ContiguousMSReader final : public MSReader {
 public:
  explicit ContiguousMSReader(ContiguousMS* contiguousMS);

  virtual ~ContiguousMSReader(){};

  size_t RowId() const override { return _currentRowId; }

  bool CurrentRowAvailable() override;

  void NextInputRow() override;

  void ReadMeta(double& u, double& v, double& w) override;

  void ReadMeta(MSProvider::MetaData& metaData) override;

  void ReadData(std::complex<float>* buffer) override;

  void ReadModel(std::complex<float>* buffer) override;

  void ReadWeights(std::complex<float>* buffer) override;

  void ReadWeights(float* buffer) override;

  void WriteImagingWeights(const float* buffer) override;

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
