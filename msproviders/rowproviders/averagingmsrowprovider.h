#ifndef AVERAGING_MS_ROW_PROVIDER_H
#define AVERAGING_MS_ROW_PROVIDER_H

#include "msrowprovider.h"

#include <stdexcept>

#include <aocommon/uvector.h>
#include <aocommon/vectorset.h>

namespace wsclean {

/**
 * This row provider averages the data in time with BDA: so input is regular data,
 * output is BDA data.
 *
 * Note that this differs from the BdaMsRowProvider, that takes BDA data as input
 * and outputs that.
 */
class AveragingMSRowProvider : public MSRowProvider {
 public:
  AveragingMSRowProvider(double nWavelengthsAveraging, const string& msPath,
                         const schaapcommon::reordering::MSSelection& selection,
                         const aocommon::VectorSet<size_t>& selectedDataDescIds, size_t fieldId,
                         const std::string& dataColumnName, const std::string& modelColumnName,
                         bool requireModel);

  virtual bool AtEnd() const final override { return _flushPosition >= _nElements; }

  virtual void NextRow() final override;

  virtual void ReadData(DataArray& data, FlagArray& flags, WeightArray& weights, double& u,
                        double& v, double& w, uint32_t& dataDescId, uint32_t& antenna1,
                        uint32_t& antenna2, uint32_t& fieldId, double& time,
                        float& interval) final override;

  virtual void ReadModel(DataArray& model) final override;

  virtual void OutputStatistics() const final override;

  virtual WriteInfo GetWriteInfo() final override {
    throw std::runtime_error(
        "Can't extrapolate baseline-dependent averaged data back to MS resolution");
  }

 private:
  class AveragingBuffer {
   public:
    AveragingBuffer()
        : _data(),
          _modelData(),
          _weights(),
          _uvw{0.0, 0.0, 0.0},
          _time(0.0),
          _unweightedTime(0.0),
          _averagedDataCount(0),
          _summedWeight(0) {}

    bool IsInitialized() const { return !_data.empty(); }

    void Initialize(size_t bufferSize, bool includeModel);

    void AddData(size_t n, const std::complex<float>* data, const bool* flags, const float* weights,
                 const double* uvw, double time);

    void AddDataAndModel(size_t n, const std::complex<float>* data,
                         const std::complex<float>* modelData, const bool* flags,
                         const float* weights, const double* uvw, double time);

    size_t AveragedDataCount() const { return _averagedDataCount; }

    void Get(size_t n, std::complex<float>* data, bool* flags, float* weights, double* uvw,
             double& time);

    void Get(size_t n, std::complex<float>* data, std::complex<float>* modelData, bool* flags,
             float* weights, double* uvw, double& time);

    void Reset(size_t n);

   private:
    aocommon::UVector<std::complex<float>> _data;
    aocommon::UVector<std::complex<float>> _modelData;
    aocommon::UVector<float> _weights;
    double _uvw[3];
    double _time;
    // In case all visibilities in an averaging buffer are flagged, it is still
    // desirable to at least return a time that is part of the observation:
    double _unweightedTime;
    size_t _averagedDataCount;
    double _summedWeight;
  };

  bool processCurrentTimestep();

  typedef std::tuple<double, double, double> Pos;

  aocommon::UVector<size_t> _averagingFactors;
  std::vector<AveragingBuffer> _buffers;
  aocommon::VectorMap<size_t> _spwIndexToDataDescId;
  aocommon::VectorMap<size_t> _dataDescIdToSpwIndex;

  size_t _nAntennae;

  DataArray _currentData, _currentModel;
  FlagArray _currentFlags;
  WeightArray _currentWeights;
  size_t _averagedDataDescId;
  size_t _averagedAntenna1Index;
  size_t _averagedAntenna2Index;
  size_t averaged_interval_ = 0;
  size_t _nElements;

  // Once the Measurement Set has completely been read, the buffer might be
  // still full. This value represents the position within that buffer at that
  // point
  size_t _flushPosition;

  // Some statistics
  size_t _averagedRowCount;
  size_t _averageFactorSum;
  size_t _rowCount;
  size_t _fieldId;
  float constant_interval_ = 0.0;
};

}  // namespace wsclean

#endif
