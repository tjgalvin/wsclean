#ifndef MSPROVIDERS_MSREADERS_MSREADER_H
#define MSPROVIDERS_MSREADERS_MSREADER_H

#include "../msprovider.h"

/**
 * The abstract MSReader class is the base class for classes that read
 * visibilities. Derived classes are usually instantiated via
 * MSProvider::MakeReader(). This class is designed such that each instance of
 * a reader provides independent read access to the underlying visibilities.
 *
 * This class maintains a reading position, that goes sequentially through the
 * data. The interface of this class is implemented in @ref ContiguousMSReader
 * and @ref PartitionedMSReader.
 */
class MSReader {
 public:
  MSReader(MSProvider* msProvider) : _msProvider(msProvider){};

  virtual ~MSReader(){};

  /**
   * This provides a unique, consecutive number that corresponds to
   * the current reading position. Note that this number does not have
   * to map directly to measurement set row indices, because unselected
   * data does not affect the RowId. @ref MSProvider::MakeIdToMSRowMapping()
   * can be used to convert this Id to a measurement row number.
   */
  virtual size_t RowId() const = 0;

  /**
   * Returns true as long as there is more data available for reading.
   */
  virtual bool CurrentRowAvailable() = 0;

  /**
   * Move the reading position to the next row.
   */
  virtual void NextInputRow() = 0;

  /**
   * @{
   * Read meta data from the current reading position.
   */
  virtual void ReadMeta(double& u, double& v, double& w) = 0;

  virtual void ReadMeta(MSProvider::MetaData& metaData) = 0;
  /** @} */

  /**
   * Read visibility data from current reading position.
   */
  virtual void ReadData(std::complex<float>* buffer) = 0;

  /**
   * Read the model visibilities from the current reading position.
   */
  virtual void ReadModel(std::complex<float>* buffer) = 0;

  virtual void ReadWeights(float* buffer) = 0;

  /**
   * Write imaging weights to the current READING position.
   * Note that despite this is a write operation, the reading position is
   * used nevertheless. This is because it is written while reading the meta
   * data inside WSClean, hence it would be inconvenient if the writing position
   * would be used.
   */
  virtual void WriteImagingWeights(const float* buffer) = 0;

 protected:
  MSProvider* _msProvider;
};

#endif
