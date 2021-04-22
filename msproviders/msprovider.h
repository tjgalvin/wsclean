#ifndef MSPROVIDER_H
#define MSPROVIDER_H

#include "synchronizedms.h"

#include <casacore/casa/Arrays/Array.h>

#include <casacore/tables/Tables/ArrayColumn.h>

#include <aocommon/polarization.h>

#include <complex>
#include <set>
#include <vector>
#include <type_traits>

namespace casacore {
class MeasurementSet;
}
class MSSelection;

/**
 * The abstract MSProvider class is the base class for classes that read and
 * write the visibilities. An MSProvider knows which rows are selected and
 * doesn't read or write to unselected rows. It provides the visibilities
 * weighted with the visibility weight and converts the visibilities to a
 * requested polarization. The @ref ContiguousMS and
 * @ref PartitionedMS classes implement the MSProvider interface.
 *
 * The class maintains two indices: one for reading and one for writing. Both
 * operations must go sequentially through the data, but reading may be at a
 * different position than writing. This for example allows reading some (meta)
 * data ahead, buffering those, calculate the corresponding visibilities and
 * write them back.
 */
class MSProvider {
 public:
  friend class MSRowProvider;
  friend class DirectMSRowProvider;

  struct MetaData {
    double uInM, vInM, wInM;
    size_t dataDescId, fieldId, antenna1, antenna2;
    double time;
  };

  virtual ~MSProvider() {}

  virtual SynchronizedMS MS() = 0;

  /**
   * The column name from which data is read.
   * Writing is done to a different column.
   */
  virtual const std::string& DataColumnName() = 0;

  /**
   * This provides a unique, consecutive number that corresponds to
   * the current reading position. Note that this number does not have
   * to map directly to measurement set row indices, because unselected
   * data does not affect the RowId. @ref MakeIdToMSRowMapping()
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
   * Move the model writing position to the next row.
   */
  virtual void NextOutputRow() = 0;

  /**
   * Reset both the reading and writing position to the first row.
   */
  virtual void Reset() = 0;

  /**
   * @{
   * Read meta data from the current reading position.
   */
  virtual void ReadMeta(double& u, double& v, double& w,
                        size_t& dataDescId) = 0;

  virtual void ReadMeta(MetaData& metaData) = 0;
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

  virtual void ReadWeights(std::complex<float>* buffer) = 0;

  /**
   * Write model visibilities to the current writing position. If add is true,
   * the provided data are add-assigned to the existing model visibilities. If
   * false, the existing model visibilities are overwritten.
   */
  virtual void WriteModel(const std::complex<float>* buffer, bool add) = 0;

  /**
   * Write imaging weights to the current READING position.
   * Note that despite this is a write operation, the reading position is
   * used nevertheless. This is because it is written while reading the meta
   * data inside WSClean, hence it would be inconvenient if the writing position
   * would be used.
   */
  virtual void WriteImagingWeights(const float* buffer) = 0;

  /**
   * Prepare the msprovider for writing. This is explicitly required before
   * writing to make it possible to image read-only sets when no writing is
   * required.
   */
  virtual void ReopenRW() = 0;

  /**
   * First MJD time of the data.
   */
  virtual double StartTime() = 0;

  /**
   * To obtain a mapping between @ref RowId() and measurement set rows.
   */
  virtual void MakeIdToMSRowMapping(std::vector<size_t>& idToMSRow) = 0;

  /**
   * Polarization that this msprovider provides.
   * Be aware that it may return 'Instrumental', which means it provides
   * all polarizations that are provided by the underlying measurement set /
   * data.
   */
  virtual aocommon::PolarizationEnum Polarization() = 0;

  /**
   * Number of channels provided by this provider. May be different from the
   * underlying measurement set if not all channels are selected.
   */
  virtual size_t NChannels() = 0;

  /**
   * Count of antennas in the underlying measurement set (irrespective of
   * whether all of them are selected). An antenna indix returned by the @ref
   * ReadMeta() methods is always less than this number.
   */
  virtual size_t NAntennas() = 0;

  /**
   * This is the number of polarizations provided by this MSProvider.
   * Note that this does not have to be equal to the nr of pol in the
   * MS, as in most cases each pol is provided by a separate msprovider.
   */
  virtual size_t NPolarizations() = 0;

  /**
   * Get a list of polarizations in the measurement set.
   * This will always list the individual polarizations, and not return
   * one of the special polarization values (like FullJones or Instrumental).
   */
  static std::vector<aocommon::PolarizationEnum> GetMSPolarizations(
      casacore::MeasurementSet& ms);

 protected:
  static void copyData(std::complex<float>* dest, size_t startChannel,
                       size_t endChannel,
                       const std::vector<aocommon::PolarizationEnum>& polsIn,
                       const casacore::Array<std::complex<float>>& data,
                       aocommon::PolarizationEnum polOut);

  template <typename NumType>
  static void copyWeights(NumType* dest, size_t startChannel, size_t endChannel,
                          const std::vector<aocommon::PolarizationEnum>& polsIn,
                          const casacore::Array<std::complex<float>>& data,
                          const casacore::Array<float>& weights,
                          const casacore::Array<bool>& flags,
                          aocommon::PolarizationEnum polOut);

  template <typename NumType>
  static bool isCFinite(const std::complex<NumType>& c) {
    return std::isfinite(c.real()) && std::isfinite(c.imag());
  }

  template <bool add>
  static void reverseCopyData(
      casacore::Array<std::complex<float>>& dest, size_t startChannel,
      size_t endChannel,
      const std::vector<aocommon::PolarizationEnum>& polsDest,
      const std::complex<float>* source, aocommon::PolarizationEnum polSource);

  static void reverseCopyWeights(
      casacore::Array<float>& dest, size_t startChannel, size_t endChannel,
      const std::vector<aocommon::PolarizationEnum>& polsDest,
      const float* source, aocommon::PolarizationEnum polSource);

  static void getRowRange(casacore::MeasurementSet& ms,
                          const MSSelection& selection, size_t& startRow,
                          size_t& endRow);

  static void getRowRangeAndIDMap(casacore::MeasurementSet& ms,
                                  const MSSelection& selection,
                                  size_t& startRow, size_t& endRow,
                                  const std::set<size_t>& dataDescIdMap,
                                  std::vector<size_t>& idToMSRow);

  static void copyRealToComplex(std::complex<float>* dest, const float* source,
                                size_t n) {
    const float* end = source + n;
    while (source != end) {
      *dest = *source;
      ++dest;
      ++source;
    }
  }

  static void initializeModelColumn(casacore::MeasurementSet& ms);

  static casacore::ArrayColumn<float> initializeImagingWeightColumn(
      casacore::MeasurementSet& ms);

  /**
   * Make an arraycolumn object for the weight spectrum column if it exists and
   * is valid. The weight spectrum column is an optional column, the weight
   * column should be used if it doesn't exist. Moreover, some measurement sets
   * have an empty or invalid sized weight spectrum column; this method only
   * returns true if the column can be used.
   */
  static bool openWeightSpectrumColumn(
      casacore::MeasurementSet& ms,
      std::unique_ptr<casacore::ROArrayColumn<float>>& weightColumn,
      const casacore::IPosition& dataColumnShape);

  static void expandScalarWeights(
      const casacore::Array<float>& weightScalarArray,
      casacore::Array<float>& weightSpectrumArray) {
    casacore::Array<float>::const_contiter src = weightScalarArray.cbegin();
    for (casacore::Array<float>::contiter i = weightSpectrumArray.cbegin();
         i != weightSpectrumArray.cend(); ++i) {
      *i = *src;
      ++src;
      if (src == weightScalarArray.cend()) src = weightScalarArray.cbegin();
    }
  }

  MSProvider() {}

 private:
  MSProvider(const MSProvider&) {}
  void operator=(const MSProvider&) {}
};

#endif
