#ifndef MSPROVIDER_H
#define MSPROVIDER_H

#include "synchronizedms.h"

#include <aocommon/multibanddata.h>
#include <aocommon/polarization.h>

#include <schaapcommon/reordering/storagemanagertype.h>

#include "../structures/msselection.h"

#include <casacore/casa/Arrays/Array.h>
#include <casacore/tables/Tables/ArrayColumn.h>

#include <complex>
#include <set>
#include <type_traits>
#include <vector>

namespace casacore {
class MeasurementSet;
}  // namespace casacore

namespace wsclean {
class MSReader;

/**
 * The abstract MSProvider class is the base class for classes that read and
 * write the visibilities. Write functionality is directly provided by classes
 * that derive from MSProvider, whereas reading functionality is provided in a
 * separate class (MSReader), which can be instantiated with @ref
 * MSProvider::MakeReader. An MSProvider knows which rows are selected and
 * doesn't write to unselected rows. This information on selected rows is also
 * passed to the @ref MSReader. MSProvider provides the visibilities weighted
 * with the visibility weight and converts the visibilities to a requested
 * polarization. The @ref ContiguousMS and
 * @ref ReorderedMs classes implement the MSProvider interface.
 *
 * The class maintains an index for the write position. The index for the
 * reading position is maintained by the closely connected @ref MSReader class.
 * Writing (and reading) goes sequentially through the data.
 *
 * An MS provider provides data for a single dataDescID, and therefore for
 * a single spectral window. Measurement sets with multiple data desc IDs
 * require creating multiple MS providers. Because of this, all rows of a
 * single MS provider have the same number of channels.
 */
class MSProvider {
 public:
  struct MetaData {
    double uInM, vInM, wInM;
    size_t fieldId, antenna1, antenna2;
    double time;
  };

  MSProvider() = default;

  virtual ~MSProvider();

  MSProvider(const MSProvider&) = delete;
  MSProvider& operator=(const MSProvider&) = delete;

  virtual SynchronizedMS MS() = 0;

  /**
   * The column name from which data is read.
   * Writing is done to a different column.
   */
  virtual const std::string& DataColumnName() = 0;

  /**
   * Move the model writing position to the next row.
   */
  virtual void NextOutputRow() = 0;

  /**
   * Reset the writing position to the first row.
   */
  virtual void ResetWritePosition() = 0;

  /**
   * Write model visibilities to the current writing position. If add is true,
   * the provided data are add-assigned to the existing model visibilities. If
   * false, the existing model visibilities are overwritten.
   */
  virtual void WriteModel(const std::complex<float>* buffer, bool add) = 0;

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
   * Be aware that it may return 'DiagonalInstrumental' or 'Instrumental',
   * which means it provides 2 or 4 polarizations that are provided by the
   * underlying measurement set / data.
   */
  virtual aocommon::PolarizationEnum Polarization() = 0;

  /**
   * The dataDescID to which this MS provider is associated. An MSProvider
   * is associated with only one dataDescId. When an MS has multiple
   * dataDescIds, it will have multiple MSProviders.
   */
  virtual size_t DataDescId() = 0;

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
   * Get the band information that this MSProvider covers. The BandData
   * includes all channels in the original data, even when not all
   * channels are selected (@sa NChannels()).
   */
  virtual const aocommon::BandData& Band() = 0;

  /**
   * Get a set of polarizations in the measurement set for a given data desc
   * id. This will always list the individual polarizations, and not return one
   * of the special polarization values (like FullJones or Instrumental).
   */
  static std::set<aocommon::PolarizationEnum> GetMSPolarizations(
      size_t data_desc_id, const casacore::MeasurementSet& ms);

  virtual std::unique_ptr<MSReader> MakeReader() = 0;

  /**
   * Reset model data in the MSProvider to zeros.
   */
  void ResetModelColumn();

  static void GetRowRange(
      casacore::MeasurementSet& ms,
      const schaapcommon::reordering::MSSelection& selection, size_t& startRow,
      size_t& endRow);

  static void GetRowRangeAndIDMap(
      casacore::MeasurementSet& ms,
      const schaapcommon::reordering::MSSelection& selection, size_t& startRow,
      size_t& endRow, const std::set<size_t>& dataDescIdMap,
      std::vector<size_t>& idToMSRow);

  static void CopyRealToComplex(std::complex<float>* dest, const float* source,
                                size_t n) {
    const float* end = source + n;
    while (source != end) {
      *dest = *source;
      ++dest;
      ++source;
    }
  }

  static void InitializeModelColumn(
      casacore::MeasurementSet& ms, const std::string& model_column_name,
      schaapcommon::reordering::StorageManagerType type);

  static casacore::ArrayColumn<float> InitializeImagingWeightColumn(
      casacore::MeasurementSet& ms);

  /**
   * Make an arraycolumn object for the weight spectrum column if it exists and
   * is valid. The weight spectrum column is an optional column, the weight
   * column should be used if it doesn't exist. Moreover, some measurement sets
   * have an empty weight spectrum column; this method only
   * returns true if the column can be used.
   */
  static bool OpenWeightSpectrumColumn(
      const casacore::MeasurementSet& ms,
      std::unique_ptr<casacore::ArrayColumn<float>>& weightColumn);

  static void ExpandScalarWeights(
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
};

}  // namespace wsclean

#endif
