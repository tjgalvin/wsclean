#ifndef MSPROVIDERS_MSWEIGHTCOLUMN_H
#define MSPROVIDERS_MSWEIGHTCOLUMN_H

#include "msrowproviderbase.h"

#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>

#include <string>

namespace casacore {
class MeasurementSet;
}  // namespace casacore

/**
 * Helper to get the Measurement set's weight.
 *
 * When available it uses the measurement set's WEIGHT_SPECTRUM column else its
 * WEIGHT column.
 */
class MsWeightColumn {
 public:
  explicit MsWeightColumn(const casacore::MeasurementSet& ms);

  /**
   * Read the weight.
   *
   * @param [out] weights Returns the result.
   * @param row The row whose information to return.
   * @param shape Sets the proper shape of @a weights when the data is
   *              retrieved from the WEIGHT column.
   */
  void ReadData(MsRowProviderBase::WeightArray& weights, size_t row,
                const casacore::IPosition& shape) const;

 private:
  std::unique_ptr<casacore::ArrayColumn<float>> weight_column_;
  enum class Type { kSpectrum, kScalar };
  Type type_{Type::kSpectrum};
};

#endif  // MSPROVIDERS_MSWEIGHTCOLUMN_H
