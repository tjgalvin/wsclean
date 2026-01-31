#ifndef WSCLEAN_MSPROVIDERS_OPERATIONS_H_
#define WSCLEAN_MSPROVIDERS_OPERATIONS_H_

#include <memory>
#include <string>

#include <aocommon/polarization.h>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/Tables/ArrayColumn.h>

#include <schaapcommon/reordering/msselection.h>
#include <schaapcommon/reordering/storagemanagertype.h>

namespace wsclean {

/**
 * Make an arraycolumn object for the weight spectrum column if it exists and
 * is valid. The weight spectrum column is an optional column, the weight
 * column should be used if it doesn't exist. Moreover, some measurement sets
 * have an empty weight spectrum column; this method only
 * returns true if the column can be used.
 */
bool OpenWeightSpectrumColumn(
    const casacore::MeasurementSet& ms,
    std::unique_ptr<casacore::ArrayColumn<float>>& weightColumn);

inline void ExpandScalarWeights(
    const casacore::Array<float>& weight_scalar_array,
    casacore::Array<float>& weight_spectrum_array) {
  casacore::Array<float>::const_contiter src = weight_scalar_array.cbegin();
  for (casacore::Array<float>::contiter i = weight_spectrum_array.cbegin();
       i != weight_spectrum_array.cend(); ++i) {
    *i = *src;
    ++src;
    if (src == weight_scalar_array.cend()) src = weight_scalar_array.cbegin();
  }
}

void GetRowRange(casacore::MeasurementSet& ms,
                 const schaapcommon::reordering::MSSelection& selection,
                 size_t& start_row, size_t& end_row);

/**
 * Get a set of polarizations in the measurement set for a given data desc
 * id. This will always list the individual polarizations, and not return one
 * of the special polarization values (like FullJones or Instrumental).
 */
std::set<aocommon::PolarizationEnum> GetMSPolarizations(
    size_t data_desc_id, const casacore::MeasurementSet& ms);

void InitializeModelColumn(casacore::MeasurementSet& ms,
                           const std::string& data_column_name,
                           const std::string& model_column_name,
                           schaapcommon::reordering::StorageManagerType type);

casacore::ArrayColumn<float> InitializeImagingWeightColumn(
    casacore::MeasurementSet& ms);

double GetMsInterval(const casacore::MeasurementSet& ms);

std::vector<std::string> GetAntennaNames(const casacore::MSAntenna& antenna);

}  // namespace wsclean

#endif
